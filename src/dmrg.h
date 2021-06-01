#ifndef CHAINDMRG_DMRG_H
#define CHAINDMRG_DMRG_H
#include <filesystem>
#include <variant>
#include "itensor/all.h"
#include "chain.h"
using namespace itensor;


// Reads and writes data coming out of DMRG
template<typename SiteSetType>
class DMRGProgress {
public:
    SiteSetType sites_;
    // Each element is the latest progress for a given state
    std::vector<int> times_swepts_;
    std::vector<int> maxdims_;
    std::vector<Real> ens_;
    std::vector<MPS> psis_;
    std::vector<MPO> Hs_;

    void SetSites(SiteSetType sites){
        sites_ = sites;
    }

    SiteSetType Sites(){
        return sites_;
    }

    // Replace the latest progress for the current state
    void Update(int times_swept, int maxdim, Real en, MPS const& psi, MPO H) {
        times_swepts_.back() = times_swept;
        maxdims_.back() = maxdim;
        ens_.back() = en;
        psis_.back() = psi;
        Hs_.back() = std::move(H);
    }

    // Next element corresponding to the NextState state
    void StateComplete() {
        times_swepts_.push_back(0);
        maxdims_.push_back(0);
        ens_.push_back(0);
        psis_.emplace_back();
        Hs_.emplace_back();
    }

    std::tuple<int,int,Real,MPS,MPO> All() {
        return std::tuple<int,int,Real,MPS,MPO>(times_swepts_.back(), maxdims_.back(), ens_.back(), psis_.back(), Hs_.back());
    }

    std::vector<MPS> PastStates() {
        std::vector<MPS> past_psis;
        past_psis.reserve((int) psis_.size() - 1);
        for (int i=0;i< (int) psis_.size() - 1; i++){
            past_psis.push_back(psis_.at(i));
        }
        return past_psis;
    }

    std::vector<MPS> States() {
        return psis_;
    }

    std::vector<Real> Energies() const {
        return ens_;
    }

    int NumPastStates() {
        return (int) psis_.size() - 1;
    }

//    int NumStates() {
//        return (int) psis_.size();
//    }

    void Write(const std::filesystem::path& p) const
    {
        std::string path = std::string(p);
        writeToFile(path + ".st", sites_);
        writeToFile(path + ".pgs", times_swepts_);
        writeToFile(path + ".md", maxdims_);
        writeToFile(path + ".en", ens_);
        writeToFile(path + ".psi", psis_);
        writeToFile(path + ".H", Hs_);

        // Uncomment to save backup files
//        writeToFile(path + ".pgs.bk", times_swepts_);
//        writeToFile(path + ".md.bk", maxdims_);
//        writeToFile(path + ".en.bk", ens_);
//        writeToFile(path + ".psi.bk", psis_);
//        writeToFile(path + ".H.bk", Hs_);
    }

    void Read(const std::filesystem::path& p)
    {
        std::string path = std::string(p);
        try{
            readFromFile(path + ".st", sites_);
            readFromFile(path + ".pgs", times_swepts_);
            readFromFile(path + ".md", maxdims_);
            readFromFile(path + ".en", ens_);
            readFromFile(path + ".psi", psis_);
            readFromFile(path + ".H", Hs_);
        }catch(std::exception const& e){
            println(e.what());
//            try{
//                readFromFile(path + ".pgs.bk", times_swepts_);
//                readFromFile(path + ".md.bk", maxdims_);
//                readFromFile(path + ".en.bk", ens_);
//                readFromFile(path + ".psi.bk", psis_);
//                readFromFile(path + ".H.bk", Hs_);
//            }catch(std::exception const& e){
                println(e.what());
                times_swepts_.clear();
                maxdims_.clear();
                ens_.clear();
                psis_.clear();
                Hs_.clear();
                StateComplete();
//            }
        }

    }
};

// Examples of SiteSetType are Golden, Haagerup, HaagerupQ
template<typename SiteSetType>
class DMRG {
    std::string site_name_; // eg. golden, haagerup
    std::string boundary_condition_; // "p" for periodic, "o" for open, "s" for sine-squared deformed
    std::string job_; // site_name_ + boundary_condition_

    // Chain parameters
    int num_sites_;
    Real theta_; // angle between K (pure identity) and J (pure rho)
    Real k_, j_, u_; // K (pure identity), J (pure rho), U (penalty)
    int charge_; // Choose total QN of the chain when conserving by specifying the charge of the central site eg. "0"

    // DMRG parameters
    int max_bond_dim;
    int gs_max_bond_dim_;
    int es_max_bond_dim_;
    double svd_cutoff_;
    double noise_;
    double ortho_weight;

    // DMRG design
    float stable_tol;
    float gs_stable_tol_;
    float es_stable_tol_; // Finish DMRG of current state if energy changes by less than stable_tol_ for num_reps_til_stable_ reps
    int num_reps_til_stable_;
    int init_num_sweeps_per_rep_;
    int num_sweeps_per_rep_;
    int num_states_; // total number of states to run DMRG

    // File I/O
    // Modify prefix_ if storing in custom location
    std::filesystem::path prefix_ = std::filesystem::current_path().parent_path();
    std::filesystem::path progress_directory_ = prefix_ / "pgs";
    std::filesystem::path ee_directory_ = prefix_ / "ee";
    std::filesystem::path en_directory_ = prefix_ / "en";
    std::string filename_; // does not include extension

    std::filesystem::path progress_path_, ee_path_, en_path_;
    SiteSetType sites_;

    DMRGProgress<SiteSetType> dmrg_progress;

    double min_en;
    int cnt;
    int times_swept;
    int bond_dim;

    // Variable declarations
    Real en;
    MPS psi;
    std::string state_name;

public:
    explicit DMRG(std::tuple<std::basic_string<char>, std::basic_string<char>, int, int, float, float, float, float, int, int, int> params){
        site_name_ = std::get<0>(params);
        boundary_condition_ = std::get<1>(params);
        num_sites_ = std::get<2>(params);
        es_max_bond_dim_= std::get<3>(params);
        svd_cutoff_= std::get<4>(params);
        es_stable_tol_= std::get<5>(params);
        theta_= std::get<6>(params);
        u_= std::get<7>(params);
        charge_= std::get<8>(params);
        num_states_= std::get<9>(params);

        sites_ = SiteSetType(num_sites_, {"ConserveQNs=", true});
        if(boundary_condition_ == "p"){
            job_ = site_name_ + " PBC";
        }else if(boundary_condition_ == "o"){
            job_ = site_name_ + " OBC";
        }else if(boundary_condition_ == "s"){
            job_ = site_name_ + " SSD";
        }else{
            job_ = site_name_ + " SSD/PBC";
        }
        // fixme: refactor so that we're passing the most basic variable type instead of a type alias

        // DMRG parameters
        gs_max_bond_dim_ = es_max_bond_dim_;
        noise_ = 0.0;
        ortho_weight = 1000.0;

        // DMRG design
        init_num_sweeps_per_rep_ = 5;
        num_sweeps_per_rep_ = 1;
        num_reps_til_stable_ = 3;
        gs_stable_tol_ = std::pow(es_stable_tol_, 1);

        // Chain parameters
        k_ = cos(theta_ * Pi);
        j_ = sin(theta_ * Pi);
        if(std::abs(k_) < 1E-5){
            k_ = 0;
        }
        if(std::abs(j_) < 1E-5){
            j_ = 0;
        }

        // File I/O
        filename_ = format("%s_%s_%d_%d_%g_%g_%g_%g_%d", site_name_, boundary_condition_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, u_, charge_);
        progress_path_ = progress_directory_ / filename_;
        ee_path_ = ee_directory_ / (filename_ + ".ee");
        en_path_ = en_directory_ / (filename_ + ".en");
    }

    bool IsSimulatingGroundState(){
        if(std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            return dmrg_progress.NumPastStates() == 0;
        }else{
            return true;
        }
    }

    // Adjust parameters depending on whether simulating ground state
    void AdjustDMRGParams(){
        if(IsSimulatingGroundState()){
            max_bond_dim = gs_max_bond_dim_;
            stable_tol = gs_stable_tol_;
        }else{
            max_bond_dim = es_max_bond_dim_;
            stable_tol = es_stable_tol_;
        }
    }

    // Reset min energy and counters to simulate new state
    void ResetDMRG(){
        min_en = 1000000;
        cnt = num_reps_til_stable_;
        times_swept = 0;
    }

    // Decide name of the current state to simulate based on number of completed states
    std::string StateName(){
        std::string state_name;
        switch(dmrg_progress.NumPastStates()){
            case 0 : state_name = "Ground state"; break;
            case 1 : state_name = "1st excited state"; break;
            case 2 : state_name = "2nd excited state"; break;
            case 3 : state_name = "3rd excited state"; break;
            default : state_name = format("%dth excited state", dmrg_progress.NumPastStates()); break;
        };
        return state_name;
    }

    // Simulate states by DMRG
    void Run() {
        if (not std::filesystem::exists(progress_directory_)){
            std::filesystem::create_directory(progress_directory_);
        }
        if (not std::filesystem::exists(ee_directory_)){
            std::filesystem::create_directory(ee_directory_);
        }
        if (not std::filesystem::exists(en_directory_)){
            std::filesystem::create_directory(en_directory_);
        }
        printf("\n> %s: num_sites_=%d max_bond_dim=%d svd_cutoff_=%g stable_tol_=%g theta_=%g k_=%g j_=%g u_=%g q=%d\n", job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, k_, j_, u_, charge_);

        // Set the parameters controlling the accuracy of the DMRG
        // calculation for each DMRG sweep.
        //
//        double min_en;
////        = 1000000;
//        //
//        // Begin the DMRG calculation
//        // for the ground state
//        //
//        int cnt;
////        = num_reps_til_stable_;
//        int times_swept;
////        = 0;
//        int bond_dim;

        // fixme : depends on ground state or not
//        int max_bond_dim = es_max_bond_dim_;
//        Real stable_tol = es_stable_tol_;

//        // Variable declarations
//        Real en;
//        MPS psi;
//        std::string state_name;

        // fixme: internal?
//        DMRGProgress<SiteSetType> dmrg_progress;

        if(std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))){
            // Job was run before
            // Read progress
            dmrg_progress.Read(progress_path_);
            sites_ = dmrg_progress.Sites();

//            switch(dmrg_progress.NumPastStates()){
//                case 0 : state_name = "Ground state"; break;
//                case 1 : state_name = "1st excited state"; break;
//                case 2 : state_name = "2nd excited state"; break;
//                case 3 : state_name = "3rd excited state"; break;
//                default : state_name = format("%dth excited state", dmrg_progress.NumPastStates()); break;
//            }
            // fixme
            state_name = StateName();

            // fixme
            if(dmrg_progress.NumPastStates() >= num_states_){
                printf("\n> Job already complete\n> Requested number of states: %d \n> Completed number of states: %d\n\n", num_states_, dmrg_progress.NumPastStates());
                return;
            }
            // printf("\n> Simulation: %s\n", state);
        }else{
            // New job
            dmrg_progress.StateComplete();
        }

        // fixme: Polymorphic based on SiteSetType
        MPO H = Hamiltonian(sites_, boundary_condition_, num_sites_, u_, k_, j_);

        // If QN conserving, create charge-neutral initial state and change center site to specified charge
        auto pre_init_state = InitState(sites_,"0");
        int center = (num_sites_ + 1) / 2;
        pre_init_state.set(center, std::to_string(charge_));
        MPS init_state = MPS(pre_init_state);

        // If using sine-squared deformed to hot start for periodic
        // hot_start keeps track of which stage (1: SSD warmup, 0: PBC warmup, -1: PBC)
        // Otherwise just -1
        int hot_start = -1;
        if(boundary_condition_ == "sp"){
            hot_start = 1;
        }

        // Loop over states
        // Each iteration includes warmup and post-warmup DMRG runs until energy is stable
        for(;;){
            if(dmrg_progress.times_swepts_.back() == 0 && hot_start != 0){
                // Not yet run
                // Begin DMRG with warmup
                state_name = StateName();
//                switch(dmrg_progress.NumPastStates()){
//                    case 0 : state_name = "Ground state"; break;
//                    case 1 : state_name = "1st excited state"; break;
//                    case 2 : state_name = "2nd excited state"; break;
//                    case 3 : state_name = "3rd excited state"; break;
//                    default : state_name = format("%dth excited state", dmrg_progress.NumPastStates()); break;
//                }
                printf("\n> Simulation: %s\n", state_name);

//                // Adjust parameters depending on whether simulating ground state
//                if(IsSimulatingGroundState()){
//                    max_bond_dim = gs_max_bond_dim_;
//                    stable_tol = gs_stable_tol_;
//                }else{
//                    max_bond_dim = es_max_bond_dim_;
//                    stable_tol = es_stable_tol_;
//                }

//                // Reset
//                min_en = 1000000;
//                cnt = num_reps_til_stable_;
//                times_swept = 0;

                AdjustDMRGParams();
                ResetDMRG();

                // DMRG Warmup
                bond_dim = 200;
                auto sweeps = Sweeps(init_num_sweeps_per_rep_);
                sweeps.maxdim() = 10,20,40,80,200;
                sweeps.cutoff() = std::max(1E-5, svd_cutoff_),std::max(1E-6, svd_cutoff_),std::max(1E-7, svd_cutoff_),std::max(1E-8, svd_cutoff_),std::max(1E-9, svd_cutoff_),std::max(1E-10, svd_cutoff_),std::max(1E-11, svd_cutoff_),std::max(1E-12, svd_cutoff_),svd_cutoff_;
                sweeps.niter() = 2;
                sweeps.noise() = noise_;
                std::tie(en,psi) = dmrg(H, dmrg_progress.PastStates(), init_state, sweeps, {"Quiet", true, "Weight=", ortho_weight});

                times_swept += init_num_sweeps_per_rep_;
                dmrg_progress.Update(times_swept, bond_dim, en, psi, H);
                dmrg_progress.SetSites(sites_);
                dmrg_progress.Write(progress_path_);
                printf("\n    > Times swept: %d\n      %s of %s\n      L=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n", dmrg_progress.times_swepts_.back(), state_name, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, k_, j_, u_, charge_);
            }else if(dmrg_progress.times_swepts_.back() == 0 && hot_start == 0) {
                // If start with sine-squared deformed to hot start for periodic
                // PBC warmup phase
//                switch(dmrg_progress.NumPastStates()){
//                    case 0 : state_name = "Ground state"; break;
//                    case 1 : state_name = "1st excited state"; break;
//                    case 2 : state_name = "2nd excited state"; break;
//                    case 3 : state_name = "3rd excited state"; break;
//                    default : state_name = format("%dth excited state", dmrg_progress.NumPastStates()); break;
//                }
                state_name = StateName();
                printf("\n> Simulation: %s\n", state_name);

//                // Adjust parameters depending on whether simulating ground state
//                if(dmrg_progress.NumPastStates() == 0){
//                    max_bond_dim = gs_max_bond_dim_;
//                    stable_tol = gs_stable_tol_;
//                }else{
//                    max_bond_dim = es_max_bond_dim_;
//                    stable_tol = es_stable_tol_;
//                }

//                // Reset
//                min_en = 1000000;
//                cnt=num_reps_til_stable_;
//                times_swept=0;

                AdjustDMRGParams();
                ResetDMRG();

                // DMRG
                // Warmup by maintaining previous warmup DMRG parameters
                bond_dim = 200;
                auto sweeps = Sweeps(init_num_sweeps_per_rep_);
                sweeps.maxdim() = bond_dim;
                sweeps.cutoff() = svd_cutoff_;
                sweeps.niter() = 2;
                sweeps.noise() = noise_;
                std::tie(en, psi) = dmrg(H, dmrg_progress.PastStates(), init_state, sweeps,
                                         {"Quiet", true, "Weight=", ortho_weight});

                times_swept += init_num_sweeps_per_rep_;
                dmrg_progress.Update(times_swept, bond_dim, en, psi, H);
                dmrg_progress.SetSites(sites_);
                dmrg_progress.Write(progress_path_);
                printf("\n    > Times swept: %d\n      %s of %s\n      L=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n",
                       dmrg_progress.times_swepts_.back(), state_name, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_,
                       gs_stable_tol_, theta_, k_, j_, u_, charge_);
            }

            // Post-warmup DMRG runs
            AdjustDMRGParams();
            for(;;){
                // Read progress
                std::tie(times_swept, bond_dim, en, psi, H) = dmrg_progress.All();

                // Dump entanglement entropy curve
                // Only for simulation of ground state
                if(IsSimulatingGroundState()){
                    DumpEE(num_sites_, CalcEE(psi, num_sites_), ee_path_);
                }

                // Gradually increase bond dimension each run until maximum
                if(bond_dim <= max_bond_dim - 200){
                    bond_dim += 200;
                }

                // DMRG
                auto sw = Sweeps(num_sweeps_per_rep_);
                sw.maxdim() = bond_dim;
                sw.cutoff() = svd_cutoff_;
                sw.niter() = 2;
                sw.noise() = noise_;
                std::tie(en,psi) = dmrg(H, dmrg_progress.PastStates(), psi, sw, {"Quiet=", true, "Weight=", ortho_weight});

                // Decide whether energy has stabilized
                if (abs(min_en-en) * num_sites_ < stable_tol){
                    cnt -= 1;
                } else {
                    cnt = num_reps_til_stable_;
                    min_en = en;
                }
                times_swept += num_sweeps_per_rep_;
                dmrg_progress.Update(times_swept, bond_dim, en, psi, H);
                // fixme
                dmrg_progress.SetSites(sites_);
                dmrg_progress.Write(progress_path_);

                printf("\n    > Times swept: %d\n      %s of %s\n      L=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n", dmrg_progress.times_swepts_.back(), state_name, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, k_, j_, u_, charge_);
                if (cnt == 0){
                    break;
                }
            }

            // Completed current state

            // fixme
            DumpEnergy(dmrg_progress.NumPastStates(), en / num_sites_, en_path_);
            printf("\n> E/num_sites_ of %s: %g\n\n", state_name, en / num_sites_);

            // If using sine-squared deformed to hot start for periodic
            // Make appropriate transition from SSD to PBC while using SSD result as initial state for PBC
            if(hot_start == 1){
                H = Hamiltonian(sites_, "p", num_sites_, u_, k_, j_);
                init_state = psi;
                dmrg_progress = DMRGProgress<SiteSetType>();
                hot_start = 0;
            }else if(hot_start == 0){
                auto pre_init_state = InitState(sites_,"0");
                int center = (num_sites_ + 1) / 2;
                pre_init_state.set(center, std::to_string(charge_));
                MPS init_state = MPS(pre_init_state);
                hot_start = -1;
            }

            // Check whether to simulate next state
            if(dmrg_progress.NumPastStates() >= num_states_){
                // Finished
                dmrg_progress.Write(progress_path_);
                break;
            }

            // Move on to next state
            dmrg_progress.StateComplete();
            // fixme
            dmrg_progress.SetSites(sites_);
            dmrg_progress.Write(progress_path_);

//            else{
//                // Move on to next state
//                dmrg_progress.NextState();
//                // fixme
//                dmrg_progress.SetSites(sites_);
//                dmrg_progress.Write(progress_path_);
//            }
        }
    }

    // Measure matrix elements of translation operator and rho defect operator
    void Analyze(){
        printf("\n%s: num_sites=%d maxmaxdim=%d svd_cutoff_=%g theta=%g K=%g J=%g U=%g\n\n", job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, theta_, k_, j_, u_);

        MPS psiT;
        MPS psiR;

//        auto f_data = HaagerupFData();
//        f_data.DumpFMathematica("/Users/yinhslin/Documents/H3FChain.num_reps_til_stable_");
//        return;

        DMRGProgress<SiteSetType> dmrg_progress;
        dmrg_progress.Read(progress_path_);

        auto states = dmrg_progress.States();
//        sites = dmrg_progress.Sites();
        int len = std::min((int) states.size(), num_states_ + 1);
        // fixme: refactor Haagerup/HaagerupQN branching
        SiteSet sites_new;
         if (site_name_ == "golden" or site_name_ == "haagerup"){
             sites_new = sites_;
             for (int i = 0; i < len; i++) {
                 auto new_psi = MPS(sites_new);
                 auto psi = states.at(i);
                 for(auto j : range1(num_sites_)) {
                     new_psi.set(j, psi(j)* delta(siteIndex(psi, j), sites_new(j)) );
                 }
                 states.at(i) = new_psi;
             }

         } else {
             sites_new = Haagerup(num_sites_);
             for (int i = 0; i < len; i++) {
                 auto states_no_qn = removeQNs(states.at(i));
                 states.at(i) = Z3FourierTransform(states_no_qn, sites_new);
             }
         }
//        auto sites_new = sites;

        auto rho_op = RhoOp(sites_new, site_name_); // changed sites to sites_new

        auto translate_op = TranslationOp(sites_new); // changed sites to sites_new
        // Act translation/rho on psi to create psiT/psiR
        std::vector<MPS> pastT;
        std::vector<MPS> pastR;

        // Dangling indices on the edges of MPO that will be contracted using an IdentityOp tensor
        // to construct a periodic operator.
        auto left_extra = Index(36, "Site");
        auto right_extra = Index(36, "Site");

        for(int i=0;i<len;i++){
            psiT = MPS(states.at(i));
            psiR = MPS(states.at(i));

            // calculate the action of translation on wavefunction
            psiT = applyMPO(translate_op, psiT);
            pastT.push_back(psiT);



            auto new_id = AugmentMPO(IdentityOp(sites_new, rho_op), left_extra, right_extra);

            // Default svd_cutoff_ for applyMPO(density matrix variant) is 1E-12. Setting to larger value will
            // speed up the program significantly at the cost of accuracy of measurements

            auto step1 = applyMPO(AugmentMPO(rho_op, left_extra, right_extra),
                                  AugmentMPS(psiR, left_extra, right_extra), {"Cutoff", 1E-3}
            );
            auto step2 = applyMPO(new_id,step1,{"Cutoff", 1E-3});
            // println(innerC(AugmentMPS(psiR, left_extra, right_extra),step2));


//            psiR = applyMPO(TranslationOp(sites, true), psiR);
//            if (site_name_ == "golden"){
//                GoldenFData f_data = GoldenFData();
//                ActLocal(psiR, f_data.RhoDefectCell(sites(1), sites(2)),1);
//            } else if (site_name_ == "haagerup"){
//                HaagerupFData f_data = HaagerupFData();
//                ActLocal(psiR, f_data.RhoDefectCell(sites(1), sites(2)),1);
//            }
//            psiR = applyMPO(TranslationOp(sites, false), psiR);
            pastR.push_back(step2);
        }

        // Create ITensors for
        // En:  Energies and
        // OpT:  Matrix elements between psi and the translated psiT
        // OpR:  Matrix elements between psi and psiR
        auto s = Index(len);
        auto sP = prime(s);
        auto En = ITensor(dag(s),sP);
        auto OpT = ITensor(dag(s),sP);
        auto OpR = ITensor(dag(s),sP);

        Real gs_en = 0;
        // dmrg_progress.Energies().at(0);

        for(int i=0;i<len;i++){
            En.set(s(i+1),sP(i+1),dmrg_progress.Energies().at(i) - gs_en);
            for(int j=0;j<len;j++){
                OpT.set(s(i+1),sP(j+1),innerC(states.at(i), pastT.at(j)));
                OpR.set(s(i+1),sP(j+1),innerC(AugmentMPS(states.at(i), left_extra, right_extra), pastR.at(j)));
            }
        }
        PrintData(En);
        PrintData(OpR);
        PrintData(OpT);


        // Diagonalize translation
        auto [UT,DT] = eigen(OpT);
        // Diagonalize rho
        auto [UR,DR] = eigen(OpR);

        PrintData(DT);
        PrintData(DR);

        // // Rotate basis for energies accordingly and create diagonal ITensor
        // En = prime(UT) * En * dag(UT);
        // std::vector<Real> diag_En;
        // diag_En.reserve(len);
        // for(int i=0;i<len;i++){
        //     diag_En.push_back(eltC(En, i+1, i+1).real());
        // }
        // En = diagITensor(diag_En, dag(s), prime(s));
        // auto spins = DT.apply([this](Cplx a){return Spin(a, num_sites_);});
        // PrintData(En);
        // PrintData(spins);
    }
};


#endif //CHAINDMRG_DMRG_H

#ifndef CHAINDMRG_DMRG_H
#define CHAINDMRG_DMRG_H
#include <filesystem>
#include <variant>
#include "itensor/all.h"
#include "chain.h"
using namespace itensor;



template<typename SiteSetType>
class DMRGProgress {
public:
    SiteSetType sites_;
    std::vector<int> times_swepts;
    std::vector<int> maxdims;
    std::vector<Real> ens;
    std::vector<MPS> psis;
    std::vector<MPO> Hs;

    void update(int times_swept, int maxdim, Real en, MPS const& psi, MPO H) {
        times_swepts.back() = times_swept;
        maxdims.back() = maxdim;
        ens.back() = en;
        psis.back() = psi;
        Hs.back() = std::move(H);
    }

    void putSites(SiteSetType sites){
        sites_ = sites;
    }

    SiteSetType Sites(){
        return sites_;
    }

    void next() {
        times_swepts.push_back(0);
        maxdims.push_back(0);
        ens.push_back(0);
        psis.emplace_back();
        Hs.emplace_back();
    }

    std::tuple<int,int,Real,MPS,MPO> get() {
        return std::tuple<int,int,Real,MPS,MPO>(times_swepts.back(), maxdims.back(), ens.back(), psis.back(), Hs.back());
    }

    std::vector<MPS> past_states() {
        std::vector<MPS> past_psis;
        past_psis.reserve((int) psis.size()-1);
        for (int i=0;i<(int) psis.size()-1;i++){
            past_psis.push_back(psis.at(i));
        }
        return past_psis;
    }

    std::vector<MPS> States() {
        return psis;
    }

    std::vector<Real> Energies() const {
        return ens;
    }

    int num_past_states() {
        return (int) psis.size()-1;
    }

    int num_states() {
        return (int) psis.size();
    }
    // fixme: string -> filesystem::path
    void write(const std::filesystem::path& p) const
    {
        std::string path = std::string(p);
        writeToFile(path + ".st", sites_);
        writeToFile(path + ".pgs", times_swepts);
        writeToFile(path + ".md", maxdims);
        writeToFile(path + ".en", ens);
        writeToFile(path + ".psi", psis);
        writeToFile(path + ".H", Hs);

        writeToFile(path + ".pgs.bk", times_swepts);
        writeToFile(path + ".md.bk", maxdims);
        writeToFile(path + ".en.bk", ens);
        writeToFile(path + ".psi.bk", psis);
        writeToFile(path + ".H.bk", Hs);
    }
    // fixme: string -> filesystem::path
    void read(const std::filesystem::path& p)
    {
        std::string path = std::string(p);
        try{
            readFromFile(path + ".st", sites_);
            readFromFile(path + ".pgs", times_swepts);
            readFromFile(path + ".md", maxdims);
            readFromFile(path + ".en", ens);
            readFromFile(path + ".psi", psis);
            readFromFile(path + ".H", Hs);
        }catch(std::exception const& e){
            try{
                readFromFile(path + ".pgs.bk", times_swepts);
                readFromFile(path + ".md.bk", maxdims);
                readFromFile(path + ".en.bk", ens);
                readFromFile(path + ".psi.bk", psis);
                readFromFile(path + ".H.bk", Hs);
            }catch(std::exception const& e){
                times_swepts.clear();
                maxdims.clear();
                ens.clear();
                psis.clear();
                Hs.clear();
                next();
            }
        }

    }
};
template<typename SiteSetType>
class DMRG {
    std::string name_;
    std::string boundary_condition;
    std::string job;
    int N;
    int gs_maxmaxdimension;
    double cutoff;
    float tolerance;
    Real theta;
    int init_num_sweeps;
    int num_sweeps;
    int maxmaxdimension;
    double noise;
    int m;
    float gs_tolerance;
    Real K, J, U;
    int charge;
    int num_past_states;
//    std::string progress_path, ee_path, en_path;
//    std::string prefix = std::filesystem::current_path();
    std::filesystem::path prefix = std::filesystem::current_path().parent_path();
    std::filesystem::path progress_path, ee_path, en_path;
    std::filesystem::path progress_directory = prefix / "pgs";
    std::filesystem::path ee_directory = prefix / "ee";
    std::filesystem::path en_directory = prefix / "en";
    std::string filename;
    SiteSetType sites;
public:
    // constructor
    explicit DMRG(std::tuple<std::basic_string<char>, std::basic_string<char>, int, int, float, float, float, float, int, int, int> params){
        name_ = std::get<0>(params);
        boundary_condition = std::get<1>(params);
        N = std::get<2>(params);
        gs_maxmaxdimension= std::get<3>(params);
        cutoff= std::get<4>(params);
        tolerance= std::get<5>(params);
        theta= std::get<6>(params);
        U= std::get<7>(params);
        charge= std::get<8>(params);
        num_past_states= std::get<9>(params);

        if(boundary_condition == "p"){
            job = name_ + " PBC";
        }else if(boundary_condition == "o"){
            job = name_ + " OBC";
        }else if(boundary_condition == "s"){
            job = name_ + " SSD";
        }else{
            job = name_ + " SSD/PBC";
        }
        // fixme: refactor so that we're passing the most basic variable type instead of a type alias
        sites = SiteSetType(N, {"ConserveQNs=", true});
        init_num_sweeps = 5;
        num_sweeps = 1;
        maxmaxdimension = gs_maxmaxdimension;
        noise = 0.0;
        m = 3;
        gs_tolerance = std::pow(tolerance, 1);
        K = cos(theta * Pi);
        J = sin(theta * Pi);
        if(std::abs(K) < 1E-5){
            K = 0;
        }
        if(std::abs(J) < 1E-5){
            J = 0;
        }
        filename = format("%s_%s_%d_%d_%g_%g_%g_%g_%d", name_, boundary_condition, N, gs_maxmaxdimension, cutoff, gs_tolerance, theta, U, charge);
        progress_path = progress_directory / filename;
        ee_path = ee_directory / (filename + ".ee");
        en_path = en_directory / (filename + ".en");
    }
    void run() {
        if (not std::filesystem::exists(progress_directory)){
            std::filesystem::create_directory(progress_directory);
        }
        if (not std::filesystem::exists(ee_directory)){
            std::filesystem::create_directory(ee_directory);
        }
        if (not std::filesystem::exists(en_directory)){
            std::filesystem::create_directory(en_directory);
        }
        printf("\n> %s: N=%d maxmaxdim=%d cutoff=%g tolerance=%g theta=%g K=%g J=%g U=%g q=%d\n",job,N,gs_maxmaxdimension,cutoff,gs_tolerance,theta,K,J,U,charge);

        // Set the parameters controlling the accuracy of the DMRG
        // calculation for each DMRG sweep.
        //
        double min_en = 1000000;
        //
        // Begin the DMRG calculation
        // for the ground state
        //
        int count=m;
        int s=0;
        int maxdim;
        int maxmaxdim = maxmaxdimension;
        Real tol = tolerance;
        Real en;
        MPS psi;
        MPS psiT;
        std::string state;

        DMRGProgress<SiteSetType> dmrg_progress;
        if(std::filesystem::exists(progress_directory / (filename + ".pgs"))){
            dmrg_progress.read(progress_path);
            sites = dmrg_progress.Sites();

            switch(dmrg_progress.num_past_states()){
                case 0 : state = "Ground state"; break;
                case 1 : state = "1st excited state"; break;
                case 2 : state = "2nd excited state"; break;
                case 3 : state = "3rd excited state"; break;
                default : state = format("%dth excited state", dmrg_progress.num_past_states()); break;
            }

            if(dmrg_progress.num_past_states() >= num_past_states){
                printf("\n> Job already completed: Computed to %s\n\n", state);
                return;
            }

            // printf("\n> Simulation: %s\n", state);
        }else{
            dmrg_progress.next();
        }

        MPO H = ConstructH(sites, boundary_condition, N, U, K, J);
//        MPS init_state = InitState(sites,"r");
        auto pre_init_state = InitState(sites,"0");
        int center = (N+1)/2;
        pre_init_state.set(center, std::to_string(charge));
        MPS init_state = MPS(pre_init_state);

        int hot_start = -1;
        if(boundary_condition == "sp"){
            hot_start = 1;
        }

        for(;;){
            if(dmrg_progress.times_swepts.back() == 0 && hot_start != 0){
                switch(dmrg_progress.num_past_states()){
                    case 0 : state = "Ground state"; break;
                    case 1 : state = "1st excited state"; break;
                    case 2 : state = "2nd excited state"; break;
                    case 3 : state = "3rd excited state"; break;
                    default : state = format("%dth excited state", dmrg_progress.num_past_states()); break;
                }
                printf("\n> Simulation: %s\n", state);

                if(dmrg_progress.num_past_states() == 0){
                    maxmaxdim = gs_maxmaxdimension;
                    tol = gs_tolerance;
                }else{
                    maxmaxdim = maxmaxdimension;
                    tol = tolerance;
                }

                min_en = 1000000;
                count=m;
                s=0;

                // warmup
                maxdim = 200;
                auto sweeps = Sweeps(init_num_sweeps);
                sweeps.maxdim() = 10,20,40,80,200;
                sweeps.cutoff() = std::max(1E-5,cutoff),std::max(1E-6,cutoff),std::max(1E-7,cutoff),std::max(1E-8,cutoff),std::max(1E-9,cutoff),std::max(1E-10,cutoff),std::max(1E-11,cutoff),std::max(1E-12,cutoff),cutoff;
                sweeps.niter() = 2;
                sweeps.noise() = noise;
                std::tie(en,psi) = dmrg(H,dmrg_progress.past_states(),init_state,sweeps,{"Quiet",true,"Weight=",1000.0});

                s += init_num_sweeps;
                dmrg_progress.update(s,maxdim,en,psi,H);
                dmrg_progress.putSites(sites);
                dmrg_progress.write(progress_path);
                printf("\n    > Times swept: %d\n      %s of %s\n      N=%d maxmaxdim=%d cutoff=%g tolerance=%g theta=%g K=%g J=%g U=%g q=%d\n",dmrg_progress.times_swepts.back(),state,job,N,gs_maxmaxdimension,cutoff,gs_tolerance,theta,K,J,U,charge);
            }else if(dmrg_progress.times_swepts.back() == 0 && hot_start == 0){
                switch(dmrg_progress.num_past_states()){
                    case 0 : state = "Ground state"; break;
                    case 1 : state = "1st excited state"; break;
                    case 2 : state = "2nd excited state"; break;
                    case 3 : state = "3rd excited state"; break;
                    default : state = format("%dth excited state", dmrg_progress.num_past_states()); break;
                }
                printf("\n> Simulation: %s\n", state);

                if(dmrg_progress.num_past_states() == 0){
                    maxmaxdim = gs_maxmaxdimension;
                    tol = gs_tolerance;
                }else{
                    maxmaxdim = maxmaxdimension;
                    tol = tolerance;
                }

                min_en = 1000000;
                count=m;
                s=0;

                // warmup
                maxdim = 200;
                auto sweeps = Sweeps(init_num_sweeps);
                sweeps.maxdim() = maxdim;
                sweeps.cutoff() = cutoff;
                sweeps.niter() = 2;
                sweeps.noise() = noise;
                std::tie(en,psi) = dmrg(H,dmrg_progress.past_states(),init_state,sweeps,{"Quiet",true,"Weight=",1000.0});

                s += init_num_sweeps;
                dmrg_progress.update(s,maxdim,en,psi,H);
                dmrg_progress.putSites(sites);
                dmrg_progress.write(progress_path);
                printf("\n    > Times swept: %d\n      %s of %s\n      N=%d maxmaxdim=%d cutoff=%g tolerance=%g theta=%g K=%g J=%g U=%g q=%d\n",dmrg_progress.times_swepts.back(),state,job,N,gs_maxmaxdimension,cutoff,gs_tolerance,theta,K,J,U,charge);
            }
            for(;;){
                std::tie(s,maxdim,en,psi,H) = dmrg_progress.get();
                if(dmrg_progress.num_past_states() == 0){
                    dumpEE(N, calcEE(psi, N), ee_path);
                }

                if(maxdim <= maxmaxdim-200){
                    maxdim += 200;
                }

                auto sw=Sweeps(num_sweeps);
                sw.maxdim()=maxdim;
                sw.cutoff()=cutoff;
                sw.niter()=2;
                sw.noise()=noise;
                std::tie(en,psi)=dmrg(H,dmrg_progress.past_states(),psi,sw,{"Quiet=",true,"Weight=",1000.0});
                if (abs(min_en-en)*N < tol){
                    count -= 1;
                } else {
                    count = m;
                    min_en = en;
                }
                s += num_sweeps;
                dmrg_progress.update(s,maxdim,en,psi,H);
                dmrg_progress.putSites(sites);
                dmrg_progress.write(progress_path);

                printf("\n    > Times swept: %d\n      %s of %s\n      N=%d maxmaxdim=%d cutoff=%g tolerance=%g theta=%g K=%g J=%g U=%g q=%d\n",dmrg_progress.times_swepts.back(),state,job,N,gs_maxmaxdimension,cutoff,gs_tolerance,theta,K,J,U,charge);
                if (count == 0){
                    break;
                }
            }

            // debug
            // print(totalQN(psi));
            print("\n");

            dumpEnergy(dmrg_progress.num_past_states(), en/N, en_path);
            printf("\n> E/N of %s: %g\n\n",state,en/N);

            if(hot_start == 1){
                H = ConstructH(sites, "p", N, U, K, J);
                init_state = psi;
                dmrg_progress = DMRGProgress<SiteSetType>();
                hot_start = 0;
            }else if(hot_start == 0){
//                init_state = InitState(sites,"r");
                auto pre_init_state = InitState(sites,"0");
                int center = (N+1)/2;
                pre_init_state.set(center, std::to_string(charge));
                MPS init_state = MPS(pre_init_state);
                hot_start = -1;
            }

            if(dmrg_progress.num_past_states() >= num_past_states){
                dmrg_progress.write(progress_path);
                break;
            }

            dmrg_progress.next();
            dmrg_progress.putSites(sites);
            dmrg_progress.write(progress_path);
        }
    }


    void analyze(){
        printf("\n%s: N=%d maxmaxdim=%d cutoff=%g theta=%g K=%g J=%g U=%g\n\n",job,N,gs_maxmaxdimension,cutoff,theta,K,J,U);

        MPS psiT;
        MPS psiR;

//        auto f_data = HaagerupFData();
//        f_data.DumpF("/Users/yinhslin/Documents/H3FChain.m");
//        return;

        DMRGProgress<SiteSetType> dmrg_progress;
        dmrg_progress.read(progress_path);

        auto states = dmrg_progress.States();
//        sites = dmrg_progress.Sites();
        int len = std::min((int) states.size(), num_past_states+1);
        // fixme: refactor Haagerup/HaagerupQN branching
        SiteSet sites_new;
         if (name_ == "golden" or name_ == "haagerup"){
             sites_new = sites;
             for (int i = 0; i < len; i++) {
                 auto new_psi = MPS(sites_new);
                 auto psi = states.at(i);
                 for(auto j : range1(N)) {
                     new_psi.set(j, psi(j)* delta(siteIndex(psi, j), sites_new(j)) );
                 }
                 states.at(i) = new_psi;
             }

         } else {
             sites_new = Haagerup(N);
             for (int i = 0; i < len; i++) {
                 auto states_no_qn = removeQNs(states.at(i));
                 states.at(i) = BasisRotate(states_no_qn, sites_new);
             }
         }
//        auto sites_new = sites;

        auto rho_op = RhoOp(sites_new, name_); // changed sites to sites_new

        auto translate_op = TranslationOp(sites_new); // changed sites to sites_new
        // Act translation/rho on psi to create psiT/psiR
        std::vector<MPS> pastT;
        std::vector<MPS> pastR;

        // Dangling indices on the edges of MPO that will be contracted using an identity tensor
        // to construct a periodic operator.
        auto left_extra = Index(36, "Site");
        auto right_extra = Index(36, "Site");

        for(int i=0;i<len;i++){
            psiT = MPS(states.at(i));
            psiR = MPS(states.at(i));

            // calculate the action of translation on wavefunction
            psiT = applyMPO(translate_op, psiT);
            pastT.push_back(psiT);



            auto new_id = augmentMPO(identity(sites_new, rho_op), left_extra, right_extra);

            // Default cutoff for applyMPO(density matrix variant) is 1E-12. Setting to larger value will
            // speed up the program significantly at the cost of accuracy of measurements

            auto step1 = applyMPO(augmentMPO(rho_op, left_extra, right_extra),
                                  augmentMPS(psiR, left_extra, right_extra),{"Cutoff", 1E-3}
            );
            auto step2 = applyMPO(new_id,step1,{"Cutoff", 1E-3});
            // println(innerC(augmentMPS(psiR, left_extra, right_extra),step2));


//            psiR = applyMPO(TranslationOp(sites, true), psiR);
//            if (name_ == "golden"){
//                GoldenFData f_data = GoldenFData();
//                ActLocal(psiR, f_data.RhoDefect(sites(1), sites(2)),1);
//            } else if (name_ == "haagerup"){
//                HaagerupFData f_data = HaagerupFData();
//                ActLocal(psiR, f_data.RhoDefect(sites(1), sites(2)),1);
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
                OpR.set(s(i+1),sP(j+1),innerC(augmentMPS(states.at(i), left_extra, right_extra), pastR.at(j)));
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
        // auto spins = DT.apply([this](Cplx a){return Spin(a, N);});
        // PrintData(En);
        // PrintData(spins);
    }
};


#endif //CHAINDMRG_DMRG_H

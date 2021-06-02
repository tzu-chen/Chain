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
    // Each element below is the latest progress for a given state
    std::vector<int> times_swepts_;
    std::vector<int> maxdims_;
    std::vector<Real> ens_;
    std::vector<MPS> psis_;
    std::vector<MPO> Hs_;

    void SetSites(SiteSetType sites) {
        sites_ = sites;
    }

    SiteSetType Sites() {
        return sites_;
    }

    // Replace the latest progress for the current state
    void Update(int times_swept, int max_dim, Real en, MPS const& psi, MPO H) {
        times_swepts_.back() = times_swept;
        maxdims_.back() = max_dim;
        ens_.back() = en;
        psis_.back() = psi;
        Hs_.back() = std::move(H);
    }

    // Extend vectors to signal completion of previous state and prepare for next state
    void NextState() {
        times_swepts_.push_back(0);
        maxdims_.push_back(0);
        ens_.push_back(0);
        psis_.emplace_back(MPS(1));
        Hs_.emplace_back(MPO(1));
    }

    std::tuple<int,int,Real,MPS,MPO> All() {
        return std::tuple<int,int,Real,MPS,MPO>(times_swepts_.back(), maxdims_.back(), ens_.back(), psis_.back(), Hs_.back());
    }

    std::vector<MPS> DoneStates() {
        std::vector<MPS> past_psis;
        past_psis.reserve((int) psis_.size() - 1);
        for (int i=0;i< (int) psis_.size() - 1; i++) {
            past_psis.push_back(psis_.at(i));
        }
        return past_psis;
    }

//    // fixme
//    std::vector<MPS> States() {
//        return psis_;
//    }

    std::vector<Real> Energies() const {
        return ens_;
    }

    int NumDoneStates() {
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
        }catch(std::exception const& e) {
//            try{
//                readFromFile(path + ".pgs.bk", times_swepts_);
//                readFromFile(path + ".md.bk", maxdims_);
//                readFromFile(path + ".en.bk", ens_);
//                readFromFile(path + ".psi.bk", psis_);
//                readFromFile(path + ".H.bk", Hs_);
//            }catch(std::exception const& e) {
                println(e.what());
                times_swepts_.clear();
                maxdims_.clear();
                ens_.clear();
                psis_.clear();
                Hs_.clear();
            NextState();
//            }
        }

    }
};

// Examples of SiteSetType are Golden, Haagerup, HaagerupQ
template<typename SiteSetType>
class DMRG {
    std::string site_type_; // eg. golden, haagerup, haagerup_q
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
    std::filesystem::path m_directory_ =  prefix_ / "m";
    std::string filename_; // does not include extension

    std::filesystem::path progress_path_, ee_path_, en_path_, m_path_;
    SiteSetType sites_;

    DMRGProgress<SiteSetType> dmrg_progress_;

    double min_en_;
    int cnt_;
    int times_swept_;
    int bond_dim_;

    // Variable declarations
    Real en_;
    MPS psi_;
    std::string state_name_;
    MPS init_state_;

public:
    explicit DMRG(std::tuple<std::basic_string<char>, std::basic_string<char>, int, int, float, float, float, float, int, int, int> params) {
        site_type_ = std::get<0>(params);
        boundary_condition_ = std::get<1>(params);
        num_sites_ = std::get<2>(params);
        es_max_bond_dim_ = std::get<3>(params);
        svd_cutoff_ = std::get<4>(params);
        es_stable_tol_ = std::get<5>(params);
        theta_ = std::get<6>(params);
        u_ = std::get<7>(params);
        charge_ = std::get<8>(params);
        num_states_ = std::get<9>(params);

        sites_ = SiteSetType(num_sites_, {"ConserveQNs=", true});
        if (boundary_condition_ == "p") {
            job_ = site_type_ + " PBC";
        } else if (boundary_condition_ == "o") {
            job_ = site_type_ + " OBC";
        } else if (boundary_condition_ == "s") {
            job_ = site_type_ + " SSD";
        } else {
            job_ = site_type_ + " SSD/PBC";
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
        if (std::abs(k_) < 1E-5) {
            k_ = 0;
        }
        if (std::abs(j_) < 1E-5) {
            j_ = 0;
        }

        // File I/O
        filename_ = format("%s_%s_%d_%d_%g_%g_%g_%g_%d", site_type_, boundary_condition_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, u_, charge_);
        progress_path_ = progress_directory_ / filename_;
        ee_path_ = ee_directory_ / (filename_ + ".ee");
        en_path_ = en_directory_ / (filename_ + ".en");
        m_path_ = m_directory_ / (filename_ + ".m");
    }

    bool IsSimulatingGroundState() {
        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            return dmrg_progress_.NumDoneStates() == 0;
        } else {
            return true;
        }
    }

    // Adjust parameters depending on whether simulating ground state
    void AdjustDMRGParams() {
        if (IsSimulatingGroundState()) {
            max_bond_dim = gs_max_bond_dim_;
            stable_tol = gs_stable_tol_;
        } else {
            max_bond_dim = es_max_bond_dim_;
            stable_tol = es_stable_tol_;
        }
    }

    // Reset min energy and counters to simulate new state
    void ResetDMRG() {
        min_en_ = 1000000;
        cnt_ = num_reps_til_stable_;
        times_swept_ = 0;
    }

    // Decide name of the current state to simulate based on number of completed states
    std::string StateName() {
        std::string state_name;
        switch(dmrg_progress_.NumDoneStates()) {
            case 0 : state_name = "Ground state"; break;
            case 1 : state_name = "1st excited state"; break;
            case 2 : state_name = "2nd excited state"; break;
            case 3 : state_name = "3rd excited state"; break;
            default : state_name = format("%dth excited state", dmrg_progress_.NumDoneStates()); break;
        }
        return state_name;
    }

    MPS DefaultInitState() {
        auto pre_init_state = InitState(sites_,"0");
        int center = (num_sites_ + 1) / 2;
        pre_init_state.set(center, std::to_string(charge_));
        MPS init_state = MPS(pre_init_state);
        return init_state;
    }

    // Simulate states by DMRG
    void Run() {
        if (not std::filesystem::exists(progress_directory_)) {
            std::filesystem::create_directory(progress_directory_);
        }
        if (not std::filesystem::exists(ee_directory_)) {
            std::filesystem::create_directory(ee_directory_);
        }
        if (not std::filesystem::exists(en_directory_)) {
            std::filesystem::create_directory(en_directory_);
        }
        printf("\n> %s: num_sites=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n", job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, k_, j_, u_, charge_);

        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            // Job was run before
            // Read progress
            dmrg_progress_.Read(progress_path_);
            sites_ = dmrg_progress_.Sites();

            // fixme
            state_name_ = StateName();

            // fixme
            if (dmrg_progress_.NumDoneStates() >= num_states_) {
                printf("\n> Job already complete\n> Requested number of states: %d \n> Completed number of states: %d\n\n", num_states_,
                       dmrg_progress_.NumDoneStates());
                return;
            }
            // printf("\n> Simulation: %s\n", state);
        } else {
            // New job
            dmrg_progress_.NextState();
        }

        // fixme: Polymorphic based on SiteSetType
        MPO H = Hamiltonian(sites_, boundary_condition_, num_sites_, u_, k_, j_);

        // If QN conserving, create charge-neutral initial state and change center site to specified charge
//        auto pre_init_state = InitState(sites_,"0");
//        int center = (num_sites_ + 1) / 2;
//        pre_init_state.set(center, std::to_string(charge_));
//        MPS init_state = MPS(pre_init_state);
        init_state_ = DefaultInitState();

        // If using sine-squared deformed to hot start for periodic
        // hot_start keeps track of which stage (1: SSD warmup, 0: PBC warmup, -1: PBC)
        // Otherwise just -1
        int hot_start = -1;
        if (boundary_condition_ == "sp") {
            hot_start = 1;
        }

        // Loop over states
        // Each iteration includes warmup and post-warmup DMRG runs until energy is stable
        for(;;) {
            if (dmrg_progress_.times_swepts_.back() == 0 && hot_start != 0) {
                // Not yet run
                // Begin DMRG with warmup
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);

                AdjustDMRGParams();
                ResetDMRG();

                // DMRG Warmup
                bond_dim_ = 200;
                auto sweeps = Sweeps(init_num_sweeps_per_rep_);
                sweeps.maxdim() = 10,20,40,80,200;
                sweeps.cutoff() = std::max(1E-5, svd_cutoff_),std::max(1E-6, svd_cutoff_),std::max(1E-7, svd_cutoff_),std::max(1E-8, svd_cutoff_),std::max(1E-9, svd_cutoff_),std::max(1E-10, svd_cutoff_),std::max(1E-11, svd_cutoff_),std::max(1E-12, svd_cutoff_),svd_cutoff_;
                sweeps.niter() = 2;
                sweeps.noise() = noise_;
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), init_state_, sweeps, {"Quiet", true, "Weight=", ortho_weight});

                times_swept_ += init_num_sweeps_per_rep_;
                dmrg_progress_.Update(times_swept_, bond_dim_, en_, psi_, H);
                dmrg_progress_.SetSites(sites_);
                dmrg_progress_.Write(progress_path_);
                printf("\n    > Times swept: %d\n      %s of %s\n      L=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n", dmrg_progress_.times_swepts_.back(), state_name_, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, k_, j_, u_, charge_);
            } else if (dmrg_progress_.times_swepts_.back() == 0 && hot_start == 0) {
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);

                AdjustDMRGParams();
                ResetDMRG();

                // DMRG
                // Warmup by maintaining previous warmup DMRG parameters
                bond_dim_ = 200;
                auto sweeps = Sweeps(init_num_sweeps_per_rep_);
                sweeps.maxdim() = bond_dim_;
                sweeps.cutoff() = svd_cutoff_;
                sweeps.niter() = 2;
                sweeps.noise() = noise_;
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), init_state_, sweeps,
                                           {"Quiet", true, "Weight=", ortho_weight});

                times_swept_ += init_num_sweeps_per_rep_;
                dmrg_progress_.Update(times_swept_, bond_dim_, en_, psi_, H);
                dmrg_progress_.SetSites(sites_);
                dmrg_progress_.Write(progress_path_);
                printf("\n    > Times swept: %d\n      %s of %s\n      L=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n",
                       dmrg_progress_.times_swepts_.back(), state_name_, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_,
                       gs_stable_tol_, theta_, k_, j_, u_, charge_);
            }

            // Post-warmup DMRG runs
            AdjustDMRGParams();
            for(;;) {
                // Read progress
                std::tie(times_swept_, bond_dim_, en_, psi_, H) = dmrg_progress_.All();

                // Dump entanglement entropy curve
                // Only for simulation of ground state
                if (IsSimulatingGroundState()) {
                    DumpEE(num_sites_, CalcEE(psi_, num_sites_), ee_path_);
                }

                // Gradually increase bond dimension each run until maximum
                if (bond_dim_ <= max_bond_dim - 200) {
                    bond_dim_ += 200;
                }

                // DMRG
                auto sw = Sweeps(num_sweeps_per_rep_);
                sw.maxdim() = bond_dim_;
                sw.cutoff() = svd_cutoff_;
                sw.niter() = 2;
                sw.noise() = noise_;
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), psi_, sw, {"Quiet=", true, "Weight=", ortho_weight});

                // Decide whether energy has stabilized
                if (abs(min_en_ - en_) * num_sites_ < stable_tol) {
                    cnt_ -= 1;
                } else {
                    cnt_ = num_reps_til_stable_;
                    min_en_ = en_;
                }
                times_swept_ += num_sweeps_per_rep_;
                dmrg_progress_.Update(times_swept_, bond_dim_, en_, psi_, H);
                // fixme
                dmrg_progress_.SetSites(sites_);
                dmrg_progress_.Write(progress_path_);

                printf("\n    > Times swept: %d\n      %s of %s\n      L=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n", dmrg_progress_.times_swepts_.back(), state_name_, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, k_, j_, u_, charge_);
                if (cnt_ == 0) {
                    break;
                }
            }

            // Completed current state

            // fixme
            DumpEnergy(dmrg_progress_.NumDoneStates(), en_ / num_sites_, en_path_);
            printf("\n> E/L of %s: %g\n\n", state_name_, en_ / num_sites_);

            // If using sine-squared deformed to hot start for periodic
            // Make appropriate transition from SSD to PBC while using SSD result as initial state for PBC
            if (hot_start == 1) {
                H = Hamiltonian(sites_, "p", num_sites_, u_, k_, j_);
                init_state_ = psi_;
                dmrg_progress_ = DMRGProgress<SiteSetType>();
                hot_start = 0;
            } else if (hot_start == 0) {
//                auto pre_init_state = InitState(sites_, "0");
//                int center = (num_sites_ + 1) / 2;
//                pre_init_state.set(center, std::to_string(charge_));
//                MPS init_state = MPS(pre_init_state);
                init_state_ = DefaultInitState();
                hot_start = -1;
            }

            dmrg_progress_.NextState(); // Even if no next state, this marks completion of current state
            // fixme
            dmrg_progress_.SetSites(sites_);
            dmrg_progress_.Write(progress_path_);

            // Check whether to simulate next state
            if (dmrg_progress_.NumDoneStates() >= num_states_) {
                break;
            }
        }
    }

    // Measure matrix elements of translation operator and rho defect operator
    void Analyze() {
        if (not std::filesystem::exists(m_directory_)) {
            std::filesystem::create_directory(m_directory_);
        }

//        printf("\n%s: num_sites=%d max_bond_dim=%d svd_cutoff=%g theta=%g K=%g J=%g U=%g\n\n", job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, theta_, k_, j_, u_);
        printf("\n> %s: num_sites=%d max_bond_dim=%d svd_cutoff=%g stable_tol=%g theta=%g K=%g J=%g U=%g Q=%d\n", job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_stable_tol_, theta_, k_, j_, u_, charge_);

        // Default svd_cutoff for applyMPO(density matrix variant) is 1E-12
        // Setting to larger value will speed up the program significantly at the cost of accuracy of measurements
        float svd_cutoff = 1E-8;

        // Variable declaration
        MPS psi_translated;
        MPS psi_acted_by_rho;
        std::vector<MPS> states_translated;
        std::vector<MPS> states_acted_by_rho;

        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            dmrg_progress_.Read(progress_path_);
        } else {
            printf("\n> No progress file available\n");
            return;
        }

        auto states = dmrg_progress_.DoneStates();
        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);

        // fixme: refactor Haagerup/HaagerupQ branching
        // Perform change of basis if appropriate in order to use F symbols
        // Rotate HaagerupQ to Haagerup
        // Golden -> Golden, Haagerup, HaagerupQ -> Haagerup
        SiteSet sites;
        if (site_type_ == "golden" or site_type_ == "haagerup") {
            sites = sites_;
            for (int i = 0; i < num_states; i++) {
                auto new_psi = MPS(sites);
                auto psi = states.at(i);
                for(auto j : range1(num_sites_)) {
                    new_psi.set(j, psi(j) * delta(siteIndex(psi, j), sites(j)) );
                }
                states.at(i) = new_psi;
            }
        } else {
            sites = Haagerup(num_sites_);
            for (int i = 0; i < num_states; i++) {
                auto states_no_qn = removeQNs(states.at(i));
                states.at(i) = Z3FourierTransform(states_no_qn, sites);
            }
        }

        auto translate_op = TranslationOp(sites); // periodic MPS
        auto rho_op = RhoOp(sites, site_type_); // dangling periodic MPS
        auto id_op = IdentityOp(sites, rho_op); // dangling periodic MPS with matching indices with rho_op

        // Dangling bond indices at the edges of MPO
        // fixme
        auto left_dangling_ind = Index(36, "Site");
        auto right_dangling_ind = Index(36, "Site");

        // fixme: Can we use OpenMP for this?
        for (int i=0; i < num_states; i++) {
            // Initialize
            psi_translated = MPS(states.at(i));
            psi_acted_by_rho = MPS(states.at(i));

            // Act by translation
            psi_translated = applyMPO(translate_op, psi_translated);
            states_translated.push_back(psi_translated);

            // Act by rho defect in two steps
            // First act by dangling rho operator
            // Then act by dangling identity operator
            psi_acted_by_rho = applyMPO(AugmentMPO(rho_op, left_dangling_ind, right_dangling_ind),
                                  AugmentMPS(psi_acted_by_rho, left_dangling_ind, right_dangling_ind),
                                  {"Cutoff", svd_cutoff}
            );
            psi_acted_by_rho = applyMPO(AugmentMPO(id_op, left_dangling_ind, right_dangling_ind),
                                        psi_acted_by_rho, {"Cutoff", svd_cutoff});
            states_acted_by_rho.push_back(psi_acted_by_rho);

            // println(innerC(AugmentMPS(psiR, left_extra, right_extra),step2));

//            psiR = applyMPO(TranslationOp(sites, true), psiR);
//            if (site_name_ == "golden") {
//                GoldenFData f_data = GoldenFData();
//                ActLocal(psiR, f_data.RhoDefectCell(sites(1), sites(2)),1);
//            } else if (site_name_ == "haagerup") {
//                HaagerupFData f_data = HaagerupFData();
//                ActLocal(psiR, f_data.RhoDefectCell(sites(1), sites(2)),1);
//            }
//            psiR = applyMPO(TranslationOp(sites, false), psiR);
        }

        // Create ITensors for
        // En: Energies
        // OpT: Matrix elements between psi and psi_translated
        // OpR: Matrix elements between psi and psi_acted_by_rho
        auto s = Index(num_states);
        auto sP = prime(s);
        auto en_matrix = ITensor(dag(s), sP);
        auto translation_matrix = ITensor(dag(s), sP);
        auto rho_matrix = ITensor(dag(s), sP);

        Real en_shift = 0;
        // Uncomment if shift by ground state energy
        // en_shift = dmrg_progress.Energies().at(0);

        for(int i=0; i < num_states; i++) {
            en_matrix.set(s(i + 1), sP(i + 1), dmrg_progress_.Energies().at(i) - en_shift);
            for(int j=0; j < num_states; j++) {
                translation_matrix.set(s(i + 1), sP(j + 1), innerC(states.at(i), states_translated.at(j)));
                rho_matrix.set(s(i + 1), sP(j + 1), innerC(AugmentMPS(states.at(i), left_dangling_ind, right_dangling_ind), states_acted_by_rho.at(j)));
            }
        }

        DumpMathematica(num_states, en_matrix, translation_matrix, rho_matrix, m_path_);

        // Diagonalize translation
        auto [UT,translation_diag] = eigen(translation_matrix);
        // Diagonalize rho
        auto [UR,rho_diag] = eigen(rho_matrix);

        printf("\n> Unordered set of energies:\n");
        PrintData(en_matrix);
        printf("\n> Unordered set of translation eigenvalues:\n");
        PrintData(translation_diag);
        printf("\n> Unordered set of rho eigenvalues:\n");
        PrintData(rho_diag);

        // // Rotate basis for energies accordingly and create diagonal ITensor
        // En = prime(UT) * En * dag(UT);
        // std::vector<Real> diag_En;
        // diag_En.reserve(len);
        // for(int i=0;i<len;i++) {
        //     diag_En.push_back(eltC(En, i+1, i+1).real());
        // }
        // En = diagITensor(diag_En, dag(s), prime(s));
        // auto spins = DT.apply([this](Cplx a) {return Spin(a, num_sites_);});
        // PrintData(En);
        // PrintData(spins);
    }
};
#endif //CHAINDMRG_DMRG_H
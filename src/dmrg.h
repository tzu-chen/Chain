#ifndef CHAINDMRG_DMRG_H
#define CHAINDMRG_DMRG_H

#include <filesystem>
#include <variant>
#include "itensor/all.h"
#include "chain.h"
#include <sstream>
#include "path.h"

using namespace itensor;


// Reads and writes data coming out of DMRG
template<typename SiteSetType>
class DMRGProgress {
public:
    SiteSetType sites_;
    // Each element below is the latest progress for a given state.
    std::vector<int> num_sweeps_vec_;
    std::vector<int> max_dims_;
    std::vector<Real> ens_;
    std::vector<MPS> psis_;
    std::vector<MPO> Hs_;
    std::vector<MPS> psis_acted_by_rho_;
    std::vector<MPS> psis_acted_by_translation_;
    std::vector<MPS> psis_acted_by_hamiltonian_;

    void SetSites(SiteSetType sites) {
        sites_ = sites;
    }

    SiteSetType Sites() {
        return sites_;
    }

    // Replace the latest progress for the current state
    void Update(int times_swept, int max_dim, Real en, MPS const& psi, MPO H) {
        num_sweeps_vec_.back() = times_swept;
        max_dims_.back() = max_dim;
        ens_.back() = en;
        psis_.back() = psi;
        Hs_.back() = std::move(H);
    }

    // Extend vectors to signal completion of previous state and prepare for next state
    void NextState() {
        num_sweeps_vec_.push_back(0);
        max_dims_.push_back(0);
        ens_.push_back(0);
        psis_.emplace_back(MPS(1));
        Hs_.emplace_back(MPO(1));
    }

    MPO H() {
        return Hs_.back();
    }

    std::tuple<int,int,Real,MPS,MPO> All() {
        return std::tuple<int,int,Real,MPS,MPO>(num_sweeps_vec_.back(), max_dims_.back(), ens_.back(), psis_.back(), Hs_.back());
    }

    std::vector<MPS> DoneStates() {
        std::vector<MPS> past_psis;
        past_psis.reserve((int) psis_.size() - 1);
        for (int i=0;i< (int) psis_.size() - 1; i++) {
            past_psis.push_back(psis_.at(i));
        }
        return past_psis;
    }

    std::vector<MPS> StatesActedByHamiltonian() {
        return psis_acted_by_hamiltonian_;
    }

    std::vector<MPS> StatesActedByTranslation() {
        return psis_acted_by_translation_;
    }

    std::vector<MPS> StatesActedByRho() {
        return psis_acted_by_rho_;
    }

    std::vector<Real> Energies() const {
        return ens_;
    }

    int NumDoneStates() {
        return (int) psis_.size() - 1;
    }

    void Write(const std::filesystem::path& p, bool backup = false) const
    {
        std::string path = std::string(p);
        if (!backup) {
            writeToFile(path + ".st", sites_);
            writeToFile(path + ".pgs", num_sweeps_vec_);
            writeToFile(path + ".md", max_dims_);
            writeToFile(path + ".en", ens_);
            writeToFile(path + ".psi", psis_);
            writeToFile(path + ".H", Hs_);
        } else {
            // Uncomment to save backup files
            writeToFile(path + ".st.bk", sites_);
            writeToFile(path + ".pgs.bk", num_sweeps_vec_);
            writeToFile(path + ".md.bk", max_dims_);
            writeToFile(path + ".en.bk", ens_);
            writeToFile(path + ".psi.bk", psis_);
            writeToFile(path + ".H.bk", Hs_);
        }
    }

    void Read(const std::filesystem::path& p, bool backup = false)
    {
        std::string path = std::string(p);
        if (!backup) {
            readFromFile(path + ".st", sites_);
            readFromFile(path + ".pgs", num_sweeps_vec_);
            readFromFile(path + ".md", max_dims_);
            readFromFile(path + ".en", ens_);
            readFromFile(path + ".psi", psis_);
            readFromFile(path + ".H", Hs_);
        } else {
            readFromFile(path + ".st.bk", sites_);
            readFromFile(path + ".pgs.bk", num_sweeps_vec_);
            readFromFile(path + ".md.bk", max_dims_);
            readFromFile(path + ".en.bk", ens_);
            readFromFile(path + ".psi.bk", psis_);
            readFromFile(path + ".H.bk", Hs_);
        }
//        catch (std::exception const& e) {
//            try {
//                readFromFile(path + ".st.bk", sites_);
//                readFromFile(path + ".pgs.bk", num_sweeps_vec_);
//                readFromFile(path + ".md.bk", max_dims_);
//                readFromFile(path + ".en.bk", ens_);
//                readFromFile(path + ".psi.bk", psis_);
//                readFromFile(path + ".H.bk", Hs_);
//            } catch (std::exception const &e) {
//                // File is corrupted, start over
//                println(e.what());
//                num_sweeps_vec_.clear();
//                max_dims_.clear();
//                ens_.clear();
//                psis_.clear();
//                Hs_.clear();
//                NextState();
//            }
//        }
    }

    void RecoverStates(const std::filesystem::path& p) {
        Read(p);

        std::vector<int> num_sweeps_vec = num_sweeps_vec_;
        std::vector<int> max_dims = max_dims_;
        std::vector<Real> ens = ens_;
        std::vector<MPS> psis = psis_;
        std::vector<MPO> Hs = Hs_;

        num_sweeps_vec_ = std::vector<int>();
        max_dims_ = std::vector<int>();
        ens_ = std::vector<Real>();
        psis_ = std::vector<MPS>();
        Hs_ = std::vector<MPO>();
        NextState();

        int num_done_states = (int) psis.size()-1;
        printf("Recover states from corrupt vector.\n");
        for (int i=0;i<num_done_states;i++) {
            innerC(psis.at(i),psis.at(i));
            printf("\rState %d/%d is OK.", i+1, num_done_states);

            num_sweeps_vec_.pop_back();
            max_dims_.pop_back();
            ens_.pop_back();
            psis_.pop_back();
            Hs_.pop_back();

            num_sweeps_vec_.push_back(num_sweeps_vec.at(i));
            max_dims_.push_back(max_dims.at(i));
            ens_.push_back(ens.at(i));
            psis_.push_back(psis.at(i));
            Hs_.push_back(Hs.at(i));
            NextState();

            Write(p);
        }
    }

    // Write and read states acted by Hamiltonian, translation, and rho.
    void WriteHamiltonian(const std::filesystem::path& p) const
    {
        std::string path = std::string(p);
        writeToFile(path + ".ham", psis_acted_by_hamiltonian_);
        writeToFile(path + ".ham.bk", psis_acted_by_hamiltonian_);
    }

    void ReadHamiltonian(const std::filesystem::path& p)
    {
        std::string path = std::string(p);
        try{
            readFromFile(path + ".ham", psis_acted_by_hamiltonian_);
        } catch (std::exception const& e) {
            try{
                readFromFile(path + ".ham.bk", psis_acted_by_hamiltonian_);
            } catch (std::exception const& e) {
                println(e.what());
                psis_acted_by_hamiltonian_.pop_back();
            }
        }
    }

    void WriteTranslation(const std::filesystem::path& p) const
    {
        std::string path = std::string(p);
        writeToFile(path + ".tra", psis_acted_by_translation_);
        writeToFile(path + ".tra.bk", psis_acted_by_translation_);
    }

    void ReadTranslation(const std::filesystem::path& p)
    {
        std::string path = std::string(p);
        try{
            readFromFile(path + ".tra", psis_acted_by_translation_);
        } catch (std::exception const& e) {
            try{
                readFromFile(path + ".tra.bk", psis_acted_by_translation_);
            } catch (std::exception const& e) {
                println(e.what());
                psis_acted_by_translation_.pop_back();
            }
        }
    }

    void WriteRho(const std::filesystem::path& p) const
    {
        std::string path = std::string(p);
        writeToFile(path + ".rho", psis_acted_by_rho_);
        writeToFile(path + ".rho.bk", psis_acted_by_rho_);
    }

    void ReadRho(const std::filesystem::path& p)
    {
        std::string path = std::string(p);
        try{
            readFromFile(path + ".rho", psis_acted_by_rho_);
        } catch (std::exception const& e) {
            try {
                readFromFile(path + ".rho.bk", psis_acted_by_rho_);
            } catch (std::exception const &e) {
                println(e.what());
                psis_acted_by_rho_.pop_back();
            }
        }
    }
};

// SiteSetType can be Golden, Haagerup, HaagerupQ
template<typename SiteSetType>
class DMRG {
    std::string site_type_; // "golden", "haagerup", "haagerupq"
    std::string boundary_condition_;
        // "p" for periodic, "o" for open, "s" for sine-squared deformed, "sp" for hot-starting "p" with "s"
    std::string job_; // site_name_ + boundary_condition_

    // Chain parameters
    int num_sites_;
    Real u_;
    std::string coupling_str_;
    std::vector<Real> couplings_;
    int charge_; // Choose total QN of the chain when conserving by specifying the charge of the central site eg. "0"

    // DMRG parameters
    int max_bond_dim;
    int gs_max_bond_dim_;
    int es_max_bond_dim_;
    double svd_cutoff_;
    double noise_;
    double orthogonal_weight;

    // DMRG design
    float precision_; // Finish DMRG of current state if energy changes by less than precision_/num_sites_ for num_reps_til_stable_ reps
    float gs_precision_;
    float es_precision_;
    int num_reps_til_stable_;
    int init_num_sweeps_per_rep_;
    int num_sweeps_per_rep_;
    int num_states_; // total number of states to run DMRG

    // File I/O
    std::filesystem::path prefix_ = kPath;
    std::filesystem::path progress_directory_ = prefix_ / "pgs";
    std::filesystem::path ee_directory_ = prefix_ / "ee";
    std::filesystem::path en_directory_ = prefix_ / "en";
    std::filesystem::path m_directory_ =  prefix_ / "m";
    std::filesystem::path an_directory_ =  prefix_ / "an";
    std::string filename_; // does not include extension

    std::filesystem::path progress_path_, ee_path_, en_path_, m_path_, an_path_;
    SiteSetType sites_;

    DMRGProgress<SiteSetType> dmrg_progress_;

    double min_en_;
    int cnt_;
    int num_sweeps_;
    int bond_dim_;

    // Variable declarations
    Real en_;
    MPS psi_;
    std::string state_name_;
    MPS init_state_;

public:
    explicit DMRG(std::tuple<std::basic_string<char>, std::basic_string<char>, int, int, float, float, std::string, float, int, int, int> params) {
        site_type_ = std::get<0>(params);
        boundary_condition_ = std::get<1>(params);
        num_sites_ = std::get<2>(params);
        es_max_bond_dim_ = std::get<3>(params);
        svd_cutoff_ = std::get<4>(params);
        es_precision_ = std::get<5>(params);
        coupling_str_ = std::get<6>(params);
        couplings_ = ParseCouplings(coupling_str_);
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
        } else if (boundary_condition_ == "d") {
            job_ = site_type_ + " DBC";
        } else if (boundary_condition_ == "sp") {
            job_ = site_type_ + " SSD/PBC";
        } else {
            printf("Invalid boundary condition");
            return;
        }
        // fixme: refactor so that we're passing the most basic variable type instead of a type alias.
        // fixme: pass site instead of SiteSet.

        // DMRG parameters
        gs_max_bond_dim_ = es_max_bond_dim_;
        noise_ = 0.0;
        orthogonal_weight = 1000.0;

        // DMRG design
        init_num_sweeps_per_rep_ = 5;
        num_sweeps_per_rep_ = 1;
        num_reps_til_stable_ = 3;
        gs_precision_ = std::pow(es_precision_, 1);

        // File I/O
        filename_ = format("%s_%s_%d_%d_%g_%g_%g_%s_%d", site_type_, boundary_condition_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_precision_, coupling_str_, u_, charge_);
        progress_path_ = progress_directory_ / filename_;
        ee_path_ = ee_directory_ / (filename_ + ".ee");
        en_path_ = en_directory_ / (filename_ + ".en");
        m_path_ = m_directory_ / (filename_ + ".m");
        an_path_ = an_directory_ / (filename_ + ".an");
    }

    // Parse coupling string by comma for coefficients of projectors.
    std::vector<Real> ParseCouplings(std::string coupling_string) {
        std::vector<Real> couplings;
        std::stringstream s_stream(coupling_string);
        while(s_stream.good()) {
            string substr;
            getline(s_stream, substr, ',');
            double val = std::stod(substr);
            if (abs(val) < 1E-5){
                val = 0;
            }
            couplings.push_back(val);
        };
        return couplings;
    }

    // Print information about the current job.
    void PrintJob(bool with_sweeps) {
        if (with_sweeps) {
            int num_sweeps = dmrg_progress_.num_sweeps_vec_.back();
            if (num_sweeps == 0) {
                printf("\n  > Sweeps #1-#%d\n    %s of %s\n    L=%d  bond_dim=%d  svd_cutoff=%g  precision=%g  penalty=%g  couplings=%s  charge=%d\n", init_num_sweeps_per_rep_, state_name_, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_precision_, u_, coupling_str_, charge_);
            } else {
                printf("\n  > Sweep #%d\n    %s of %s\n    L=%d  bond_dim=%d  svd_cutoff=%g  precision=%g  penalty=%g  couplings=%s  charge=%d\n", num_sweeps+1, state_name_, job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_precision_, u_, coupling_str_, charge_);
            }
        } else {
            printf("\n> %s:  L=%d  bond_dim=%d  svd_cutoff=%g  precision=%g  penalty=%g  couplings=%s  charge=%d\n", job_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_precision_, u_, coupling_str_, charge_);
        }
    }

    // Some actions, such as computing entanglement entropy curve, are to be done only for first state.
    // Better DMRG parameters might also be desired for first state.
    bool IsSimulatingFirstState() {
        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            return dmrg_progress_.NumDoneStates() == 0;
        } else {
            return true;
        }
    }

    // Adjust parameters depending on whether simulating first state.
    void AdjustDMRGParams() {
        if (IsSimulatingFirstState()) {
            max_bond_dim = gs_max_bond_dim_;
            precision_ = gs_precision_;
        } else {
            max_bond_dim = es_max_bond_dim_;
            precision_ = es_precision_;
        }
    }

    // Reset min energy and counters to simulate new state.
    void ResetDMRG() {
        min_en_ = 1000000;
        cnt_ = num_reps_til_stable_;
        num_sweeps_ = 0;
    }

    // Decide name of the current state to simulate based on number of completed states.
    // Use 1st, 2nd instead of ground state, 1st excited state because for charged sector first state is not special.
    std::string StateName() {
        std::string state_name;
        switch(dmrg_progress_.NumDoneStates()) {
            case 0 : state_name = "1st state"; break;
            case 1 : state_name = "2nd state"; break;
            case 2 : state_name = "3rd state"; break;
            default : state_name = format("%dth state", dmrg_progress_.NumDoneStates()+1); break;
        }
        return state_name;
    }

    // If QN conserving, create charge-neutral initial state and change center site to specified charge.
    MPS DefaultInitState() {
        auto pre_init_state = InitState(sites_,"0");
        int center = (num_sites_+1)/2;
        pre_init_state.set(center, std::to_string(charge_));
//        MPS init_state = MPS(pre_init_state);
//        return init_state;
        return randomMPS(pre_init_state);
//        return randomMPS(sites_);
    }

    // Cannot catch EXC_BAD_ACCESS on iOS
    // https://stackoverflow.com/questions/16202029/is-there-a-way-to-catch-or-handle-exc-bad-access
    void Simulate(bool analyze = false) {
        try {
            SimulateRaw(analyze);
        } catch (std::exception const& e) {
            println(e.what());
            // Replace progress files with backup progress files
            // Put in try-catch because .bk files may not exist
            try {
                printf("Replace progress files with backup progress files.\n");
                dmrg_progress_.Read(progress_path_, true);
                dmrg_progress_.Write(progress_path_);
            } catch (std::exception const& e) {
            }
            try {
                SimulateRaw(analyze);
            } catch (std::exception const& e) {
                dmrg_progress_.RecoverStates(progress_path_);
                Simulate(analyze);
            }
        }
    }

    // Simulate states by DMRG.
    void SimulateRaw(bool analyze = false) {
        if (not std::filesystem::exists(progress_directory_)) {
            std::filesystem::create_directory(progress_directory_);
        }
        if (not std::filesystem::exists(ee_directory_)) {
            std::filesystem::create_directory(ee_directory_);
        }
        if (not std::filesystem::exists(en_directory_)) {
            std::filesystem::create_directory(en_directory_);
        }
        PrintJob(false);

        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            // Job was run before, read progress
            dmrg_progress_.Read(progress_path_);
            sites_ = dmrg_progress_.Sites();
            if (dmrg_progress_.NumDoneStates() >= num_states_) {
                printf("\n> Job already complete\n> Requested number of states: %d \n> Completed number of states: %d\n\n", num_states_,
                       dmrg_progress_.NumDoneStates());
                // Assumes that for charged sectors there is no need to measure rho, valid for Fibonacci and Haagerup.
                if (analyze) {
                    if (boundary_condition_ == "p") {
                        if (charge_ == 0) {
                            Analyze();
                        } else {
                            AnalyzeWithoutRho();
                        }
                    } else {
                        Energies();
                    }
                }
                return;
            }
        } else {
            // Job not run before, new job.
            dmrg_progress_.NextState();
            dmrg_progress_.SetSites(sites_);
            dmrg_progress_.Write(progress_path_);
        }

        MPO H = sites_.Hamiltonian(boundary_condition_, num_sites_, u_, couplings_);
        init_state_ = DefaultInitState();

        // Only meaningful if using sine-squared deformed to hot start for periodic.
        // hot_start keeps track of which stage (1: SSD warmup, 0: PBC warmup, -1: PBC)
        // Otherwise just -1.
        int hot_start = -1;
        if (boundary_condition_ == "sp") {
            hot_start = 1;
        }

        // Loop over states, each iteration includes warmup and post-warmup DMRG runs until energy is stable.
        for (;;) {
            // TODO: refactor the multiple appearances of AdjustDMRGParams to here?
//            AdjustDMRGParams();
            if (dmrg_progress_.num_sweeps_vec_.back() == 0 && hot_start != 0) {
                // Did not begin to simulate, begin DMRG with warmup.
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);
                AdjustDMRGParams();
                ResetDMRG();
                PrintJob(true);
                // DMRG Warmup.
                bond_dim_ = 200;
                auto sweeps = Sweeps(init_num_sweeps_per_rep_);
                sweeps.maxdim() = 10,20,40,80,200;
                sweeps.cutoff() = std::max(1E-5, svd_cutoff_),std::max(1E-6, svd_cutoff_),std::max(1E-7, svd_cutoff_),std::max(1E-8, svd_cutoff_),std::max(1E-9, svd_cutoff_),std::max(1E-10, svd_cutoff_),std::max(1E-11, svd_cutoff_),std::max(1E-12, svd_cutoff_),svd_cutoff_;
                sweeps.niter() = 2;
                sweeps.noise() = noise_;
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), init_state_, sweeps, {"Quiet", true, "Weight=", orthogonal_weight});
                num_sweeps_ += init_num_sweeps_per_rep_;
                dmrg_progress_.Update(num_sweeps_, bond_dim_, en_, psi_, H);
                dmrg_progress_.Write(progress_path_);
            } else if (dmrg_progress_.num_sweeps_vec_.back() == 0 && hot_start == 0) {
                // Only meaningful if using sine-squared deformed to hot start for periodic.
                // Switching from SSD to PBC so warmup again.
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);
                AdjustDMRGParams();
                ResetDMRG();
                PrintJob(true);
                // DMRG warmup by maintaining previous warmup DMRG parameters.
                bond_dim_ = 200;
                auto sweeps = Sweeps(init_num_sweeps_per_rep_);
                sweeps.maxdim() = bond_dim_;
                sweeps.cutoff() = svd_cutoff_;
                sweeps.niter() = 2;
                sweeps.noise() = noise_;
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), init_state_, sweeps,
                                           {"Quiet", true, "Weight=", orthogonal_weight});
                num_sweeps_ += init_num_sweeps_per_rep_;
                dmrg_progress_.Update(num_sweeps_, bond_dim_, en_, psi_, H);
                dmrg_progress_.Write(progress_path_);
            } else {
                // Warmup previously completed, no more warmup.
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);
            }

            // Regular (post-warmup) DMRG runs.
            AdjustDMRGParams();
            for (;;) {
                // Read progress.
                std::tie(num_sweeps_, bond_dim_, en_, psi_, H) = dmrg_progress_.All();
                // Dump entanglement entropy curve, only for simulation of first state
                if (IsSimulatingFirstState()) {
                    DumpEE(num_sites_, CalcEE(psi_, num_sites_), ee_path_);
                }
                // Gradually increase bond dimension each run until maximum.
                if (bond_dim_ <= max_bond_dim - 200) {
                    bond_dim_ += 200;
                }
                PrintJob(true);
                // DMRG.
                auto sw = Sweeps(num_sweeps_per_rep_);
                sw.maxdim() = bond_dim_;
                sw.cutoff() = svd_cutoff_;
                sw.niter() = 2;
                sw.noise() = noise_;
                // TODO
//                std::tie(en_, psi_) = dmrg(H, psi_, sw, {"Quiet=", true, "Weight=", orthogonal_weight});
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), psi_, sw, {"Quiet=", true, "Weight=", orthogonal_weight});
                // Decide whether energy has stabilized.
                if (abs(min_en_ - en_) * num_sites_ < precision_) {
                    cnt_ -= 1;
                } else {
                    cnt_ = num_reps_til_stable_;
                    min_en_ = en_;
                }
                num_sweeps_ += num_sweeps_per_rep_;
                dmrg_progress_.Update(num_sweeps_, bond_dim_, en_, psi_, H);
                dmrg_progress_.Write(progress_path_);
                if (cnt_ == 0) {
                    break;
                }
            }

            // Only meaningful if using sine-squared deformed to hot start for periodic.
            // Make appropriate transition from SSD to PBC while using SSD result as initial state for PBC.
            if (hot_start == 1) {
                H = sites_.Hamiltonian("p", num_sites_, u_, couplings_);
                init_state_ = psi_;
                dmrg_progress_ = DMRGProgress<SiteSetType>();
                hot_start = 0;
            } else if (hot_start == 0) {
                init_state_ = DefaultInitState();
                hot_start = -1;
            }

            // Even if no next state, this marks completion of current state.
            dmrg_progress_.NextState();
            dmrg_progress_.Write(progress_path_);
            // Write backup progress files
            if (dmrg_progress_.NumDoneStates() % 1 == 0) {
                dmrg_progress_.Write(progress_path_, true);
            }

            // Completed current state.
            DumpEnergy(dmrg_progress_.NumDoneStates()+1, en_, en_path_);
            printf("\n> Energy of %s: %g\n\n", state_name_, en_);

            // Assumes that for charged sectors there is no need to measure rho, valid for Fibonacci and Haagerup.
            if (analyze) {
                if (boundary_condition_ == "p") {
                    if (charge_ == 0) {
                        Analyze();
                    } else {
                        AnalyzeWithoutRho();
                    }
                } else {
                    Energies();
                }
            }

            // Check whether to simulate next state.
            if (dmrg_progress_.NumDoneStates() >= num_states_) {
                break;
            }
        }
    }

//    DEPRECIATED
//    The following methods decide filters out states with bad translation or rho eigenvalues.
//    Decided it is better to just output the eigenvalues regardless.
//
//    // Check if the translation eigenvalue is a num_sites_ root of unity to within numerical precision
//    // If close to the nth power of the generating root of unity, make it n and return filtered vector
//    std::vector<Cplx> FilterTranslation(std::vector<Cplx> translation_eigenvalues) {
//        std::vector<Cplx> result;
//        for (int i=0; i<translation_eigenvalues.size(); i++) {
//            Cplx tmp = log(translation_eigenvalues.at(i)) / (2 * Pi * 1_i) * num_sites_;
//            if (abs(tmp.imag()) < 1E-1 && abs(tmp.real()-round(tmp.real())) < 1E-1) {
//                if (round(tmp.real()) == -0) {
//                    result.push_back( 0 );
//                } else if (round(tmp.real()) == - num_sites_/2) {
//                    result.push_back( num_sites_/2 );
//                } else {
//                    result.push_back( round(tmp.real()) );
//                }
//            }
//        }
//        return result;
//    }
//
//    // Check if the rho eigenvalue is one of the possibilities to within numerical precision
//    // If close to a possibility p, make it p and return filtered vector
//    std::vector<Cplx> FilterRho(std::vector<Cplx> rho_eigenvalues, std::vector<Cplx> possibilities) {
//        std::vector<Cplx> result;
//        for (int i=0; i<rho_eigenvalues.size(); i++) {
//            Cplx tmp = rho_eigenvalues.at(i);
////            println(tmp);
////            if (abs(tmp.imag()) < 1E-3) {
//                for (int j=0; j<possibilities.size(); j++) {
//                    if (abs(tmp - possibilities.at(j)) < 1E-1) {
//                        result.push_back(possibilities.at(j));
//                    }
//                }
////            }
//        }
//        return result;
//    }
//
//    std::vector<Cplx> OrderedAppend(std::vector<Cplx> old_vector, std::vector<Cplx> new_vector) {
//        std::vector<Cplx> old_vector_copy = old_vector;
//        std::vector<Cplx> new_ones;
//        for (int i=0; i<new_vector.size(); i++) {
//            bool is_contained = false;
//            for (int j=0; j<old_vector_copy.size(); j++) {
//                if (new_vector.at(i) == old_vector_copy.at(j)) {
//                    old_vector_copy.erase(old_vector_copy.begin() + j);
//                    is_contained = true;
//                    break;
//                }
//            }
//            if (!is_contained) {
//                new_ones.push_back(new_vector.at(i));
//            }
//        }
//        for (int i=0; i<new_ones.size(); i++) {
//            old_vector.push_back(new_ones.at(i));
//        }
//        return old_vector;
//    }

    // Measure matrix elements of translation operator and rho defect operator
    ITensor Measure(std::string observable) {
        // SVD cutoff for applying MPOs, default is 1E-12.
        // Setting to coarser value will speed up the program at the cost of accuracy of measurements.
        float svd_cutoff = 1E-5;

        // Variable declaration.
        MPS psi_acted;
        std::vector<MPS> states_acted;

        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            dmrg_progress_.Read(progress_path_);
        } else {
            printf("\n> No progress file available\n");
            return ITensor();
        }

        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);
        auto s = Index(num_states);
        auto sP = prime(s);
        auto matrix = ITensor(dag(s), sP);
        // fixme
//        if (observable == "energy") {
//            std::vector<Real> eigenvalues;
//            for (int i=0;i<num_states;i++) {
//                eigenvalues.push_back(dmrg_progress_.Energies().at(i));
//                matrix.set(i+1,i+1,dmrg_progress_.Energies().at(i));
//            }
////            printf("\n> %s eigenvalues:\n", observable);
////            PrintVector(eigenvalues);
////            printf("\n\n");
//            DumpMathematicaSingle(observable, num_states, matrix, m_path_);
            std::filesystem::remove(en_path_);
            for (int i=0;i<dmrg_progress_.NumDoneStates();i++) {
                DumpEnergy(i+1, dmrg_progress_.Energies().at(i), en_path_);
            }
            // fixme
        {
            auto states = dmrg_progress_.DoneStates();

            // fixme: refactor Haagerup/HaagerupQ branching.
            // Perform change of basis on HaagerupQ state to turn into Haagerup state in order to use F-symbols.
            // Golden -> Golden, Haagerup, HaagerupQ -> Haagerup.
            SiteSet sites;
            if (site_type_ == "golden") {
                sites = Golden(num_sites_);
            } else {
                sites = Haagerup(num_sites_);
            }

            if (observable != "energy") {
                if (site_type_ == "golden" or site_type_ == "haagerup") {
                    for (int i = 0; i < num_states; i++) {
                        auto new_psi = MPS(sites);
                        auto psi = states.at(i);
                        for (auto j : range1(num_sites_)) {
                            new_psi.set(j, psi(j) * delta(siteIndex(psi, j), sites(j)));
                        }
                        states.at(i) = new_psi;
                    }
                } else {
                    for (int i = 0; i < num_states; i++) {
                        auto states_no_qn = removeQNs(states.at(i));
                        states.at(i) = Z3FourierTransform(states_no_qn, sites);
                    }
                }
            }

            if (observable == "energy") {
                // TODO: Compute the Hamiltonian matrix for simulated states.
                // Due to simulation error, there will be small off-diagonal elements.
                // Idea is to diagonalize the Hamiltonian matrix to reduce simulation errors.
                // Below are a few attempts, by applyMPO, inner, and davidson.

                // applyMPO: density matrix method too heavy, fit method too inaccurate.
//                if (std::filesystem::exists(progress_directory_ / (filename_ + ".ham"))) {
//                    dmrg_progress_.ReadHamiltonian(progress_path_);
//                    for (int i=0; i<std::min(num_states, (int) dmrg_progress_.StatesActedByHamiltonian().size()); i++) {
//                        states_acted.push_back(dmrg_progress_.StatesActedByHamiltonian().at(i));
//                    }
//                }
//                MPO H;
//                println("\nAct Hamiltonian...");
//                for (int i = states_acted.size(); i < num_states; i++) {
//                    psi_acted = MPS(states.at(i));
//                    H = dmrg_progress_.Hs_.at(i);
//                    printf("\r%d/%d", i+1, num_states);
//                    psi_acted = applyMPO(H, psi_acted, {"Method", "Fit", "Cutoff", svd_cutoff, "Nsweep", 1});
//                    states_acted.push_back(psi_acted);
//                    dmrg_progress_.psis_acted_by_hamiltonian_.push_back(psi_acted);
//                    dmrg_progress_.WriteHamiltonian(progress_path_);
//                }
//                printf("\r...finished.\n");

                // fixme: inner:
//                H = dmrg_progress_.Hs_.at(0);
//                auto psi = MPS(states.at(0));
//                println(innerC(psi, H, psi));

                // davidson: doesn't work for more than one state.
//                println("davidson");
//                H = dmrg_progress_.Hs_.at(0);
//                LocalMPO PH(H, {});
//                LocalMPO PH1(H, {});
//                auto psi = MPS(states.at(0));
//                auto psi1 = MPS(states.at(1));
//                int b = (int) num_sites_/2;
//                b = 1;
//                PH.position(b,psi);
//                PH.position(b,psi1);
//                PH1.position(b,psi1);
//                auto phi = psi(b)*psi(b+1);
//                auto phi1 = psi1(b)*psi1(b+1);
//                std::vector<ITensor> phis = {phi, phi1};
//                double energy = davidson(PH,phi,{});
//                double energy1 = davidson(PH1,phi1,{});
//                println(energy);
//                println(energy1);
//                PrintVector(davidson(PH,phis,{}));
//                println("");
//                println("end");
            }

            if (observable == "translation") {
                if (std::filesystem::exists(progress_directory_ / (filename_ + ".tra"))) {
                    dmrg_progress_.ReadTranslation(progress_path_);
                    for (int i=0; i<std::min(num_states, (int) dmrg_progress_.StatesActedByTranslation().size()); i++) {
                        states_acted.push_back(dmrg_progress_.StatesActedByTranslation().at(i));
                    }
                }
                // Act translation by applyMPO, slower than consecutive swaps.
//                auto translation_op = TranslationOp(sites);
                println("\nAct translation...");
                for (int i = states_acted.size(); i < num_states; i++) {
                    psi_acted = MPS(states.at(i));
                    printf("\r%d/%d", i+1, num_states);
//                    psi_acted = applyMPO(translation_op, psi_acted, {"Method", "Fit", "Cutoff", svd_cutoff, "Verbose", true});
                    for (int j=1; j<num_sites_; j++) {
                        Swap(psi_acted, j);
                    }
                    states_acted.push_back(psi_acted);
                    dmrg_progress_.psis_acted_by_translation_.push_back(psi_acted);
                    dmrg_progress_.WriteTranslation(progress_path_);
                }
                printf("\r...finished.\n");
            }

            int d = dim(sites(1));
            auto left_dangling_ind = Index(d*d, "Site");
            auto right_dangling_ind = Index(d*d, "Site");

            if (observable == "rho") {
                // TODO: Zipper algorithm, currently slower and more inaccurate than applyMPO fit.

                // test ActThree using older rho operator specific to L=3
//                auto dummy_FF = HaagerupFData();
//                auto rho = dummy_FF.RhoDefect3(sites(1), sites(2), sites(3));
//                auto dummy_psi0 = MPS(states.at(0));
//                auto dummy_psi = MPS(states.at(0));
//                ActThree(dummy_psi, rho, 1);
//                println(innerC(dummy_psi, dummy_psi0));
//                for(int i=1;i<=3;i++){
//                    rho *= dummy_psi0(i) * prime(dummy_psi0(i)).conj();
//                }
//                PrintData(rho);

//                // zipper start
//                for (int k=0;k<1;k++) {
//                    auto start = std::chrono::high_resolution_clock::now();
//
//                    auto dummy_l = Index(6, "Site");
//                    auto dummy_ll = Index(6, "Site");
//                    auto zipperMPS0 = MPS(ZipperAugmentMPS(MPS(states.at(k)), dummy_l, dummy_ll));
//                    auto zipperMPS = MPS(zipperMPS0);
//                    auto dummy_F = HaagerupFData();
//                    auto FGate = dummy_F.ZipperAugmentGate(siteIndex(zipperMPS, 1),
//                                                           siteIndex(zipperMPS, 2),
//                                                           siteIndex(zipperMPS, 3));
//                    auto ns = num_sites_ + 2;
//                    ActThree(zipperMPS0, dummy_F.ZipperDeltaGate(siteIndex(zipperMPS,1),
//                                                           siteIndex(zipperMPS,2),
//                                                           siteIndex(zipperMPS,3)),
//                             1
//                    );
//                    ActThree(zipperMPS, dummy_F.ZipperDeltaGate(siteIndex(zipperMPS,1),
//                                                                 siteIndex(zipperMPS,2),
//                                                                 siteIndex(zipperMPS,3)),
//                             1
//                    );
//                    ActThree(zipperMPS, FGate, 1);
//                    for (int i = 2; i < ns - 1; i++) {
//                        ActThree(zipperMPS, dummy_F.ZipperGate(siteIndex(zipperMPS, i),
//                                                               siteIndex(zipperMPS, i + 1),
//                                                               siteIndex(zipperMPS, i + 2)),
//                                 i
//                        );
//                    }
//                    for (int j = 1; j < ns; j++) {
//                        localSwap(zipperMPS, j);
//                    }
//                    ActThree(zipperMPS, dummy_F.ZipperGate(siteIndex(zipperMPS, ns - 2),
//                                                           siteIndex(zipperMPS, ns - 1),
//                                                           siteIndex(zipperMPS, ns)),
//                             ns - 2
//
//                    );
//                    for (int j = 1; j < ns; j++) {
//                        localSwap(zipperMPS, j);
//                    }
//                    ActThree(zipperMPS, dummy_F.ZipperReductionGate(siteIndex(zipperMPS, ns - 2),
//                                                                    siteIndex(zipperMPS, ns - 1),
//                                                                    siteIndex(zipperMPS, ns)),
//                             ns - 2
//                             );
//                    for (int j = ns - 1; j > 0; j--) {
//                        localSwap(zipperMPS, j);
//                    }
//                    for (int j = ns - 1; j > 0; j--) {
//                        localSwap(zipperMPS, j);
//                    }
//                    for (int j=ns-1; j >0; j--) {
//                        localSwap(zipperMPS, j);
//                    }
//                    println(innerC(zipperMPS, zipperMPS0));
//
//                    auto stop = std::chrono::high_resolution_clock::now();
//                    printf("\nTime spent: %gs\n", std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000);
//                }
//                // zipper end

                if (std::filesystem::exists(progress_directory_ / (filename_ + ".rho"))) {
                    dmrg_progress_.ReadRho(progress_path_);
                    for (int i=0; i<std::min(num_states, (int) dmrg_progress_.StatesActedByRho().size()); i++) {
                        states_acted.push_back(dmrg_progress_.StatesActedByRho().at(i));
                    }
                }

                // Dangling MPOs.
//                auto rho_op = RhoOp(sites, site_type_);
//                auto id_op = IdentityOp(sites, rho_op); // Matching indices with rho_op.
                println("\nAct rho...");

                // TODO: Measure rho by davidson?

                // Use nmultMPO to contract rho_op and id_op to create periodic MPO.
//                println("nmultMPO");
//                psi_acted = MPS(states.at(0));
//                id_op.prime("Site");
//                // Can verify that ordering of id_op and rho_op doesn't matter.
//                auto rho_op_nmult = nmultMPO(id_op, rho_op);
//                println(innerC(psi_acted, rho_op_nmult, psi_acted));
//                println("end");

                // TODO: Make applyMPO fit method work if first act identity MPO, density matrix method works but slow.
//                rho_op = RhoOp(sites, site_type_);
//                id_op = IdentityOp(sites, rho_op);
//                println("inner");
//                auto psi0 = MPS(states.at(0));
//                psi0 = AugmentMPS(psi0, left_dangling_ind, right_dangling_ind);
//                psi_acted = MPS(states.at(0));
//                psi_acted = AugmentMPS(psi_acted, left_dangling_ind, right_dangling_ind);
////                fitApplyMPOImpl(psi0,AugmentMPO(id_op, left_dangling_ind, right_dangling_ind),psi_acted,{"Method", "Fit", "Cutoff", svd_cutoff, "Nsweep", 1});
//                psi_acted = applyMPO(AugmentMPO(id_op, left_dangling_ind, right_dangling_ind),
//                                     psi_acted, {"Method", "Fit", "Cutoff", svd_cutoff, "Nsweep", 1});
//                println("checkpoint 1");
//                println(innerC(psi0, psi_acted));
//                println("checkpoint 2");
//                println(innerC(psi_acted, AugmentMPO(rho_op, left_dangling_ind, right_dangling_ind),
//                               psi0));
//                println("end");

                // TODO: Try to apply identity MPO by direct construction from the state. Failed.
//                auto psi0 = MPS(states.at(0));
//                auto psi_acted = MPS(states.at(0));
//                auto psi1 = ApplyIdentityOp(psi0, d*d, left_dangling_ind, right_dangling_ind);
//                psi0 = AugmentMPS(psi0, left_dangling_ind, right_dangling_ind);
//                psi_acted = AugmentMPS(psi_acted, left_dangling_ind, right_dangling_ind);
//                psi_acted = applyMPO(AugmentMPO(id_op, left_dangling_ind, right_dangling_ind), psi_acted);
//                PrintData(psi1);
//                PrintData(psi_acted);
//                println(innerC(psi0, AugmentMPO(rho_op, left_dangling_ind, right_dangling_ind),
//                               psi1));
//                println(innerC(psi0, AugmentMPO(rho_op, left_dangling_ind, right_dangling_ind),
//                               psi_acted));

                MPO rho_op;
                MPO id_op;
//                if (boundary_condition_ == "p") {
                rho_op = RhoOp(sites, site_type_);
                id_op = IdentityOp(sites, rho_op);
                for (int i = states_acted.size(); i < num_states; i++) {
//                    auto start = std::chrono::high_resolution_clock::now();
                    // Initialize
                    psi_acted = MPS(states.at(i));
                    // Act rho in two steps, first by dangling rho operator, then by dangling identity operator.
                    printf("\r%d/%d", i + 1, num_states);
                    psi_acted = applyMPO(AugmentMPO(rho_op, left_dangling_ind, right_dangling_ind),
                                         AugmentMPS(psi_acted, left_dangling_ind, right_dangling_ind),
                                         {"Method", "Fit", "Cutoff", svd_cutoff, "Nsweep", 2}
                    );
                    psi_acted = applyMPO(AugmentMPO(id_op, left_dangling_ind, right_dangling_ind),
                                         psi_acted, {"Method", "Fit", "Cutoff", svd_cutoff, "Nsweep", 1});
                    states_acted.push_back(psi_acted);
                    dmrg_progress_.psis_acted_by_rho_.push_back(psi_acted);
                    dmrg_progress_.WriteRho(progress_path_);

//                    auto stop = std::chrono::high_resolution_clock::now();
//                    printf("\nTime spent: %gs\n", std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000);
                }
                printf("\r...finished.\n");
//                } else {
//                    rho_op = OpenRhoOp(sites, site_type_);
//                    for (int i = states_acted.size(); i < num_states; i++) {
//                        psi_acted = MPS(states.at(i));
//                        printf("\r%d/%d", i + 1, num_states);
//                        psi_acted = applyMPO(rho_op, psi_acted);
//                        states_acted.push_back(psi_acted);
//                        dmrg_progress_.psis_acted_by_rho_.push_back(psi_acted);
//                        dmrg_progress_.WriteRho(progress_path_);
//                    }
//                    printf("\r...finished.\n");
//                }
            }

            Real en_shift = 0;
            // Uncomment if shift by ground state energy.
            // en_shift = dmrg_progress.Energies().at(0);

            for (int i = 0; i < num_states; i++) {
                // Use energies obtained from DMRG (by davidson)
                if (observable == "energy") {
                    matrix.set(s(i+1), sP(i+1), dmrg_progress_.Energies().at(i) - en_shift);
                }

                for (int j = 0; j < num_states; j++) {
                    // Uncomment to use energies obtained by diagonalizing Hamiltonian matrix.
//                    if (observable == "energy" || observable == "translation") {
                    if (observable == "translation") {
                        matrix.set(s(i+1), sP(j+1), innerC(states.at(i), states_acted.at(j)));
                    }
                    if (observable == "rho") {
//                    if (boundary_condition_ == "p") {
                        matrix.set(s(i + 1), sP(j + 1),
                                   innerC(AugmentMPS(states.at(i), left_dangling_ind, right_dangling_ind),
                                          states_acted.at(j)));
//                        } else {
//                            matrix.set(s(i + 1), sP(j + 1), innerC(states.at(i),states_acted.at(j)));
//                        }
                    }
                }
            }
            DumpMathematicaSingle(observable, num_states, matrix, m_path_);

//            DEPRECIATED
//            // Analysis of translation and rho eigenvalues assuming that the states are simulated well enough.
//            // Otherwise, better work with the Mathematica .m file created above.
//            std::vector<Cplx> eigenvalues;
//            for (int i = 0; i < num_states; i++) {
//                auto submatrix = ITensor(dag(s), sP);
//                for (int j = 1; j <= i+1; j++) {
//                    for (int k = 1; k <= i+1; k++) {
//                        submatrix.set(s(j), sP(k), eltC(matrix, j, k));
//                        submatrix.set(s(j), sP(k), eltC(matrix, j, k));
//                    }
//                }
//                auto[U, sub_eigenvalues] = eigen(submatrix);
//                std::vector<Cplx> eigenvalues_tmp;
//                for (int j = 1; j <= i+1; j++) {
//                    eigenvalues_tmp.push_back(eltC(sub_eigenvalues, j, j));
//                }
//                if (observable == "translation") {
//                    eigenvalues_tmp = FilterTranslation(eigenvalues_tmp);
//                }
//                if (observable == "rho") {
//                    std::vector<Cplx> rho_possibilities;
//                    if (site_type_ == "golden") {
//                        rho_possibilities = {(1 + sqrt(5)) / 2, (1 - sqrt(5)) / 2};
//                    } else {
//                        rho_possibilities = {(3 + sqrt(13)) / 2, (3 - sqrt(13)) / 2, 1, -1};
//                    }
//                    eigenvalues_tmp = FilterRho(eigenvalues_tmp, rho_possibilities);
//                }
//                eigenvalues = OrderedAppend(eigenvalues, eigenvalues_tmp);
//            }
//            printf("> %s eigenvalues:\n", observable);
//            PrintVector(eigenvalues);
//            printf("\n\n");
        }
        return matrix;
    }

    // Mode 1
    void Analyze() {
        if (not std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            printf("\n> No progress file available\n");
            return;
        }
        if (not std::filesystem::exists(m_directory_)) {
            std::filesystem::create_directory(m_directory_);
        }
        if (not std::filesystem::exists(an_directory_)) {
            std::filesystem::create_directory(an_directory_);
        }
        PrintJob(false);
        std::filesystem::remove(m_path_);
        auto energy_matrix = Measure("energy");
        auto translation_matrix = Measure("translation");
        auto rho_matrix = Measure("rho");

        // Uniformize the ITensor Index structure
        auto s = findInds(energy_matrix, {})(1);
        auto s_translation = findInds(translation_matrix, {})(1).noPrime();
        auto s_rho = findInds(rho_matrix, {})(1).noPrime();
        translation_matrix = translation_matrix * delta(s, s_translation) * delta(prime(s), prime(s_translation));
        rho_matrix = rho_matrix * delta(s, s_rho) * delta(prime(s), prime(s_rho));

        ITensor m;
        if (boundary_condition_ == "p") {
            m = energy_matrix/num_sites_ + translation_matrix + rho_matrix;
        } else {
            m = energy_matrix/num_sites_;
        }
        // TODO: Diagonalization does not work well for haagerup if we do not include rho here.
        auto[U,D] = eigen(m);
        U.conj();
        energy_matrix = dag(U)*energy_matrix*prime(U);
        translation_matrix = dag(U)*translation_matrix*prime(U);
        rho_matrix = dag(U)*rho_matrix*prime(U);
        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);
        std::vector<std::vector<Real>> result;
        auto precision = 1000000.0;
        for (int i=1;i<=num_states;i++) {
            auto momentum = std::round((log(eltC(translation_matrix,i,i))/(2*Pi*1_i) * num_sites_).real()*precision)/precision;
            if (momentum == -num_sites_/2) { momentum = num_sites_/2;}
            auto rho = std::round(eltC(rho_matrix,i,i).real()*precision)/precision;
            if (rho == -0) { rho = 0;}
            result.push_back( std::vector<Real>({ eltC(energy_matrix,i,i).real(), momentum, rho }) );
        }
        struct {
            bool operator()(std::vector<Real> a, std::vector<Real> b) const { return a.at(0) < b.at(0); }
        } compare;
        std::sort(result.begin(), result.end(), compare);
        printf("\nSpectrum:  {energy, momentum, rho}\n");
        PrintMatrix(result);
        DumpMatrix(result, an_path_);
    }

    // Mode 2
    void AnalyzeWithoutRho() {
        if (not std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            printf("\n> No progress file available\n");
            return;
        }
        if (not std::filesystem::exists(m_directory_)) {
            std::filesystem::create_directory(m_directory_);
        }
        if (not std::filesystem::exists(an_directory_)) {
            std::filesystem::create_directory(an_directory_);
        }
        PrintJob(false);
        std::filesystem::remove(m_path_);
        auto energy_matrix = Measure("energy");
        auto translation_matrix = Measure("translation");

        // Uniformize the ITensor Index structure
        auto s = findInds(energy_matrix, {})(1);
        auto s_translation = findInds(translation_matrix, {})(1).noPrime();
        translation_matrix = translation_matrix * delta(s, s_translation) * delta(prime(s), prime(s_translation));

        auto m = energy_matrix/num_sites_ + translation_matrix;
        auto[U,D] = eigen(m);
        U.conj();
        energy_matrix = dag(U)*energy_matrix*prime(U);
        translation_matrix = dag(U)*translation_matrix*prime(U);
        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);
        std::vector<std::vector<Real>> result;
        int precision = 100000.0;
        for (int i=1;i<=num_states;i++) {
            auto momentum = std::round((log(eltC(translation_matrix,i,i))/(2*Pi*1_i) * num_sites_).real()*precision)/precision;
            if (momentum == -num_sites_/2) { momentum = num_sites_/2;}
            result.push_back( std::vector<Real>({ eltC(energy_matrix,i,i).real(), momentum }) );
        }
        struct {
            bool operator()(std::vector<Real> a, std::vector<Real> b) const { return a.at(0) < b.at(0); }
        } compare;
        std::sort(result.begin(), result.end(), compare);
        printf("\nSpectrum:  {energy, momentum}\n");
        PrintMatrix(result);
        DumpMatrix(result, an_path_);
    }

    // Mode 3
    void Energies() {
        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            dmrg_progress_.Read(progress_path_);
        } else {
            printf("\n> No progress file available\n");
            return;
        }
        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);
        std::vector<Real> energies;
        for (int i=0; i<num_states; i++) {
            energies.push_back(dmrg_progress_.Energies().at(i));
        }
        printf("{{%s},", coupling_str_);
        PrintVector(energies);
        print("},\n");

        // Still create .an file for convenience
        auto energy_matrix = Measure("energy");
        std::vector<std::vector<Real>> result;
        for (int i=1;i<=num_states;i++) {
            result.push_back( std::vector<Real>({ eltC(energy_matrix,i,i).real() }) );
        }
        DumpMatrix(result, an_path_);
    }

    // Mode 4
    void NormalizeEnergies() {
        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            dmrg_progress_.Read(progress_path_);
        } else {
            printf("\n> No progress file available\n");
            return;
        }
        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);
        std::vector<Real> energies;
        for (int i=0; i<num_states; i++) {
            energies.push_back(dmrg_progress_.Energies().at(i));
        }
        if (energies.size() > 2) {
            std::sort(energies.begin(), energies.end());
            std::vector<Real> result;
            result.push_back(0);
            result.push_back(1);
            auto gs_en = energies.at(0);
            auto gap = energies.at(1) - gs_en;
            for (int i=2; i<energies.size(); i++) {
                result.push_back((energies.at(i)-gs_en)/gap);
            }
            printf("{{%s},%g,", coupling_str_, gap);
//            printf("\n> Normalized energies:\n");
            PrintVector(result);
            print("},\n");
        } else {
            println("Fewer than 3 states available");
        }
    }

    // Mode 6
    void Repair() {
        dmrg_progress_.RecoverStates(progress_path_);
    }

    void Test() {
        MPO H = sites_.Hamiltonian(boundary_condition_, num_sites_, u_, couplings_);
        H = removeQNs(H);
        ITensor T = H(1);
        std::vector<Index> sites = {sites_(1)};
        for (int i=2;i<=num_sites_;i++) {
            T *= H(i);
            sites.push_back(sites_(i));
        }
//        sites = removeQNs(sites);
        auto [U,S,V] = svd(expHermitian(T,-1),IndexSet(sites),{"MaxDim=",num_states_});
        std::vector<Real> energies;
        int num_states = std::min(num_states_, (int) dim(inds(S)(1)));
        for (int i=1;i<=num_states;i++) {
            energies.push_back(-log(elt(S,i,i)));
        }
        PrintVector(energies);
    }
};

#endif
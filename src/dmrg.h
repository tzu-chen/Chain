#ifndef CHAINDMRG_DMRG_H
#define CHAINDMRG_DMRG_H

#include <filesystem>
#include <variant>
#include "itensor/all.h"
#include "chain.h"
#include <sstream>
#include "path.h"

using namespace itensor;

//#ifdef ITENSOR_USE_OMP
//#include "/usr/local/opt/libomp/include/omp.h"
//#endif


// Reads and writes data coming out of DMRG
template<typename SiteSetType>
class DMRGProgress {
public:
    SiteSetType sites_;
    // Each element below is the latest progress for a given state
    std::vector<int> num_sweeps_vec_;
    std::vector<int> max_dims_;
    std::vector<Real> ens_;
    std::vector<MPS> psis_;
    std::vector<MPO> Hs_;
    std::vector<MPS> psis_acted_by_rho_;

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

    std::vector<MPS> StatesActedByRho() {
        return psis_acted_by_rho_;
    }

    std::vector<Real> Energies() const {
        return ens_;
    }

    int NumDoneStates() {
        return (int) psis_.size() - 1;
    }

    void Write(const std::filesystem::path& p) const
    {
        std::string path = std::string(p);
        writeToFile(path + ".st", sites_);
        writeToFile(path + ".pgs", num_sweeps_vec_);
        writeToFile(path + ".md", max_dims_);
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
            readFromFile(path + ".pgs", num_sweeps_vec_);
            readFromFile(path + ".md", max_dims_);
            readFromFile(path + ".en", ens_);
            readFromFile(path + ".psi", psis_);
            readFromFile(path + ".H", Hs_);
        } catch (std::exception const& e) {
//            try{
//                readFromFile(path + ".pgs.bk", times_swepts_);
//                readFromFile(path + ".md.bk", maxdims_);
//                readFromFile(path + ".en.bk", ens_);
//                readFromFile(path + ".psi.bk", psis_);
//                readFromFile(path + ".H.bk", Hs_);
//            } catch (std::exception const& e) {
                println(e.what());
                num_sweeps_vec_.clear();
                max_dims_.clear();
                ens_.clear();
                psis_.clear();
                Hs_.clear();
            NextState();
//            }
        }

    }

    void WriteRho(const std::filesystem::path& p) const
    {
        std::string path = std::string(p);
        writeToFile(path + ".rho", psis_acted_by_rho_);
    }

    void ReadRho(const std::filesystem::path& p)
    {
        std::string path = std::string(p);
        try{
            readFromFile(path + ".rho", psis_acted_by_rho_);
        } catch (std::exception const& e) {
            println(e.what());
            psis_acted_by_rho_.pop_back();
        }

    }
};

// Examples of SiteSetType are Golden, Haagerup, HaagerupQ
template<typename SiteSetType>
class DMRG {
    std::string site_type_; // eg. golden, haagerup, haagerupq
    std::string boundary_condition_; // "p" for periodic, "o" for open, "s" for sine-squared deformed
    std::string job_; // site_name_ + boundary_condition_

    // Chain parameters
    int num_sites_;
    Real theta_; // angle between K (identity) and J (rho)
    Real phi_;
//    Real k_, j_, u_, m_; // K (identity), J (rho), M (a * rho), U (penalty)
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
    double ortho_weight;

    // DMRG design
    float precision_; // Finish DMRG of current state if energy changes by less than precision_/num_sites_ for num_reps_til_stable_ reps
    float gs_precision_;
    float es_precision_;
    int num_reps_til_stable_;
    int init_num_sweeps_per_rep_;
    int num_sweeps_per_rep_;
    int num_states_; // total number of states to run DMRG

    // File I/O
    // Modify prefix_ if storing in custom location
    std::filesystem::path prefix_ = kPath;
//    std::filesystem::current_path().parent_path();
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
//        for (int i=0;i<couplings_.size();i++) {
//            println(couplings_.at(i));
//        }
//        PrintData(couplings_);
//        theta_ = std::get<6>(params);
//        phi_ = std::get<7>(params);
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
        // fixme: refactor so that we're passing the most basic variable type instead of a type alias
        // fixme: pass site instead of SiteSet

        // DMRG parameters
        gs_max_bond_dim_ = es_max_bond_dim_;
        noise_ = 0.0;
        ortho_weight = 1000.0;

        // DMRG design
        init_num_sweeps_per_rep_ = 5;
        num_sweeps_per_rep_ = 1;
        num_reps_til_stable_ = 3;
        gs_precision_ = std::pow(es_precision_, 1);

        // Chain parameters
//        k_ = cos(theta_ * Pi);
//        j_ = sin(theta_ * Pi) * cos(phi_ * Pi);
//        m_ = sin(theta_ * Pi) * sin(phi_ * Pi);
//        if (std::abs(k_) < 1E-5) {
//            k_ = 0;
//        }
//        if (std::abs(j_) < 1E-5) {
//            j_ = 0;
//        }
//        if (std::abs(m_) < 1E-5) {
//            m_ = 0;
//        }

        // File I/O
        filename_ = format("%s_%s_%d_%d_%g_%g_%g_%s_%d", site_type_, boundary_condition_, num_sites_, gs_max_bond_dim_, svd_cutoff_, gs_precision_, coupling_str_, u_, charge_);
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

    // Adjust parameters depending on whether simulating ground state
    void AdjustDMRGParams() {
        if (IsSimulatingGroundState()) {
            max_bond_dim = gs_max_bond_dim_;
            precision_ = gs_precision_;
        } else {
            max_bond_dim = es_max_bond_dim_;
            precision_ = es_precision_;
        }
    }

    // Reset min energy and counters to simulate new state
    void ResetDMRG() {
        min_en_ = 1000000;
        cnt_ = num_reps_til_stable_;
        num_sweeps_ = 0;
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

    // If QN conserving, create charge-neutral initial state and change center site to specified charge
    MPS DefaultInitState() {
        auto pre_init_state = InitState(sites_,"0");
        int center = (num_sites_ + 1) / 2;
        pre_init_state.set(center, std::to_string(charge_));
        MPS init_state = MPS(pre_init_state);
        return init_state;
    }

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

    // Simulate states by DMRG
    void Simulate() {
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
            // Job was run before
            // Read progress
            dmrg_progress_.Read(progress_path_);
            sites_ = dmrg_progress_.Sites();
            if (dmrg_progress_.NumDoneStates() >= num_states_) {
                printf("\n> Job already complete\n> Requested number of states: %d \n> Completed number of states: %d\n\n", num_states_,
                       dmrg_progress_.NumDoneStates());
                return;
            }
        } else {
            // New job
            dmrg_progress_.NextState();
            dmrg_progress_.SetSites(sites_);
            dmrg_progress_.Write(progress_path_);

        }

        MPO H = sites_.Hamiltonian(boundary_condition_, num_sites_, u_, couplings_);
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
        for (;;) {
            if (dmrg_progress_.num_sweeps_vec_.back() == 0 && hot_start != 0) {
                // Not yet run
                // Begin DMRG with warmup
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);

                AdjustDMRGParams();
                ResetDMRG();

                PrintJob(true);
                // DMRG Warmup
                bond_dim_ = 200;
                auto sweeps = Sweeps(init_num_sweeps_per_rep_);
                sweeps.maxdim() = 10,20,40,80,200;
                sweeps.cutoff() = std::max(1E-5, svd_cutoff_),std::max(1E-6, svd_cutoff_),std::max(1E-7, svd_cutoff_),std::max(1E-8, svd_cutoff_),std::max(1E-9, svd_cutoff_),std::max(1E-10, svd_cutoff_),std::max(1E-11, svd_cutoff_),std::max(1E-12, svd_cutoff_),svd_cutoff_;
                sweeps.niter() = 2;
                sweeps.noise() = noise_;
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), init_state_, sweeps, {"Quiet", true, "Weight=", ortho_weight});
                num_sweeps_ += init_num_sweeps_per_rep_;
                dmrg_progress_.Update(num_sweeps_, bond_dim_, en_, psi_, H);
                dmrg_progress_.Write(progress_path_);
            } else if (dmrg_progress_.num_sweeps_vec_.back() == 0 && hot_start == 0) {
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);

                AdjustDMRGParams();
                ResetDMRG();

                PrintJob(true);
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
                num_sweeps_ += init_num_sweeps_per_rep_;
                dmrg_progress_.Update(num_sweeps_, bond_dim_, en_, psi_, H);
                dmrg_progress_.Write(progress_path_);
            } else {
                state_name_ = StateName();
                printf("\n> Simulation: %s\n", state_name_);
            }

            // Post-warmup DMRG runs
            AdjustDMRGParams();
            for (;;) {
                // Read progress
                std::tie(num_sweeps_, bond_dim_, en_, psi_, H) = dmrg_progress_.All();

                // Dump entanglement entropy curve
                // Only for simulation of ground state
                if (IsSimulatingGroundState()) {
                    DumpEE(num_sites_, CalcEE(psi_, num_sites_), ee_path_);
                }
                // Gradually increase bond dimension each run until maximum
                if (bond_dim_ <= max_bond_dim - 200) {
                    bond_dim_ += 200;
                }
                PrintJob(true);
                // DMRG
                auto sw = Sweeps(num_sweeps_per_rep_);
                sw.maxdim() = bond_dim_;
                sw.cutoff() = svd_cutoff_;
                sw.niter() = 2;
                sw.noise() = noise_;
                std::tie(en_, psi_) = dmrg(H, dmrg_progress_.DoneStates(), psi_, sw, {"Quiet=", true, "Weight=", ortho_weight});
                // Decide whether energy has stabilized
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

            // Completed current state
            DumpEnergy(dmrg_progress_.NumDoneStates(), en_ / num_sites_, en_path_);
            printf("\n> E/L of %s: %g\n\n", state_name_, en_ / num_sites_);

            // If using sine-squared deformed to hot start for periodic
            // Make appropriate transition from SSD to PBC while using SSD result as initial state for PBC
            if (hot_start == 1) {
                H = sites_.Hamiltonian("p", num_sites_, u_, couplings_);
                init_state_ = psi_;
                dmrg_progress_ = DMRGProgress<SiteSetType>();
                hot_start = 0;
            } else if (hot_start == 0) {
                init_state_ = DefaultInitState();
                hot_start = -1;
            }

            dmrg_progress_.NextState(); // Even if no next state, this marks completion of current state
            dmrg_progress_.Write(progress_path_);

            // Check whether to simulate next state
            if (dmrg_progress_.NumDoneStates() >= num_states_) {
                break;
            }
        }
    }

    // Print the real parts of std::vector<Cplx> object into a Mathematica-compatible array.
    void PrintVector(std::vector<Cplx> vector) {
        print("{");
        if (vector.size() > 0) {
            print(vector.at(0).real());
            for (int i=1; i<vector.size(); i++) {
                print(",");
                print(vector.at(i).real());
            }
        }
        print("}");
    }

    // Print std::vector<Real> object into a Mathematica-compatible array.
    void PrintVector(std::vector<Real> vector) {
        print("{");
        if (vector.size() > 0) {
            print(vector.at(0));
            for (int i = 1; i < vector.size(); i++) {
                print(",");
                print(vector.at(i));
            }
        }
        print("}");
    }

    // Check if the translation eigenvalue is a num_sites_ root of unity to within numerical precision
    // If close to the nth power of the generating root of unity, make it n and return filtered vector
    std::vector<Cplx> FilterTranslation(std::vector<Cplx> translation_eigenvalues) {
        std::vector<Cplx> result;
        for (int i=0; i<translation_eigenvalues.size(); i++) {
//            Cplx tmp = pow(translation_eigenvalues.at(i), num_sites_);
//            if (abs(tmp.imag()) < 1E-3 && abs(tmp.real()-1) < 1E-3) {
//                result.push_back(tmp);
//            }
            Cplx tmp = log(translation_eigenvalues.at(i)) / (2 * Pi * 1_i) * num_sites_;
            if (abs(tmp.imag()) < 1E-1 && abs(tmp.real()-round(tmp.real())) < 1E-1) {
                if (round(tmp.real()) == -0) {
                    result.push_back( 0 );
                } else if (round(tmp.real()) == - num_sites_/2) {
                    result.push_back( num_sites_/2 );
                } else {
                    result.push_back( round(tmp.real()) );
                }
            }
        }
        return result;
    }

    // Check if the rho eigenvalue is one of the possibilities to within numerical precision
    // If close to a possibility p, make it p and return filtered vector
    std::vector<Cplx> FilterRho(std::vector<Cplx> rho_eigenvalues, std::vector<Cplx> possibilities) {
        std::vector<Cplx> result;
        for (int i=0; i<rho_eigenvalues.size(); i++) {
            Cplx tmp = rho_eigenvalues.at(i);
            println(tmp);
//            if (abs(tmp.imag()) < 1E-3) {
                for (int j=0; j<possibilities.size(); j++) {
                    if (abs(tmp - possibilities.at(j)) < 1E-1) {
                        result.push_back(possibilities.at(j));
                    }
                }
//            }
        }
        return result;
    }

    std::vector<Cplx> OrderedAppend(std::vector<Cplx> old_vector, std::vector<Cplx> new_vector) {
        std::vector<Cplx> old_vector_copy = old_vector;
        std::vector<Cplx> new_ones;
        for (int i=0; i<new_vector.size(); i++) {
            bool is_contained = false;
            for (int j=0; j<old_vector_copy.size(); j++) {
                if (new_vector.at(i) == old_vector_copy.at(j)) {
                    old_vector_copy.erase(old_vector_copy.begin() + j);
                    is_contained = true;
                    break;
                }
            }
            if (!is_contained) {
                new_ones.push_back(new_vector.at(i));
            }
        }
        for (int i=0; i<new_ones.size(); i++) {
            old_vector.push_back(new_ones.at(i));
        }
        return old_vector;
    }

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
    }

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

//    // Measure matrix elements of translation operator and rho defect operator
//    void Analyze() {
//        if (not std::filesystem::exists(m_directory_)) {
//            std::filesystem::create_directory(m_directory_);
//        }
//
//        PrintJob(false);
//
//        // Default svd_cutoff for applyMPO(density matrix variant) is 1E-12
//        // Setting to larger value will speed up the program significantly at the cost of accuracy of measurements
//        float svd_cutoff = 1E-8;
//
//        // Variable declaration
//        MPS psi_translated;
//        MPS psi_acted_by_rho;
//        std::vector<MPS> states_translated;
//        std::vector<MPS> states_acted_by_rho;
//
//        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
//            dmrg_progress_.Read(progress_path_);
//        } else {
//            printf("\n> No progress file available\n");
//            return;
//        }
//
//        auto states = dmrg_progress_.DoneStates();
//        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);
//
//        // fixme: refactor Haagerup/HaagerupQ branching
//        // Perform change of basis if appropriate in order to use F symbols
//        // Rotate HaagerupQ to Haagerup
//        // Golden -> Golden, Haagerup, HaagerupQ -> Haagerup
//        SiteSet sites;
//        if (site_type_ == "golden" or site_type_ == "haagerup") {
//            sites = sites_;
//            for (int i = 0; i < num_states; i++) {
//                auto new_psi = MPS(sites);
//                auto psi = states.at(i);
//                for (auto j : range1(num_sites_)) {
//                    new_psi.set(j, psi(j) * delta(siteIndex(psi, j), sites(j)) );
//                }
//                states.at(i) = new_psi;
//            }
//        } else {
//            sites = Haagerup(num_sites_);
//            for (int i = 0; i < num_states; i++) {
//                auto states_no_qn = removeQNs(states.at(i));
//                states.at(i) = Z3FourierTransform(states_no_qn, sites);
//            }
//        }
//
//        auto translate_op = TranslationOp(sites); // periodic MPS
//        auto rho_op = RhoOp(sites, site_type_); // dangling periodic MPS
//        auto id_op = IdentityOp(sites, rho_op); // dangling periodic MPS with matching indices with rho_op
//
//        // Dangling bond indices at the edges of MPO
//        // fixme: 36?
//        auto left_dangling_ind = Index(36, "Site");
//        auto right_dangling_ind = Index(36, "Site");
//
//        // fixme: Can we use OpenMP for this?
//        for (int i=0;i<num_states;i++) {
//            states_translated.push_back(MPS());
//            states_acted_by_rho.push_back(MPS());
//        }
//
//        {
////            auto np = omp_get_num_threads();
////            println(np);
////            #pragma omp parallel for default(shared) private(i, psi_translated, psi_acted_by_rho)
//            for (int i = 0; i < num_states; i++) {
//                // Initialize
//                psi_translated = MPS(states.at(i));
//                psi_acted_by_rho = MPS(states.at(i));
//
//                // Act by translation
//                psi_translated = applyMPO(translate_op, psi_translated);
////            states_translated.push_back(psi_translated);
//                states_translated.at(i) = psi_translated;
//
//                // Act by rho defect in two steps
//                // First act by dangling rho operator
//                // Then act by dangling identity operator
//                psi_acted_by_rho = applyMPO(AugmentMPO(rho_op, left_dangling_ind, right_dangling_ind),
//                                            AugmentMPS(psi_acted_by_rho, left_dangling_ind, right_dangling_ind),
//                                            {"Cutoff", svd_cutoff}
//                );
//                psi_acted_by_rho = applyMPO(AugmentMPO(id_op, left_dangling_ind, right_dangling_ind),
//                                            psi_acted_by_rho, {"Cutoff", svd_cutoff});
////            states_acted_by_rho.push_back(psi_acted_by_rho);
//                states_acted_by_rho.at(i) = psi_acted_by_rho;
//
//                // println(innerC(AugmentMPS(psiR, left_extra, right_extra),step2));
//
////            psiR = applyMPO(TranslationOp(sites, true), psiR);
////            if (site_name_ == "golden") {
////                GoldenFData f_data = GoldenFData();
////                ActLocal(psiR, f_data.RhoDefectCell(sites(1), sites(2)),1);
////            } else if (site_name_ == "haagerup") {
////                HaagerupFData f_data = HaagerupFData();
////                ActLocal(psiR, f_data.RhoDefectCell(sites(1), sites(2)),1);
////            }
////            psiR = applyMPO(TranslationOp(sites, false), psiR);
//            }
//        }
//
//        // Create ITensors for
//        // En: Energies
//        // OpT: Matrix elements between psi and psi_translated
//        // OpR: Matrix elements between psi and psi_acted_by_rho
//        auto s = Index(num_states);
//        auto sP = prime(s);
//        auto en_matrix = ITensor(dag(s), sP);
//        auto translation_matrix = ITensor(dag(s), sP);
//        auto rho_matrix = ITensor(dag(s), sP);
//
//        Real en_shift = 0;
//        // Uncomment if shift by ground state energy
//        // en_shift = dmrg_progress.Energies().at(0);
//
//        for (int i=0; i < num_states; i++) {
//            en_matrix.set(s(i+1), sP(i+1), dmrg_progress_.Energies().at(i) - en_shift);
//            for (int j=0; j < num_states; j++) {
//                translation_matrix.set(s(i+1), sP(j+1), innerC(states.at(i), states_translated.at(j)));
//                rho_matrix.set(s(i+1), sP(j+1), innerC(AugmentMPS(states.at(i), left_dangling_ind, right_dangling_ind), states_acted_by_rho.at(j)));
//            }
//        }
//        DumpMathematicaAll(num_states, en_matrix, translation_matrix, rho_matrix, m_path_);
//
//        // Diagonalize translation
//        auto [UT,translation_diag] = eigen(translation_matrix);
//        // Diagonalize rho
//        auto [UR,rho_diag] = eigen(rho_matrix);
//
////        printf("\n> Ordered set of energies:\n");
////        PrintData(en_matrix);
////        printf("\n> Unordered set of translation eigenvalues:\n");
////        PrintData(translation_diag);
////        printf("\n> Unordered set of rho eigenvalues:\n");
////        PrintData(rho_diag);
//
//        // Analysis of translation and rho eigenvalues assuming that the states are simulated well enough
//        // Otherwise, better work with the Mathematica .m file created above
//        // fixme: make method of sites?
//        std::vector<Cplx> rho_possibilities;
//        if (site_type_ == "golden") {
//            rho_possibilities = { (1+sqrt(5))/2, (1-sqrt(5))/2 };
//        } else {
//            rho_possibilities = { (3+sqrt(13))/2, (3-sqrt(13))/2, 1, -1 };
//        }
//        std::vector<Cplx> translation_eigenvalues;
//        std::vector<Cplx> rho_eigenvalues;
//        for (int i=0; i<num_states; i++) {
//            auto translation_submatrix = ITensor(dag(s), sP);
//            auto rho_submatrix = ITensor(dag(s), sP);
//            for (int j=1; j<=i+1; j++) {
//                for (int k=1; k<=i+1; k++) {
//                    translation_submatrix.set(s(j), sP(k), eltC(translation_matrix, j, k));
//                    rho_submatrix.set(s(j), sP(k), eltC(rho_matrix, j, k));
//                }
//            }
//            auto [UT,translation_sub_eigenvalues] = eigen(translation_submatrix);
//            auto [UR,rho_sub_eigenvalues] = eigen(rho_submatrix);
//            std::vector<Cplx> translation_eigenvalues_tmp;
//            std::vector<Cplx> rho_eigenvalues_tmp;
//            for (int j=1; j<=i+1; j++) {
//                translation_eigenvalues_tmp.push_back(eltC(translation_sub_eigenvalues, j, j));
//                rho_eigenvalues_tmp.push_back(eltC(rho_sub_eigenvalues, j, j));
//            }
//            translation_eigenvalues_tmp = FilterTranslation(translation_eigenvalues_tmp);
//            rho_eigenvalues_tmp = FilterRho(rho_eigenvalues_tmp, rho_possibilities);
//            translation_eigenvalues = OrderedAppend(translation_eigenvalues, translation_eigenvalues_tmp);
//            rho_eigenvalues = OrderedAppend(rho_eigenvalues, rho_eigenvalues_tmp);
//        }
//
////        printf("{{%s},", coupling_str_);
////        PrintVector(rho_eigenvalues);
////        print("},\n");
//
//        printf("\n> Energies:\n");
//        std::vector<Real> energies = dmrg_progress_.Energies();
//        energies.pop_back();
//        PrintVector(energies);
//        printf("\n\n> Translation eigenvalues:\n");
//        PrintVector(translation_eigenvalues);
//        printf("\n\n> rho eigenvalues:\n");
//        PrintVector(rho_eigenvalues);
//        printf("\n\n");
//    }

    // Measure matrix elements of translation operator and rho defect operator
    void Measure(std::string observable) {
        // Default svd_cutoff for applyMPO(density matrix variant) is 1E-12
        // Setting to larger value will speed up the program significantly at the cost of accuracy of measurements
        float svd_cutoff = 1E-5;

        // Variable declaration
        MPS psi_acted;
        std::vector<MPS> states_acted;

        if (std::filesystem::exists(progress_directory_ / (filename_ + ".pgs"))) {
            dmrg_progress_.Read(progress_path_);
        } else {
            printf("\n> No progress file available\n");
            return;
        }

        int num_states = std::min(dmrg_progress_.NumDoneStates(), num_states_);

        if (observable == "energy") {
            std::vector<Real> eigenvalues;
            auto s = Index(num_states);
            auto sP = prime(s);
            auto matrix = ITensor(dag(s), sP);
            for (int i=0;i<num_states;i++) {
                eigenvalues.push_back(dmrg_progress_.Energies().at(i));
                matrix.set(i+1,i+1,dmrg_progress_.Energies().at(i));
            }
            printf("\n> %s eigenvalues:\n", observable);
            PrintVector(eigenvalues);
            printf("\n\n");
            DumpMathematicaSingle(observable, num_states, matrix, m_path_);
        } else {
            auto states = dmrg_progress_.DoneStates();

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
                    for (auto j : range1(num_sites_)) {
                        new_psi.set(j, psi(j) * delta(siteIndex(psi, j), sites(j)));
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

//            for (int i = 0; i < num_states; i++) {
//                states_acted.push_back(MPS());
//            }

            if (observable == "translation") {
                auto translation_op = TranslationOp(sites); // periodic MPS
                for (int i = 0; i < num_states; i++) {
                    psi_acted = MPS(states.at(i));
                    println("Checkpoint 1");
//                    psi_acted = applyMPO(translation_op, psi_acted, {"Method", "Fit", "Cutoff", svd_cutoff, "Verbose", true});
                    for (int j=1; j<num_sites_; j++) {
                        Swap(psi_acted, sites, j);
                    }
                    println("Checkpoint 2");
//                    states_acted.at(i) = psi_acted;
                    states_acted.push_back(psi_acted);
                }
            }

            auto left_dangling_ind = Index(36, "Site");
            auto right_dangling_ind = Index(36, "Site");

            if (observable == "rho") {
                if (std::filesystem::exists(progress_directory_ / (filename_ + ".rho"))) {
                    dmrg_progress_.ReadRho(progress_path_);
                    for (int i=0; i<std::min(num_states, (int) dmrg_progress_.StatesActedByRho().size()); i++) {
                        states_acted.push_back(dmrg_progress_.StatesActedByRho().at(i));
                    }
                }
                auto rho_op = RhoOp(sites, site_type_); // dangling periodic MPS
                auto id_op = IdentityOp(sites, rho_op); // dangling periodic MPS with matching indices with rho_op
                for (int i = states_acted.size(); i < num_states; i++) {
                    // Initialize
                    psi_acted = MPS(states.at(i));
                    // fixme
//                    auto left_dangling_ind = Index(6, "Site");
//                    auto left_left_dangling_ind = Index(6, "Site");
//                    auto augmented_psi = AugmentMPSZipper(psi_acted, left_dangling_ind, left_left_dangling_ind);
//                    PrintData(innerC(augmented_psi, augmented_psi));
//                    return;
                    // Act by rho defect in two steps
                    // First act by dangling rho operator
                    // Then act by dangling identity operator
                    println("Checkpoint 1");
                    psi_acted = applyMPO(AugmentMPO(rho_op, left_dangling_ind, right_dangling_ind),
                                         AugmentMPS(psi_acted, left_dangling_ind, right_dangling_ind),
                                         {"Method", "Fit", "Cutoff", svd_cutoff, "Nsweep", 2, "Verbose", true}
                    );
//                    psi_acted = applyMPO(AugmentMPO(id_op, left_dangling_ind, right_dangling_ind),
//                                         AugmentMPS(psi_acted, left_dangling_ind, right_dangling_ind),
//                                         {"Cutoff", svd_cutoff, "Verbose", true}
//                    );
                    println("Checkpoint 2");
                    psi_acted = applyMPO(AugmentMPO(id_op, left_dangling_ind, right_dangling_ind),
                                         psi_acted, {"Method", "Fit", "Cutoff", svd_cutoff, "Nsweep", 1, "Verbose", true});
                    println("Checkpoint 3");
//                    states_acted.at(i) = psi_acted;
                    states_acted.push_back(psi_acted);

                    dmrg_progress_.psis_acted_by_rho_.push_back(psi_acted);
                    dmrg_progress_.WriteRho(progress_path_);
                }
            }

//            println("Acted operator");

            // Create ITensor
            auto s = Index(num_states);
            auto sP = prime(s);
            auto matrix = ITensor(dag(s), sP);

            Real en_shift = 0;
            // Uncomment if shift by ground state energy
            // en_shift = dmrg_progress.Energies().at(0);

            for (int i = 0; i < num_states; i++) {
                if (observable == "energy") {
                    matrix.set(s(i+1), sP(i+1), dmrg_progress_.Energies().at(i) - en_shift);
                }
                for (int j = 0; j < num_states; j++) {
                    if (observable == "translation") {
                        matrix.set(s(i+1), sP(j+1), innerC(states.at(i), states_acted.at(j)));
                    }
                    if (observable == "rho") {
                        matrix.set(s(i+1), sP(j+1),
                                   innerC(AugmentMPS(states.at(i), left_dangling_ind, right_dangling_ind),
                                          states_acted.at(j)));
                    }
                }
            }
            DumpMathematicaSingle(observable, num_states, matrix, m_path_);
//            println("Computed inner product");

//            // Diagonalize
//            auto [U,diag] = eigen(matrix);
//            printf("\n> %s eigenvalues (unordered):\n", observable);
//            PrintData(diag);

            // Analysis of translation and rho eigenvalues assuming that the states are simulated well enough
            // Otherwise, better work with the Mathematica .m file created above

            std::vector<Cplx> eigenvalues;
            for (int i = 0; i < num_states; i++) {
                auto submatrix = ITensor(dag(s), sP);
                for (int j = 1; j <= i+1; j++) {
                    for (int k = 1; k <= i+1; k++) {
                        submatrix.set(s(j), sP(k), eltC(matrix, j, k));
                        submatrix.set(s(j), sP(k), eltC(matrix, j, k));
                    }
                }
                auto[U, sub_eigenvalues] = eigen(submatrix);
                std::vector<Cplx> eigenvalues_tmp;
                for (int j = 1; j <= i+1; j++) {
                    eigenvalues_tmp.push_back(eltC(sub_eigenvalues, j, j));
                }
                if (observable == "translation") {
                    eigenvalues_tmp = FilterTranslation(eigenvalues_tmp);
                }
                if (observable == "rho") {
                    std::vector<Cplx> rho_possibilities;
                    if (site_type_ == "golden") {
                        rho_possibilities = {(1 + sqrt(5)) / 2, (1 - sqrt(5)) / 2};
                    } else {
                        rho_possibilities = {(3 + sqrt(13)) / 2, (3 - sqrt(13)) / 2, 1, -1};
                    }
                    eigenvalues_tmp = FilterRho(eigenvalues_tmp, rho_possibilities);
                }
                eigenvalues = OrderedAppend(eigenvalues, eigenvalues_tmp);
            }
            printf("> %s eigenvalues:\n", observable);
            PrintVector(eigenvalues);
            printf("\n\n");
        }

//        printf("{{%s},", coupling_str_);
//        PrintVector(rho_eigenvalues);
//        print("},\n");

    }

    void AnalyzeNoRho() {
        if (not std::filesystem::exists(m_directory_)) {
            std::filesystem::create_directory(m_directory_);
        }
        PrintJob(false);
        CleanMathematica(m_path_);
        Measure("energy");
        Measure("translation");
    }

    void Analyze() {
        if (not std::filesystem::exists(m_directory_)) {
            std::filesystem::create_directory(m_directory_);
        }
        PrintJob(false);
        CleanMathematica(m_path_);
        Measure("energy");
        Measure("translation");
        Measure("rho");
    }
};

#endif
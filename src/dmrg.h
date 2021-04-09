#ifndef CHAINDMRG_DMRG_H
#define CHAINDMRG_DMRG_H
#include "itensor/all.h"
#include "chain.h"
using namespace itensor;

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
    int num_past_states;
//    std::string progress_path, ee_path, en_path;
//    std::string prefix = std::filesystem::current_path();
    std::filesystem::path prefix = std::filesystem::path("/home/tzuchen/CLionProjects/ChainDMRG");
    std::filesystem::path progress_path, ee_path, en_path;
    std::filesystem::path progress_directory = prefix / "pgs";
    std::filesystem::path ee_directory = prefix / "ee";
    std::filesystem::path en_directory = prefix / "en";
    std::string filename;
    SiteSetType sites;
public:
    // constructor
    explicit DMRG(char** argv){
        name_ = std::string(argv[1]);
        boundary_condition = std::string(argv[2]);
        if(boundary_condition == "p"){
            job = name_ + " PBC";
        }else if(boundary_condition == "o"){
            job = name_ + " OBC";
        }else if(boundary_condition == "s"){
            job = name_ + " SSD";
        }else{
            job = name_ + " SSD/PBC";
        }
        N = std::stoi(argv[3]);
        // fixme: refactor so that we're passing the most basic variable type instead of a type alias
        sites = SiteSetType(N, {"ConserveQNs=", true});
        gs_maxmaxdimension = std::stoi(argv[4]);
        cutoff = std::stof(argv[5]);
        tolerance = std::stof(argv[6]);
        theta = std::stof(argv[7]);
        init_num_sweeps = 5;
        num_sweeps = 1;
        maxmaxdimension = gs_maxmaxdimension;
        noise = 0.0;
        m = 3;
        gs_tolerance = std::pow(tolerance, 2);
        K = cos(theta * Pi);
        J = sin(theta * Pi);
        if(std::abs(K) < 1E-5){
            K = 0;
        }
        if(std::abs(J) < 1E-5){
            J = 0;
        }
        U = std::stof(argv[8]);
        num_past_states = std::stoi(argv[9]);
        filename = format("%s_%s_%d_%g_%g_%g", name_, boundary_condition, N, cutoff, theta, U);
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
        printf("\n> %s: N=%d maxmaxdim=%d cutoff=%g theta=%g K=%g J=%g U=%g\n",job,N,gs_maxmaxdimension,cutoff,theta,K,J,U);

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
        MPS init_state = InitState(sites,"r");
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
                std::tie(en,psi) = dmrg(H,dmrg_progress.past_states(),init_state,sweeps,{"Quiet",true,"Weight=",100.0});

                s += init_num_sweeps;
                dmrg_progress.update(s,maxdim,en,psi,H);
                dmrg_progress.putSites(sites);
                dmrg_progress.write(progress_path);
                printf("\n    > Times swept: %d\n      %s of %s\n      N=%d maxdim=%d cutoff=%g theta=%g K=%g J=%g U=%g\n",dmrg_progress.times_swepts.back(),state,job,N,maxdim,cutoff,theta,K,J,U);
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
                std::tie(en,psi) = dmrg(H,dmrg_progress.past_states(),init_state,sweeps,{"Quiet",true,"Weight=",100.0});

                s += init_num_sweeps;
                dmrg_progress.update(s,maxdim,en,psi,H);
                dmrg_progress.putSites(sites);
                dmrg_progress.write(progress_path);
                printf("\n    > Times swept: %d\n      %s of %s\n      N=%d maxdim=%d cutoff=%g theta=%g K=%g J=%g U=%g\n",dmrg_progress.times_swepts.back(),state,job,N,maxdim,cutoff,theta,K,J,U);
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
                std::tie(en,psi)=dmrg(H,dmrg_progress.past_states(),psi,sw,{"Quiet=",true,"Weight=",20.0});
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

                printf("\n    > Times swept: %d\n      %s of %s\n      N=%d maxdim=%d cutoff=%g theta=%g K=%g J=%g U=%g\n",dmrg_progress.times_swepts.back(),state,job,N,maxdim,cutoff,theta,K,J,U);
                if (count == 0){
                    break;
                }
            }

            dumpEnergy(dmrg_progress.num_past_states(), en/N, en_path);
            printf("\n> E/N of %s: %g\n\n",state,en/N);

            if(hot_start == 1){
                H = ConstructH(sites, "p", N, U, K, J);
                init_state = psi;
                dmrg_progress = DMRGProgress<SiteSetType>();
                hot_start = 0;
            }else if(hot_start == 0){
                init_state = InitState(sites,"r");
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

        DMRGProgress<SiteSetType> dmrg_progress;
        dmrg_progress.read(progress_path);

        auto states = dmrg_progress.States();
        sites = dmrg_progress.Sites();
        int len = std::min((int) states.size(), num_past_states+1);

        // Act translation/rho on psi to create psiT/psiR
        std::vector<MPS> pastT;
        std::vector<MPS> pastR;
        GoldenFData f_data = GoldenFData();
        for(int i=0;i<len;i++){
            psiT = MPS(states.at(i));
            psiR = MPS(states.at(i));

            // calculate the action of translation on wavefunction
//            ActGlobal(psiT, sites, &FData::SwapITensor); // Old method
            psiT = applyMPO(TranslationOp(sites), psiT);
            pastT.push_back(psiT);

            // calculate the action of rho on wavefunction
            psiR = applyMPO(RhoOp(sites), psiR);
            psiR = applyMPO(TranslationOp(sites, true), psiR);
            ActLocal(psiR, f_data.RhoDefect(sites(1), sites(2)),1);
            psiR = applyMPO(TranslationOp(sites, false), psiR);
            pastR.push_back(psiR);
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

        Real gs_en = dmrg_progress.Energies().at(0);

        for(int i=0;i<len;i++){
            En.set(s(i+1),sP(i+1),dmrg_progress.Energies().at(i) - gs_en);
            for(int j=0;j<len;j++){
                OpT.set(s(i+1),sP(j+1),Chop(innerC(states.at(i), pastT.at(j))));
                OpR.set(s(i+1),sP(j+1),Chop(innerC(states.at(i), pastR.at(j))));
            }
        }
        PrintData(En);
        PrintData(OpR);
        PrintData(OpT);


        // Diagonalize translation
        auto [UT,DT] = eigen(OpT);
        // Rotate basis for energies accordingly and create diagonal ITensor
        En = prime(UT) * En * dag(UT);
        std::vector<Real> diag_En;
        diag_En.reserve(len);
        for(int i=0;i<len;i++){
            diag_En.push_back(eltC(En, i+1, i+1).real());
        }
        En = diagITensor(diag_En, dag(s), prime(s));
        auto spins = DT.apply(Spin);
        PrintData(En);
        PrintData(spins);
    }
};
#endif //CHAINDMRG_DMRG_H

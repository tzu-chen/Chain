#include "itensor/all.h"
#include "analysis.cpp"

using namespace itensor;

std::vector<double> calcEE(MPS psi, int N){
    std::vector<double> SvNs;
    for(auto b=1;b<N;b++){
        psi.position(b);

        //SVD this wavefunction to get the spectrum
        //of density-matrix eigenvalues
        auto l = leftLinkIndex(psi,b);
        auto s = siteIndex(psi,b);
        auto [U,S,V] = svd(psi(b),{l,s});
        auto u = commonIndex(U,S);

        //Apply von Neumann formula
        //to the squares of the singular values
        Real SvN = 0.;
        for(auto n : range1(dim(u)))
            {
            auto Sn = elt(S,n,n);
            auto p = sqr(Sn);
            if(p > 1E-12) SvN += -p*log(p);
            }
        SvNs.push_back(SvN);
    }
    return SvNs;
}

void dumpEE(int N, std::vector<double> SvNs, std::string const& filename){
    FILE * eefile;
    eefile = fopen(filename.c_str(), "a");
    fprintf(eefile, "{\n");
    for(auto b=1;b<N;b++){
        fprintf(eefile, "{%d,  %.10f},\n",b,SvNs[b-1]);
    }
    fprintf(eefile, "},\n");
    fclose(eefile);
}

void dumpEnergy(int state, Real en, std::string const& filename){
    FILE * file;
    file = fopen(filename.c_str(), "a");
    fprintf(file, "{");
    fprintf(file, "%d",state);
    fprintf(file, ",");
    fprintf(file, "%.10f",en);
    fprintf(file, "},\n");
    fclose(file);
}

// std::string name(typename SiteSetType) {
//     if(SiteSetType == Golden){
//         return "Golden";
//     }else{
//         return "Haagerup";
//     }
// }

template<typename SiteSetType>
class DMRG {
public:
    int run(char **argv) {
        std::string name_ = std::string(argv[1]);
        std::string boundary_condition = std::string(argv[2]);
        std::string job;
        if(boundary_condition == "p"){
            job = name_ + " PBC";
        }else if(boundary_condition == "o"){
            job = name_ + " OBC";
        }else if(boundary_condition == "s"){
            job = name_ + " SSD";
        }else{
            job = name_ + " SSD/PBC";
        }
        int N = std::stoi(argv[3]);
        int init_num_sweeps = 5;
        int num_sweeps = 1;
        int gs_maxmaxdimension = std::stoi(argv[4]);
        int maxmaxdimension = gs_maxmaxdimension;
        double cutoff = std::stof(argv[5]);
        double noise = 0.0;
        int m = 3;
        float tolerance = std::stof(argv[6]);
        float gs_tolerance = std::pow(tolerance,2);
        Real theta = std::stof(argv[7]);
        Real K = cos(theta*Pi);
        Real J = sin(theta*Pi);
        if(std::abs(K) < 1E-5){
            K = 0;
        }
        if(std::abs(J) < 1E-5){
            J = 0;
        }
        Real U=std::stof(argv[8]);
        int num_past_states = std::stoi(argv[9]);

        std::string const& progress_path = format("%s/pgs/%s_%s_%d_%g_%g_%g",prefix,name_,boundary_condition,N,cutoff,theta,U);
        std::string const& ee_path = format("%s/ee/%s_%s_%d_%g_%g_%g.ee",prefix,name_,boundary_condition,N,cutoff,theta,U);
        std::string const& en_path = format("%s/en/%s_%s_%d_%g_%g_%g.en",prefix,name_,boundary_condition,N,cutoff,theta,U);

        printf("\n> %s: N=%d maxmaxdim=%d cutoff=%g theta=%g K=%g J=%g U=%g\n",job,N,gs_maxmaxdimension,cutoff,theta,K,J,U);

        SiteSetType sites = SiteSetType(N,{"ConserveQNs=",true});

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

        if(fileExists(progress_path + ".pgs")){
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
                return 0;
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
        return 0;
    }
};

int main(int argc, char** argv){
    if(argv[10] != nullptr){
        return analyze(argc,argv);
    }

    if(std::string(argv[1]) == "golden"){
        return DMRG<Golden>().run(argv);
    }else if(std::string(argv[1]) == "haagerup"){
        return DMRG<Haagerup>().run(argv);
    }
}
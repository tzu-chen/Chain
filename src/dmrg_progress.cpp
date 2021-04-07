#include <utility>

#include "itensor/all.h"

using namespace itensor;

// const std::string prefix = "/n/holyscratch01/yin_lab/Users/yhlin";
const std::string prefix = "/home/tzuchen/CLionProjects/ChainDMRG";

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

    void write(const std::string& path) const
    {
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

    void read(const std::string& path)
    {
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
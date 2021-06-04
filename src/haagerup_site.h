#include "itensor/all.h"
#include "f_data.h"

using namespace itensor;

class HaagerupSite;
using Haagerup=BasicSiteSet<HaagerupSite>;
// 1    1
// 2    a
// 3    b=a^2
// 4    r=\rho
// 5    ar
// 6    br


using  std::sqrt;
class HaagerupSite
{
    const double zetainv=(sqrt(13)-3)/2;
    const double sqrtzetainv=sqrt(zetainv);
    const double x=(2-sqrt(13))/3;
    const double z=(1+sqrt(13))/6;
    const double y1=(5-sqrt(13)-sqrt(6*(1+sqrt(13))))/12;
    const double y2=(5-sqrt(13)+sqrt(6*(1+sqrt(13))))/12;
    Index s;
    HaagerupFData haagerup_f_data;
public:
    explicit HaagerupSite(Index const& I);
    explicit HaagerupSite(Args const& args = Args::global());

    Index
    index() const;

    IndexVal
    state(std::string const& state);
    
    ITensor
    proj(int i) const;

    ITensor
    FF(int projector, int x, int y) const;

//    ITensor
//    FF(int i) const;
//
//    ITensor
//    Frr(int i) const;
//
//    ITensor
//    Frar(int i) const;
//
//    ITensor
//    Frbr(int i) const;
    
    ITensor
    op(std::string const& opname,Args const& args = Args::global()) const;
};



// Hamiltonian appropriate for the boundary condition, number of sites, and couplings
// Polymorphic function specified by SiteSetType = Haagerup
MPO Hamiltonian(Haagerup sites, std::string boundary_condition, int num_sites, Real U, Real K, Real J, Real M);



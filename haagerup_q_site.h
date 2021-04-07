#ifndef  HAAGERUP_Q_SITE_H
#define HAAGERUP_Q_SITE_H

#include "itensor/all.h"

using namespace itensor;

class HaagerupSite;
using Haagerup=BasicSiteSet<HaagerupSite>;
// 1    1
// 2    a
// 3    b=a^2
// 4    r=\rho
// 5    ar
// 6    br
class HaagerupSite
{
    Index s;
public:
    explicit HaagerupSite(Index const& I);
    explicit HaagerupSite(Args const& args = Args::global());
    Index
    index() const;

    IndexVal
    state(std::string const& state);
    
    ITensor
    m(int i) const;
    
    ITensor
    q(int i) const;

    ITensor
    qr(int i, int j) const;
    
    ITensor
    op(std::string const& opname,Args const& args = Args::global()) const;
};

MPO ConstructH(Haagerup sites, std::string boundary_condition, int N, Real U, Real K, Real J);

#endif
#ifndef GOLDEN_SITE_H
#define GOLDEN_SITE_H

#include <utility>
#include "itensor/all.h"
#include "f_data.h"

using namespace itensor;

class GoldenSite;
using Golden=BasicSiteSet<GoldenSite>;
// 1    1
// 2    t
class GoldenSite
{
    Index s;
public:
    GoldenFData golden_f_data = GoldenFData();
    const double kPhiInv = golden_f_data.kPhiInv;
    const double kSqrtPhiInv = golden_f_data.kSqrtPhiInv;

    explicit GoldenSite(Index  I);
    explicit GoldenSite(Args const& args = Args::global());
    Index
    index() const;

    IndexVal
    state(std::string const& state);
    
    ITensor
    proj(int i) const;
    
    ITensor
    FF() const;
    
    ITensor
    op(std::string const& opname,Args const& args = Args::global()) const;
};

// Returns Hamiltonian appropriate for the boundary condition, number of sites, and couplings
MPO ConstructH(const Golden& sites, const std::string& boundary_condition, int N, Real U, Real K, Real J);
#endif
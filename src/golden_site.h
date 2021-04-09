#ifndef GOLDEN_SITE_H
#define GOLDEN_SITE_H

#include "itensor/all.h"
#include "f_data.h"

using namespace itensor;

// 1    1
// 2    t
class GoldenSite
{
public:
    explicit GoldenSite(Index  I);
    explicit GoldenSite(Args const& args = Args::global());
    [[nodiscard]] Index
    index() const;

    IndexVal
    state(std::string const& state);
    
    [[nodiscard]] ITensor
    proj(int i) const;
    
    [[nodiscard]] ITensor
    FF() const;
    
    [[nodiscard]] ITensor
    op(std::string const& opname,Args const& args = Args::global()) const;

private:
    Index s;
    GoldenFData golden_f_data;
    const double kPhiInv = golden_f_data.kPhiInv;
    const double kSqrtPhiInv = golden_f_data.kSqrtPhiInv;
};
using Golden=BasicSiteSet<GoldenSite>;

// Returns Hamiltonian appropriate for the boundary condition, number of sites, and couplings
MPO ConstructH(const Golden& sites, const std::string& boundary_condition, int N, Real U, Real K, Real J);

#endif //GOLDEN_SITE_H
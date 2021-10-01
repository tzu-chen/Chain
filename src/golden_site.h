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

    Index index() const;

    IndexVal state(std::string const& state);

    ITensor proj(int i) const;
    
    ITensor FF() const;
    
    ITensor op(std::string const& opname,Args const& args = Args::global()) const;

private:
    Index s_;
    GoldenFData golden_f_data_;
    const double phi_inv_ = golden_f_data_.phi_inv_;
    const double sqrt_phi_inv_ = golden_f_data_.sqrt_phi_inv_;
};

class Golden: public BasicSiteSet<GoldenSite> {
public:
    Golden() : BasicSiteSet<GoldenSite>() {};
    explicit Golden(int N, Args const& args = Args::global()) : BasicSiteSet<GoldenSite>(N, args) {};

    // Hamiltonian appropriate for the boundary condition, number of sites, and couplings.
    MPO Hamiltonian(const std::string& boundary_condition, int num_sites, Real U, std::vector<Real> couplings);
    std::vector<MPO> Hprojs(const std::string& boundary_condition, int num_sites);
};

#endif
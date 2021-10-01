#ifndef HAAGERUP_Q_SITE_H
#define HAAGERUP_Q_SITE_H

#include "itensor/all.h"

using namespace itensor;

// 1    1
// 2    a
// 3    b=a^2
// 4    r=\rho
// 5    ar
// 6    br
class HaagerupQSite
{
    Index s_;
public:
    explicit HaagerupQSite(Index  I);
    explicit HaagerupQSite(Args const& args = Args::global());

    Index index() const;

    IndexVal state(std::string const& state);

    ITensor id() const;

    ITensor m(int i) const;
    
    ITensor q(int i) const;

    ITensor qr(int i, int j) const;

    ITensor qar(int i, int j) const;

    ITensor qbr(int i, int j) const;

    ITensor q_all(Cplx c_3, Cplx c_5, Cplx c_6, Cplx c_7, int i, int j) const;

    ITensor op(std::string const& opname,Args const& args = Args::global()) const;
};

class HaagerupQ: public BasicSiteSet<HaagerupQSite> {
public:
    HaagerupQ() : BasicSiteSet<HaagerupQSite>() {};
    explicit HaagerupQ(int N, Args const& args = Args::global()) : BasicSiteSet<HaagerupQSite>(N, args) {};

    // Hamiltonian appropriate for the boundary condition, number of sites, and couplings.
    MPO Hamiltonian(const std::string& boundary_condition, int num_sites, Real U, std::vector<Real> couplings);
    std::vector<MPO> Hprojs(const std::string& boundary_condition, int num_sites);
};

#endif
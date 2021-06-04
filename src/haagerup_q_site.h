#ifndef HAAGERUP_Q_SITE_H
#define HAAGERUP_Q_SITE_H

#include "itensor/all.h"

using namespace itensor;

class HaagerupQSite;
using HaagerupQ=BasicSiteSet<HaagerupQSite>;
// 1    1
// 2    a
// 3    b=a^2
// 4    r=\rho
// 5    ar
// 6    br

class HaagerupQSite
{
    Index s;
public:
    explicit HaagerupQSite(Index  I);
    explicit HaagerupQSite(Args const& args = Args::global());
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
    qar(int i, int j) const;
    
    ITensor
    op(std::string const& opname,Args const& args = Args::global()) const;
};

// Hamiltonian appropriate for the boundary condition, number of sites, and couplings
// Polymorphic function specified by SiteSetType = HaagerupQ
MPO Hamiltonian(const HaagerupQ& sites, const std::string& boundary_condition, int num_sites, Real U, Real K, Real J, Real M);

#endif
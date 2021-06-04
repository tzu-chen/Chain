#include "golden_site.h"

inline int mod(int x,int N) {
    if (x>N)
        return x-N;
    return x;
}

GoldenSite::GoldenSite(Index I) : s(std::move(I)) {}

GoldenSite::GoldenSite(const Args &args) {
    auto tags = TagSet("Site,Golden");
    s = Index(2,tags);
}

ITensor GoldenSite::op(const string &opname, const Args &args) const {
    if (opname=="n1") {
        return proj(1);
    } else if (opname=="nt") {
        return proj(2);
    } else if (opname=="FF") {
        return FF();
    }
    throw ITError("Operator name "+opname+" not recognized");
}

ITensor GoldenSite::FF() const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    Op.set(s(1),sP(1),kPhiInv*kPhiInv);
    Op.set(s(1),sP(2),kSqrtPhiInv*kPhiInv);
    Op.set(s(2),sP(1),kSqrtPhiInv*kPhiInv);
    Op.set(s(2),sP(2),kPhiInv);
    return Op;
}

ITensor GoldenSite::proj(int i) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    Op.set(s(i),sP(i),1);
    return Op;
}

IndexVal GoldenSite::state(const string &state) {
    if (state=="1") {
        return s(1);
    }
//    else if (state=="r") {
//        return s(2);
//    }
    else {
        return s(2);
    }
    throw ITError("State "+state+" not recognized");
}

Index GoldenSite::index() const { return s; }

MPO Hamiltonian(const Golden& sites, const std::string& boundary_condition, int num_sites, Real U, Real K, Real J, Real M) {
    auto ampo = AutoMPO(sites);
    int L = num_sites;
    if (boundary_condition != "p") {
        L = num_sites - 1;
    }
    if (boundary_condition == "sp") {
        L = num_sites - 2;
    }
    // set up excluded pairs
    Real U_j;
    if (U != 0) {
        for(int j = 1; j <= L; ++j) {
            if (boundary_condition == "s" || boundary_condition == "sp") {
                // Uj = u_ * std::pow(sin(Pi*(j-0.5)/(num_sites_-1)),2);
                // Uj = u_ * std::pow(sin(Pi*(j+(j%2)-0.5)/(L+1)),2);
                U_j = U * std::pow(sin(Pi * (j) / (L + 1)), 2);
            } else {
                U_j = U;
            }

            ampo += U_j,"n1",j,"n1",mod(j+1, num_sites);
        }
    }

    L = num_sites;
    if (boundary_condition != "p") {
        L = num_sites - 2;
    }
    if (boundary_condition == "sp") {
        L = num_sites - 3;
    }
    // projectors
    Real K_j;
    if (K != 0) {
        for(int j = 1; j <= L; ++j) {
            if (boundary_condition == "s" || boundary_condition == "sp") {
                // Kj = k_ * std::pow(sin(Pi*(j)/(num_sites_-1)),2);
                // Kj = k_ * std::pow(sin(Pi*(j+0.25)/(num_sites_-0.5)),2);
                // Kj = k_ * std::pow(sin(Pi*(j+(j%2)-0.5+0.5)/(L+2)),2);
                K_j = K * std::pow(sin(Pi * (j + 0.5) / (L + 2)), 2);
            } else {
                K_j = K;
            }

            ampo += K_j,"n1",j,"nt",mod(j+1, num_sites),"n1",mod(j+2, num_sites);
            ampo += K_j,"nt",j,"FF",mod(j+1, num_sites),"nt",mod(j+2, num_sites);
        }
    }

    auto H = toMPO(ampo);

    return H;
}
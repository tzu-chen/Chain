#include "haagerup_site.h"

inline int mod(int x,int N) {
    if (x>N)
        return x-N;
    return x;
}

HaagerupSite::HaagerupSite(const Index &I) : s(I) {}
HaagerupSite::HaagerupSite(const Args &args) {
    auto tags=TagSet("Site,Haagerup");
    s=Index(6,tags);
}
Index HaagerupSite::index() const { return s; }

IndexVal HaagerupSite::state(const string &state) {
    if (state=="r") {
        return s(4);
    } else {
        return s(4);
    }
//    throw ITError("State "+state+" not recognized");
}

ITensor HaagerupSite::proj(int i) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    Op.set(s(i),sP(i),1);
    return Op;
}

ITensor HaagerupSite::FF(int projector, int left, int right) const {
    auto sP=prime(s);
    auto OpL=ITensor(dag(s));
    auto OpR=ITensor(sP);
    for (int i=1;i<=6;i++) {
        OpL.set(s(i),haagerup_f_data.FSymbol(left, 4, 4, right, i, projector));
        OpR.set(sP(i),haagerup_f_data.FSymbol(left, 4, 4, right, i, projector));
    }
    return OpL * OpR;
}


ITensor HaagerupSite::op(const string &opname, const Args &args) const {
    if (opname=="id") {
        return toDense(delta(s, prime(s)));
    } else if (opname=="n1") {
        return proj(1);
    } else if (opname=="na") {
        return proj(2);
    } else if (opname=="nb") {
        return proj(3);
    } else if (opname=="nr") {
        return proj(4);
    } else if (opname=="nar") {
        return proj(5);
    } else if (opname=="nbr") {
        return proj(6);
    } else if (opname=="FF1") {
        return FF(1,4,4);
    } else if (opname=="FFa") {
        return FF(1,5,5);
    } else if (opname=="FFb") {
        return FF(1,6,6);
    } else if (opname=="Frr") {
        return FF(4,4,4);
    } else if (opname=="Farar") {
        return FF(4,5,5);
    } else if (opname=="Fbrbr") {
        return FF(4,6,6);
    } else if (opname=="Frar") {
        return FF(4,4,5);
    } else if (opname=="Farbr") {
        return FF(4,5,6);
    } else if (opname=="Fbrr") {
        return FF(4,6,4);
    } else if (opname=="Frbr") {
        return FF(4,4,6);
    } else if (opname=="Farr") {
        return FF(4,5,4);
    } else if (opname=="Fbrar") {
        return FF(4,6,5);
    } else if (opname=="FFar_1") {
        return FF(5,1,4);
    } else if (opname=="FFar_a") {
        return FF(5,2,5);
    } else if (opname=="FFar_b") {
        return FF(5,3,6);
    } else if (opname=="Far_rr") {
        return FF(5,4,4);
    } else if (opname=="Far_arar") {
        return FF(5,5,5);
    } else if (opname=="Far_brbr") {
        return FF(5,6,6);
    } else if (opname=="Far_rar") {
        return FF(5,4,5);
    } else if (opname=="Far_arbr") {
        return FF(5,5,6);
    } else if (opname=="Far_brr") {
        return FF(5,6,4);
    } else if (opname=="Far_rbr") {
        return FF(5,4,6);
    } else if (opname=="Far_arr") {
        return FF(5,5,4);
    } else if (opname=="Far_brar") {
        return FF(5,6,5);
    }
    throw ITError("Operator name "+opname+" not recognized");
}

MPO Haagerup::Hamiltonian(const std::string& boundary_condition, int num_sites, Real U, std::vector<Real> couplings) {
    Real K = couplings.at(0);
    Real J = couplings.at(1);
    Real M = couplings.at(2);

    auto mpo = AutoMPO(*this);

    int L = num_sites;
    if (boundary_condition != "p") {
        L = num_sites - 1;
    }
    if (boundary_condition == "sp") {
        L = num_sites - 2;
    }
    // set up excluded pairs
    if (U!=0) {
        for (int j = 1; j <= L; ++j) {
            mpo += U,"n1",j,"n1",mod(j+1, num_sites);
            mpo += U,"n1",j,"na",mod(j+1, num_sites);
            mpo += U,"n1",j,"nb",mod(j+1, num_sites);
            mpo += U,"n1",j,"nar",mod(j+1, num_sites);
            mpo += U,"n1",j,"nbr",mod(j+1, num_sites);

            mpo += U,"na",j,"n1",mod(j+1, num_sites);
            mpo += U,"na",j,"na",mod(j+1, num_sites);
            mpo += U,"na",j,"nb",mod(j+1, num_sites);
            mpo += U,"na",j,"nr",mod(j+1, num_sites);
            mpo += U,"na",j,"nbr",mod(j+1, num_sites);

            mpo += U,"nb",j,"n1",mod(j+1, num_sites);
            mpo += U,"nb",j,"na",mod(j+1, num_sites);
            mpo += U,"nb",j,"nb",mod(j+1, num_sites);
            mpo += U,"nb",j,"nr",mod(j+1, num_sites);
            mpo += U,"nb",j,"nar",mod(j+1, num_sites);

            mpo += U,"nr",j,"na",mod(j+1, num_sites);
            mpo += U,"nr",j,"nb",mod(j+1, num_sites);

            mpo += U,"nar",j,"n1",mod(j+1, num_sites);
            mpo += U,"nar",j,"nb",mod(j+1, num_sites);

            mpo += U,"nbr",j,"n1",mod(j+1, num_sites);
            mpo += U,"nbr",j,"na",mod(j+1, num_sites);
        }
        if (boundary_condition == "d") {
            Real UL = num_sites * U;
            mpo += UL,"id",1;
            mpo += UL,"id",num_sites;
            mpo += -UL,"nr",1;
            mpo += -UL,"nr",num_sites;
        }
    }

    L = num_sites;
    if (boundary_condition != "p") {
        L = num_sites - 2;
    }
    if (boundary_condition == "sp") {
        L = num_sites - 3;
    }
    // Projectors
    // Identity
    if (K!=0) {
        for (int j = 1; j <= L; ++j) {
            mpo += K,"n1",j,"nr",mod(j+1, num_sites),"n1",mod(j+2, num_sites);
            mpo += K,"na",j,"nar",mod(j+1, num_sites),"na",mod(j+2, num_sites);
            mpo += K,"nb",j,"nbr",mod(j+1, num_sites),"nb",mod(j+2, num_sites);
            mpo += K,"nr",j,"FF1",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += K,"nar",j,"FFa",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += K,"nbr",j,"FFb",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
        }
    }

    // rho
    if (J!=0) {
        for (int j = 1; j <= L; ++j) {
            mpo += J,"n1",j,"nr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += J,"nr",j,"nr",mod(j+1, num_sites),"n1",mod(j+2, num_sites);
            mpo += J,"na",j,"nar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += J,"nar",j,"nar",mod(j+1, num_sites),"na",mod(j+2, num_sites);
            mpo += J,"nb",j,"nbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            mpo += J,"nbr",j,"nbr",mod(j+1, num_sites),"nb",mod(j+2, num_sites);

            mpo += J,"nr",j,"Frr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += J,"nr",j,"Frar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += J,"nr",j,"Frbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            mpo += J,"nar",j,"Farr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += J,"nar",j,"Farar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += J,"nar",j,"Farbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            mpo += J,"nbr",j,"Fbrr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += J,"nbr",j,"Fbrar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += J,"nbr",j,"Fbrbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
        }
    }

    // a * rho
    if (M!=0) {
        for (int j = 1; j <= L; ++j) {
            mpo += M,"n1",j,"nr",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += M,"nar",j,"nr",mod(j+1, num_sites),"n1",mod(j+2, num_sites);
            mpo += M,"na",j,"nar",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            mpo += M,"nbr",j,"nar",mod(j+1, num_sites),"na",mod(j+2, num_sites);
            mpo += M,"nb",j,"nbr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += M,"nr",j,"nbr",mod(j+1, num_sites),"nb",mod(j+2, num_sites);

            mpo += M,"nr",j,"Far_rr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += M,"nr",j,"Far_rar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += M,"nr",j,"Far_rbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            mpo += M,"nar",j,"Far_arr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += M,"nar",j,"Far_arar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += M,"nar",j,"Far_arbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            mpo += M,"nbr",j,"Far_brr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            mpo += M,"nbr",j,"Far_brar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            mpo += M,"nbr",j,"Far_brbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
        }
    }

    auto H = toMPO(mpo);

    return H;
}
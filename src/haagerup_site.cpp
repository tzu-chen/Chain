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
    throw ITError("State "+state+" not recognized");
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
//    for (int i=1;i<=6;i++) {
//        println(haagerup_f_data.FSymbol(left, 4, 4, right, i, projector));
//    }
//    println(sqrtzetainv);
//    println(x);
//    println(y1);
//    println(y2);
    for (int i=1;i<=6;i++) {
        OpL.set(s(i),haagerup_f_data.FSymbol(left, 4, 4, right, i, projector));
        OpR.set(sP(i),haagerup_f_data.FSymbol(left, 4, 4, right, i, projector));
    }
//    PrintData(OpL);
    return OpL * OpR;
}

//ITensor HaagerupSite::FF(int i) const {
//    auto sP=prime(s);
//    auto Op=ITensor(dag(s),sP);
//    Op.set(s(i),sP(i),zetainv*zetainv);
//    for(int j : {4,5,6}) {
//        Op.set(s(i),sP(j),zetainv*sqrtzetainv);
//        Op.set(s(j),sP(i),zetainv*sqrtzetainv);
//    }
//    for(int j: {4,5,6}) {
//        for(int k: {4,5,6}) {
//            Op.set(s(j),sP(k),zetainv);
//        }
//    }
//    return Op;
//}
//
//ITensor HaagerupSite::Frr(int i) const {
//    auto sP=prime(s);
//    auto OpL=ITensor(dag(s));
//    auto OpR=ITensor(sP);
//    OpL.set(s(i),sqrtzetainv);
//    OpL.set(s(3+i),x);
//    OpL.set(s(3+mod(i+1,3)),y1);
//    OpL.set(s(3+mod(i+2,3)),y2);
//    OpR.set(sP(i),sqrtzetainv);
//    OpR.set(sP(3+i),x);
//    OpR.set(sP(3+mod(i+1,3)),y1);
//    OpR.set(sP(3+mod(i+2,3)),y2);
//    return OpL*OpR;
//}
//
//ITensor HaagerupSite::Frar(int i) const {
//    auto sP=prime(s);
//    auto OpL=ITensor(dag(s));
//    auto OpR=ITensor(sP);
//    OpL.set(s(3+i),y1);
//    OpL.set(s(3+mod(i+1,3)),y2);
//    OpL.set(s(3+mod(i+2,3)),z);
//    OpR.set(sP(3+i),y1);
//    OpR.set(sP(3+mod(i+1,3)),y2);
//    OpR.set(sP(3+mod(i+2,3)),z);
//    return OpL*OpR;
//}
//
//ITensor HaagerupSite::Frbr(int i) const {
//    auto sP=prime(s);
//    auto OpL=ITensor(dag(s));
//    auto OpR=ITensor(sP);
//    OpL.set(s(3+i),y2);
//    OpL.set(s(3+mod(i+1,3)),z);
//    OpL.set(s(3+mod(i+2,3)),y1);
//    OpR.set(sP(3+i),y2);
//    OpR.set(sP(3+mod(i+1,3)),z);
//    OpR.set(sP(3+mod(i+2,3)),y1);
//    return OpL*OpR;
//}
//
//ITensor HaagerupSite::op(const string &opname, const Args &args) const {
//    if (opname=="n1") {
//        return proj(1);
//    }else if (opname=="na") {
//        return proj(2);
//    }else if (opname=="nb") {
//        return proj(3);
//    }else if (opname=="nr") {
//        return proj(4);
//    }else if (opname=="nar") {
//        return proj(5);
//    }else if (opname=="nbr") {
//        return proj(6);
//    }else if (opname=="FF1") {
//        return FF(1);
//    }else if (opname=="FFa") {
//        return FF(2);
//    }else if (opname=="FFb") {
//        return FF(3);
//    }else if (opname=="Frr") {
//        return Frr(1);
//    }else if (opname=="Farar") {
//        return Frr(2);
//    }else if (opname=="Fbrbr") {
//        return Frr(3);
//    }else if (opname=="Frar") {
//        return Frar(1);
//    }else if (opname=="Farbr") {
//        return Frar(2);
//    }else if (opname=="Fbrr") {
//        return Frar(3);
//    }else if (opname=="Frbr") {
//        return Frbr(1);
//    }else if (opname=="Farr") {
//        return Frbr(2);
//    }else if (opname=="Fbrar") {
//        return Frbr(3);
//    }
//    throw ITError("Operator name "+opname+" not recognized");
//}

ITensor HaagerupSite::op(const string &opname, const Args &args) const {
    if (opname=="n1") {
        return proj(1);
    }else if (opname=="na") {
        return proj(2);
    }else if (opname=="nb") {
        return proj(3);
    }else if (opname=="nr") {
        return proj(4);
    }else if (opname=="nar") {
        return proj(5);
    }else if (opname=="nbr") {
        return proj(6);
    }else if (opname=="FF1") {
        return FF(4,1,4);
    }else if (opname=="FFa") {
        return FF(4,2,5);
    }else if (opname=="FFb") {
        return FF(4,3,6);
    }else if (opname=="Frr") {
        return FF(4,4,4);
    }else if (opname=="Farar") {
        return FF(4,5,5);
    }else if (opname=="Fbrbr") {
        return FF(4,6,6);
    }else if (opname=="Frar") {
        return FF(4,4,5);
    }else if (opname=="Farbr") {
        return FF(4,5,6);
    }else if (opname=="Fbrr") {
        return FF(4,6,4);
    }else if (opname=="Frbr") {
        return FF(4,4,6);
    }else if (opname=="Farr") {
        return FF(4,5,4);
    }else if (opname=="Fbrar") {
        return FF(4,6,5);
    }else if (opname=="FFar_1") {
        return FF(5,1,4);
    }else if (opname=="FFar_a") {
        return FF(5,2,5);
    }else if (opname=="FFar_b") {
        return FF(5,3,6);
    }else if (opname=="Far_rr") {
        return FF(5,4,4);
    }else if (opname=="Far_arar") {
        return FF(5,5,5);
    }else if (opname=="Far_brbr") {
        return FF(5,6,6);
    }else if (opname=="Far_rar") {
        return FF(5,4,5);
    }else if (opname=="Far_arbr") {
        return FF(5,5,6);
    }else if (opname=="Far_brr") {
        return FF(5,6,4);
    }else if (opname=="Far_rbr") {
        return FF(5,4,6);
    }else if (opname=="Far_arr") {
        return FF(5,5,4);
    }else if (opname=="Far_brar") {
        return FF(5,6,5);
    }
    throw ITError("Operator name "+opname+" not recognized");
}

MPO Hamiltonian(Haagerup sites, std::string boundary_condition, int num_sites, Real U, Real K, Real J, Real M) {
    auto ampo = AutoMPO(sites);
    int L = num_sites;
    if (boundary_condition != "p") {
        L = num_sites - 1;
    }
    if (boundary_condition == "sp") {
        L = num_sites - 2;
    }
    // set up excluded pairs
    if (U!=0) {
        for(int j = 1; j <= L; ++j) {
            ampo += U,"n1",j,"n1",mod(j+1, num_sites);
            ampo += U,"n1",j,"na",mod(j+1, num_sites);
            ampo += U,"n1",j,"nb",mod(j+1, num_sites);
            ampo += U,"n1",j,"nar",mod(j+1, num_sites);
            ampo += U,"n1",j,"nbr",mod(j+1, num_sites);

            ampo += U,"na",j,"n1",mod(j+1, num_sites);
            ampo += U,"na",j,"na",mod(j+1, num_sites);
            ampo += U,"na",j,"nb",mod(j+1, num_sites);
            ampo += U,"na",j,"nr",mod(j+1, num_sites);
            ampo += U,"na",j,"nbr",mod(j+1, num_sites);

            ampo += U,"nb",j,"n1",mod(j+1, num_sites);
            ampo += U,"nb",j,"na",mod(j+1, num_sites);
            ampo += U,"nb",j,"nb",mod(j+1, num_sites);
            ampo += U,"nb",j,"nr",mod(j+1, num_sites);
            ampo += U,"nb",j,"nar",mod(j+1, num_sites);

            ampo += U,"nr",j,"na",mod(j+1, num_sites);
            ampo += U,"nr",j,"nb",mod(j+1, num_sites);

            ampo += U,"nar",j,"n1",mod(j+1, num_sites);
            ampo += U,"nar",j,"nb",mod(j+1, num_sites);

            ampo += U,"nbr",j,"n1",mod(j+1, num_sites);
            ampo += U,"nbr",j,"na",mod(j+1, num_sites);
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
        for(int j = 1; j <= L; ++j) {
            ampo += K,"n1",j,"nr",mod(j+1, num_sites),"n1",mod(j+2, num_sites);
            ampo += K,"na",j,"nar",mod(j+1, num_sites),"na",mod(j+2, num_sites);
            ampo += K,"nb",j,"nbr",mod(j+1, num_sites),"nb",mod(j+2, num_sites);
            ampo += K,"nr",j,"FF1",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += K,"nar",j,"FFa",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += K,"nbr",j,"FFb",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
        }
    }

    // rho
    if (J!=0) {
        for(int j = 1; j <= L; ++j) {
            ampo += J,"n1",j,"nr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += J,"nr",j,"nr",mod(j+1, num_sites),"n1",mod(j+2, num_sites);
            ampo += J,"na",j,"nar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += J,"nar",j,"nar",mod(j+1, num_sites),"na",mod(j+2, num_sites);
            ampo += J,"nb",j,"nbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            ampo += J,"nbr",j,"nbr",mod(j+1, num_sites),"nb",mod(j+2, num_sites);

            ampo += J,"nr",j,"Frr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += J,"nr",j,"Frar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += J,"nr",j,"Frbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            ampo += J,"nar",j,"Farr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += J,"nar",j,"Farar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += J,"nar",j,"Farbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            ampo += J,"nbr",j,"Fbrr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += J,"nbr",j,"Fbrar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += J,"nbr",j,"Fbrbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
        }
    }

    // a * rho
    if (M!=0) {
        for(int j = 1; j <= L; ++j) {
            ampo += M,"n1",j,"nr",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += M,"nar",j,"nr",mod(j+1, num_sites),"n1",mod(j+2, num_sites);
            ampo += M,"na",j,"nar",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            ampo += M,"nbr",j,"nar",mod(j+1, num_sites),"na",mod(j+2, num_sites);
            ampo += M,"nb",j,"nbr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += M,"nr",j,"nbr",mod(j+1, num_sites),"nb",mod(j+2, num_sites);

            ampo += M,"nr",j,"Far_rr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += M,"nr",j,"Far_rar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += M,"nr",j,"Far_rbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            ampo += M,"nar",j,"Far_arr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += M,"nar",j,"Far_arar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += M,"nar",j,"Far_arbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
            ampo += M,"nbr",j,"Far_brr",mod(j+1, num_sites),"nr",mod(j+2, num_sites);
            ampo += M,"nbr",j,"Far_brar",mod(j+1, num_sites),"nar",mod(j+2, num_sites);
            ampo += M,"nbr",j,"Far_brbr",mod(j+1, num_sites),"nbr",mod(j+2, num_sites);
        }
    }

    auto H = toMPO(ampo);

    return H;
}
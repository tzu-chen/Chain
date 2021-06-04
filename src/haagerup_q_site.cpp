#include "haagerup_q_site.h"

#include <utility>

inline int mod(int x,int N) {
    if (x>N)
        return x-N;
    return x;
}

HaagerupQSite::HaagerupQSite(Index I) : s(std::move(I)) {}

HaagerupQSite::HaagerupQSite(const Args &args) {
    auto ts=TagSet("Site,Haagerup");
    if (args.getBool("ConserveQNs",true))
    {
        s = Index(QN({"P",0,3}),1,
                  QN({"P",1,3}),1,
                  QN({"P",2,3}),1,
                  QN({"P",0,3}),1,
                  QN({"P",1,3}),1,
                  QN({"P",2,3}),1,Out,ts);
    }
    else
    {
        s = Index(6,ts);
    }    }

Index HaagerupQSite::index() const { return s; }

IndexVal HaagerupQSite::state(const string &state) {
//    if (state=="0") {
//        return s(4);
//    } else

    if (state=="1") {
        return s(5);
    } else if (state=="2") {
        return s(6);
    } else {
        return s(4);
    }
//    throw ITError("State "+state+" not recognized");
}

ITensor HaagerupQSite::m(int i) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    if (i==1) {
        Op.set(s(1),sP(1),1.0/3);
        Op.set(s(2),sP(2),1.0/3);
        Op.set(s(3),sP(3),1.0/3);
    } else if (i==2) {
        Op.set(s(1),sP(3),1.0/3);
        Op.set(s(2),sP(1),1.0/3);
        Op.set(s(3),sP(2),1.0/3);
    } else if (i==3) {
        Op.set(s(1),sP(2),1.0/3);
        Op.set(s(2),sP(3),1.0/3);
        Op.set(s(3),sP(1),1.0/3);
    } else if (i==4) {
        Op.set(s(4),sP(4),1.0/3);
        Op.set(s(5),sP(5),1.0/3);
        Op.set(s(6),sP(6),1.0/3);
    } else if (i==5) {
        Op.set(s(4),sP(6),1.0/3);
        Op.set(s(5),sP(4),1.0/3);
        Op.set(s(6),sP(5),1.0/3);
    } else if (i==6) {
        Op.set(s(4),sP(5),1.0/3);
        Op.set(s(5),sP(6),1.0/3);
        Op.set(s(6),sP(4),1.0/3);
    }
    return Op;
}

ITensor HaagerupQSite::q(int i) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    if (i==1) {
        Op.set(s(1),sP(1),0.030557695601338686774);
        Op.set(s(2),sP(2),0.030557695601338686774);
        Op.set(s(3),sP(3),0.030557695601338686774);
        Op.set(s(4),sP(4),0.90832691319598393968);
        Op.set(s(1),sP(4),0.16660245292295808691);
        Op.set(s(4),sP(1),0.16660245292295808691);
    } else if (i==2) {
        Op.set(s(1),sP(3),0.030557695601338686774);
        Op.set(s(2),sP(1),0.030557695601338686774);
        Op.set(s(3),sP(2),0.030557695601338686774);
        Op.set(s(2),sP(4),0.16660245292295808691);
        Op.set(s(4),sP(3),0.16660245292295808691);
    } else if (i==3) {
        Op.set(s(1),sP(2),0.030557695601338686774);
        Op.set(s(2),sP(3),0.030557695601338686774);
        Op.set(s(3),sP(1),0.030557695601338686774);
        Op.set(s(3),sP(4),0.16660245292295808691);
        Op.set(s(4),sP(2),0.16660245292295808691);
    }
    return Op;
}

ITensor HaagerupQSite::q_all(Cplx c_3, Cplx c_5, Cplx c_6, Cplx c_7, int i, int j) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    Real c_1 = 0.033641737525777182951;
    Real c_2 = -0.018511383658106454101;
    Real c_4 = 0.23240812075600178448;
    Real c_8 = -0.10092521257733154885;
    Cplx c_3_bar = conj(c_3);
    Cplx c_5_bar = conj(c_5);
    Cplx c_6_bar = conj(c_6);
    Cplx c_7_bar = conj(c_7);
    if (i==1) {
        if (j==1) {
            Op.set(s(1),sP(1),c_1);
            Op.set(s(1),sP(4),c_2);
            Op.set(s(2),sP(2),c_1);
            Op.set(s(2),sP(5),c_3);
            Op.set(s(3),sP(3),c_1);
            Op.set(s(3),sP(6),c_3_bar);

            Op.set(s(4),sP(1),c_2);
            Op.set(s(4),sP(4),c_4);
            Op.set(s(5),sP(2),c_3_bar);
            Op.set(s(5),sP(5),1.0/3);
            Op.set(s(6),sP(3),c_3);
            Op.set(s(6),sP(6),1.0/3);
        } else if (j==2) {
            Op.set(s(1),sP(3),c_1);
            Op.set(s(1),sP(6),c_3_bar);
            Op.set(s(2),sP(1),c_1);
            Op.set(s(2),sP(4),c_2);
            Op.set(s(3),sP(2),c_1);
            Op.set(s(3),sP(5),c_3);

            Op.set(s(4),sP(3),c_2);
            Op.set(s(4),sP(6),c_5);
            Op.set(s(5),sP(1),c_3_bar);
            Op.set(s(5),sP(4),c_5);
            Op.set(s(6),sP(2),c_3);
            Op.set(s(6),sP(5),c_6);
        } else if (j==3) {
            Op.set(s(1),sP(2),c_1);
            Op.set(s(1),sP(5),c_3);
            Op.set(s(2),sP(3),c_1);
            Op.set(s(2),sP(6),c_3_bar);
            Op.set(s(3),sP(1),c_1);
            Op.set(s(3),sP(4),c_2);

            Op.set(s(4),sP(2),c_2);
            Op.set(s(4),sP(5),c_5_bar);
            Op.set(s(5),sP(3),c_3_bar);
            Op.set(s(5),sP(6),c_6_bar);
            Op.set(s(6),sP(1),c_3);
            Op.set(s(6),sP(4),c_5_bar);
        }
    } else if (i==2) {
        if (j==1) {
            Op.set(s(1),sP(3),c_1);
            Op.set(s(1),sP(6),c_3_bar);
            Op.set(s(2),sP(1),c_1);
            Op.set(s(2),sP(4),c_2);
            Op.set(s(3),sP(2),c_1);
            Op.set(s(3),sP(5),c_3);

            Op.set(s(4),sP(3),c_2);
            Op.set(s(4),sP(6),c_5);
            Op.set(s(5),sP(1),c_3_bar);
            Op.set(s(5),sP(4),c_5);
            Op.set(s(6),sP(2),c_3);
            Op.set(s(6),sP(5),c_6);
        } else if (j==2) {
            Op.set(s(1),sP(2),c_1);
            Op.set(s(1),sP(5),c_3);
            Op.set(s(2),sP(3),c_1);
            Op.set(s(2),sP(6),c_3_bar);
            Op.set(s(3),sP(1),c_1);
            Op.set(s(3),sP(4),c_2);

            Op.set(s(4),sP(2),c_2);
            Op.set(s(4),sP(5),c_7);
            Op.set(s(5),sP(3),c_3_bar);
            Op.set(s(5),sP(6),c_5_bar);
            Op.set(s(6),sP(1),c_3);
            Op.set(s(6),sP(4),c_7);
        } else if (j==3) {
            Op.set(s(1),sP(1),c_1);
            Op.set(s(1),sP(4),c_2);
            Op.set(s(2),sP(2),c_1);
            Op.set(s(2),sP(5),c_3);
            Op.set(s(3),sP(3),c_1);
            Op.set(s(3),sP(6),c_3_bar);

            Op.set(s(4),sP(1),c_2);
            Op.set(s(4),sP(4),c_8);
            Op.set(s(5),sP(2),c_3_bar);
            Op.set(s(6),sP(3),c_3);
        }
    } else if (i==3) {
        if (j==1) {
            Op.set(s(1),sP(2),c_1);
            Op.set(s(1),sP(5),c_3);
            Op.set(s(2),sP(3),c_1);
            Op.set(s(2),sP(6),c_3_bar);
            Op.set(s(3),sP(1),c_1);
            Op.set(s(3),sP(4),c_2);

            Op.set(s(4),sP(2),c_2);
            Op.set(s(4),sP(5),c_5_bar);
            Op.set(s(5),sP(3),c_3_bar);
            Op.set(s(5),sP(6),c_6_bar);
            Op.set(s(6),sP(1),c_3);
            Op.set(s(6),sP(4),c_5_bar);
        } else if (j==2) {
            Op.set(s(1),sP(1),c_1);
            Op.set(s(1),sP(4),c_2);
            Op.set(s(2),sP(2),c_1);
            Op.set(s(2),sP(5),c_3);
            Op.set(s(3),sP(3),c_1);
            Op.set(s(3),sP(6),c_3_bar);

            Op.set(s(4),sP(1),c_2);
            Op.set(s(4),sP(4),c_8);
            Op.set(s(5),sP(2),c_3_bar);
            Op.set(s(6),sP(3),c_3);
        } else if (j==3) {
            Op.set(s(1),sP(3),c_1);
            Op.set(s(1),sP(6),c_3_bar);
            Op.set(s(2),sP(1),c_1);
            Op.set(s(2),sP(4),c_2);
            Op.set(s(3),sP(2),c_1);
            Op.set(s(3),sP(5),c_3);

            Op.set(s(4),sP(3),c_2);
            Op.set(s(4),sP(6),c_7_bar);
            Op.set(s(5),sP(1),c_3_bar);
            Op.set(s(5),sP(4),c_7_bar);
            Op.set(s(6),sP(2),c_3);
            Op.set(s(6),sP(5),c_5);
        }
    }
    return Op;
}

ITensor HaagerupQSite::qr(int i, int j) const {
    Cplx c_3 = -0.039825165312405310998+0.046388867673581487842 * 1_i;
    Cplx c_5 = -0.050462606288665774427+0.109830493882197205744 * 1_i;
    Cplx c_6 = -0.21966098776439441149 * 1_i;
    Cplx c_7 = 0.16666666666666666667+0.14308449170979940122 * 1_i;
    Real c_8 = -0.10092521257733154885;
    return q_all(c_3, c_5, c_6, c_7, i, j);
}

ITensor HaagerupQSite::qar(int i, int j) const {
    Cplx c_3 = -0.020261355201913645545-0.057684038707248573062 * 1_i;
    Real c_4 = 0.23240812075600178448;
    Cplx c_5 = 0.120347300956507040838-0.011213347953941672653 * 1_i;
    Cplx c_6 = -0.19023199562434830725+0.10983049388219720574 * 1_i;
    Cplx c_7 = -0.20724813804160352388+0.07279532144250674052 * 1_i;
    return q_all(c_3, c_5, c_6, c_7, i, j);
}

ITensor HaagerupQSite::qbr(int i, int j) const {
    Cplx c_3 = 0.060086520514318956543+0.011295171033667085220 * 1_i;
    Cplx c_5 = -0.069884694667841266411-0.098617145928255533092 * 1_i;
    Cplx c_6 = 0.19023199562434830725+0.10983049388219720574 * 1_i;
    Cplx c_7 = 0.04058147137493685721-0.21587981315230614174 * 1_i;
    return q_all(c_3, c_5, c_6, c_7, i, j);
}

ITensor HaagerupQSite::op(const string &opname, const Args &args) const {
    if (opname=="m1") {
        return m(1);
    } else if (opname=="m2") {
        return m(2);
    } else if (opname=="m3") {
        return m(3);
    } else if (opname=="m4") {
        return m(4);
    } else if (opname=="m5") {
        return m(5);
    } else if (opname=="m6") {
        return m(6);
    } else if (opname=="q1") {
        return q(1);
    } else if (opname=="q2") {
        return q(2);
    } else if (opname=="q3") {
        return q(3);
    } else if (opname=="qr11") {
        return qr(1,1);
    } else if (opname=="qr12") {
        return qr(1,2);
    } else if (opname=="qr13") {
        return qr(1,3);
    } else if (opname=="qr21") {
        return qr(2,1);
    } else if (opname=="qr22") {
        return qr(2,2);
    } else if (opname=="qr23") {
        return qr(2,3);
    } else if (opname=="qr31") {
        return qr(3,1);
    } else if (opname=="qr32") {
        return qr(3,2);
    } else if (opname=="qr33") {
        return qr(3,3);
    } else if (opname=="qar11") {
        return qar(1,1);
    } else if (opname=="qar12") {
        return qar(1,2);
    } else if (opname=="qar13") {
        return qar(1,3);
    } else if (opname=="qar21") {
        return qar(2,1);
    } else if (opname=="qar22") {
        return qar(2,2);
    } else if (opname=="qar23") {
        return qar(2,3);
    } else if (opname=="qar31") {
        return qar(3,1);
    } else if (opname=="qar32") {
        return qar(3,2);
    } else if (opname=="qar33") {
        return qar(3,3);
    } else if (opname=="qbr11") {
        return qbr(1,1);
    } else if (opname=="qbr12") {
        return qbr(1,2);
    } else if (opname=="qbr13") {
        return qbr(1,3);
    } else if (opname=="qbr21") {
        return qbr(2,1);
    } else if (opname=="qbr22") {
        return qbr(2,2);
    } else if (opname=="qbr23") {
        return qbr(2,3);
    } else if (opname=="qbr31") {
        return qbr(3,1);
    } else if (opname=="qbr32") {
        return qbr(3,2);
    } else if (opname=="qbr33") {
        return qbr(3,3);
    }
    throw ITError("Operator name "+opname+" not recognized");
}

MPO Hamiltonian(const HaagerupQ& sites, const std::string& boundary_condition, int num_sites, Real U, Real K, Real J, Real M) {
    auto ampo = AutoMPO(sites);

    // set up excluded pairs
    int L = num_sites;
    if (boundary_condition != "p") {
        L = num_sites - 1;
    }
    if (boundary_condition == "sp") {
        L = num_sites - 2;
    }
    Real Uj;
    if (U!=0) {
        for(int j = 1; j <= L; ++j) {
            if (boundary_condition == "s" || boundary_condition == "sp") {
                Uj = U * std::pow(sin(Pi*(j)/(L+1)),2);
            } else {
                Uj = U;
            }

            ampo += 9*Uj,"m1",j,"m1",mod(j+1, num_sites);
            ampo += 6*Uj,"m1",j,"m4",mod(j+1, num_sites);
            ampo += -3*Uj,"m2",j,"m6",mod(j+1, num_sites);
            ampo += -3*Uj,"m3",j,"m5",mod(j+1, num_sites);
            ampo += 6*Uj,"m4",j,"m1",mod(j+1, num_sites);
            ampo += -3*Uj,"m5",j,"m3",mod(j+1, num_sites);
            ampo += -3*Uj,"m6",j,"m2",mod(j+1, num_sites);
        }
    }

//    if (boundary_condition_ == "o") {
//        ampo += 6*Uj,"m1",1;
//        ampo += 6*Uj,"m1",num_sites_;
//    }

    // Projectors
    L = num_sites;
    if (boundary_condition != "p") {
        L = num_sites - 2;
    }
    if (boundary_condition == "sp") {
        L = num_sites - 3;
    }

    // Identity
    Real Kj;
    if (K!=0) {
        for(int j = 1; j <= L; ++j) {
            if (boundary_condition == "s" || boundary_condition == "sp") {
                Kj = K * std::pow(sin(Pi*(j+0.5)/(L+2)),2);
            } else {
                Kj = K;
            }

            ampo += 3*Kj,"m1",j,"m4",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*Kj,"m1",j,"m5",mod(j+1, num_sites),"m3",mod(j+2, num_sites);
            ampo += 3*Kj,"m1",j,"m6",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*Kj,"m2",j,"m4",mod(j+1, num_sites),"m3",mod(j+2, num_sites);
            ampo += 3*Kj,"m2",j,"m5",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*Kj,"m2",j,"m6",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*Kj,"m3",j,"m4",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*Kj,"m3",j,"m5",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*Kj,"m3",j,"m6",mod(j+1, num_sites),"m3",mod(j+2, num_sites);

            ampo += 3*Kj,"m4",j,"q1",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*Kj,"m4",j,"q2",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 3*Kj,"m4",j,"q3",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Kj,"m5",j,"q1",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 3*Kj,"m5",j,"q2",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Kj,"m5",j,"q3",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*Kj,"m6",j,"q1",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Kj,"m6",j,"q2",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*Kj,"m6",j,"q3",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
        }
    }

//    if (boundary_condition_ == "o") {
//        ampo += 6*Uj,"m1",1;
//        ampo += 6*Uj,"m1",num_sites_;
//    }

    // rho
    Real Jj;
    if (J!=0) {
        for(int j = 1; j <= L; ++j) {
            if (boundary_condition == "s" || boundary_condition == "sp") {
                Jj = J * std::pow(sin(Pi*(j+0.5)/(L+2)),2);
            } else {
                Jj = J;
            }

            ampo += 3*Jj,"m1",j,"m4",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*Jj,"m1",j,"m5",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 3*Jj,"m1",j,"m6",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Jj,"m2",j,"m4",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 3*Jj,"m2",j,"m5",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Jj,"m2",j,"m6",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*Jj,"m3",j,"m4",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Jj,"m3",j,"m5",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*Jj,"m3",j,"m6",mod(j+1, num_sites),"m6",mod(j+2, num_sites);

            ampo += 3*Jj,"m4",j,"m4",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*Jj,"m4",j,"m5",mod(j+1, num_sites),"m3",mod(j+2, num_sites);
            ampo += 3*Jj,"m4",j,"m6",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*Jj,"m5",j,"m4",mod(j+1, num_sites),"m3",mod(j+2, num_sites);
            ampo += 3*Jj,"m5",j,"m5",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*Jj,"m5",j,"m6",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*Jj,"m6",j,"m4",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*Jj,"m6",j,"m5",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*Jj,"m6",j,"m6",mod(j+1, num_sites),"m3",mod(j+2, num_sites);

            ampo += 9*Jj,"m4",j,"qr11",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 9*Jj,"m4",j,"qr12",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 9*Jj,"m4",j,"qr13",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 9*Jj,"m5",j,"qr31",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 9*Jj,"m5",j,"qr32",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 9*Jj,"m5",j,"qr33",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 9*Jj,"m6",j,"qr21",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 9*Jj,"m6",j,"qr22",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 9*Jj,"m6",j,"qr23",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
        }
    }

    // a * rho
    Real Mj;
    if (M!=0) {
        for(int j = 1; j <= L; ++j) {
            if (boundary_condition == "s" || boundary_condition == "sp") {
                Mj = M * std::pow(sin(Pi*(j+0.5)/(L+2)),2);
            } else {
                Mj = M;
            }

            Cplx omega = -0.5+sqrt(3)/2 * 1_i;
            Cplx omega_bar = -0.5-sqrt(3)/2 * 1_i;

            ampo += 3*Mj,"m1",j,"m4",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*omega_bar*Mj,"m1",j,"m5",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 3*omega*Mj,"m1",j,"m6",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*omega_bar*Mj,"m2",j,"m4",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 3*omega*Mj,"m2",j,"m5",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Mj,"m2",j,"m6",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*omega*Mj,"m3",j,"m4",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 3*Mj,"m3",j,"m5",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 3*omega_bar*Mj,"m3",j,"m6",mod(j+1, num_sites),"m6",mod(j+2, num_sites);

            ampo += 3*Mj,"m4",j,"m4",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*Mj,"m4",j,"m5",mod(j+1, num_sites),"m3",mod(j+2, num_sites);
            ampo += 3*Mj,"m4",j,"m6",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*omega*Mj,"m5",j,"m4",mod(j+1, num_sites),"m3",mod(j+2, num_sites);
            ampo += 3*omega*Mj,"m5",j,"m5",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*omega*Mj,"m5",j,"m6",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*omega_bar*Mj,"m6",j,"m4",mod(j+1, num_sites),"m2",mod(j+2, num_sites);
            ampo += 3*omega_bar*Mj,"m6",j,"m5",mod(j+1, num_sites),"m1",mod(j+2, num_sites);
            ampo += 3*omega_bar*Mj,"m6",j,"m6",mod(j+1, num_sites),"m3",mod(j+2, num_sites);

            ampo += 9*Mj,"m4",j,"qar11",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 9*Mj,"m4",j,"qar12",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 9*Mj,"m4",j,"qar13",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 9*Mj,"m5",j,"qar31",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 9*Mj,"m5",j,"qar32",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 9*Mj,"m5",j,"qar33",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
            ampo += 9*Mj,"m6",j,"qar21",mod(j+1, num_sites),"m4",mod(j+2, num_sites);
            ampo += 9*Mj,"m6",j,"qar22",mod(j+1, num_sites),"m6",mod(j+2, num_sites);
            ampo += 9*Mj,"m6",j,"qar23",mod(j+1, num_sites),"m5",mod(j+2, num_sites);
        }
    }

    auto H = toMPO(ampo);
    return H;
}
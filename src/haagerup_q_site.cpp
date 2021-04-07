#include "haagerup_q_site.h"

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}

HaagerupSite::HaagerupSite(const Index &I) : s(I) {}

HaagerupSite::HaagerupSite(const Args &args) {
    auto ts=TagSet("Site,Haagerup");
    if(args.getBool("ConserveQNs",true))
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

Index HaagerupSite::index() const { return s; }

IndexVal HaagerupSite::state(const string &state) {
    if(state=="0"){
        return s(4);
    }else if(state=="1"){
        return s(5);
    }else if(state=="2"){
        return s(6);
    }else{
        return s(4);
    }
    throw ITError("State "+state+" not recognized");
}

ITensor HaagerupSite::m(int i) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    if(i==1){
        Op.set(s(1),sP(1),1.0/3);
        Op.set(s(2),sP(2),1.0/3);
        Op.set(s(3),sP(3),1.0/3);
    }else if(i==2){
        Op.set(s(1),sP(3),1.0/3);
        Op.set(s(2),sP(1),1.0/3);
        Op.set(s(3),sP(2),1.0/3);
    }else if(i==3){
        Op.set(s(1),sP(2),1.0/3);
        Op.set(s(2),sP(3),1.0/3);
        Op.set(s(3),sP(1),1.0/3);
    }else if(i==4){
        Op.set(s(4),sP(4),1.0/3);
        Op.set(s(5),sP(5),1.0/3);
        Op.set(s(6),sP(6),1.0/3);
    }else if(i==5){
        Op.set(s(4),sP(6),1.0/3);
        Op.set(s(5),sP(4),1.0/3);
        Op.set(s(6),sP(5),1.0/3);
    }else if(i==6){
        Op.set(s(4),sP(5),1.0/3);
        Op.set(s(5),sP(6),1.0/3);
        Op.set(s(6),sP(4),1.0/3);
    }
    return Op;
}

ITensor HaagerupSite::q(int i) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    if(i==1){
        Op.set(s(1),sP(1),0.030557695601338686774);
        Op.set(s(2),sP(2),0.030557695601338686774);
        Op.set(s(3),sP(3),0.030557695601338686774);
        Op.set(s(4),sP(4),0.90832691319598393968);
        Op.set(s(1),sP(4),0.16660245292295808691);
        Op.set(s(4),sP(1),0.16660245292295808691);
    }else if(i==2){
        Op.set(s(1),sP(3),0.030557695601338686774);
        Op.set(s(2),sP(1),0.030557695601338686774);
        Op.set(s(3),sP(2),0.030557695601338686774);
        Op.set(s(2),sP(4),0.16660245292295808691);
        Op.set(s(4),sP(3),0.16660245292295808691);
    }else if(i==3){
        Op.set(s(1),sP(2),0.030557695601338686774);
        Op.set(s(2),sP(3),0.030557695601338686774);
        Op.set(s(3),sP(1),0.030557695601338686774);
        Op.set(s(3),sP(4),0.16660245292295808691);
        Op.set(s(4),sP(2),0.16660245292295808691);
    }
    return Op;
}

ITensor HaagerupSite::qr(int i, int j) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    if(i==1){
        if(j==1){
            Op.set(s(1),sP(1),0.033641737525777182951);
            Op.set(s(1),sP(4),-0.018511383658106454101);
            Op.set(s(2),sP(2),0.033641737525777182951);
            Op.set(s(2),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(3),sP(3),0.033641737525777182951);
            Op.set(s(3),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);

            Op.set(s(4),sP(1),-0.018511383658106454101);
            Op.set(s(4),sP(4),0.23240812075600178448);
            Op.set(s(5),sP(2),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(5),sP(5),1.0/3);
            Op.set(s(6),sP(3),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(6),1.0/3);
        }else if(j==2){
            Op.set(s(1),sP(3),0.033641737525777182951);
            Op.set(s(1),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(2),sP(1),0.033641737525777182951);
            Op.set(s(2),sP(4),-0.018511383658106454101);
            Op.set(s(3),sP(2),0.033641737525777182951);
            Op.set(s(3),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);

            Op.set(s(4),sP(3),-0.018511383658106454101);
            Op.set(s(4),sP(6),-0.050462606288665774427+0.109830493882197205744 * 1_i);
            Op.set(s(5),sP(1),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(5),sP(4),-0.050462606288665774427+0.109830493882197205744 * 1_i);
            Op.set(s(6),sP(2),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(5),-0.21966098776439441149 * 1_i);
        }else if(j==3){
            Op.set(s(1),sP(2),0.033641737525777182951);
            Op.set(s(1),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(2),sP(3),0.033641737525777182951);
            Op.set(s(2),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(3),sP(1),0.033641737525777182951);
            Op.set(s(3),sP(4),-0.018511383658106454101);

            Op.set(s(4),sP(2),-0.018511383658106454101);
            Op.set(s(4),sP(5),-0.050462606288665774427-0.109830493882197205744 * 1_i);
            Op.set(s(5),sP(3),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(5),sP(6),0.21966098776439441149 * 1_i);
            Op.set(s(6),sP(1),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(4),-0.050462606288665774427-0.109830493882197205744 * 1_i);
        }
    }else if(i==2){
        if(j==1){
            Op.set(s(1),sP(3),0.033641737525777182951);
            Op.set(s(1),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(2),sP(1),0.033641737525777182951);
            Op.set(s(2),sP(4),-0.018511383658106454101);
            Op.set(s(3),sP(2),0.033641737525777182951);
            Op.set(s(3),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);

            Op.set(s(4),sP(3),-0.018511383658106454101);
            Op.set(s(4),sP(6),-0.050462606288665774427+0.109830493882197205744 * 1_i);
            Op.set(s(5),sP(1),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(5),sP(4),-0.050462606288665774427+0.109830493882197205744 * 1_i);
            Op.set(s(6),sP(2),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(5),-0.21966098776439441149 * 1_i);
        }else if(j==2){
            Op.set(s(1),sP(2),0.033641737525777182951);
            Op.set(s(1),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(2),sP(3),0.033641737525777182951);
            Op.set(s(2),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(3),sP(1),0.033641737525777182951);
            Op.set(s(3),sP(4),-0.018511383658106454101);

            Op.set(s(4),sP(2),-0.018511383658106454101);
            Op.set(s(4),sP(5),0.16666666666666666667+0.14308449170979940122 * 1_i);
            Op.set(s(5),sP(3),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(5),sP(6),-0.050462606288665774427-0.109830493882197205744 * 1_i);
            Op.set(s(6),sP(1),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(4),0.16666666666666666667+0.14308449170979940122 * 1_i);
        }else if(j==3){
            Op.set(s(1),sP(1),0.033641737525777182951);
            Op.set(s(1),sP(4),-0.018511383658106454101);
            Op.set(s(2),sP(2),0.033641737525777182951);
            Op.set(s(2),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(3),sP(3),0.033641737525777182951);
            Op.set(s(3),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);

            Op.set(s(4),sP(1),-0.018511383658106454101);
            Op.set(s(4),sP(4),-0.10092521257733154885);
            Op.set(s(5),sP(2),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(3),-0.039825165312405310998+0.046388867673581487842 * 1_i);
        }
    }else if(i==3){
        if(j==1){
            Op.set(s(1),sP(2),0.033641737525777182951);
            Op.set(s(1),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(2),sP(3),0.033641737525777182951);
            Op.set(s(2),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(3),sP(1),0.033641737525777182951);
            Op.set(s(3),sP(4),-0.018511383658106454101);

            Op.set(s(4),sP(2),-0.018511383658106454101);
            Op.set(s(4),sP(5),-0.050462606288665774427-0.109830493882197205744 * 1_i);
            Op.set(s(5),sP(3),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(5),sP(6),0.21966098776439441149 * 1_i);
            Op.set(s(6),sP(1),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(4),-0.050462606288665774427-0.109830493882197205744 * 1_i);
        }else if(j==2){
            Op.set(s(1),sP(1),0.033641737525777182951);
            Op.set(s(1),sP(4),-0.018511383658106454101);
            Op.set(s(2),sP(2),0.033641737525777182951);
            Op.set(s(2),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(3),sP(3),0.033641737525777182951);
            Op.set(s(3),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);

            Op.set(s(4),sP(1),-0.018511383658106454101);
            Op.set(s(4),sP(4),-0.10092521257733154885);
            Op.set(s(5),sP(2),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(3),-0.039825165312405310998+0.046388867673581487842 * 1_i);
        }else if(j==3){
            Op.set(s(1),sP(3),0.033641737525777182951);
            Op.set(s(1),sP(6),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(2),sP(1),0.033641737525777182951);
            Op.set(s(2),sP(4),-0.018511383658106454101);
            Op.set(s(3),sP(2),0.033641737525777182951);
            Op.set(s(3),sP(5),-0.039825165312405310998+0.046388867673581487842 * 1_i);

            Op.set(s(4),sP(3),-0.018511383658106454101);
            Op.set(s(4),sP(6),0.16666666666666666667-0.14308449170979940122 * 1_i);
            Op.set(s(5),sP(1),-0.039825165312405310998-0.046388867673581487842 * 1_i);
            Op.set(s(5),sP(4),0.16666666666666666667-0.14308449170979940122 * 1_i);
            Op.set(s(6),sP(2),-0.039825165312405310998+0.046388867673581487842 * 1_i);
            Op.set(s(6),sP(5),-0.050462606288665774427+0.109830493882197205744 * 1_i);
        }
    }
    return Op;
}

ITensor HaagerupSite::op(const string &opname, const Args &args) const {
    if(opname=="m1"){
        return m(1);
    }else if(opname=="m2"){
        return m(2);
    }else if(opname=="m3"){
        return m(3);
    }else if(opname=="m4"){
        return m(4);
    }else if(opname=="m5"){
        return m(5);
    }else if(opname=="m6"){
        return m(6);
    }else if(opname=="q1"){
        return q(1);
    }else if(opname=="q2"){
        return q(2);
    }else if(opname=="q3"){
        return q(3);
    }else if(opname=="qr11"){
        return qr(1,1);
    }else if(opname=="qr12"){
        return qr(1,2);
    }else if(opname=="qr13"){
        return qr(1,3);
    }else if(opname=="qr21"){
        return qr(2,1);
    }else if(opname=="qr22"){
        return qr(2,2);
    }else if(opname=="qr23"){
        return qr(2,3);
    }else if(opname=="qr31"){
        return qr(3,1);
    }else if(opname=="qr32"){
        return qr(3,2);
    }else if(opname=="qr33"){
        return qr(3,3);
    }
    throw ITError("Operator name "+opname+" not recognized");
}
MPO ConstructH(Haagerup sites, std::string boundary_condition, int N, Real U, Real K, Real J) {
    auto ampo = AutoMPO(sites);

    int L = N;
    if(boundary_condition != "p"){
        L = N-1;
    }
    if(boundary_condition == "sp"){
        L = N-2;
    }
    // set up excluded pairs
    Real Uj;
    if(U!=0){
        for(int j = 1; j <= L; ++j){
            if(boundary_condition == "s" || boundary_condition == "sp"){
                Uj = U * std::pow(sin(Pi*(j)/(L+1)),2);
            }else{
                Uj = U;
            }

            ampo += 9*Uj,"m1",j,"m1",mod(j+1,N);
            ampo += 6*Uj,"m1",j,"m4",mod(j+1,N);
            ampo += -3*Uj,"m2",j,"m6",mod(j+1,N);
            ampo += -3*Uj,"m3",j,"m5",mod(j+1,N);
            ampo += 6*Uj,"m4",j,"m1",mod(j+1,N);
            ampo += -3*Uj,"m5",j,"m3",mod(j+1,N);
            ampo += -3*Uj,"m6",j,"m2",mod(j+1,N);
        }
    }

    L = N;
    if(boundary_condition != "p"){
        L = N-2;
    }
    if(boundary_condition == "sp"){
        L = N-3;
    }
    // projectors
    Real Kj;
    if(K!=0){
        for(int j = 1; j <= L; ++j){
            if(boundary_condition == "s" || boundary_condition == "sp"){
                Kj = K * std::pow(sin(Pi*(j+0.5)/(L+2)),2);
            }else{
                Kj = K;
            }

            ampo += 3*Kj,"m1",j,"m4",mod(j+1,N),"m1",mod(j+2,N);
            ampo += 3*Kj,"m1",j,"m5",mod(j+1,N),"m3",mod(j+2,N);
            ampo += 3*Kj,"m1",j,"m6",mod(j+1,N),"m2",mod(j+2,N);
            ampo += 3*Kj,"m2",j,"m4",mod(j+1,N),"m3",mod(j+2,N);
            ampo += 3*Kj,"m2",j,"m5",mod(j+1,N),"m2",mod(j+2,N);
            ampo += 3*Kj,"m2",j,"m6",mod(j+1,N),"m1",mod(j+2,N);
            ampo += 3*Kj,"m3",j,"m4",mod(j+1,N),"m2",mod(j+2,N);
            ampo += 3*Kj,"m3",j,"m5",mod(j+1,N),"m1",mod(j+2,N);
            ampo += 3*Kj,"m3",j,"m6",mod(j+1,N),"m3",mod(j+2,N);

            ampo += 3*Kj,"m4",j,"q1",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 3*Kj,"m4",j,"q2",mod(j+1,N),"m6",mod(j+2,N);
            ampo += 3*Kj,"m4",j,"q3",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 3*Kj,"m5",j,"q1",mod(j+1,N),"m6",mod(j+2,N);
            ampo += 3*Kj,"m5",j,"q2",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 3*Kj,"m5",j,"q3",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 3*Kj,"m6",j,"q1",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 3*Kj,"m6",j,"q2",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 3*Kj,"m6",j,"q3",mod(j+1,N),"m6",mod(j+2,N);
        }
    }

    Real Jj;
    if(J!=0){
        for(int j = 1; j <= L; ++j){
            if(boundary_condition == "s" || boundary_condition == "sp"){
                Jj = J * std::pow(sin(Pi*(j+0.5)/(L+2)),2);
            }else{
                Jj = J;
            }

            ampo += 3*Jj,"m1",j,"m4",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 3*Jj,"m1",j,"m5",mod(j+1,N),"m6",mod(j+2,N);
            ampo += 3*Jj,"m1",j,"m6",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 3*Jj,"m2",j,"m4",mod(j+1,N),"m6",mod(j+2,N);
            ampo += 3*Jj,"m2",j,"m5",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 3*Jj,"m2",j,"m6",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 3*Jj,"m3",j,"m4",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 3*Jj,"m3",j,"m5",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 3*Jj,"m3",j,"m6",mod(j+1,N),"m6",mod(j+2,N);

            ampo += 3*Jj,"m4",j,"m4",mod(j+1,N),"m1",mod(j+2,N);
            ampo += 3*Jj,"m4",j,"m5",mod(j+1,N),"m3",mod(j+2,N);
            ampo += 3*Jj,"m4",j,"m6",mod(j+1,N),"m2",mod(j+2,N);
            ampo += 3*Jj,"m5",j,"m4",mod(j+1,N),"m3",mod(j+2,N);
            ampo += 3*Jj,"m5",j,"m5",mod(j+1,N),"m2",mod(j+2,N);
            ampo += 3*Jj,"m5",j,"m6",mod(j+1,N),"m1",mod(j+2,N);
            ampo += 3*Jj,"m6",j,"m4",mod(j+1,N),"m2",mod(j+2,N);
            ampo += 3*Jj,"m6",j,"m5",mod(j+1,N),"m1",mod(j+2,N);
            ampo += 3*Jj,"m6",j,"m6",mod(j+1,N),"m3",mod(j+2,N);

            ampo += 9*Jj,"m4",j,"qr11",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 9*Jj,"m4",j,"qr12",mod(j+1,N),"m6",mod(j+2,N);
            ampo += 9*Jj,"m4",j,"qr13",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 9*Jj,"m5",j,"qr31",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 9*Jj,"m5",j,"qr32",mod(j+1,N),"m6",mod(j+2,N);
            ampo += 9*Jj,"m5",j,"qr33",mod(j+1,N),"m5",mod(j+2,N);
            ampo += 9*Jj,"m6",j,"qr21",mod(j+1,N),"m4",mod(j+2,N);
            ampo += 9*Jj,"m6",j,"qr22",mod(j+1,N),"m6",mod(j+2,N);
            ampo += 9*Jj,"m6",j,"qr23",mod(j+1,N),"m5",mod(j+2,N);
        }
    }

    auto H = toMPO(ampo);
    return H;
}
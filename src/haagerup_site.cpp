#include "haagerup_site.h"

inline int mod(int x,int N){
    if(x>N)
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
    if(state=="r"){
        return s(4);
    }else{
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

ITensor HaagerupSite::FF(int i) const {
    auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    Op.set(s(i),sP(i),zetainv*zetainv);
    for(int j : {4,5,6}){
        Op.set(s(i),sP(j),zetainv*sqrtzetainv);
        Op.set(s(j),sP(i),zetainv*sqrtzetainv);
    }
    for(int j: {4,5,6}){
        for(int k: {4,5,6}){
            Op.set(s(j),sP(k),zetainv);
        }
    }
    return Op;
}

ITensor HaagerupSite::Frr(int i) const {
    auto sP=prime(s);
    auto OpL=ITensor(dag(s));
    auto OpR=ITensor(sP);
    OpL.set(s(i),sqrtzetainv);
    OpL.set(s(3+i),x);
    OpL.set(s(3+mod(i+1,3)),y1);
    OpL.set(s(3+mod(i+2,3)),y2);
    OpR.set(sP(i),sqrtzetainv);
    OpR.set(sP(3+i),x);
    OpR.set(sP(3+mod(i+1,3)),y1);
    OpR.set(sP(3+mod(i+2,3)),y2);
    return OpL*OpR;
}

ITensor HaagerupSite::Frar(int i) const {
    auto sP=prime(s);
    auto OpL=ITensor(dag(s));
    auto OpR=ITensor(sP);
    OpL.set(s(3+i),y1);
    OpL.set(s(3+mod(i+1,3)),y2);
    OpL.set(s(3+mod(i+2,3)),z);
    OpR.set(sP(3+i),y1);
    OpR.set(sP(3+mod(i+1,3)),y2);
    OpR.set(sP(3+mod(i+2,3)),z);
    return OpL*OpR;
}

ITensor HaagerupSite::Frbr(int i) const {
    auto sP=prime(s);
    auto OpL=ITensor(dag(s));
    auto OpR=ITensor(sP);
    OpL.set(s(3+i),y2);
    OpL.set(s(3+mod(i+1,3)),z);
    OpL.set(s(3+mod(i+2,3)),y1);
    OpR.set(sP(3+i),y2);
    OpR.set(sP(3+mod(i+1,3)),z);
    OpR.set(sP(3+mod(i+2,3)),y1);
    return OpL*OpR;
}

ITensor HaagerupSite::op(const string &opname, const Args &args) const {
    if(opname=="n1"){
        return proj(1);
    }else if(opname=="na"){
        return proj(2);
    }else if(opname=="nb"){
        return proj(3);
    }else if(opname=="nr"){
        return proj(4);
    }else if(opname=="nar"){
        return proj(5);
    }else if(opname=="nbr"){
        return proj(6);
    }else if(opname=="FF1"){
        return FF(1);
    }else if(opname=="FFa"){
        return FF(2);
    }else if(opname=="FFb"){
        return FF(3);
    }else if(opname=="Frr"){
        return Frr(1);
    }else if(opname=="Farar"){
        return Frr(2);
    }else if(opname=="Fbrbr"){
        return Frr(3);
    }else if(opname=="Frar"){
        return Frar(1);
    }else if(opname=="Farbr"){
        return Frar(2);
    }else if(opname=="Fbrr"){
        return Frar(3);
    }else if(opname=="Frbr"){
        return Frbr(1);
    }else if(opname=="Farr"){
        return Frbr(2);
    }else if(opname=="Fbrar"){
        return Frbr(3);
    }
    throw ITError("Operator name "+opname+" not recognized");
}
// fixme: (maybe) make ConstructH a method of the site specific class?
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
    if(U!=0){
        for(int j = 1; j <= L; ++j){
            ampo += U,"n1",j,"n1",mod(j+1,N);
            ampo += U,"n1",j,"na",mod(j+1,N);
            ampo += U,"n1",j,"nb",mod(j+1,N);
            ampo += U,"n1",j,"nar",mod(j+1,N);
            ampo += U,"n1",j,"nbr",mod(j+1,N);

            ampo += U,"na",j,"n1",mod(j+1,N);
            ampo += U,"na",j,"na",mod(j+1,N);
            ampo += U,"na",j,"nb",mod(j+1,N);
            ampo += U,"na",j,"nr",mod(j+1,N);
            ampo += U,"na",j,"nbr",mod(j+1,N);

            ampo += U,"nb",j,"n1",mod(j+1,N);
            ampo += U,"nb",j,"na",mod(j+1,N);
            ampo += U,"nb",j,"nb",mod(j+1,N);
            ampo += U,"nb",j,"nr",mod(j+1,N);
            ampo += U,"nb",j,"nar",mod(j+1,N);

            ampo += U,"nr",j,"na",mod(j+1,N);
            ampo += U,"nr",j,"nb",mod(j+1,N);

            ampo += U,"nar",j,"n1",mod(j+1,N);
            ampo += U,"nar",j,"nb",mod(j+1,N);

            ampo += U,"nbr",j,"n1",mod(j+1,N);
            ampo += U,"nbr",j,"na",mod(j+1,N);
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
    if(K!=0){
        for(int j = 1; j <= L; ++j){
            ampo += K,"n1",j,"nr",mod(j+1,N),"n1",mod(j+2,N);
            ampo += K,"na",j,"nar",mod(j+1,N),"na",mod(j+2,N);
            ampo += K,"nb",j,"nbr",mod(j+1,N),"nb",mod(j+2,N);
            ampo += K,"nr",j,"FF1",mod(j+1,N),"nr",mod(j+2,N);
            ampo += K,"nar",j,"FFa",mod(j+1,N),"nar",mod(j+2,N);
            ampo += K,"nbr",j,"FFb",mod(j+1,N),"nbr",mod(j+2,N);
        }
    }

    if(J!=0){
        for(int j = 1; j <= L; ++j){
            ampo += J,"n1",j,"nr",mod(j+1,N),"nr",mod(j+2,N);
            ampo += J,"nr",j,"nr",mod(j+1,N),"n1",mod(j+2,N);
            ampo += J,"na",j,"nar",mod(j+1,N),"nar",mod(j+2,N);
            ampo += J,"nar",j,"nar",mod(j+1,N),"na",mod(j+2,N);
            ampo += J,"nb",j,"nbr",mod(j+1,N),"nbr",mod(j+2,N);
            ampo += J,"nbr",j,"nbr",mod(j+1,N),"nb",mod(j+2,N);
            ampo += J,"nr",j,"Frr",mod(j+1,N),"nr",mod(j+2,N);
            ampo += J,"nr",j,"Frar",mod(j+1,N),"nar",mod(j+2,N);
            ampo += J,"nr",j,"Frbr",mod(j+1,N),"nbr",mod(j+2,N);
            ampo += J,"nar",j,"Farr",mod(j+1,N),"nr",mod(j+2,N);
            ampo += J,"nar",j,"Farar",mod(j+1,N),"nar",mod(j+2,N);
            ampo += J,"nar",j,"Farbr",mod(j+1,N),"nbr",mod(j+2,N);
            ampo += J,"nbr",j,"Fbrr",mod(j+1,N),"nr",mod(j+2,N);
            ampo += J,"nbr",j,"Fbrar",mod(j+1,N),"nar",mod(j+2,N);
            ampo += J,"nbr",j,"Fbrbr",mod(j+1,N),"nbr",mod(j+2,N);
        }
    }

    auto H = toMPO(ampo);

    return H;
}
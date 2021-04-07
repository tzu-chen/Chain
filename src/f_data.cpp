#include "f_data.h"

bool FData::IsInvertible(int i) const {
    if(i<=kNu){
        return true;
    }else{
        return false;
    }
}

int FData::Dual(int i) const {
    if(i<=kNu){
        return 1+((1-i)%kNu);
    }else{
        return i;
    }
}

std::set<int> FData::Fusion(int a, int b) {
    std::set<int> ans;
    if(a<=kNu && b<=kNu){
        ans.insert(1+((a+b-2)%kNu));
    }else if(a<=kNu && b>kNu){
        ans.insert(kRho+((a+b-2)%kNu));
    }else if(a>kNu && b<=kNu){
        return Fusion(1+((kRho-b)%kNu),a);
    }else{
        ans.insert(1+((a-b)%kNu));
        for(int i=kRho;i<=kRk;i++){
            ans.insert(i);
        }
    }
    return ans;
}

bool FData::HasFusion(int i, int j, int k) {
    auto fused = Fusion(i,j);
    return fused.find(k) != fused.end();
}

ITensor FData::HasFusionITensor(const SiteSet &sites, int i1, int i2, int i3) {
    auto s1 = sites(i1);
    auto s2 = sites(i2);
    auto s3 = sites(i3);
    auto Op = ITensor(s1, s2, s3);
    for(int i=1;i<=kRk;i++){
        for(int j=1;j<=kRk;j++){
            for(int k=1;k<=kRk;k++){
                if(HasFusion(i,j,k)!=0){
                    Op.set( s1(i),s2(j),s3(k), HasFusion(i,j,k) );
                }
            }
        }
    }
    return Op;
}

double FData::D(int i) const {
    if(i<=kNu){
        return 1;
    }else{
        return kPhi;
    }
}

FData::FData(int nu) {
    kNu = nu;
    kRho = kNu+1;
    kRk = 2 * kNu;
    kPhi = (sqrt(4+kNu*kNu)+kNu)/2;
    kPhiInv = 1/kPhi;
    kSqrtPhiInv = sqrt(kPhiInv);
}

int FData::add(int i, int j) const {
    return kRho + (i+j-1%kNu);
}

double FData::FSymbolPattern(int i, int j, int k, int l, int m, int n) {
    if(!(HasFusion(i,j,m) && HasFusion(k,Dual(l),Dual(m)) && HasFusion(Dual(l),i,Dual(n)) && HasFusion(j,k,n))){
        return 0;
    }
    if(IsInvertible(i) || IsInvertible(j) || IsInvertible(k) || IsInvertible(l)){
        return 1;
    }
    if(IsInvertible(m) && IsInvertible(n)){
        return kPhiInv;
    }
    if(IsInvertible(m) || IsInvertible(n)){
        return kSqrtPhiInv;
    }
    if(i!=kRho){
        return FSymbol(kRho, j, add(k, i-kRho), l, m, n);
    }
    if(j!=kRho){
        return FSymbol(kRho, kRho, k, add(l, j-kRho), m, n);
    }
    if(k!=kRho){
        return FSymbol(kRho, kRho, kRho, add(l, kRho-k), m, add(n, kRho-k));
    }
    if(m!=kRho){
        return FSymbol(kRho, kRho, kRho, l, kRho, add(n, m-kRho));
    }
    // Should never get here
    // Throw exception
    return 0;
}

ITensor FData::FSymbolITensor(const SiteSet &sites, int i1, int i2, int i3) {
    auto s1 = sites(i1);
    auto s2 = sites(i2);
    auto s3 = sites(i3);
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);

    for(int i=1;i<=kRk;i++){
        for(int j=1;j<=kRk;j++){
            for(int k=1;k<=kRk;k++){
                for(int l=1;l<=kRk;l++){
                    for(int m=1;m<=kRk;m++){
                        for(int n=1;n<=kRk;n++){
                            auto f = FSymbol(i,j,k,l,m,n);
                            if(f!=0){
                                Op.set( s1(i),s1P(j),s2(k),s2P(l),s3(m),s3P(n), f);
                            }
                        }
                    }
                }
            }
        }
    }
    return Op;
}

ITensor FData::RhoDefect(const Index &s1, const Index &s2) {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P);
    for(int i1=1;i1<=kRk;i1++){
        for(int j1=1;j1<=kRk;j1++){
            for(int i2=1;i2<=kRk;i2++){
                for(int j2=1;j2<=kRk;j2++){
                    auto f = FSymbol(kRho,i1,kRho,j2,j1,i2) * sqrt(D(i1)/D(j1));
                    if(f!=0){
                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2), f);
                    }
                }
            }
        }
    }
    return Op;
}

ITensor FData::SwapITensor(const Index &s1, const Index &s2) {
    // auto s1 = sites(i1);
    // auto s2 = sites(i2);
    auto a = ITensor(dag(s1),prime(s2));
    auto b = ITensor(dag(s2),prime(s1));
    for(auto j : range1(s1))
    {
        a.set(dag(s1)(j),prime(s2)(j),1.);
        b.set(dag(s2)(j),prime(s1)(j),1.);
    }
    return a*b;
}

ITensor FData::RhoDefect2(const Index &s1, const Index &s2) {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P);
    for(int i1=1;i1<=kRk;i1++){
        for(int j1=1;j1<=kRk;j1++){
            for(int i2=1;i2<=kRk;i2++){
                for(int j2=1;j2<=kRk;j2++){
                    auto f = FSymbol(kRho,i1,kRho,j2,j1,i2) * sqrt(D(i1)/D(j1)) * FSymbol(kRho,i2,kRho,j1,j2,i1) * sqrt(D(i2)/D(j2));
                    if(f!=0){
                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2), f);
                    }
                }
            }
        }
    }
    return Op;
}

ITensor FData::RhoDefect3(const Index &s1, const Index &s2, const Index &s3) {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
    for(int i1=1;i1<=kRk;i1++){
        for(int j1=1;j1<=kRk;j1++){
            for(int i2=1;i2<=kRk;i2++){
                for(int j2=1;j2<=kRk;j2++){
                    for(int i3=1;i3<=kRk;i3++){
                        for(int j3=1;j3<=kRk;j3++){
                            auto f = FSymbol(kRho,i1,kRho,j2,j1,i2) * sqrt(D(i1)/D(j1)) * FSymbol(kRho,i2,kRho,j3,j2,i3) * sqrt(D(i2)/D(j2)) * FSymbol(kRho,i3,kRho,j1,j3,i1) * sqrt(D(i3)/D(j3));
                            if(f!=0){
                                Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2),s3(i3),s3P(j3), f);
                            }
                        }
                    }
                }
            }
        }
    }
    return Op;
}

ITensor FData::RhoDefect4(const Index &s1, const Index &s2, const Index &s3, const Index &s4) {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto s4P = prime(s4);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P, dag(s4), s4P);
    for(int i1=1;i1<=kRk;i1++){
        for(int j1=1;j1<=kRk;j1++){
            for(int i2=1;i2<=kRk;i2++){
                for(int j2=1;j2<=kRk;j2++){
                    for(int i3=1;i3<=kRk;i3++){
                        for(int j3=1;j3<=kRk;j3++){
                            for(int i4=1;i4<=kRk;i4++){
                                for(int j4=1;j4<=kRk;j4++){
                                    auto f = FSymbol(kRho,i1,kRho,j2,j1,i2) * sqrt(D(i1)/D(j1)) * FSymbol(kRho,i2,kRho,j3,j2,i3) * sqrt(D(i2)/D(j2)) * FSymbol(kRho,i3,kRho,j4,j3,i4) * sqrt(D(i3)/D(j3)) * FSymbol(kRho,i4,kRho,j1,j4,i1) * sqrt(D(i4)/D(j4));
                                    if(f!=0){
                                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2),s3(i3),s3P(j3),s4(i4),s4P(j4), f);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return Op;
}

ITensor
FData::RhoDefect6(const Index &s1, const Index &s2, const Index &s3, const Index &s4, const Index &s5, const Index &s6) {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto s4P = prime(s4);
    auto s5P = prime(s5);
    auto s6P = prime(s6);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P, dag(s4), s4P, dag(s5), s5P, dag(s6), s6P);
    for(int i1=1;i1<=kRk;i1++){
        for(int j1=1;j1<=kRk;j1++){
            for(int i2=1;i2<=kRk;i2++){
                for(int j2=1;j2<=kRk;j2++){
                    for(int i3=1;i3<=kRk;i3++){
                        for(int j3=1;j3<=kRk;j3++){
                            for(int i4=1;i4<=kRk;i4++){
                                for(int j4=1;j4<=kRk;j4++){
                                    for(int i5=1;i5<=kRk;i5++){
                                        for(int j5=1;j5<=kRk;j5++){
                                            for(int i6=1;i6<=kRk;i6++){
                                                for(int j6=1;j6<=kRk;j6++){
                                                    auto f = FSymbol(kRho,i1,kRho,j2,j1,i2) * sqrt(D(i1)/D(j1)) * FSymbol(kRho,i2,kRho,j3,j2,i3) * sqrt(D(i2)/D(j2)) * FSymbol(kRho,i3,kRho,j4,j3,i4) * sqrt(D(i3)/D(j3)) * FSymbol(kRho,i4,kRho,j5,j4,i5) * sqrt(D(i4)/D(j4)) * FSymbol(kRho,i5,kRho,j6,j5,i6) * sqrt(D(i5)/D(j5)) * FSymbol(kRho,i6,kRho,j1,j6,i1) * sqrt(D(i6)/D(j6));
                                                    if(f!=0){
                                                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2),s3(i3),s3P(j3),s4(i4),s4P(j4),s5(i5),s5P(j5),s6(i6),s6P(j6), f);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return Op;
}

GoldenFData::GoldenFData() : FData(1) {
}

double GoldenFData::FSymbol(int i, int j, int k, int l, int m, int n) {
    if(i==kRho && j==kRho && k==kRho && m==kRho && !(IsInvertible(l)) && !(IsInvertible(n))){
        return -kPhiInv;
    }
    return FSymbolPattern(i, j, k, l, m, n);
}

HaagerupFData::HaagerupFData() : FData(3) {
    x = (2-sqrt(13))/3;
    y1 = (5-sqrt(13)-sqrt(6*(1+sqrt(13))))/12;
    y2 = (5-sqrt(13)+sqrt(6*(1+sqrt(13))))/12;
    z = (1+sqrt(13))/6;
}

double HaagerupFData::FSymbol(int i, int j, int k, int l, int m, int n) {
    if(i==kRho && j==kRho && k==kRho && m==kRho && !(IsInvertible(l)) && !(IsInvertible(n))){
        switch(l+n){
            case 8 : return x;
            case 9 : return y1;
            case 10: return y2;
            case 11: return z;
            case 12: return y1;
        }
    }
    return FSymbolPattern(i, j, k, l, m, n);
}

#include "f_data.h"


// TODO: Allow fusion categories beyond Haagerup-Izumi.

FData::FData(int nu) {
    nu_ = nu;
    rho_ = nu_ + 1;
    rk_ = 2 * nu_;
    phi_ = (sqrt(4 + nu_ * nu_) + nu_) / 2;
    phi_inv_ = 1 / phi_;
    sqrt_phi_inv_ = sqrt(phi_inv_);
}

bool FData::IsInvertible(int i) const {
    if (i <= nu_) {
        return true;
    } else {
        return false;
    }
}

int FData::Dual(int i) const {
    if (i <= nu_) {
        return 1+((nu_ + 1 - i) % nu_);
    } else {
        return i;
    }
}

std::set<int> FData::Fusion(int a, int b) const {
    std::set<int> ans;
    if (a <= nu_ && b <= nu_) {
        ans.insert(1+((a+b-2) % nu_));
    } else if (a <= nu_ && b > nu_) {
        ans.insert(rho_ + ((a + b - 2) % nu_));
    } else if (a > nu_ && b <= nu_) {
        return Fusion(1+((rho_ - b) % nu_), a);
    } else {
        ans.insert(1+((nu_ + a - b) % nu_));
        for (int i=rho_; i <= rk_; i++) {
            ans.insert(i);
        }
    }
    return ans;
}

bool FData::HasFusion(int i, int j, int k) const {
    auto fused = Fusion(i,j);
    return fused.find(k) != fused.end();
}

double FData::QD(int i) const {
    if (i <= nu_) {
        return 1;
    } else {
        return phi_;
    }
}

int FData::Add(int i, int j) const {
    return rho_ + ((i + j - 1) % nu_);
}

double FData::FSymbolPattern(int i, int j, int k, int l, int m, int n) const {
    if (!(HasFusion(i,j,m) && HasFusion(k,Dual(l),Dual(m)) && HasFusion(Dual(l),i,Dual(n)) && HasFusion(j,k,n))) {
        return 0;
    }
    if (IsInvertible(i) || IsInvertible(j) || IsInvertible(k) || IsInvertible(l)) {
        return 1;
    }
    if (IsInvertible(m) && IsInvertible(n)) {
        return phi_inv_;
    }
    if (IsInvertible(m) || IsInvertible(n)) {
        return sqrt_phi_inv_;
    }
    if (i != rho_) {
        return FSymbol(rho_, j, Add(k, i - rho_), l, m, n);
    }
    if (j != rho_) {
        return FSymbol(rho_, rho_, k, Add(l, j - rho_), m, n);
    }
    if (k != rho_) {
        return FSymbol(rho_, rho_, rho_, Add(l, rho_ - k), m, Add(n, rho_ - k));
    }
    if (m != rho_) {
        return FSymbol(rho_, rho_, rho_, l, rho_, Add(n, m - rho_));
    }
    throw ITError("FSymbolPattern failed.");
    return 0;
}

ITensor FData::RhoDefectCell(const Index &s1, const Index &s2) {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P);
    for (int i1=1; i1 <= rk_; i1++) {
        for (int j1=1; j1 <= rk_; j1++) {
            for (int i2=1; i2 <= rk_; i2++) {
                for (int j2=1; j2 <= rk_; j2++) {
                    auto f = FSymbol(rho_, i1, rho_, j2, j1, i2) * sqrt(QD(i1) / QD(i2));
                    if (f!=0) {
                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2), f);
                    }
                }
            }
        }
    }
    return Op;
}

double FData::FSymbol(int i, int j, int k, int l, int m, int n) const { return 0;}

//double FData::TetrahedralSymbol(int i, int j, int k, int l, int m, int n) const {
//    return FSymbol(i, j, k, l, m, n);
////    / sqrt( QD(m) * QD(n) );
//}

void FData::DumpFMathematica(const std::string filename) {
    FILE * file;
    file = fopen(filename.c_str(),"w");
    for (int i=1; i <= rk_; i++) {
        for (int j=1; j <= rk_; j++) {
            for (int k=1; k <= rk_; k++) {
                for (int l=1; l <= rk_; l++) {
                    for (int m=1; m <= rk_; m++) {
                        for (int n=1; n <= rk_; n++) {
                            double f = FSymbol(i,j,k,l,m,n);
                            fprintf(file, "FFF[%d,%d,%d,%d][%d,%d] = %f;\n",i,j,k,l,m,n,f);
                        }
                    }
                }
            }
        }
    }
    fclose(file);
}

ITensor FData::ZipperDeltaGate(const Index &s1, const Index &s2, const Index &s3) const {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
    for (int i1=1; i1 <= rk_; i1++) {
        for (int j1=1; j1 <= rk_; j1++) {
            for (int i2=1; i2 <= rk_; i2++) {
                for (int j2=1; j2 <= rk_; j2++) {
                    for (int i3=1; i3 <= rk_; i3++) {
                        for (int j3 = 1; j3 <= rk_; j3++) {
                            if (i3==j1 && i3==j2 && i3==j3) {
                                Op.set(s1(i1), s1P(j1), s2(i2), s2P(j2), s3(i3), s3P(j3), 1);
                            }
                        }
                    }
                }
            }
        }
    }
    return Op;
}

ITensor FData::ZipperAugmentGate(const Index &s1, const Index &s2, const Index &s3) const {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
    for (int i1=1; i1 <= rk_; i1++) {
        for (int j1=1; j1 <= rk_; j1++) {
            for (int i2=1; i2 <= rk_; i2++) {
                for (int j2=1; j2 <= rk_; j2++) {
                    for (int i3=1; i3 <= rk_; i3++) {
                        for (int j3 = 1; j3 <= rk_; j3++) {
                            if (i3==j1 && i3==j3) {
                                auto f = FSymbol(rho_, rho_, i3, i3, 1, j2);
                                if (f != 0) {
                                    Op.set(s1(i1), s1P(j1), s2(i2), s2P(j2), s3(i3), s3P(j3), f);
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

ITensor FData::ZipperGate(const Index &s1, const Index &s2, const Index &s3) const {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
    for (int i1=1; i1 <= rk_; i1++) {
        for (int j1=1; j1 <= rk_; j1++) {
            for (int i2=1; i2 <= rk_; i2++) {
                for (int j2=1; j2 <= rk_; j2++) {
                    for (int i3=1; i3 <= rk_; i3++) {
                        for (int j3 = 1; j3 <= rk_; j3++) {
                            if (i1==j1 && i3==j3) {
                                auto f = FSymbol(rho_, i1, rho_, i3, i2, j2);
                                if (f != 0) {
                                    Op.set(s1(i1), s1P(j1), s2(i2), s2P(j2), s3(i3), s3P(j3), f);
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

ITensor FData::ZipperReductionGate(const Index &s1, const Index &s2, const Index &s3) const {
    auto s1P = prime(s1);
    auto s2P = prime(s2);
    auto s3P = prime(s3);
    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
    for (int i1=1; i1 <= rk_; i1++) {
        for (int j1=1; j1 <= rk_; j1++) {
            for (int i2=1; i2 <= rk_; i2++) {
                for (int j2=1; j2 <= rk_; j2++) {
                    for (int i3=1; i3 <= rk_; i3++) {
                        for (int j3 = 1; j3 <= rk_; j3++) {
                            if (i1==i3 && i1==j1 && i1==j2 && i1==j3) {
                                if (HasFusion(rho_, i1, i2) && HasFusion(rho_, i2, i1))
                                {
                                    auto f = sqrt(QD(i2) * QD(rho_) / QD(i1));
                                    Op.set(s1(i1), s1P(j1), s2(i2), s2P(j2), s3(i3), s3P(j3), f);
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

double GoldenFData::FSymbol(int i, int j, int k, int l, int m, int n) const {
    if (i == rho_ && j == rho_ && k == rho_ && m == rho_ && !(IsInvertible(l)) && !(IsInvertible(n))) {
        return -phi_inv_;
    }
    return FSymbolPattern(i, j, k, l, m, n);
}

HaagerupFData::HaagerupFData() : FData(3) {
    x_ = (2 - sqrt(13)) / 3;
    y1_ = (5 - sqrt(13) - sqrt(6 * (1 + sqrt(13)))) / 12;
    y2_ = (5 - sqrt(13) + sqrt(6 * (1 + sqrt(13)))) / 12;
    z_ = (1 + sqrt(13)) / 6;
}

double HaagerupFData::FSymbol(int i, int j, int k, int l, int m, int n) const {
    if (i == rho_ && j == rho_ && k == rho_ && m == rho_ && !(IsInvertible(l)) && !(IsInvertible(n))) {
        switch(l+n) {
            case 8 : return x_;
            case 9 : return y1_;
            case 10: return y2_;
            case 11: return z_;
            case 12: return y1_;
        }
    }
    return FSymbolPattern(i, j, k, l, m, n);
}

// DEPRECIATED
//
//ITensor FData::HasFusionITensor(const SiteSet &sites, int i1, int i2, int i3) const {
//    auto s1 = sites(i1);
//    auto s2 = sites(i2);
//    auto s3 = sites(i3);
//    auto Op = ITensor(s1, s2, s3);
//    for (int i=1; i <= kRk_; i++) {
//        for (int j=1; j <= kRk_; j++) {
//            for (int k=1; k <= kRk_; k++) {
//                if (HasFusion(i,j,k)!=0) {
//                    Op.set( s1(i),s2(j),s3(k), HasFusion(i,j,k) );
//                }
//            }
//        }
//    }
//    return Op;
//}
//
//ITensor FData::FSymbolITensor(const SiteSet &sites, int i1, int i2, int i3) {
//    auto s1 = sites(i1);
//    auto s2 = sites(i2);
//    auto s3 = sites(i3);
//    auto s1P = prime(s1);
//    auto s2P = prime(s2);
//    auto s3P = prime(s3);
//    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
//
//    for (int i=1; i <= kRk_; i++) {
//        for (int j=1; j <= kRk_; j++) {
//            for (int k=1; k <= kRk_; k++) {
//                for (int l=1; l <= kRk_; l++) {
//                    for (int m=1; m <= kRk_; m++) {
//                        for (int n=1; n <= kRk_; n++) {
//                            auto f = FSymbol(i,j,k,l,m,n);
//                            if (f!=0) {
//                                Op.set( s1(i),s1P(j),s2(k),s2P(l),s3(m),s3P(n), f);
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return Op;
//}
//
//ITensor FData::SwapITensor(const Index &s1, const Index &s2) {
//    auto a = ITensor(dag(s1),prime(s2));
//    auto b = ITensor(dag(s2),prime(s1));
//    for (auto j : range1(s1))
//    {
//        a.set(dag(s1)(j),prime(s2)(j),1.);
//        b.set(dag(s2)(j),prime(s1)(j),1.);
//    }
//    return a*b;
//}
//
//ITensor FData::RhoDefect2(const Index &s1, const Index &s2) {
//    auto s1P = prime(s1);
//    auto s2P = prime(s2);
//    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P);
//    for (int i1=1; i1 <= kRk_; i1++) {
//        for (int j1=1; j1 <= kRk_; j1++) {
//            for (int i2=1; i2 <= kRk_; i2++) {
//                for (int j2=1; j2 <= kRk_; j2++) {
//                    auto f = FSymbol(kRho_, i1, kRho_, j2, j1, i2) * sqrt(QD(i1) / QD(j1)) * FSymbol(kRho_, i2, kRho_, j1, j2, i1) *
//                             sqrt(QD(i2) /
//                                  QD(
//                                                  j2));
//                    if (f!=0) {
//                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2), f);
//                    }
//                }
//            }
//        }
//    }
//    return Op;
//}
//

//ITensor FData::RhoDefect3(const Index &s1, const Index &s2, const Index &s3){
//    auto s1P = prime(s1);
//    auto s2P = prime(s2);
//    auto s3P = prime(s3);
//    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
//    for(int i1=1; i1 <= rk_; i1++){
//        for(int j1=1; j1 <= rk_; j1++){
//            for(int i2=1; i2 <= rk_; i2++){
//                for(int j2=1; j2 <= rk_; j2++){
//                    for(int i3=1; i3 <= rk_; i3++){
//                        for(int j3=1; j3 <= rk_; j3++){
//                            // auto f = FSymbol(kRho,i1,kRho,Dual(j2),j1,i2) * sqrt(D(i1)/D(i2)) * FSymbol(kRho,i2,kRho,Dual(j3),j2,i3) * sqrt(D(i2)/D(i3)) * FSymbol(kRho,i3,kRho,Dual(j1),j3,i1) * sqrt(D(i3)/D(i1));
//                            auto f = FSymbol(rho_, i1, rho_, j2, j1, i2) * sqrt(QD(i1) / QD(i2)) * FSymbol(rho_, i2, rho_, j3, j2, i3) * sqrt(QD(i2) / QD(i3)) * FSymbol(rho_, i3, rho_, j1, j3, i1) * sqrt(QD(i3) / QD(i1));
//                            // auto f = FSymbol(kRho,i1,kRho,j2,Dual(j1),Dual(i2)) * sqrt(D(i1)/D(i2)) * FSymbol(kRho,i2,kRho,j3,Dual(j2),Dual(i3)) * sqrt(D(i2)/D(i3)) * FSymbol(kRho,i3,kRho,j1,Dual(j3),Dual(i1)) * sqrt(D(i3)/D(i1));
//                            if(f!=0){
//                                Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2),s3(i3),s3P(j3), f);
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return Op;
//}

//ITensor FData::RhoDefect3(const Index &s1, const Index &s2, const Index &s3) {
//    auto s1P = prime(s1);
//    auto s2P = prime(s2);
//    auto s3P = prime(s3);
//    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P);
//    for (int i1=1; i1 <= kRk_; i1++) {
//        for (int j1=1; j1 <= kRk_; j1++) {
//            for (int i2=1; i2 <= kRk_; i2++) {
//                for (int j2=1; j2 <= kRk_; j2++) {
//                    for (int i3=1; i3 <= kRk_; i3++) {
//                        for (int j3=1; j3 <= kRk_; j3++) {
//                            auto f = FSymbol(kRho_, i1, kRho_, j2, j1, i2) * sqrt(QD(i1) / QD(j1)) * FSymbol(kRho_, i2, kRho_, j3, j2, i3) *
//                                     sqrt(QD(i2) / QD(j2)) * FSymbol(kRho_, i3, kRho_, j1, j3, i1) *
//                                     sqrt(QD(i3) / QD(j3));
//                            if (f!=0) {
//                                Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2),s3(i3),s3P(j3), f);
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return Op;
//}

//
//ITensor FData::RhoDefect4(const Index &s1, const Index &s2, const Index &s3, const Index &s4) {
//    auto s1P = prime(s1);
//    auto s2P = prime(s2);
//    auto s3P = prime(s3);
//    auto s4P = prime(s4);
//    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P, dag(s4), s4P);
//    for (int i1=1; i1 <= kRk_; i1++) {
//        for (int j1=1; j1 <= kRk_; j1++) {
//            for (int i2=1; i2 <= kRk_; i2++) {
//                for (int j2=1; j2 <= kRk_; j2++) {
//                    for (int i3=1; i3 <= kRk_; i3++) {
//                        for (int j3=1; j3 <= kRk_; j3++) {
//                            for (int i4=1; i4 <= kRk_; i4++) {
//                                for (int j4=1; j4 <= kRk_; j4++) {
//                                    auto f = FSymbol(kRho_, i1, kRho_, j2, j1, i2) * sqrt(
//                                            QD(i1) / QD(j1)) * FSymbol(kRho_, i2, kRho_, j3, j2, i3) *
//                                             sqrt(QD(i2) /
//                                                  QD(
//                                                                  j2)) * FSymbol(kRho_, i3, kRho_, j4, j3, i4) *
//                                             sqrt(QD(i3) /
//                                                  QD(
//                                                                  j3)) * FSymbol(kRho_, i4, kRho_, j1, j4, i1) *
//                                             sqrt(QD(i4) /
//                                                  QD(
//                                                                  j4));
//                                    if (f!=0) {
//                                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2),s3(i3),s3P(j3),s4(i4),s4P(j4), f);
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return Op;
//}
//
//ITensor
//FData::RhoDefect6(const Index &s1, const Index &s2, const Index &s3, const Index &s4, const Index &s5, const Index &s6) {
//    auto s1P = prime(s1);
//    auto s2P = prime(s2);
//    auto s3P = prime(s3);
//    auto s4P = prime(s4);
//    auto s5P = prime(s5);
//    auto s6P = prime(s6);
//    auto Op = ITensor(dag(s1), s1P, dag(s2), s2P, dag(s3), s3P, dag(s4), s4P, dag(s5), s5P, dag(s6), s6P);
//    for (int i1=1; i1 <= kRk_; i1++) {
//        for (int j1=1; j1 <= kRk_; j1++) {
//            for (int i2=1; i2 <= kRk_; i2++) {
//                for (int j2=1; j2 <= kRk_; j2++) {
//                    for (int i3=1; i3 <= kRk_; i3++) {
//                        for (int j3=1; j3 <= kRk_; j3++) {
//                            for (int i4=1; i4 <= kRk_; i4++) {
//                                for (int j4=1; j4 <= kRk_; j4++) {
//                                    for (int i5=1; i5 <= kRk_; i5++) {
//                                        for (int j5=1; j5 <= kRk_; j5++) {
//                                            for (int i6=1; i6 <= kRk_; i6++) {
//                                                for (int j6=1; j6 <= kRk_; j6++) {
//                                                    auto f = FSymbol(kRho_, i1, kRho_, j2, j1, i2) * sqrt(QD(i1) /
//                                                                                                          QD(
//                                                                                                                   j1)) * FSymbol(kRho_, i2, kRho_, j3, j2, i3) *
//                                                             sqrt(QD(i2) /
//                                                                  QD(
//                                                                                  j2)) * FSymbol(kRho_, i3, kRho_, j4, j3, i4) *
//                                                             sqrt(QD(i3) /
//                                                                  QD(
//                                                                                  j3)) * FSymbol(kRho_, i4, kRho_, j5, j4, i5) *
//                                                             sqrt(QD(i4) /
//                                                                  QD(
//                                                                                  j4)) * FSymbol(kRho_, i5, kRho_, j6, j5, i6) *
//                                                             sqrt(QD(i5) /
//                                                                  QD(
//                                                                                  j5)) * FSymbol(kRho_, i6, kRho_, j1, j6, i1) *
//                                                             sqrt(QD(i6) /
//                                                                  QD(
//                                                                                  j6));
//                                                    if (f!=0) {
//                                                        Op.set(s1(i1),s1P(j1),s2(i2),s2P(j2),s3(i3),s3P(j3),s4(i4),s4P(j4),s5(i5),s5P(j5),s6(i6),s6P(j6), f);
//                                                    }
//                                                }
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return Op;
//}
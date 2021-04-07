#ifndef F_DATA_H
#define F_DATA_H


#include "itensor/all.h"
using namespace itensor;

class FData {
public:
    int kNu;
    int kRho;
    int kRk;
    double kPhi;
    double kPhiInv;
    double kSqrtPhiInv;

    explicit FData(int nu);

    // Returns true if invertible
    bool IsInvertible(int i) const;

    int Dual(int i) const;

    std::set<int> Fusion(int a, int b);

    bool HasFusion(int i, int j, int k);

    // Package F symbols into a single ITensor
    ITensor HasFusionITensor(const SiteSet& sites, int i1, int i2, int i3);

    // Quantum dimension
    double D(int i) const;

    inline int add(int i, int j) const;

    virtual double FSymbol(int i, int j, int k, int l, int m, int n) { return 0; }

    double FSymbolPattern(int i, int j, int k, int l, int m, int n);

    // Package F symbols into a single ITensor
    ITensor FSymbolITensor(const SiteSet& sites, int i1, int i2, int i3);

    ITensor RhoDefect(const Index& s1, const Index& s2);
    ITensor SwapITensor(const Index& s1, const Index& s2);
    
    ITensor RhoDefect2(const Index& s1, const Index& s2);

    ITensor RhoDefect3(const Index& s1, const Index& s2, const Index& s3);

    ITensor RhoDefect4(const Index& s1, const Index& s2, const Index& s3, const Index& s4);

    ITensor RhoDefect6(const Index& s1, const Index& s2, const Index& s3, const Index& s4, const Index& s5, const Index& s6);
};

class GoldenFData : public FData {
public:
    GoldenFData();

    // F symbol for the Fibonacci category
    double FSymbol(int i, int j, int k, int l, int m, int n) override;
};

class HaagerupFData : public FData {
public:
    double x;
    double y1;
    double y2;
    double z;

    HaagerupFData();

    // F symbol for the Fibonacci category
    double FSymbol(int i, int j, int k, int l, int m, int n) override;
};
#endif
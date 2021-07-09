#ifndef F_DATA_H
#define F_DATA_H


#include "itensor/all.h"
using namespace itensor;

// F-symbol data for the Z_nu Haagerup-Izumi categories.
class FData {
public:
    int nu_;
    int rho_;
    int rk_;
    double phi_;
    double phi_inv_;
    double sqrt_phi_inv_;

    explicit FData(int nu);

    bool IsInvertible(int i) const;

    int Dual(int i) const;

    std::set<int> Fusion(int a, int b) const;

    bool HasFusion(int i, int j, int k) const;

    // Quantum dimension.
    double QD(int i) const;

    int Add(int i, int j) const;

    // ITensor at each site in the periodic MPO form of the rho defect.
    // If T is the output of this method, then RhoOp = -T-T-T-...-T-.
    ITensor RhoDefectCell(const Index& s1, const Index& s2);

    // Relate any nontrivial F-symbol to one of the form F^{rho,rho,rho}_{*}(rho,*).
    // Valid for transparent Haagerup-Izumi.
    double FSymbolPattern(int i, int j, int k, int l, int m, int n) const;

    virtual double FSymbol(int i, int j, int k, int l, int m, int n) const;

    // Save the F-symbols to a Mathematica .m file
    void DumpFMathematica(const string filename);

    // Gates for the zipper algorithm (so far too slow)
    ITensor ZipperDeltaGate(const Index &s1, const Index &s2, const Index &s3) const;

    ITensor ZipperAugmentGate(const Index &s1, const Index &s2, const Index &s3) const;

    ITensor ZipperGate(const Index &s1, const Index &s2, const Index &s3) const;

    ITensor ZipperReductionGate(const Index &s1, const Index &s2, const Index &s3) const;

//    UNUSED
//
//    // Package fusion coefficients into a single ITensor
//    ITensor HasFusionITensor(const SiteSet& sites, int i1, int i2, int i3) const;
//
//    // Package F symbols into a single ITensor
//    ITensor FSymbolITensor(const SiteSet& sites, int i1, int i2, int i3);
//
//    ITensor SwapITensor(const Index& s1, const Index& s2);
//
//    ITensor RhoDefect2(const Index& s1, const Index& s2);
//
//    ITensor RhoDefect3(const Index& s1, const Index& s2, const Index& s3);
//
//    ITensor RhoDefect4(const Index& s1, const Index& s2, const Index& s3, const Index& s4);
//
//    ITensor RhoDefect6(const Index& s1, const Index& s2, const Index& s3, const Index& s4, const Index& s5, const Index& s6);

//    double TetrahedralSymbol(int i, int j, int k, int l, int m, int n) const;
};


// F-symbols for specific Haagerup-Izumi categories.
//
// Fibonacci
class GoldenFData : public FData {
public:
    GoldenFData();

    // F-symbol for the Fibonacci category.
    double FSymbol(int i, int j, int k, int l, int m, int n) const override;
};

// Haagerup H3
class HaagerupFData : public FData {
public:
    double x_;
    double y1_;
    double y2_;
    double z_;

    HaagerupFData();

    // F-symbol for the Haagerup H3 category.
    double FSymbol(int i, int j, int k, int l, int m, int n) const override;
};

#endif
#ifndef CHAINDMRG_CHAIN_H
#define CHAINDMRG_CHAIN_H
#include "itensor/all.h"
#include "golden_site.h"
#include "haagerup_site.h"
#include "haagerup_q_site.h"

using namespace itensor;

// MPO that translates by one site, inv controls direction.
// Slower than using Swap.
MPO TranslationOp(const SiteSet& sites, bool inv=false);

// 3-index Kronecker delta.
ITensor Delta3ITensor(const Index& s1, const Index& s2, const Index& s3);

//MPO OpenRhoOp(const SiteSet& sites, const std::string& site_type);

// Dangling form of the periodic MPO for the rho defect.
// If contracted with IdentityOp forms the periodic MPO for the rho defect.
MPO RhoOp(const SiteSet& sites, const std::string& site_type);

// Dangling form of a trivial periodic MPO, with site and dangling indices matching the input dangling MPO op.
// If contracted with the input dangling MPO op forms the periodic MPO for op.
MPO IdentityOp(const SiteSet& sites, MPO const& op);

// Recast a dangling MPS/MPO as an open MPS/MPO with two extra sites.
MPS AugmentMPS(MPS const& original_psi, Index const& sl, Index const& sr);
MPO AugmentMPO(MPO const& original_mpo, Index const& sl, Index const& sr);
//MPS ApplyIdentityOp(MPS const& original_psi, int bond_dim, Index const& sl, Index const& sr);

// Matrix for transforming Z3 charge_ basis (HaagerupQ) to standard basis Haagerup.
ITensor Z3FourierMatrix(Index const& s, Index const& sP);

// Perform the transform on a HaagerupQ MPS and return a Haagerup MPS.
MPS Z3FourierTransform(MPS const& psi, SiteSet const& sites_new);

// Swap sites b and b+1.
void Swap(MPS &psi, int b);

// Copy of the applyMPO function of ITensor using the density matrix method.
// For customization and examining bottleneck.
MPS mydensityMatrixApplyMPOImpl(MPO const& K, MPS const& psi, Args args);

// TODO: Zipper algorithm.

MPS ZipperAugmentMPS(MPS const& original_psi, Index const& sl, Index const& sll);

class ThreeSiteGate {
public:
    ThreeSiteGate(SiteSet const& sites, int i1, int i2, int i3);
    ThreeSiteGate(SiteSet const& sites, int i1, int i2, int i3, ITensor gate);
    operator const ITensor&() const { return gate_; }
    int i1() const { return i1_; }
    int i2() const { return i2_; }
    int i3() const { return i3_; }
    ITensor const& gate() const { return gate_; }
private:
    int i1_, i2_, i3_;
    ITensor gate_;
    void makeSwapGate(SiteSet const& sites);
};

ITensor inline
operator*(ThreeSiteGate const& G, ITensor T) { T *= G.gate(); return T; }

ITensor inline
operator*(ITensor T, ThreeSiteGate const& G) { T *= G.gate(); return T; }

void ActThreeSiteGate(MPS & psi, ITensor const& G, int b);

// DEPRECIATED
//
//typedef ITensor (FData::*TwoSiteGate) (const Index& s1, const Index& s2);
//void ActLocal(MPS& psi, const ITensor& G, int b);
//void ActGlobal(MPS& psi, const SiteSet& sites, TwoSiteGate gate, const std::string& sitetype);
//MPO RhoOp2(const SiteSet &sites);
//MPO NewRhoOp(const SiteSet &sites);
//void SetRho(const SiteSet &sites, MPO& m, int min, int max, bool boundary);
//MPO RhoOp_old(const SiteSet &sites);

#endif
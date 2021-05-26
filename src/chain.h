#ifndef CHAINDMRG_CHAIN_H
#define CHAINDMRG_CHAIN_H
#include "itensor/all.h"
#include "golden_site.h"
#include "haagerup_site.h"
#include "haagerup_q_site.h"
using namespace itensor;

typedef ITensor (FData::*TwoSiteGate)(const Index& s1, const Index& s2);
void Swap(MPS& psi, const SiteSet& sites, int b);
MPO TranslationOp(const SiteSet& sites, bool inv=false);
// MPO RhoOp(const SiteSet& sites, const std::string& sitetype_);
MPO RhoOp(const SiteSet &sites, const std::string& sitetype_);
MPO RhoOp2(const SiteSet &sites);
MPO NewRhoOp(const SiteSet &sites);
void SetRho(const SiteSet &sites, MPO& m, int min, int max, bool boundary);
ITensor Delta3ITensor(const Index& s1, const Index& s2, const Index& s3);
void ActLocal(MPS& psi, const ITensor& G, int b);
MPS
mydensityMatrixApplyMPOImpl(MPO const& K,
                            MPS const& psi,
                            Args args);
MPS augmentMPS(MPS const& psi, Index const& sl, Index const& sr);
MPO augmentMPO(MPO const& K, Index const& sl, Index const& sr);
MPO identity(const SiteSet & sites, MPO const& m);
MPO RhoOp_old(const SiteSet &sites);
void ActGlobal(MPS& psi, const SiteSet& sites, TwoSiteGate gate);
ITensor BasisRotation(Index const& s, Index const& sP);
MPS BasisRotate(MPS const& psi, SiteSet const& sites_new);
#endif //CHAINDMRG_CHAIN_H

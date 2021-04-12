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
MPO RhoOp(const SiteSet& sites, const std::string& sitetype_);
MPO NewRhoOp(const SiteSet &sites);
void ActLocal(MPS& psi, const ITensor& G, int b);

void ActGlobal(MPS& psi, const SiteSet& sites, TwoSiteGate gate);

#endif //CHAINDMRG_CHAIN_H

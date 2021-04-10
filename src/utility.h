#ifndef CHAINDMRG_UTILITY_H
#define CHAINDMRG_UTILITY_H

#include <filesystem>
#include "itensor/all.h"
using namespace itensor;
Real Spin(Cplx num, int NN);
Cplx Chop(Cplx num);
std::vector<double> calcEE(MPS psi, int N);
void dumpEE(int N, std::vector<double> SvNs, const std::filesystem::path& p);
void dumpEnergy(int state, Real en, const std::filesystem::path& p);

#endif //CHAINDMRG_UTILITY_H

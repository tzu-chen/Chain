#ifndef CHAINDMRG_UTILITY_H
#define CHAINDMRG_UTILITY_H

#include <filesystem>
#include "itensor/all.h"
using namespace itensor;
Real Spin(Cplx num, int NN);
Cplx Chop(Cplx num);
std::vector<double> CalcEE(MPS psi, int N);
void DumpEE(int N, std::vector<double> SvNs, const std::filesystem::path& p);
void DumpEnergy(int state, Real en, const std::filesystem::path& p);
void DumpMeasurement(const std::string& name, int len, const ITensor& i_tensor, const std::filesystem::path &p);
void DumpMathematica(int len, const ITensor& En, const ITensor& OpT, const std::filesystem::path &p);

#endif //CHAINDMRG_UTILITY_H

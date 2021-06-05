#ifndef CHAIN_DMRG_UTILITY_H
#define CHAIN_DMRG_UTILITY_H

#include <filesystem>
#include "itensor/all.h"
using namespace itensor;

// Calculate entanglement entropy curve
std::vector<double> CalcEE(MPS psi, int num_sites);

// Save entanglement entropy curve to file
void DumpEE(int num_sites, std::vector<double> SvNs, const std::filesystem::path& p);

// Save energies to file
void DumpEnergy(int state_order, Real en, const std::filesystem::path& p);

// Save a particular measurement stored as an ITensor object in Mathematica format to file
void DumpMathematicaSingle(const std::string& name, int len, const ITensor& tensor, const std::filesystem::path &p);

// Save measurements of energies, translation matrix elements, and rho defect matrix elements
// stored as ITensor objects in Mathematica format to file
void DumpMathematicaAll(int len, const ITensor& En, const ITensor& OpT, const ITensor& OpR, const std::filesystem::path& p);

// UNUSED
//Real Spin(Cplx num, int NN);
//Cplx Chop(Cplx num);

#endif
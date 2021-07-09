#ifndef CHAIN_DMRG_UTILITY_H
#define CHAIN_DMRG_UTILITY_H

#include <filesystem>
#include "itensor/all.h"
using namespace itensor;

// Calculate entanglement entropy curve.
std::vector<double> CalcEE(MPS psi, int num_sites);

// Append entanglement entropy curve to path (extension should be ee).
void DumpEE(int num_sites, std::vector<double> SvNs, const std::filesystem::path& p);

// Append energy to path (extension should be en).
void DumpEnergy(int state_order, Real en, const std::filesystem::path& p);

// Append tensor in Mathematica matrix format to path (extension should be .m).
void DumpMathematicaSingle(const std::string& name, int len, const ITensor& tensor, const std::filesystem::path &p);

// Append std::vector<std::vector<Real>> object as a Mathematica format matrix to file (extension should be .an).
void DumpMatrix(std::vector<std::vector<Real>> matrix, const std::filesystem::path &p);

// Print the real parts of std::vector<Cplx> object as a Mathematica-compatible array.
void PrintVector(std::vector<Cplx> vector);

// Print std::vector<Real> object as a Mathematica-compatible array.
void PrintVector(std::vector<Real> vector);

// Print std::vector<std::vector<Real>> object as a Mathematica-compatible matrix.
void PrintMatrix(std::vector<std::vector<Real>> matrix);

// DEPRECIATED
// Save measurements of energies, translation matrix elements, and rho defect matrix elements
// stored as ITensor objects in Mathematica format to file
//void DumpMathematicaAll(int len, const ITensor& En, const ITensor& OpT, const ITensor& OpR, const std::filesystem::path& p);
//Real Spin(Cplx num, int NN);
//Cplx Chop(Cplx num);

#endif
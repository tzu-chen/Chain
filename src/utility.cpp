#include "utility.h"


void CleanFile(const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * file;
    file = fopen(filename.c_str(), "w");
    fclose(file);
}

std::vector<double> CalcEE(MPS psi, int num_sites) {
    std::vector<double> SvNs;
    for (auto b=1; b < num_sites; b++) {
        psi.position(b);
        // SVD this wave function to get the spectrum
        // of density-matrix eigenvalues
        auto l = leftLinkIndex(psi,b);
        auto s = siteIndex(psi,b);
        auto [U,S,V] = svd(psi(b),{l,s});
        auto u = commonIndex(U,S);
        // Apply von Neumann formula
        // to the squares of the singular values
        Real SvN = 0.;
        for (auto n : range1(dim(u)))
        {
            auto Sn = elt(S,n,n);
            auto p = sqr(Sn);
            if (p > 1E-12) SvN += -p*log(p);
        }
        SvNs.push_back(SvN);
    }
    return SvNs;
}

void DumpEE(int num_sites, std::vector<double> SvNs, const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * ee_file;
    ee_file = fopen(filename.c_str(), "a");
    fprintf(ee_file, "{\n");
    for (auto b=1; b < num_sites; b++) {
        fprintf(ee_file, "{%d,  %.10f},\n", b, SvNs[b - 1]);
    }
    fprintf(ee_file, "},\n");
    fclose(ee_file);
}

void DumpEnergy(int state_order, Real en, const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * file;
    file = fopen(filename.c_str(), "a");
    fprintf(file, "{");
    fprintf(file, "%d", state_order);
    fprintf(file, ",");
    fprintf(file, "%.10f",en);
    fprintf(file, "},\n");
    fclose(file);
}

void DumpMathematicaSingle(const string& name, int len, const ITensor& tensor, const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * file;
    file = fopen(filename.c_str(),"a");
    fprintf(file, "%s = {\n", name.c_str());
    for (int i=1;i<=len;i++) {
        fprintf(file, "{");
        for (int j=1;j<=len;j++) {
            fprintf(file, "%.10f",eltC(tensor, i, j).real());
            fprintf(file, " + I * ");
            fprintf(file, "%.10f",eltC(tensor, i, j).imag());
            if (j!=len) {
                fprintf(file, ",  ");
            }
        }
        if (i!=len) {
            fprintf(file, "},\n");
        } else {
            fprintf(file, "}\n};\n\n");
        }
    }
    fclose(file);
}

// Print the real parts of std::vector<Cplx> object into a Mathematica-compatible array.
void PrintVector(std::vector<Cplx> vector) {
    print("{");
    if (vector.size() > 0) {
        print(vector.at(0).real());
        for (int i=1; i<vector.size(); i++) {
            print(",");
            print(vector.at(i).real());
        }
    }
    print("}");
}

// Print std::vector<Real> object into a Mathematica-compatible array.
void PrintVector(std::vector<Real> vector) {
    print("{");
    if (vector.size() > 0) {
        print(vector.at(0));
        for (int i = 1; i < vector.size(); i++) {
            print(",");
            print(vector.at(i));
        }
    }
    print("}");
}

// Print std::vector<std::vector<Real>> object into a Mathematica-compatible array.
void PrintMatrix(std::vector<std::vector<Real>> matrix) {
    print("{\n");
    if (matrix.size() > 0) {
        PrintVector(matrix.at(0));
        for (int i=1; i<matrix.size(); i++) {
            print(",\n");
            PrintVector(matrix.at(i));
        }
    }
    print("\n}\n");
}

// Print std::vector<std::vector<Real>> object into a Mathematica-compatible array.
void DumpMatrix(std::vector<std::vector<Real>> matrix, const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * file;
    file = fopen(filename.c_str(),"w");
    fprintf(file, "{\n");
    if (matrix.size() > 0) {
        auto vector = matrix.at(0);
        fprintf(file, "{");
        if (vector.size() > 0) {
            fprintf(file, "%.10f", vector.at(0));
            for (int i = 1; i < vector.size(); i++) {
                fprintf(file, ",");
                fprintf(file, "%.10f", vector.at(i));
            }
        }
        fprintf(file, "}");
        for (int i=1; i<matrix.size(); i++) {
            fprintf(file, ",\n");
            vector = matrix.at(i);
            fprintf(file, "{");
            if (vector.size() > 0) {
                fprintf(file, "%.10f", vector.at(0));
                for (int i = 1; i < vector.size(); i++) {
                    fprintf(file, ",");
                    fprintf(file, "%.10f", vector.at(i));
                }
            }
            fprintf(file, "}");
        }
    }
    fprintf(file, "\n}\n");
    fclose(file);
}

//void CleanMathematica(const std::filesystem::path& p) {
//    std::string filename = std::string(p);
//    // Delete original content
//    FILE * file;
//    file = fopen(filename.c_str(),"w");
//    fclose(file);
//}

//void DumpMathematicaAll(int len, const ITensor& En, const ITensor& OpT, const ITensor& OpR, const std::filesystem::path& p) {
//    std::string filename = std::string(p);
//    // Delete original content
//    FILE * file;
//    file = fopen(filename.c_str(),"w");
//    fclose(file);
//    DumpMathematicaSingle("En", len, En, filename);
//    DumpMathematicaSingle("OpT", len, OpT, filename);
//    DumpMathematicaSingle("OpR", len, OpR, filename);
//}

//// UNUSED
//Real Spin(Cplx num, int NN) {
//    Real spin = log(num).imag()/(2*Pi)*NN;
//    if (spin>NN/2) {
//        spin -= NN;
//    }
//    return spin;
//}
//
//// UNUSED
//Cplx Chop(Cplx num) {
//    Real r = num.real();
//    Real i = num.imag();
//    if (abs(r)<1E-3) {
//        r=0;
//    }
//    if (abs(i)<1E-3) {
//        i=0;
//    }
//    return r + i * 1_i;
//}
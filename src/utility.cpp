#include "utility.h"


std::vector<double> CalcEE(MPS psi, int N) {
    std::vector<double> SvNs;
    for(auto b=1;b<N;b++){
        psi.position(b);

        //SVD this wavefunction to All the spectrum
        //of density-matrix eigenvalues
        auto l = leftLinkIndex(psi,b);
        auto s = siteIndex(psi,b);
        auto [U,S,V] = svd(psi(b),{l,s});
        auto u = commonIndex(U,S);

        //Apply von Neumann formula
        //to the squares of the singular values
        Real SvN = 0.;
        for(auto n : range1(dim(u)))
        {
            auto Sn = elt(S,n,n);
            auto p = sqr(Sn);
            if(p > 1E-12) SvN += -p*log(p);
        }
        SvNs.push_back(SvN);
    }
    return SvNs;
}


void DumpEE(int N, std::vector<double> SvNs, const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * eefile;
    eefile = fopen(filename.c_str(), "a");
    fprintf(eefile, "{\n");
    for(auto b=1;b<N;b++){
        fprintf(eefile, "{%d,  %.10f},\n",b,SvNs[b-1]);
    }
    fprintf(eefile, "},\n");
    fclose(eefile);
}

void DumpEnergy(int state, Real en, const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * file;
    file = fopen(filename.c_str(), "a");
    fprintf(file, "{");
    fprintf(file, "%d",state);
    fprintf(file, ",");
    fprintf(file, "%.10f",en);
    fprintf(file, "},\n");
    fclose(file);
}

void DumpMeasurement(const string &name, int len, const ITensor &i_tensor, const std::filesystem::path &p) {
    std::string filename = std::string(p);
    FILE * file;
    file = fopen(filename.c_str(),"a");
    fprintf(file, "%s = {\n", name.c_str());
    for(int i=1;i<=len;i++){
        fprintf(file, "{");
        for(int j=1;j<=len;j++){
            fprintf(file, "%.10f",eltC(i_tensor, i, j).real());
            fprintf(file, " + I * ");
            fprintf(file, "%.10f",eltC(i_tensor, i, j).imag());
            if(j!=len){
                fprintf(file, ",  ");
            }
        }
        if(i!=len){
            fprintf(file, "},\n");
        }else{
            fprintf(file, "}\n};\n\n");
        }

    }
    fclose(file);
}

void DumpMathematica(int len, const ITensor& En, const ITensor& OpT, const ITensor& OpR, const std::filesystem::path& p) {
    std::string filename = std::string(p);
    FILE * file;
    file = fopen(filename.c_str(),"w");
    fclose(file);
    DumpMeasurement("En", len, En, filename);
    DumpMeasurement("OpT", len, OpT, filename);
    DumpMeasurement("OpR", len, OpR, filename);
}

//// UNUSED
//Real Spin(Cplx num, int NN) {
//    Real spin = log(num).imag()/(2*Pi)*NN;
//    if(spin>NN/2){
//        spin -= NN;
//    }
//    return spin;
//}
//
//// UNUSED
//Cplx Chop(Cplx num) {
//    Real r = num.real();
//    Real i = num.imag();
//    if(abs(r)<1E-3){
//        r=0;
//    }
//    if(abs(i)<1E-3){
//        i=0;
//    }
//    return r + i * 1_i;
//}
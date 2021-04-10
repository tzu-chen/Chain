#include "utility.h"

Real Spin(Cplx num, int NN) {
    Real spin = log(num).imag()/(2*Pi)*NN;
    if(spin>NN/2){
        spin -= NN;
    }
    return spin;
}



Cplx Chop(Cplx num) {
    Real r = num.real();
    Real i = num.imag();
    if(abs(r)<1E-3){
        r=0;
    }
    if(abs(i)<1E-3){
        i=0;
    }
    return r + i * 1_i;
}




std::vector<double> calcEE(MPS psi, int N) {
    std::vector<double> SvNs;
    for(auto b=1;b<N;b++){
        psi.position(b);

        //SVD this wavefunction to get the spectrum
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

// fixme: string -> filesystem::path


void dumpEE(int N, std::vector<double> SvNs, const std::filesystem::path &p) {
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

// fixme: string -> filesystem::path


void dumpEnergy(int state, Real en, const std::filesystem::path &p) {
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

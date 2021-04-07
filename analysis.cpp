//Need to change Golden to SiteSetType
#include "itensor/all.h"
#include "dmrg_progress.cpp"
#include "golden_site.h"
#include "haagerup_q_site.h"
#include <functional>
#include <eigen3/Eigen/Dense>
// using namespace Eigen;

void Swap(MPS& psi, const SiteSet& sites, int b){
    // Store original tags
    psi.position(b);
    auto tag = tags(rightLinkIndex(psi,1));

    auto G = BondGate(sites, b, b+1);
    auto wf = psi(b) * psi(b+1);
    wf *= G;
    wf.noPrime();

    auto [U,S,V] = svd(wf,inds(psi(b)));
    // auto [U,S,V] = svd(wf,inds(psi(b)),{"Cutoff=",1E-8});
    U.replaceTags(TagSet("U,Link,0"), tag);
    S.replaceTags(TagSet("U,Link,0"), tag);

    psi.set(b,U);
    psi.set(b+1,S*V);
}

MPO TranslationOp(const SiteSet& sites){
    int N = sites.length();
    auto G = std::vector<ITensor>(N);
    for(auto j : range1(N-1))
    {
        G[j] = BondGate(sites, j, j+1);
    }
    auto A = std::vector<ITensor>(N);
    auto B = std::vector<ITensor>(N);
    for(auto j : range1(N-1))
    {
        auto [Aj,Bj] = factor(G[j],{sites(j),prime(sites(j))});
        A[j] = Aj;
        B[j] = Bj;
    }
    auto t = std::vector<Index>(N+1);
    for(auto j : range1(N))
    {
        t[j] = sim(sites(j));
    }
    for(auto j : range1(2,N-1))
    {
        B[j-1] *= delta(prime(sites(j)),t[j]);
        A[j] *= delta(sites(j),t[j]);
    }
    auto Rho = MPO(N);
    Rho.set(1,A[1]);
    for(auto j : range1(2,N-1))
    {
        Rho.set(j,B[j-1]*A[j]);
    }
    Rho.set(N,B[N-1]);
    return Rho;
}

void Translate(MPS& psi, const SiteSet& sites){
    int len = psi.length();
    for(int i=1;i<len;i++){
        Swap(psi, sites, i);
    }
}

void ActLocal(MPS& psi, const ITensor& G, int b){
    // Store original tags
    psi.position(b);
    auto tag = tags(rightLinkIndex(psi,1));

    auto wf = psi(b) * psi(b+1);
    wf *= G;
    wf.noPrime();

    auto [U,S,V] = svd(wf,inds(psi(b)));
    U.replaceTags(TagSet("U,Link,0"), tag);
    S.replaceTags(TagSet("U,Link,0"), tag);

    psi.set(b,U);
    psi.set(b+1,S*V);
}
GoldenFData f_data = GoldenFData();
typedef ITensor (FData::*TwoSiteGate)(const Index& s1, const Index& s2);
void ActGlobal(MPS& psi, const SiteSet& sites, TwoSiteGate gate){
    int len = psi.length();
    for(int b=1;b<len;b++){
        ActLocal(psi, std::invoke(gate, f_data, sites(b), sites(b+1)), b);
    }
}
void ActGlobal2(MPS& psi, const SiteSet& sites, TwoSiteGate gate){
    int len = psi.length();
    for(int b=1;b<len;b++){
        ActLocal(psi, std::invoke(gate, f_data, sites(b), sites(b+1)), b);
    }
    ActLocal(psi, std::invoke(gate, f_data, sites(len), sites(1)), len);
}

int NN = 0;
Real Spin(Cplx num){
    Real spin = log(num).imag()/(2*Pi)*NN;
    if(spin>NN/2){
        spin -= NN;
    }
    return spin;
}

Cplx Chop(Cplx num){
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

// Real Chop(Real num){
//     if(abs(num)<1E-3){
//         return 0;
//     }else{
//         return num;
//     }
// }

int analyze(int argc,char** argv){
    std::string name_ = std::string(argv[1]);
    std::string boundary_condition = std::string(argv[2]);
    std::string job;
    if(boundary_condition == "p"){
        job = name_ + " PBC";
    }else if(boundary_condition == "o"){
        job = name_ + " OBC";
    }else if(boundary_condition == "s"){
        job = name_ + " SSD";
    }else{
        job = name_ + " SSD/PBC";
    }
    int N = std::stoi(argv[3]);
    NN = N;
    int gs_maxmaxdimension = std::stoi(argv[4]);
    double cutoff = std::stof(argv[5]);
    Real theta = std::stof(argv[7]);
    Real K = cos(theta*Pi);
    Real J = sin(theta*Pi);
    if(std::abs(K) < 1E-5){
        K = 0;
    }
    if(std::abs(J) < 1E-5){
        J = 0;
    }
    Real U=std::stof(argv[8]);
    int num_past_states = std::stoi(argv[9]);

    std::string progress_path = format("%s/pgs/%s_%s_%d_%g_%g_%g",prefix,name_,boundary_condition,N,cutoff,theta,U);
    std::string ee_path = format("%s/ee/%s_%s_%d_%g_%g_%g.ee",prefix,name_,boundary_condition,N,cutoff,theta,U);
    std::string en_path = format("%s/en/%s_%s_%d_%g_%g_%g.en",prefix,name_,boundary_condition,N,cutoff,theta,U);

    printf("\n%s: N=%d maxmaxdim=%d cutoff=%g theta=%g K=%g J=%g U=%g\n\n",job,N,gs_maxmaxdimension,cutoff,theta,K,J,U);

    MPS psiT;
    MPS psiR;

    DMRGProgress<Golden> dmrg_progress;
    dmrg_progress.read(progress_path);

    // auto sites = dmrg_progress.Sites();
    // auto sites = Golden(N,{"ConserveQNs=",true});
    auto states = dmrg_progress.States();
    auto sites = dmrg_progress.Sites();
    // int len = (int) states.size();
    int len = std::min((int) states.size(), num_past_states+1);

    // ITensor rd = [](Index s1, Index s2){
    //     return & f_data.RhoDefect(s1, s2);
    // };

    // auto rd = & f_data.RhoDefect;

    // Act translation on psi to create psiT
    std::vector<MPS> pastT;
    std::vector<MPS> pastR;
    for(int i=0;i<len;i++){
        psiT = MPS(states.at(i));
        psiR = MPS(states.at(i));
        ActGlobal(psiT, sites, &FData::SwapITensor);
        // Translate(psiT, sites);
        // ActGlobal(psiR, sites, &FData::RhoDefect);
        pastT.push_back(psiT);
        pastR.push_back(psiR);
    }

    // Eigen::MatrixXcf en(len,len);
    // Eigen::MatrixXcf shift(len,len);
    // for(int i=0;i<len;i++){
    //     for(int j=0;j<len;j++){
    //         if(i==j){
    //             en(i,i) = dmrg_progress.Energies().at(i);
    //         }else{
    //             en(i,j) = 0;
    //         }
    //         shift(i,j) = Chop(innerC(states.at(i), pastT.at(j)));
    //     }
    // }
    // Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces(shift);
    // PrintData((ces.eigenvectors().adjoint() * en * ces.eigenvectors()));

    // Create ITensors for
    // En:  Energies and
    // OpT:  Matrix elements between psi and the translated psiT
    auto s = Index(len);
    auto sP = prime(s);
    auto En = ITensor(dag(s),sP);
    auto OpT = ITensor(dag(s),sP);
    auto OpT2 = ITensor(dag(s),sP);
    auto OpR = ITensor(dag(s),sP);

    Real gs_en = dmrg_progress.Energies().at(0);

    for(int i=0;i<len;i++){
        En.set(s(i+1),sP(i+1),dmrg_progress.Energies().at(i) - gs_en);
        for(int j=0;j<len;j++){
            OpT.set(s(i+1),sP(j+1),Chop(innerC(states.at(i), pastT.at(j))));
            //OpR.set(s(i+1),sP(j+1),Chop(innerC(states.at(i), pastR.at(j))));
        }
    }
    PrintData(En);
    // PrintData(OpR);
    PrintData(OpT);
    println(innerC(states.at(1),applyMPO(TranslationOp(sites), states.at(1))));
    // Diagonalize translation
    auto [UT,DT] = eigen(OpT);
   
    // Rotate basis for energies accordingly and create diagonal ITensor
    En = prime(UT) * En * dag(UT);
    std::vector<Real> diag_En;
    for(int i=0;i<len;i++){
        diag_En.push_back(eltC(En, i+1, i+1).real());
    }
    En = diagITensor(diag_En, dag(s), prime(s));
    auto spins = DT.apply(Spin);
    PrintData(En);
    PrintData(spins);
    

    // auto [UEn,DEn] = eigen(En);
    // spins = prime(UEn) * spins * dag(UEn);

    // std::vector<Real> diag_spins;
    // for(int i=0;i<len;i++){
    //     diag_spins.push_back(eltC(spins, i+1, i+1).real());
    // }
    // spins = diagITensor(diag_spins, dag(s), prime(s));

    // PrintData(DEn);
    // PrintData(spins);

    
    // DT = diagITensor(diag_T, dag(s), prime(s));

    // PrintData(DEn);
    // PrintData(DT.apply(Spin));
    // PrintData(En);


    // Measure rho symmetry

    // int pos = (N+1)/2;
    // auto f_data = GoldenFData();

    // // // PrintData(f_data.FSymbolITensor(sites, 1, 2, 3));
    // // // return 0;

    // for(int state=0;state<len;state++){
    //     printf("\nState %d\n", state);

    //     auto psi = states.at(state);

    //     auto s1 = siteIndex(psi,pos);
    //     auto s2 = siteIndex(psi,pos+1);

    //     ITensor reduced_psi;

    //     // // Length 2
    //     // auto reduced_psi = psi(pos) * psi(pos+1) * delta(leftLinkIndex(psi,pos), rightLinkIndex(psi,pos+1));
    //     // PrintData(reduced_psi * prime(reduced_psi.conj()) * f_data.RhoDefect2(s1, s2));

    //     // Length 3
    //     if(N==3){
    //         auto s3 = siteIndex(psi,pos-1);
    //         auto rho_defect_3 = f_data.RhoDefect3(s3, s1, s2);
    //         for(int i=pos-1;i<=pos+1;i++){
    //             rho_defect_3 *= psi(i) * prime(psi(i)).conj();
    //         }
    //         PrintData(rho_defect_3);
    //     }

    //     // Length 4
    //     if(N==4){
    //         auto s3 = siteIndex(psi,pos-1);
    //         auto s4 = siteIndex(psi,pos+2);
    //         auto rho_defect_4 = f_data.RhoDefect4(s3, s1, s2, s4);
    //         for(int i=pos-1;i<=pos+2;i++){
    //             rho_defect_4 *= psi(i) * prime(psi(i)).conj();
    //         }
    //         PrintData(rho_defect_4);
    //     }

    //     // Length 6
    //     if(N==6){
    //         auto s3 = siteIndex(psi,pos-1);
    //         auto s4 = siteIndex(psi,pos+2);
    //         auto s5 = siteIndex(psi,pos-2);
    //         auto s6 = siteIndex(psi,pos+3);
    //         auto rho_defect_6 = f_data.RhoDefect6(s5, s3, s1, s2, s4, s6);
    //         for(int i=pos-2;i<=pos+3;i++){
    //             rho_defect_6 *= psi(i) * prime(psi(i)).conj();
    //         }
    //         PrintData(rho_defect_6);
    //     }

    //     // reduced_psi = psi(pos-1) * psi(pos) * psi(pos+1) * psi(pos+2);
    //     // auto T = reduced_psi * RhoDefect4(s3, s1, s2, s4);
    //     // T *= prime(reduced_psi);
    //     // auto [U,S,V] = svd(T,{leftLinkIndex(psi, pos-1),prime(leftLinkIndex(psi, pos-1))},{"MaxDim",5});
    //     // PrintData(S);

    //     // reduced_psi = psi(pos-2) * psi(pos-1) * psi(pos) * psi(pos+1) * psi(pos+2) * psi(pos+3);
    //     // auto T = reduced_psi * RhoDefect6(s5, s3, s1, s2, s4, s6);
    //     // T *= prime(reduced_psi.conj());
    //     // auto [U,S,V] = svd(T,{leftLinkIndex(psi, pos-2),prime(leftLinkIndex(psi, pos-2))},{"MaxDim",5});
    //     // PrintData(S);
    // }

    return 0;
}
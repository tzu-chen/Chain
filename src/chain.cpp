#include "chain.h"
#include "itensor/all.h"
#include <functional>

using namespace itensor;





void Swap(MPS &psi, const SiteSet &sites, int b) {
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



MPO TranslationOp(const SiteSet &sites, bool inv) {
    int N = sites.length();
    auto G = std::vector<ITensor>(N);
    for(auto j : range1(N-1))
    {
        G[j] = BondGate(sites, j, j + 1);
    }
    auto A = std::vector<ITensor>(N);
    auto B = std::vector<ITensor>(N);
    if (not inv) {
        for (auto j : range1(N - 1)) {
            auto[Aj, Bj] = factor(G[j], {sites(j), prime(sites(j))});
            A[j] = Aj;
            B[j] = Bj;
        }
    } else {
        for (auto j : range1(N - 1)) {
            auto[Aj, Bj] = factor(G[N-j], {sites(N-j), prime(sites(N-j))});
            A[N-j] = Aj;
            B[N-j] = Bj;
        }
    }
    auto t = std::vector<Index>(N+1);
    for(auto j : range1(N))
    {
        t[j] = sim(sites(j));
    }
    for(auto j : range1(2,N-1))
    {
        // fixme: investigate why this fails for HaagerupQ
        // update 1: some operations involving delta is not currently defined for sparse storage. see ticket #176
        A[j] *= delta(sites(j),t[j]);
        B[j-1] *= delta(prime(sites(j)),t[j]);
    }
    auto m = MPO(N);
    m.set(1, A[1]);
    for(auto j : range1(2,N-1))
    {
        m.set(j, B[j - 1] * A[j]);
    }
    m.set(N, B[N - 1]);
    return m;
}

MPO RhoOp_old(const SiteSet &sites) {
    // , const std::string& sitetype_) {
    int N = sites.length();

    auto G = std::vector<ITensor>(N+1);
    // if (sitetype_ == "golden"){
    //     auto f = GoldenFData();
    //     for (auto j : range1(N - 1)) {
    //         G[j] = f.RhoDefect(sites(j), sites(j + 1));
    //     }
    //     G[N] = GoldenFData().RhoDefect(sites(N), sites(1));
    // } else if (sitetype_ == "haagerup"){
    auto f = HaagerupFData();
    for (auto j : range1(N - 1)) {
        G[j] = f.RhoDefect(sites(j), sites(j + 1));
    }
    G[N] = f.RhoDefect(sites(N), sites(1));
    // }

    auto A = std::vector<ITensor>(N+1);
    auto B = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        auto [Aj,Bj] = factor(G[j],{sites(j),prime(sites(j))}); // Fixme: Can refactor
        // ,{"MaxDim", 10, "Cutoff", 1E-3});
        Aj.prime("Site").prime("Site");
        Bj.prime("Site").prime("Site").prime("Site").prime("Site");
        A[j] = Aj;
        B[j] = Bj;
    }
    auto long_link = commonIndex(A[N], B[N]);
    auto BNP = ITensor(B[N]);
    BNP.prime("Link");
    auto D12 = std::vector<ITensor>(N+1);
    auto D21 = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        D12[j] = Delta3ITensor(dag(sites(j)), prime(sites(j),2), prime(sites(j),4) );
        D21[j] = Delta3ITensor(prime(sites(j)), dag(prime(sites(j), 3)), dag(prime(sites(j),5)) );
    }

    auto m = MPO(N);
    m.set(1, BNP * A[1] * D12[1] * D21[1]);
    for(auto j : range1(2,N))
    {
        m.set(j, B[j-1] * A[j] * D12[j] * D21[j]);
    }

    auto l = sim(leftLinkIndex(m, 2));
    auto r = sim(rightLinkIndex(m, 2));
    auto s = sim(sites(2));
    auto dummy_i = Index(1);
    auto dummy_j = Index(1);
    auto dummy = ITensor(dummy_i, dummy_j);
    dummy.set(dummy_i(1), dummy_j(1), 1);
    auto identity_itensor = dummy * delta(l, r) * delta(s, prime(s)) * dummy;

    auto identity_mpo = MPO(N);
    identity_mpo.set(1, identity_itensor * delta(l, prime(long_link)) * delta(r, rightLinkIndex(m, 1)) * delta(s, sites(1)) * delta(prime(s), prime(sites(1))) );
    identity_mpo.set(N, identity_itensor * delta(l, leftLinkIndex(m, N)) * delta(r, long_link) * delta(s, sites(N)) * delta(prime(s), prime(sites(N))) );
    for(auto j : range1(2,N-1))
    {
        identity_mpo.set(j, identity_itensor * delta(l, leftLinkIndex(m, j)) * delta(r, rightLinkIndex(m, j)) * delta(s, sites(j)) * delta(prime(s), prime(sites(j))) );
    }
    identity_mpo = prime(identity_mpo);
    identity_mpo.set(1, identity_mpo(1) * delta(prime(long_link), prime(long_link, 2)));
    identity_mpo.set(N, identity_mpo(N) * delta(long_link, prime(long_link)));

    m = nmultMPO(m, identity_mpo, {"MaxDim", 800,"Cutoff", 1E-10});

    return m;
}

MPO RhoOp(const SiteSet &sites
     , const std::string& sitetype_) {
    int N = sites.length();

    auto G = std::vector<ITensor>(N+1);
     if (sitetype_ == "golden"){
         auto f = GoldenFData();
         for (auto j : range1(N - 1)) {
             G[j] = f.RhoDefect(sites(j), sites(j + 1));
         }
         G[N] = GoldenFData().RhoDefect(sites(N), sites(1));
     } else if (sitetype_ == "haagerup" || sitetype_ == "haagerupQN"){
        auto f = HaagerupFData();
        for (auto j : range1(N - 1)) {
            G[j] = f.RhoDefect(sites(j), sites(j + 1));
        }
        G[N] = f.RhoDefect(sites(N), sites(1));
     }

    auto A = std::vector<ITensor>(N+1);
    auto B = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        auto [Aj,Bj] = factor(G[j],{sites(j),prime(sites(j))}); // Fixme: Can refactor
        // ,{"MaxDim", 10, "Cutoff", 1E-3});
        Aj.prime("Site").prime("Site");
        Bj.prime("Site").prime("Site").prime("Site").prime("Site");
        A[j] = Aj;
        B[j] = Bj;
    }
    auto long_link = commonIndex(A[N], B[N]);
    auto BNP = ITensor(B[N]);
    BNP.prime("Link");
    auto D12 = std::vector<ITensor>(N+1);
    auto D21 = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        D12[j] = Delta3ITensor(dag(sites(j)), prime(sites(j),2), prime(sites(j),4) );
        D21[j] = Delta3ITensor(prime(sites(j)), dag(prime(sites(j), 3)), dag(prime(sites(j),5)) );
    }

    auto m = MPO(N);
    m.set(1, BNP * A[1] * D12[1] * D21[1]);
    for(auto j : range1(2,N))
    {
        m.set(j, B[j-1] * A[j] * D12[j] * D21[j]);
    }


//
//    m = nmultMPO(m, identity_mpo, {"MaxDim", 800,"Cutoff", 1E-10});

    return m;
}

MPO identity(const SiteSet & sites, MPO const& m){
    int N = sites.length();

    auto l = sim(leftLinkIndex(m, 2));
    auto r = sim(rightLinkIndex(m, 2));
    auto s = sim(sites(2));


    auto dummy_i = Index(1);
    auto dummy_j = Index(1);
    auto dummy = ITensor(dummy_i, dummy_j);
    dummy.set(dummy_i(1), dummy_j(1), 1);
    auto identity_itensor = dummy * delta(l, r) * delta(s, prime(s)) * dummy;

    auto long_link = uniqueIndex(m(1), m(2), "Link");
    auto identity_mpo = MPO(N);
    identity_mpo.set(1, identity_itensor * delta(l, prime(long_link)) * delta(r, rightLinkIndex(m, 1)) * delta(s, sites(1)) * delta(prime(s), prime(sites(1))) );
    identity_mpo.set(N, identity_itensor * delta(l, leftLinkIndex(m, N)) * delta(r, long_link) * delta(s, sites(N)) * delta(prime(s), prime(sites(N))) );
    for(auto j : range1(2,N-1))
    {
        identity_mpo.set(j, identity_itensor * delta(l, leftLinkIndex(m, j)) * delta(r, rightLinkIndex(m, j)) * delta(s, sites(j)) * delta(prime(s), prime(sites(j))) );
    }

    identity_mpo.set(1, identity_mpo(1) * delta(prime(long_link), prime(long_link, 2)));
    identity_mpo.set(N, identity_mpo(N) * delta(long_link, prime(long_link)));
    return  identity_mpo;
}

MPS
augmentMPS(MPS const& psi, Index const& sl, Index const& sr){
    auto res = MPS(length(psi)+2);
    auto extrasiteindexleft = Index(1, "l=0,Link");
    auto extral = ITensor(sl, extrasiteindexleft);
    for (auto i:range1(dim(sl))){
        extral.set(i, 1,1./dim(sl));
    }
    res.set(1, extral);
    auto vl = ITensor(extrasiteindexleft);
    vl.set(1, 1);
    res.set(2, vl*psi(1));


    for (auto i:range1(length(psi)-2)){
        res.set(i+2, psi(i+1));
    }

    auto extrasiteindexright = Index(1, "l=inf,Link");
    auto extrar = ITensor(sr, extrasiteindexright);
    for (auto i:range1(dim(sr))){
        extrar.set(i, 1,1./dim(sr));
    }
    res.set(length(psi)+2, extrar);
    auto vr = ITensor(extrasiteindexright);
    vr.set(1, 1);
    res.set(length(psi)+1, vr*psi(length(psi)));

    return res;
}

MPO augmentMPO(MPO const& K, Index const& sl, Index const& sr){
    auto res = MPO(length(K)+2);

    auto vl = ITensor(sl);
    vl.fill(1);
    auto indl = uniqueIndex(K(1), K(2), "Link");
    auto Tl = delta(prime(sl,1), indl);
    res.set(1, vl*Tl);


    for (auto i:range1(length(K)+1)){
        res.set(i+1, K(i));
    }

    auto vr = ITensor(sr);
    vr.fill(1);
    auto indr = uniqueIndex(K(length(K)), K(length(K)-1), "Link");
    auto Tr = delta(prime(sr,1), indr);
    res.set(length(K)+2, vr*Tr);

    return res;
}

// Currently unused. Convenient for modifying applyMPO func without re-compiling itensor.

MPS
mydensityMatrixApplyMPOImpl(MPO const& K,
                          MPS const& psi,
                          Args args)
{
    if( args.defined("Maxm") )
    {
        if( args.defined("MaxDim") )
        {
            Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
        else
        {
            Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
            args.add("MaxDim",args.getInt("Maxm"));
        }
    }

    auto cutoff = args.getReal("Cutoff",1E-13);
    auto dargs = Args{"Cutoff",cutoff};
    auto maxdim_set = args.defined("MaxDim");
    if(maxdim_set) dargs.add("MaxDim",args.getInt("MaxDim"));
    dargs.add("RespectDegenerate",args.getBool("RespectDegenerate",true));
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);

    auto N = length(psi);

    for( auto n : range1(N) )
    {
        if( commonIndex(psi(n),K(n)) != siteIndex(psi,n) )
            Error("MPS and MPO have different site indices in applyMPO method 'DensityMatrix'");
    }

    auto rand_plev = 14741;

    auto res = psi;

    //Set up conjugate psi and K
    auto psic = psi;
    auto Kc = K;
    //TODO: use sim(linkInds), sim(siteInds)
    psic.dag().prime(rand_plev);
    Kc.dag().prime(rand_plev);

    // Make sure the original and conjugates match

    for(auto j : range1(N-1)) {
        Kc.ref(j).prime(-rand_plev, uniqueSiteIndex(Kc, psic, j));
    }

    //Build environment tensors from the left
    if(verbose) print("Building environment tensors...");
    auto E = std::vector<ITensor>(N+1);
    E[1] = psi(1)*K(1)*Kc(1)*psic(1);
    for(int j = 2; j < N; ++j)
    {
        E[j] = E[j-1]*psi(j)*K(j)*Kc(j)*psic(j);
    }
    if(verbose) println("done");

    //O is the representation of the product of K*psi in the new MPS basis
    auto O = psi(N)*K(N);

    auto rho = E[N-1] * O * dag(prime(O,rand_plev));

    ITensor U,D;
    auto ts = tags(linkIndex(psi,N-1));
    auto spec = diagPosSemiDef(rho,U,D,{dargs,"Tags=",ts});
    if(verbose) printfln("  j=%02d truncerr=%.2E dim=%d",N-1,spec.truncerr(),dim(commonIndex(U,D)));

    res.ref(N) = dag(U);

    O = O*U*psi(N-1)*K(N-1);

    for(int j = N-1; j > 1; --j)
    {
        if(not maxdim_set)
        {
            //Infer maxdim from bond dim of original MPS
            //times bond dim of MPO
            //i.e. upper bound on order of rho
            auto cip = commonIndex(psi(j),E[j-1]);
            auto ciw = commonIndex(K(j),E[j-1]);
            auto maxdim = (cip) ? dim(cip) : 1l;
            maxdim *= (ciw) ? dim(ciw) : 1l;
            dargs.add("MaxDim",maxdim);
        }
        rho = E[j-1] * O * dag(prime(O,rand_plev));
        ts = tags(linkIndex(psi,j-1));
        auto spec = diagPosSemiDef(rho,U,D,{dargs,"Tags=",ts});
        O = O*U*psi(j-1)*K(j-1);
        res.ref(j) = dag(U);
        if(verbose) printfln("  j=%02d truncerr=%.2E dim=%d",j,spec.truncerr(),dim(commonIndex(U,D)));
    }

    if(normalize) O /= norm(O);
    res.ref(1) = O;
    res.leftLim(0);
    res.rightLim(2);

    return res;
}


// MPO RhoOp(const SiteSet &sites, const std::string& sitetype_) {
//     int N = sites.length();
//     auto G = std::vector<ITensor>(N);
//     // fixme: modify to allow for gates which are member functions of generic site type
//     if (sitetype_ == "golden"){
//         auto f = GoldenFData();
//         for (auto j : range1(N - 1)) {
//             G[j] = f.RhoDefect(sites(j), sites(j + 1));
//         }
//     } else if (sitetype_ == "haagerup"){
//         auto f = HaagerupFData();
//         for (auto j : range1(N - 1)) {
//             G[j] = f.RhoDefect(sites(j), sites(j + 1));
//         }
//     }

//     auto A = std::vector<ITensor>(N);
//     auto B = std::vector<ITensor>(N);
//     for(auto j : range1(N-1))
//     {
//         auto [Aj,Bj] = factor(G[j],{sites(j),prime(sites(j))},{"Truncate", false});
//         A[j] = Aj;
//         B[j] = Bj;
//     }
//     auto t = std::vector<Index>(N+1);
//     for(auto j : range1(N))
//     {
//         t[j] = sim(sites(j));
//     }
//     for(auto j : range1(2,N-1))
//     {
//         B[j-1] *= delta(prime(sites(j)),t[j]);
//         A[j] *= delta(sites(j),t[j]);
//     }
//     auto m = MPO(N);
//     m.set(1, A[1]);
//     for(auto j : range1(2,N-1))
//     {
//         m.set(j, B[j - 1] * A[j]);
//     }
//     m.set(N, B[N - 1]);
//     return m;
// }



MPO RhoOp2(const SiteSet &sites) {
    int N = sites.length();
    auto m = MPO(N);
    SetRho(sites, m, 1, N, true);
    return m;
}

MPO NewRhoOp(const SiteSet &sites) {
    int N = sites.length();
    auto aux_sites = Golden(2*N);
    auto G = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        G[j] = GoldenFData().RhoDefect(aux_sites(2*j-1), aux_sites(2*j));
    }
    auto A = std::vector<ITensor>(N+1);
    auto B = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        auto [Aj,Bj] = factor(G[j],{aux_sites(2*j-1),prime(aux_sites(2*j-1))});
        // Aj.prime("Site").prime("Site");
        // Bj.prime("Site").prime("Site").prime("Site").prime("Site");
        A[j] = Aj;
        B[j] = Bj;
    }
    auto m = MPO(2*N);
    for(auto j : range1(N))
    {
        m.set(2*j-1, A[j]);
        m.set(2*j, B[j]);
    }

    m = nmultMPO(prime(TranslationOp(aux_sites, false)), m);
    m.mapPrime(1,0);
    m.mapPrime(2,1);

    return m;
}




ITensor Delta3ITensor(const Index &s1, const Index &s2, const Index &s3) {
    auto a = ITensor(s1,s2,s3);
    for(auto j : range1(s1))
    {
        a.set(s1(j),s2(j),s3(j),1.);
    }
    return a;
}

void SetRho(const SiteSet &sites, MPO &m, int min, int max, bool boundary) {
    int N = max;
    auto G = std::vector<ITensor>(N+1);
    for(auto j : range1(N-1))
    {
        G[j] = GoldenFData().RhoDefect(sites(j), sites(j+1));
    }
    G[N] = GoldenFData().RhoDefect(sites(N), sites(1));

    auto A = std::vector<ITensor>(N+1);
    auto B = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        auto [Aj,Bj] = factor(G[j],{sites(j),prime(sites(j))});
        Aj.prime("Site").prime("Site");
        Bj.prime("Site").prime("Site").prime("Site").prime("Site");
        A[j] = Aj;
        B[j] = Bj;
    }
    auto D12 = std::vector<ITensor>(N+1);
    auto D21 = std::vector<ITensor>(N+1);
    for(auto j : range1(N))
    {
        D12[j] = Delta3ITensor(dag(sites(j)), prime(sites(j),2), prime(sites(j),4) );
        D21[j] = Delta3ITensor(prime(sites(j)), dag(prime(sites(j), 3)), dag(prime(sites(j),5)) );
    }
    // auto t = std::vector<Index>(N+1);
    // for(auto j : range1(N))
    // {
    //     t[j] = sim(sites(j));
    // }
    // for(auto j : range1(2,N-1))
    // {
    //     B[j-1] *= delta(prime(sites(j)),t[j]);
    //     A[j] *= delta(sites(j),t[j]);
    // }
    // auto m = MPO(N);
    if(boundary){
        // m.set(1, B[N] * A[1] * D12[1] * D21[1]);
        // m.set(N, B[N-1] * A[N] * D12[N] * D21[N]);
        m.set(min, A[min] * delta(dag(sites(min)), prime(sites(min), 2)) * delta(dag(prime(sites(min))), prime(sites(min), 3)) );
        m.set(N, B[N-1] * delta(dag(sites(N)), prime(sites(N), 4)) * delta(dag(prime(sites(N))), prime(sites(N), 5)) );
    }else{
        auto ind_min_old = commonIndex(B[min], A[min]);
        auto ind_max_old = commonIndex(B[max-1], A[max-1]);
        auto ind_min_new = rightLinkIndex(m, min);
        auto ind_max_new = leftLinkIndex(m, max);
        B[min] *= delta(ind_min_old, ind_min_new);
        A[max-1] *= delta(ind_max_old, ind_max_new);
    }
    for(auto j : range1(min+1,N-1))
    {
        // print(j);
        // print("\n");
        m.set(j, B[j-1] * A[j] * D12[j] * D21[j]);
    }
}

void ActLocal(MPS &psi, const ITensor &G, int b) {
    // Store original tags
    psi.position(b);
    auto tag = tags(rightLinkIndex(psi,1));

    auto wf = psi(b) * psi(b+1);
    wf *= G;
    wf.noPrime();

    auto [U,S,V] = svd(wf,inds(psi(b)),{"Truncate", false});
    U.replaceTags(TagSet("U,Link,0"), tag);
    S.replaceTags(TagSet("U,Link,0"), tag);

    psi.set(b,U);
    psi.set(b+1,S*V);
}


// fixme: modify to allow for gates which are member functions of generic site type
void ActGlobal(MPS &psi, const SiteSet &sites, TwoSiteGate gate, const std::string& sitetype_) {
    int len = psi.length();
    if (sitetype_ == "golden"){
        auto f = GoldenFData();
        for (int b = 1; b < len; b++) {
            ActLocal(psi, std::invoke(gate, f, sites(b), sites(b + 1)), b);
        }
    } else if (sitetype_ == "haagerup") {
        auto f = HaagerupFData();
        for (int b = 1; b < len; b++) {
            ActLocal(psi, std::invoke(gate, f, sites(b), sites(b + 1)), b);
        }
    }
}
ITensor BasisRotation(Index const& s, Index const& sP) {
    // auto sP=prime(s);
    auto Op=ITensor(dag(s),sP);
    // Op.set(s(1), sP(1), 3);
    auto norm = 1/sqrt(3);
    auto omega = -0.5+sqrt(3)/2 * 1_i;
    auto omega_bar = -0.5-sqrt(3)/2 * 1_i;
    for (int i=1;i<=3;i++){
        Op.set(s(i), sP(1), norm);
        Op.set(s(1), sP(i), norm);
        Op.set(s(i+3), sP(4), norm);
        Op.set(s(4), sP(i+3), norm);
    }
    Op.set(s(2), sP(2), norm * omega);
    Op.set(s(3), sP(3), norm * omega);
    Op.set(s(5), sP(5), norm * omega);
    Op.set(s(6), sP(6), norm * omega);

    Op.set(s(2), sP(3), norm * omega_bar);
    Op.set(s(3), sP(2), norm * omega_bar);
    Op.set(s(5), sP(6), norm * omega_bar);
    Op.set(s(6), sP(5), norm * omega_bar);

    return Op;
}

MPS BasisRotate(MPS const& psi, SiteSet const& sites_new){
    int N = length(psi);
    auto new_psi = MPS(sites_new);
    for(auto j : range1(N))
    {
//        new_psi.set(j, psi(j) * delta(siteIndex(psi, j), sites_new(j)));
        new_psi.set(j, psi(j) * BasisRotation(siteIndex(psi, j), sites_new(j)));
    }
    return new_psi;
}
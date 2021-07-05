#include "chain.h"
#include "itensor/all.h"
#include <functional>
#include <vector>
#include <algorithm>

using namespace itensor;



MPO TranslationOp(const SiteSet& sites, bool inv) {
    int N = sites.length();
    auto G = std::vector<ITensor>(N);
    for (auto j : range1(N-1))
    {
        G[j] = BondGate(sites, j, j+1);
    }
    auto A = std::vector<ITensor>(N);
    auto B = std::vector<ITensor>(N);
    if (not inv) {
        for (auto j : range1(N-1)) {
            auto[Aj, Bj] = factor(G[j], {sites(j), prime(sites(j))}, {"MaxDim", 100});
            A[j] = Aj;
            B[j] = Bj;
        }
    } else {
        for (auto j : range1(N-1)) {
            auto[Aj, Bj] = factor(G[N-j], {sites(N-j), prime(sites(N-j))}, {"MaxDim", 100});
            A[N-j] = Aj;
            B[N-j] = Bj;
        }
    }
    auto t = std::vector<Index>(N+1);
    for (auto j : range1(N))
    {
        t[j] = sim(sites(j));
    }
    for (auto j : range1(2,N-1))
    {
        A[j] *= delta(sites(j),t[j]);
        B[j-1] *= delta(prime(sites(j)),t[j]);
    }
    auto m = MPO(N);
    m.set(1, A[1]);
    for (auto j : range1(2,N-1))
    {
        m.set(j, B[j-1] * A[j]);
    }
    m.set(N, B[N-1]);
    return m;
}

MPO IdentityOp(const SiteSet& sites, MPO const& op) {
    int N = sites.length();

    auto l = sim(leftLinkIndex(op, 2));
    auto r = sim(rightLinkIndex(op, 2));
    auto s = sim(sites(2));

//    auto dummy_i = Index(1);
//    auto dummy_j = Index(1);
//    auto dummy = ITensor(dummy_i, dummy_j);
//    dummy.set(dummy_i(1), dummy_j(1), 1);
//    auto identity_tensor = dummy * delta(l, r) * delta(s, prime(s)) * dummy;
//    See https://github.com/ITensor/ITensor/issues/176.
    auto identity_tensor = toDense(delta(l, r) ) * toDense( delta(s, prime(s)) );

    auto long_link = uniqueIndex(op(1), op(2), "Link");
    auto identity_mpo = MPO(N);
    identity_mpo.set(1, identity_tensor * delta(l, prime(long_link)) * delta(r, rightLinkIndex(op, 1)) * delta(s, sites(1)) * delta(prime(s), prime(sites(1))) );
    identity_mpo.set(N, identity_tensor * delta(l, leftLinkIndex(op, N)) * delta(r, long_link) * delta(s, sites(N)) * delta(prime(s), prime(sites(N))) );
    for (auto j : range1(2,N-1))
    {
        identity_mpo.set(j, identity_tensor * delta(l, leftLinkIndex(op, j)) * delta(r, rightLinkIndex(op, j)) * delta(s, sites(j)) * delta(prime(s), prime(sites(j))) );
    }

    identity_mpo.set(1, identity_mpo(1) * delta(prime(long_link), prime(long_link, 2)));
    identity_mpo.set(N, identity_mpo(N) * delta(long_link, prime(long_link)));
    return  identity_mpo;
}

ITensor Delta3ITensor(const Index& s1, const Index& s2, const Index& s3) {
    auto a = ITensor(s1,s2,s3);
    for (auto j : range1(s1))
    {
        a.set(s1(j),s2(j),s3(j),1.);
    }
    return a;
}

MPO RhoOp(const SiteSet& sites, const std::string& site_type) {
    int N = sites.length();

    auto G = std::vector<ITensor>(N+1);
    if (site_type == "golden") {
        auto f = GoldenFData();
        for (auto j : range1(N-1)) {
            G[j] = f.RhoDefectCell(sites(j), sites(j+1));
        }
        G[N] = GoldenFData().RhoDefectCell(sites(N), sites(1));
    } else if (site_type == "haagerup" || site_type == "haagerupq") {
        auto f = HaagerupFData();
        for (auto j : range1(N-1)) {
            G[j] = f.RhoDefectCell(sites(j), sites(j+1));
        }
        G[N] = f.RhoDefectCell(sites(N), sites(1));
    }
    auto A = std::vector<ITensor>(N+1);
    auto B = std::vector<ITensor>(N+1);
    for (auto j : range1(N))
    {
        auto [Aj,Bj] = factor(G[j],{sites(j),prime(sites(j))}); // fixme: Can refactor
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
    for (auto j : range1(N))
    {
        D12[j] = Delta3ITensor(dag(sites(j)), prime(sites(j),2), prime(sites(j),4) );
        D21[j] = Delta3ITensor(prime(sites(j)), dag(prime(sites(j), 3)), dag(prime(sites(j),5)) );
    }

    auto m = MPO(N);
    m.set(1, BNP * A[1] * D12[1] * D21[1]);
    for (auto j : range1(2,N))
    {
        m.set(j, B[j-1] * A[j] * D12[j] * D21[j]);
    }

    return m;
}

MPS AugmentMPS(MPS const& original_psi, Index const& sl, Index const& sr) {
    auto augmented_psi = MPS(length(original_psi) + 2);
    auto extra_site_index_left = Index(1, "l=0,Link");
    auto extra_l = ITensor(sl, extra_site_index_left);
    for (auto i:range1(dim(sl))) {
        extra_l.set(i, 1, 1. / dim(sl));
    }
    augmented_psi.set(1, extra_l);
    auto vl = ITensor(extra_site_index_left);
    vl.set(1, 1);
    augmented_psi.set(2, vl * original_psi(1));

    for (auto i:range1(length(original_psi) - 2)) {
        augmented_psi.set(i+2, original_psi(i+1));
    }

    auto extra_site_index_right = Index(1, "l=inf,Link");
    auto extra_r = ITensor(sr, extra_site_index_right);
    for (auto i:range1(dim(sr))) {
        extra_r.set(i, 1, 1. / dim(sr));
    }
    augmented_psi.set(length(original_psi) + 2, extra_r);
    auto vr = ITensor(extra_site_index_right);
    vr.set(1, 1);
    augmented_psi.set(length(original_psi) + 1, vr * original_psi(length(original_psi)));

    return augmented_psi;
}

MPO AugmentMPO(MPO const& original_mpo, Index const& sl, Index const& sr) {
    auto augmented_mpo = MPO(length(original_mpo) + 2);

    auto vl = ITensor(sl);
    vl.fill(1);
    auto index_left = uniqueIndex(original_mpo(1), original_mpo(2), "Link");
    auto Tl = delta(prime(sl,1), index_left);
    augmented_mpo.set(1, vl * Tl);

    for (auto i:range1(length(original_mpo) + 1)) {
        augmented_mpo.set(i+1, original_mpo(i));
    }

    auto vr = ITensor(sr);
    vr.fill(1);
    auto index_right = uniqueIndex(original_mpo(length(original_mpo)), original_mpo(length(original_mpo) - 1), "Link");
    auto Tr = delta(prime(sr,1), index_right);
    augmented_mpo.set(length(original_mpo) + 2, vr * Tr);

    return augmented_mpo;
}

ITensor Z3FourierMatrix(Index const& s, Index const& sP) {
    auto op = ITensor(dag(s), sP);
    auto norm = 1/sqrt(3);
    auto omega = -0.5+sqrt(3)/2 * 1_i;
    auto omega_bar = -0.5-sqrt(3)/2 * 1_i;

    for (int i=1;i<=3;i++) {
        op.set(s(i), sP(1), norm);
        op.set(s(1), sP(i), norm);
        op.set(s(i+3), sP(4), norm);
        op.set(s(4), sP(i+3), norm);
    }

    op.set(s(2), sP(2), norm * omega);
    op.set(s(3), sP(3), norm * omega);
    op.set(s(5), sP(5), norm * omega);
    op.set(s(6), sP(6), norm * omega);

    op.set(s(2), sP(3), norm * omega_bar);
    op.set(s(3), sP(2), norm * omega_bar);
    op.set(s(5), sP(6), norm * omega_bar);
    op.set(s(6), sP(5), norm * omega_bar);

    return op;
}

MPS Z3FourierTransform(MPS const& psi, SiteSet const& sites_new) {
    int N = length(psi);
    auto new_psi = MPS(sites_new);
    for (auto j : range1(N))
    {
        new_psi.set(j, psi(j) * Z3FourierMatrix(siteIndex(psi, j), sites_new(j)));
    }
    return new_psi;
}

// fixme
MPS ZipperAugmentMPS(MPS const& original_psi, Index const& sl, Index const& sll) {
    auto augmented_psi = MPS(length(original_psi) + 2);
    auto extra_bond_index_left = Index(1, "l=0,Link");
    auto extra_l = ITensor(sl, extra_bond_index_left);
    auto extra_bond_index_left_left = Index(1, "l=-1,Link");
    auto extra_ll = ITensor(sll, extra_bond_index_left_left);
    for (auto i:range1(dim(sl))) {
        extra_l.set(i, 1, 1. / dim(sl));
    }
    for (auto i:range1(dim(sll))) {
        extra_ll.set(i, 1, 1. / dim(sll));
    }

    augmented_psi.set(1, extra_ll);
    auto vll = ITensor(extra_bond_index_left_left);
    vll.set(1, 1);

    augmented_psi.set(2, vll * extra_l);
    auto vl = ITensor(extra_bond_index_left);
    vl.set(1, 1);

    augmented_psi.set(3, vl * original_psi(1));

    for (auto i:range1(length(original_psi) - 1)) {
        augmented_psi.set(i+3, original_psi(i+1));
    }

    return augmented_psi;
}



// fixme
//ITensor TetrahedralGate(const SiteSet& sites, const std::string& site_type, int b, int i, int j, int k, int l, int m, int n) {
//    auto G = ITensor(sites(b), sites(b+1), sites(b+2), prime(sites(b)), prime(sites(b+1)), prime(sites(b+2)));
//    if (site_type == "golden") {
//        auto f = GoldenFData();
//        for (auto j : range1(N-1)) {
//            G[j] = f.RhoDefectCell(sites(j), sites(j+1));
//        }
//        G[N] = GoldenFData().RhoDefectCell(sites(N), sites(1));
//    } else if (site_type == "haagerup" || site_type == "haagerupq") {
//        auto f = HaagerupFData();
//        for (auto j : range1(N-1)) {
//            G[j] = f.RhoDefectCell(sites(j), sites(j+1));
//        }
//        G[N] = f.RhoDefectCell(sites(N), sites(1));
//    }
//    auto A = std::vector<ITensor>(N+1);
//    auto B = std::vector<ITensor>(N+1);
//    for (auto j : range1(N))
//    {
//        auto [Aj,Bj] = factor(G[j],{sites(j),prime(sites(j))}); // fixme: Can refactor
//        Aj.prime("Site").prime("Site");
//        Bj.prime("Site").prime("Site").prime("Site").prime("Site");
//        A[j] = Aj;
//        B[j] = Bj;
//    }
//    auto long_link = commonIndex(A[N], B[N]);
//    auto BNP = ITensor(B[N]);
//    BNP.prime("Link");
//    auto D12 = std::vector<ITensor>(N+1);
//    auto D21 = std::vector<ITensor>(N+1);
//    for (auto j : range1(N))
//    {
//        D12[j] = Delta3ITensor(dag(sites(j)), prime(sites(j),2), prime(sites(j),4) );
//        D21[j] = Delta3ITensor(prime(sites(j)), dag(prime(sites(j), 3)), dag(prime(sites(j),5)) );
//    }
//
//    auto m = MPO(N);
//    m.set(1, BNP * A[1] * D12[1] * D21[1]);
//    for (auto j : range1(2,N))
//    {
//        m.set(j, B[j-1] * A[j] * D12[j] * D21[j]);
//    }
//
//    return m;
//}

//void InitialZip(MPS &psi, const SiteSet &sites) {
//    // Store original tags
//    psi.position(1);
//    auto tag = tags(rightLinkIndex(psi,1));
//
//    auto G = BondGate(sites, b, b+1);
//    auto wf = psi(b) * psi(b+1);
//    wf *= G;
//    wf.noPrime();
//
////    auto [U,S,V] = svd(wf,inds(psi(b)));
//    auto [U,S,V] = svd(wf,inds(psi(b)),{"Cutoff=",1E-5});
//    U.replaceTags(TagSet("u_,Link,0"), tag);
//    S.replaceTags(TagSet("u_,Link,0"), tag);
//
//    psi.set(b,U);
//    psi.set(b+1,S*V);
//}

// UNUSED

MPS
mydensityMatrixApplyMPOImpl(MPO const& K, MPS const& psi, Args args)
{
    if ( args.defined("Maxm") )
    {
        if ( args.defined("MaxDim") )
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
    if (maxdim_set) dargs.add("MaxDim",args.getInt("MaxDim"));
    dargs.add("RespectDegenerate",args.getBool("RespectDegenerate",true));
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);

    auto N = length(psi);

    for ( auto n : range1(N) )
    {
        if ( commonIndex(psi(n),K(n)) != siteIndex(psi,n) )
            Error("MPS and MPO have different site indices in applyMPO method 'DensityMatrix'");
    }

    auto rand_plev = 14741;

    auto res = psi;

    //Set up conjugate psi and k_
    auto psic = psi;
    auto Kc = K;
    psic.dag().prime(rand_plev);
    Kc.dag().prime(rand_plev);

    // Make sure the original and conjugates match

    for (auto j : range1(N-1)) {
        Kc.ref(j).prime(-rand_plev, uniqueSiteIndex(Kc, psic, j));
    }

    //Build environment tensors from the left
    if (verbose) print("Building environment tensors...");
    auto E = std::vector<ITensor>(N+1);
    for (int i=1;i<=psi.length();i++) {
        PrintData(linkInds(psi, i));
    }
    E[1] = psi(1)*K(1)*Kc(1)*psic(1);
    for (int j = 2; j < N; ++j)
    {
        // fixme: size of E[j] goes like MPSBondDim^2 * MPOBondDim^2 and MPSBondDim is typically large when j is in the middle of the chain
        E[j] = E[j-1]*psi(j)*K(j)*Kc(j)*psic(j);
    }
    if (verbose) println("done");

    //O is the representation of the product of k_*psi in the new MPS basis
    auto O = psi(N)*K(N);

    auto rho = E[N-1] * O * dag(prime(O,rand_plev));

    ITensor U,D;
    auto ts = tags(linkIndex(psi,N-1));
    auto spec = diagPosSemiDef(rho,U,D,{dargs,"Tags=",ts});
    if (verbose) printfln("  j=%02d truncerr=%.2E dim=%d",N-1,spec.truncerr(),dim(commonIndex(U,D)));

    res.ref(N) = dag(U);

    O = O*U*psi(N-1)*K(N-1);

    for (int j = N-1; j > 1; --j)
    {
        if (not maxdim_set)
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
        if (verbose) printfln("  j=%02d truncerr=%.2E dim=%d",j,spec.truncerr(),dim(commonIndex(U,D)));
    }

    if (normalize) O /= norm(O);
    res.ref(1) = O;
    res.leftLim(0);
    res.rightLim(2);

    return res;
}
//
void Swap(MPS &psi, const SiteSet &sites, int b) {
    // Store original tags
    psi.position(b);
    auto tag = tags(rightLinkIndex(psi,1));

    auto G = BondGate(sites, b, b+1);
    auto wf = psi(b) * psi(b+1);
    wf *= G;
    wf.noPrime();

//    auto [U,S,V] = svd(wf,inds(psi(b)));
    auto [U,S,V] = svd(wf,inds(psi(b)),{"Cutoff=",1E-5});
    // fixme: U vs u_
    U.replaceTags(TagSet("U,Link,0"), tag);
    S.replaceTags(TagSet("U,Link,0"), tag);

    psi.set(b,U);
    psi.set(b+1,S*V);
}

void localSwap(MPS &psi, int b){
    psi.position(b);
    auto tag = tags(rightLinkIndex(psi,1));
    auto s1 = siteIndex(psi, b);
    auto s2 = siteIndex(psi, b+1);
    auto x = ITensor(dag(s1),prime(s2));
    auto y = ITensor(dag(s2),prime(s1));
    for(auto j : range1(s1))
    {
        x.set(dag(s1)(j),prime(s2)(j),1.);
        y.set(dag(s2)(j),prime(s1)(j),1.);
    }
    auto wf = psi(b) * psi(b+1) * x * y;
    wf.noPrime();
    auto [U,S,V] = svd(wf,inds(psi(b)),{"Cutoff=",1E-5});

    U.replaceTags(TagSet("U,Link,0"), tag);
    S.replaceTags(TagSet("U,Link,0"), tag);

    psi.set(b,U);
    psi.set(b+1,S*V);
}

void ActThree(MPS & psi, ITensor const& G, int b){
    psi.position(b);
    auto tag1 = tags(rightLinkIndex(psi, 1));
    auto tag2 = tags(rightLinkIndex(psi, 2));
    auto wf = psi(b) * psi(b+1) * psi(b+2);
    wf = wf * G;
    wf.noPrime();

    auto [U1,S1,V1] = svd(wf, inds(psi(b)));
    U1.replaceTags(TagSet("U,Link,0"), tag1);
    S1.replaceTags(TagSet("U,Link,0"), tag1);
    auto [U2,S2,V2] = svd(S1*V1, IndexSet(siteIndex(psi,b+1), commonIndex(U1, S1)));
    U2.replaceTags(TagSet("U,Link,0"), tag2);
    S2.replaceTags(TagSet("U,Link,0"), tag2);
    psi.set(b, U1);
    psi.set(b+1, U2);
    psi.set(b+2, S2*V2);
}
//
//void ActLocal(MPS &psi, const ITensor &G, int b) {
//    // Store original tags
//    psi.position(b);
//    auto tag = tags(rightLinkIndex(psi,1));
//
//    auto wf = psi(b) * psi(b+1);
//    wf *= G;
//    wf.noPrime();
//
//    auto [U,S,V] = svd(wf,inds(psi(b)),{"Truncate", false});
//    U.replaceTags(TagSet("u_,Link,0"), tag);
//    S.replaceTags(TagSet("u_,Link,0"), tag);
//
//    psi.set(b,U);
//    psi.set(b+1,S*V);
//}
//
//
//// fixme: modify to allow for gates which are member functions of generic site type
//void ActGlobal(MPS &psi, const SiteSet &sites, TwoSiteGate gate, const std::string& sitetype_) {
//    int len = psi.length();
//    if (sitetype_ == "golden") {
//        auto f = GoldenFData();
//        for (int b = 1; b < len; b++) {
//            ActLocal(psi, std::invoke(gate, f, sites(b), sites(b + 1)), b);
//        }
//    } else if (sitetype_ == "haagerup") {
//        auto f = HaagerupFData();
//        for (int b = 1; b < len; b++) {
//            ActLocal(psi, std::invoke(gate, f, sites(b), sites(b + 1)), b);
//        }
//    }
//}
ThreeSiteGate::ThreeSiteGate(const SiteSet &sites, int i1, int i2, int i3) {
    std::vector<int> inds = {i1, i2, i3};
    std::sort (inds.begin(), inds.end());
    i1_ = inds[0];
    i2_ = inds[1];
    i3_ = inds[2];
    makeSwapGate(sites);
}

void ThreeSiteGate::makeSwapGate(const SiteSet &sites) {
    auto s1 = sites(i1_);
    auto s2 = sites(i2_);
    auto s3 = sites(i3_);
    auto a = ITensor(dag(s1),prime(s3));
    auto b = ITensor(dag(s2),prime(s1));
    auto c = ITensor(dag(s3),prime(s2));
    for(auto j : range1(s1))
    {
        a.set(dag(s1)(j),prime(s3)(j),1.);
        b.set(dag(s2)(j),prime(s1)(j),1.);
        c.set(dag(s3)(j),prime(s2)(j),1.);
    }
    gate_ = a*b*c;
}

ThreeSiteGate::ThreeSiteGate(const SiteSet &sites, int i1, int i2, int i3, ITensor gate) {
    std::vector<int> inds = {i1, i2, i3};
    std::sort (inds.begin(), inds.end());
    i1_ = inds[0];
    i2_ = inds[1];
    i3_ = inds[2];
    gate_ = gate;
}
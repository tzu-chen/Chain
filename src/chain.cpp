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


MPO RhoOp(const SiteSet &sites) {
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

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
        // fixme: investigate why this fails for Haagerup
        B[j-1] *= delta(prime(sites(j)),t[j]);
        A[j] *= delta(sites(j),t[j]);
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
    int N = sites.length();
    auto G = std::vector<ITensor>(N);
    // fixme: modify to allow for gates which are member functions of generic site type
    for (auto j : range1(N - 1)) {
        G[j] = GoldenFData().RhoDefect(sites(j), sites(j + 1));
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
    auto m = MPO(N);
    m.set(1, A[1]);
    for(auto j : range1(2,N-1))
    {
        m.set(j, B[j - 1] * A[j]);
    }
    m.set(N, B[N - 1]);
    return m;
}



void ActLocal(MPS &psi, const ITensor &G, int b) {
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


// fixme: modify to allow for gates which are member functions of generic site type
void ActGlobal(MPS &psi, const SiteSet &sites, TwoSiteGate gate) {
    int len = psi.length();
    for (int b = 1; b < len; b++) {
        ActLocal(psi, std::invoke(gate, GoldenFData(), sites(b), sites(b + 1)), b);
    }

}

module Model

using ITensors
using ITensorMPS: AutoMPO, MPO
import ITensors: op, space, OpName, SiteType, @SiteType_str

struct Anyon end

export AnyonModel, AnyonSite, AnyonChain, FSymbol, site, hamiltonian,
       fibonacci_model, haagerup_model

"""
    AnyonModel(rank; fsymbols, qdims)

Create a generic fusion category model. `rank` is the number of simple
anyons. `fsymbols` is a 6-index array `fsymbols[a,b,c,d,e,f]` giving the
`F`-symbols and should have dimensions `(rank,rank,rank,rank,rank,rank)`.
`qdims` is a vector of quantum dimensions of length `rank`.
"""
struct AnyonModel
    rank::Int
    fsymbols::Array{ComplexF64,6}
    qdims::Vector{Float64}
end

const active_model = Ref{AnyonModel}()

function AnyonModel(rank::Int; fsymbols=zeros(ComplexF64,rank,rank,rank,rank,rank,rank), qdims=ones(Float64,rank))
    size(fsymbols) == (rank,rank,rank,rank,rank,rank) ||
        error("fsymbols must have size (rank,rank,rank,rank,rank,rank)")
    length(qdims) == rank || error("qdims must have length rank")
    return AnyonModel(rank, ComplexF64.(fsymbols), Float64.(qdims))
end

"""Return the F-symbol F^{abc}_{def}."""
FSymbol(model::AnyonModel, a,b,c,d,e,f) = model.fsymbols[a,b,c,d,e,f]

"""Return an Index for site `n`."""
site(model::AnyonModel, n::Integer) = Index(model.rank, "Anyon,Site,n=$n")

"""
    struct AnyonSite

Represents a lattice site for a given `AnyonModel`.
"""
struct AnyonSite
    model::AnyonModel
    s::Index
end

space(::SiteType"Anyon"; kwargs...) = active_model[].rank

function op(o::OpName, ::SiteType"Anyon"; kwargs...)
    m = active_model[]
    name = String(ITensors.name(o))
    if name == "id"
        return Matrix{Float64}(I, m.rank, m.rank)
    elseif startswith(name, "n")
        i = parse(Int, last(split(name, "n")))
        M = zeros(Float64, m.rank, m.rank)
        M[i,i] = 1.0
        return M
    elseif name == "FF"
        F = zeros(ComplexF64, m.rank, m.rank)
        ρ = div(m.rank, 2) + 1
        for i in 1:m.rank, j in 1:m.rank
            F[i,j] = FSymbol(m, ρ, ρ, ρ, ρ, i, 1) *
                     FSymbol(m, ρ, ρ, ρ, ρ, 1, j)
        end
        return F
    elseif startswith(name, "FF_")
        parts = split(name[4:end], "_")
        p = parse(Int, parts[1])
        l = parse(Int, parts[2])
        r = parse(Int, parts[3])
        F = zeros(ComplexF64, m.rank, m.rank)
        for i in 1:m.rank
            val = FSymbol(m, l, div(m.rank, 2) + 1, div(m.rank, 2) + 1, r, i, p)
            F[i,i] = val * val
        end
        return F
    else
        return nothing
    end
end

AnyonSite(model::AnyonModel, n::Integer) =
    AnyonSite(model, site(model, n))

"""Projector onto anyon type `i`."""
function proj(site::AnyonSite, i::Integer)
    s = site.s
    sp = prime(s)
    op = ITensor(dag(s), sp)
    op[s(i), sp(i)] = 1.0
    return op
end

"""
    FF(site::AnyonSite)

Return the F-move operator associated with the `AnyonModel` of `site`.
It is constructed entirely from the `F`-symbols stored in the model.

Currently this assumes the nontrivial anyon is labelled `rank(model)`
and the trivial anyon is `1`, which is appropriate for the golden
chain.
"""
function FF(site::AnyonSite)
    m = site.model
    s = site.s
    sp = prime(s)
    ρ = div(m.rank, 2) + 1
    op = ITensor(dag(s), sp)
    for i in 1:m.rank, j in 1:m.rank
        op[s(i), sp(j)] = FSymbol(m, ρ, ρ, ρ, ρ, i, 1) *
                          FSymbol(m, ρ, ρ, ρ, ρ, 1, j)
    end
    return op
end

function FF(site::AnyonSite, projector::Int, left::Int, right::Int)
    m = site.model
    s = site.s
    sp = prime(s)
    op = ITensor(dag(s), sp)
    for i in 1:m.rank
        val = FSymbol(m, left, div(m.rank, 2) + 1, div(m.rank, 2) + 1, right, i, projector)
        if val != 0
            op[s(i), sp(i)] = val * val
        end
    end
    return op
end

function op(site::AnyonSite, name::String)
    s = site.s
    if name == "id"
        return ITensor(delta(s, prime(s)))
    elseif startswith(name, "n")
        i = parse(Int, name[2:end])
        return proj(site, i)
    elseif name == "FF"
        return FF(site)
    elseif startswith(name, "FF_")
        parts = split(name[4:end], "_")
        p = parse(Int, parts[1])
        l = parse(Int, parts[2])
        r = parse(Int, parts[3])
        return FF(site, p, l, r)
    else
        error("Operator $name not recognized")
    end
end

struct AnyonChain{M<:AnyonModel} <: AbstractVector{AnyonSite}
    model::M
    sites::Vector{AnyonSite}
end

AnyonChain(model::AnyonModel, N::Int) = begin
    active_model[] = model
    AnyonChain(model, [AnyonSite(model, n) for n in 1:N])
end

Base.length(c::AnyonChain) = length(c.sites)
Base.getindex(c::AnyonChain, i::Int) = c.sites[i]

"""Generic three-body Hamiltonian used for the golden chain."""
function hamiltonian(c::AnyonChain; boundary::String="p", couplings=[1.0])
    L = length(c)
    ampo = AutoMPO()
    Lproj = boundary == "p" ? L : L - 2
    if boundary == "sp"
        Lproj = L - 3
    end

    if c.model.rank == 2
        K = couplings[1]
        for j in 1:Lproj
            jp1 = mod(j, L) + 1
            jp2 = mod(j + 1, L) + 1
            ampo += K, "n1", j, "n2", jp1, "n1", jp2
            ampo += K, "n2", j, "FF", jp1, "n2", jp2
        end
    elseif c.model.rank == 6
        K = couplings[1]; J = couplings[2]; M = couplings[3]
        for j in 1:Lproj
            jp1 = mod(j, L) + 1
            jp2 = mod(j + 1, L) + 1
            if K != 0
                ampo += K, "n1", j, "n4", jp1, "n1", jp2
                ampo += K, "n2", j, "n5", jp1, "n2", jp2
                ampo += K, "n3", j, "n6", jp1, "n3", jp2
                ampo += K, "n4", j, "FF_1_4_4", jp1, "n4", jp2
                ampo += K, "n5", j, "FF_1_5_5", jp1, "n5", jp2
                ampo += K, "n6", j, "FF_1_6_6", jp1, "n6", jp2
            end
            if J != 0
                ampo += J, "n1", j, "n4", jp1, "n4", jp2
                ampo += J, "n4", j, "n4", jp1, "n1", jp2
                ampo += J, "n2", j, "n5", jp1, "n5", jp2
                ampo += J, "n5", j, "n5", jp1, "n2", jp2
                ampo += J, "n3", j, "n6", jp1, "n6", jp2
                ampo += J, "n6", j, "n6", jp1, "n3", jp2

                ampo += J, "n4", j, "FF_4_4_4", jp1, "n4", jp2
                ampo += J, "n4", j, "FF_4_4_5", jp1, "n5", jp2
                ampo += J, "n4", j, "FF_4_4_6", jp1, "n6", jp2
                ampo += J, "n5", j, "FF_4_5_4", jp1, "n4", jp2
                ampo += J, "n5", j, "FF_4_5_5", jp1, "n5", jp2
                ampo += J, "n5", j, "FF_4_5_6", jp1, "n6", jp2
                ampo += J, "n6", j, "FF_4_6_4", jp1, "n4", jp2
                ampo += J, "n6", j, "FF_4_6_5", jp1, "n5", jp2
                ampo += J, "n6", j, "FF_4_6_6", jp1, "n6", jp2
            end
            if M != 0
                ampo += M, "n1", j, "n4", jp1, "n5", jp2
                ampo += M, "n5", j, "n4", jp1, "n1", jp2
                ampo += M, "n2", j, "n5", jp1, "n6", jp2
                ampo += M, "n6", j, "n5", jp1, "n2", jp2
                ampo += M, "n3", j, "n6", jp1, "n4", jp2
                ampo += M, "n4", j, "n6", jp1, "n3", jp2

                ampo += M, "n4", j, "FF_5_4_4", jp1, "n4", jp2
                ampo += M, "n4", j, "FF_5_4_5", jp1, "n5", jp2
                ampo += M, "n4", j, "FF_5_4_6", jp1, "n6", jp2
                ampo += M, "n5", j, "FF_5_5_4", jp1, "n4", jp2
                ampo += M, "n5", j, "FF_5_5_5", jp1, "n5", jp2
                ampo += M, "n5", j, "FF_5_5_6", jp1, "n6", jp2
                ampo += M, "n6", j, "FF_5_6_4", jp1, "n4", jp2
                ampo += M, "n6", j, "FF_5_6_5", jp1, "n5", jp2
                ampo += M, "n6", j, "FF_5_6_6", jp1, "n6", jp2
            end
        end
    else
        error("hamiltonian not implemented for model with rank $(c.model.rank)")
    end
    return MPO(ampo, [s.s for s in c.sites])
end

"""Return a predefined Fibonacci anyon model."""
function fibonacci_model()
    phi = (sqrt(5)+1)/2
    phi_inv = 1/phi
    sqrt_phi_inv = sqrt(phi_inv)
    fs = Array{ComplexF64,6}(undef, 2,2,2,2,2,2)
    dual(i) = i
    isinv(i) = i == 1
    function add(i,j)
        2
    end
    function fusion(a,b)
        if a == 1 && b == 1
            return (1,)
        elseif a == 1 && b == 2
            return (2,)
        elseif a == 2 && b == 1
            return (2,)
        else
            return (1,2)
        end
    end
    function hasfusion(i,j,k)
        k in fusion(i,j)
    end
    function fsym(i,j,k,l,m,n)
        if i == 2 && j == 2 && k == 2 && m == 2 && !isinv(l) && !isinv(n)
            return -phi_inv
        end
        if !(hasfusion(i,j,m) && hasfusion(k,dual(l),dual(m)) &&
              hasfusion(dual(l), i, dual(n)) && hasfusion(j,k,n))
            return 0.0
        end
        if isinv(i) || isinv(j) || isinv(k) || isinv(l)
            return 1.0
        elseif isinv(m) && isinv(n)
            return phi_inv
        elseif isinv(m) || isinv(n)
            return sqrt_phi_inv
        else
            return -phi_inv
        end
    end
    for a in 1:2, b in 1:2, c in 1:2, d in 1:2, e in 1:2, f in 1:2
        fs[a,b,c,d,e,f] = fsym(a,b,c,d,e,f)
    end
    AnyonModel(2; fsymbols=fs, qdims=[1.0, phi])
end

"""Return the Haagerup anyon model."""
function haagerup_model()
    nu = 3
    rho = nu + 1
    rk = 2 * nu
    phi = (sqrt(4 + nu^2) + nu) / 2
    phi_inv = 1 / phi
    sqrt_phi_inv = sqrt(phi_inv)

    x = (2 - sqrt(13)) / 3
    y1 = (5 - sqrt(13) - sqrt(6 * (1 + sqrt(13)))) / 12
    y2 = (5 - sqrt(13) + sqrt(6 * (1 + sqrt(13)))) / 12
    z = (1 + sqrt(13)) / 6

    fs = Array{ComplexF64,6}(undef, rk, rk, rk, rk, rk, rk)

    isinv(i) = i <= nu
    dual(i) = i <= nu ? 1 + ((nu + 1 - i) % nu) : i
    add(i, j) = rho + ((i + j - 1) % nu)

    function fusion(a, b)
        if a <= nu && b <= nu
            return (1 + ((a + b - 2) % nu),)
        elseif a <= nu && b > nu
            return (rho + ((a + b - 2) % nu),)
        elseif a > nu && b <= nu
            return fusion(1 + ((rho - b) % nu), a)
        else
            return vcat((1 + ((nu + a - b) % nu)), collect(rho:rk))
        end
    end

    hasfusion(i, j, k) = k in fusion(i, j)

    F = Ref{Function}()

    function fsym_pattern(i, j, k, l, m, n)
        if !(hasfusion(i, j, m) && hasfusion(k, dual(l), dual(m)) &&
              hasfusion(dual(l), i, dual(n)) && hasfusion(j, k, n))
            return 0.0
        elseif isinv(i) || isinv(j) || isinv(k) || isinv(l)
            return 1.0
        elseif isinv(m) && isinv(n)
            return phi_inv
        elseif isinv(m) || isinv(n)
            return sqrt_phi_inv
        elseif i != rho
            return F[](rho, j, add(k, i - rho), l, m, n)
        elseif j != rho
            return F[](rho, rho, k, add(l, j - rho), m, n)
        elseif k != rho
            return F[](rho, rho, rho, add(l, rho - k), m, add(n, rho - k))
        elseif m != rho
            return F[](rho, rho, rho, l, rho, add(n, m - rho))
        else
            error("FSymbolPattern failed")
        end
    end

    F[] = function (i, j, k, l, m, n)
        if i == rho && j == rho && k == rho && m == rho && !isinv(l) && !isinv(n)
            s = l + n
            if s == 8
                return x
            elseif s == 9
                return y1
            elseif s == 10
                return y2
            elseif s == 11
                return z
            elseif s == 12
                return y1
            end
        end
        return fsym_pattern(i, j, k, l, m, n)
    end

    for a in 1:rk, b in 1:rk, c in 1:rk, d in 1:rk, e in 1:rk, f in 1:rk
        fs[a, b, c, d, e, f] = F[](a, b, c, d, e, f)
    end

    qdims = [ones(Float64, nu); fill(phi, rk - nu)]
    AnyonModel(rk; fsymbols=fs, qdims=qdims)
end

end # module

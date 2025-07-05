module GoldenModel

using ITensors
using ITensorMPS

import ITensors: op, space, OpName, SiteType

struct GoldenSite
    s::Index
end

GoldenSite() = GoldenSite(Index(2, "Site,Golden"))

function op(site::GoldenSite, name::String)
    s = site.s
    if name == "id"
        return ITensor(delta(s, prime(s)))
    elseif name == "n1"
        return proj(site, 1)
    elseif name == "nt"
        return proj(site, 2)
    elseif name == "FF"
        return FF(site)
    else
        error("Operator $name not recognized")
    end
end

function proj(site::GoldenSite, i::Int)
    s = site.s
    sp = prime(s)
    op = ITensor(dag(s), sp)
    op[s(i), sp(i)] = 1.0
    return op
end

const phi = (sqrt(5)+1)/2
const phi_inv = 1/phi
const sqrt_phi_inv = sqrt(phi_inv)

function FF(site::GoldenSite)
    s = site.s
    sp = prime(s)
    Op = ITensor(dag(s), sp)
    Op[s(1), sp(1)] = phi_inv^2
    Op[s(1), sp(2)] = sqrt_phi_inv*phi_inv
    Op[s(2), sp(1)] = sqrt_phi_inv*phi_inv
    Op[s(2), sp(2)] = phi_inv
    return Op
end

# Define site type interface for ITensors

space(::SiteType"Golden"; kwargs...) = 2

function op(o::OpName, ::SiteType"Golden"; kwargs...)
    name = String(ITensors.name(o))
    if name == "id"
        return [1.0 0.0; 0.0 1.0]
    elseif name == "n1"
        return [1.0 0.0; 0.0 0.0]
    elseif name == "nt"
        return [0.0 0.0; 0.0 1.0]
    elseif name == "FF"
        return [phi_inv^2 sqrt_phi_inv*phi_inv; sqrt_phi_inv*phi_inv phi_inv]
    else
        return nothing
    end
end

struct Golden <: AbstractVector{GoldenSite}
    sites::Vector{GoldenSite}
end

function Golden(N::Int)
    Golden([GoldenSite() for _ in 1:N])
end

Base.length(g::Golden) = length(g.sites)
Base.getindex(g::Golden, i::Int) = g.sites[i]

function hamiltonian(g::Golden; boundary::String="p", U::Real=0.0, couplings=[1.0])
    L = length(g)
    K = couplings[1]
    mpo = AutoMPO()
    if boundary != "p"
        L = L - 1
    end
    if boundary == "sp"
        L = L - 1
    end
    for j in 1:L
        mpo += K, "n1", j, "nt", mod(j+1-1,length(g))+1, "n1", mod(j+2-1,length(g))+1
        mpo += K, "nt", j, "FF", mod(j+1-1,length(g))+1, "nt", mod(j+2-1,length(g))+1
    end
    return MPO(mpo, [s.s for s in g.sites])
end

end

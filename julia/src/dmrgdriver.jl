module DMRGDriver

using ITensors
using ITensorMPS: random_mps, Sweeps, maxdim!, cutoff!, dmrg
using JLD2
using AMDGPU
using Adapt
using ..Model

export run_dmrg

function run_dmrg(;L::Int, maxdim::Int=100, cutoff::Float64=1e-8, sweeps::Int=5,
                  boundary::String="p", U::Real=0.0,
                  couplings::Vector{<:Real}=[1.0],
                  out::String="dmrg.jld2", model::AnyonModel=fibonacci_model(),
                  amdgpu::Bool=false)
    sites = AnyonChain(model, L)
    ampo = hamiltonian(sites; boundary=boundary, U=U, couplings=couplings)
    H = ampo
    psi0 = random_mps([s.s for s in sites.sites])
    if amdgpu
        H = AMDGPU.roc(H)
        psi0 = AMDGPU.roc(psi0)
    end
    sweepset = Sweeps(sweeps)
    maxdim!(sweepset, maxdim)
    cutoff!(sweepset, cutoff)
    energy, psi = dmrg(H, psi0, sweepset; outputlevel=0)
    if amdgpu
        psi = Adapt.adapt(Array, psi)
    end
    @save out energy psi
    return energy
end

end

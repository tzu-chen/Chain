module DMRGDriver

using ITensors
using ITensorMPS: random_mps, Sweeps, maxdim!, cutoff!, dmrg
using JLD2
using ..GoldenModel

export run_dmrg

function run_dmrg(;L::Int, maxdim::Int=100, cutoff::Float64=1e-8, sweeps::Int=5,
                   boundary::String="p", out::String="dmrg.jld2")
    sites = GoldenModel.Golden(L)
    ampo = GoldenModel.hamiltonian(sites; boundary=boundary)
    H = ampo
    psi0 = random_mps([s.s for s in sites.sites])
    sweepset = Sweeps(sweeps)
    maxdim!(sweepset, maxdim)
    cutoff!(sweepset, cutoff)
    energy, psi = dmrg(H, psi0, sweepset, silent=true)
    @save out energy psi
    return energy
end

end

#!/usr/bin/env julia
using Chain
using Chain.Model
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--length", "-l"
            help = "number of sites"
            arg_type = Int
        "--maxdim", "-d"
            help = "maximum bond dimension"
            arg_type = Int
            default = 100
        "--sweeps", "-s"
            help = "number of sweeps"
            arg_type = Int
            default = 5
        "--boundary", "-b"
            help = "boundary condition (p/o/s/sp)"
            arg_type = String
            default = "p"
        "--couplings", "-j"
            help = "comma-separated list of coupling coefficients"
            arg_type = String
            default = "1.0"
        "--penalty", "-u"
            help = "penalty coefficient"
            arg_type = Float64
            default = 0.0
        "--model", "-m"
            help = "anyon model (golden/haagerup)"
            arg_type = String
            default = "golden"
        "--nstates", "-n"
            help = "number of states to target"
            arg_type = Int
            default = 1
        "--weight", "-w"
            help = "orthogonality weight"
            arg_type = Float64
            default = 1000.0
        "--out", "-o"
            help = "output JLD2 file"
            arg_type = String
            default = "dmrg.jld2"
        "--amdgpu", "-g"
            help = "use AMDGPU backend"
            action = :store_true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    L = args["length"] === nothing ? args["l"] : args["length"]
    maxdim = args["maxdim"] === nothing ? args["d"] : args["maxdim"]
    sweeps = args["sweeps"] === nothing ? args["s"] : args["sweeps"]
    boundary = args["boundary"] === nothing ? args["b"] : args["boundary"]
    couplings_str = args["couplings"] === nothing ? args["j"] : args["couplings"]
    couplings = parse.(Float64, split(couplings_str, ","))
    U = args["penalty"] === nothing ? args["u"] : args["penalty"]
    modelname = args["model"] === nothing ? args["m"] : args["model"]
    model = lowercase(modelname) == "haagerup" ? Model.haagerup_model() : Model.fibonacci_model()
    nstates = args["nstates"] === nothing ? args["n"] : args["nstates"]
    weight = args["weight"] === nothing ? args["w"] : args["weight"]
    out = args["out"] === nothing ? args["o"] : args["out"]
    use_amdgpu = get(args, "amdgpu", false)
    run_dmrg(L=L, maxdim=maxdim, sweeps=sweeps,
             boundary=boundary, couplings=couplings, U=U, out=out,
             model=model, amdgpu=use_amdgpu,
             nstates=nstates, weight=weight)
end

main()

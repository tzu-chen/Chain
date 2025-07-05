#!/usr/bin/env julia
using Chain
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--length" "-l"
            help = "number of sites"
            arg_type = Int
        "--maxdim" "-d"
            help = "maximum bond dimension"
            arg_type = Int
            default = 100
        "--sweeps" "-s"
            help = "number of sweeps"
            arg_type = Int
            default = 5
        "--boundary" "-b"
            help = "boundary condition (p/o/s/sp)"
            arg_type = String
            default = "p"
        "--out" "-o"
            help = "output JLD2 file"
            arg_type = String
            default = "dmrg.jld2"
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    L = args["length"] === nothing ? args["l"] : args["length"]
    maxdim = args["maxdim"] === nothing ? args["d"] : args["maxdim"]
    sweeps = args["sweeps"] === nothing ? args["s"] : args["sweeps"]
    boundary = args["boundary"] === nothing ? args["b"] : args["boundary"]
    out = args["out"] === nothing ? args["o"] : args["out"]
    run_dmrg(L=L, maxdim=maxdim, sweeps=sweeps,
             boundary=boundary, out=out)
end

main()

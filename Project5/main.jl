using LinearAlgebra
using Random
# using Pkg
# Pkg.add("NPZ")
# Pkg.add("ProgressMeter")
# Pkg.add("SparseArrays")
# Pkg.add("IterativeSolvers")
using NPZ
using ProgressMeter
using SparseArrays
using IterativeSolvers


include("Crank_Nicolson.jl")
include("sparse_matrices.jl")


function problem7_no_slit()
    # Free parameters
    T = 0.008
    sy = 0.05
    v0 = 0
    slits = 0

    # The 7 other system parameters should not be changed
    args = [0.005, 2.5e-5, T, 0.25, 0.5, 200, 0, 0.05, sy, v0, slits]
    name = "p7_no_slit"
    simulate(args, name)
end

function problem7_double_slit()
    # Free parameters
    T = 0.002
    sy = 0.10
    v0 = 1e10
    slits = 2

    # The 7 other system parameters should not be changed
    args = [0.005, 2.5e-5, T, 0.25, 0.5, 200, 0, 0.05, sy, v0, slits]
    name = "p7_double_slit"
    simulate(args, name)
end

function problem8()
    # Free parameters
    T = 0.002
    sy = 0.20
    v0 = 1e10
    slits = 2

    realimag = true

    # The 7 other system parameters should not be changed
    args = [0.005, 2.5e-5, T, 0.25, 0.5, 200, 0, 0.05, sy, v0, slits]
    name = "p8"
    simulate(args, name, realimag)
end

function problem9_one_slit()
    # Free parameters
    T = 0.002
    sy = 0.20
    v0 = 1e10
    slits = 1

    # The 7 other system parameters should not be changed
    args = [0.005, 2.5e-5, T, 0.25, 0.5, 200, 0, 0.05, sy, v0, slits]
    name = "p9_one_slit"
    simulate(args, name)
end

function problem9_two_slit()
    # Free parameters
    T = 0.002
    sy = 0.20
    v0 = 1e10
    slits = 2

    # The 7 other system parameters should not be changed
    args = [0.005, 2.5e-5, T, 0.25, 0.5, 200, 0, 0.05, sy, v0, slits]
    name = "p9_two_slit"
    simulate(args, name)
end

function problem9_three_slit()
    # Free parameters
    T = 0.002
    sy = 0.20
    v0 = 1e10
    slits = 3

    # The 7 other system parameters should not be changed
    args = [0.005, 2.5e-5, T, 0.25, 0.5, 200, 0, 0.05, sy, v0, slits]
    name = "p9_three_slit"
    simulate(args, name)
end

if abspath(PROGRAM_FILE) == @__FILE__
    # potential is not made, only this one works
    # problem7_no_slit()
    # problem7_double_slit()
    problem8()
    # problem9_one_slit()
    # problem9_two_slit()
    # problem9_three_slit()
end

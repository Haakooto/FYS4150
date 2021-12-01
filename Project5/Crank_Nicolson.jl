using LinearAlgebra
using Random
using SparseArrays
# using Pkg
# Pkg.add("NPZ")
using NPZ

include("sparse_matrices.jl")

function make_AB(M, h, dt, V)
    m = M - 2
    a = ones(ComplexF64, m*m)
    b = ones(ComplexF64, m*m)
    r = 1im * dt / 2  # defined without the h^2
    for k in 1:m*m
        ch = 4r / h^2 + r * V[k]
        a[k] += ch
        b[k] -= ch
    end
    return SparseMat(a, -r), SparseMat(b, r)
end

function make_V(M, h, v0, slits=2)
    m = M - 2
    V = zeros(m, m)
    return V
end

function Gauss_wave_packet(M, h, xc, yc, px, py, sx, sy)
    x = [i * h for i in 0:M-1]
    X, Y = meshgrid(x, x)

    t1x = -((X .- xc).^2) ./ (2 * sx^2)
    t1y = -((Y .- yc).^2) ./ (2 * sy^2)
    t2x = 1im * px .* (X .- xc)
    t2y = 1im * py .* (Y .- yc)

    u = exp.(t1x .+ t1y .+ t2x .+ t2y)

    u /= sqrt(sum(abs2.(u)))  # normalize u
    return u
end

function initialize_system(M, dt, h, v0, P)
    V = make_V(M, h, v0)
    u = Gauss_wave_packet(M, h, P...)
    A, B = make_AB(M, h, dt, V)
    return A, B, u
end

function main()
    h = 0.005
    dt = 2.5e-5
    T = 0.008

    xc = 0.25
    sx = 0.005
    px = 200
    yc = 0.5
    sy = 0.05
    py = 0
    v0 = 0

    gauss_params = [xc, yc, px, py, sx, sy]

    M = Int(1 / h)
    nT = Int(T / dt)

    A, B, u = initialize_system(M, dt, h, v0, gauss_params)
    u = vec(u[2:end-1, 2:end-1])

    probs = zeros(Float64, nT, (M - 2)^2 + 1)
    probs[:, 1] = collect(range(0, T, nT))
    probs[1, 2:end] = abs2.(u)

    for i in 2:nT
        b = B * u
        u, _ = SOR(A, b, initial_guess=u, omega=1, tol=1e-13)
        probs[i, 2:end] = abs2.(u)
    end
    println(sum(abs2.(u)))
    npzwrite("evolve3.npz", probs)
end

main()

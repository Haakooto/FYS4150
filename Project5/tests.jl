using Random
using LinearAlgebra
using NPZ
include("sparse_matrices.jl")
include("Crank_Nicolson.jl")


function test_multiplication()
    #=
    verify that multiplaction is done correctly with the sparse matrix
    =#
    n = 4^2  # must be a square number
    r = 1im / 2
    d = ones(n) .+ 4r  # zero potential
    A_sp = SparseMat(d, r)

    x = randn(ComplexF64, n)  # random vector

    Ax_sp = A_sp * x
    Ax_ds = A_ds * x

    r = sum(Ax_sp .- Ax_ds)
    @assert real(r) < 1e-14
    @assert imag(r) < 1e-14

    println("Multiplication test passed")

end

function test_SOR()
    #=
    Verify that the successive over relaxation for
    solving matrix problems is implemented correctly,
    using the sparse matricies
    =#
    rng = MersenneTwister(14)

    n = 30 ^ 2
    r = 1im / 2
    d = ones(n) .+ 4r

    A_sp = SparseMat(d, r)

    b = randn(rng, ComplexF64, n)

    x_SOR, i = SOR(A_sp, b, omega=1)
    println("SOR used ", i, " iterations")
    x_INV = inv(A_ds) * b

    r = sum(x_SOR .- x_INV)
    @assert real(r) < 1e-14
    @assert imag(r) < 1e-14
    println("SOR test passed")
end

function determine_best_omega()
    rng = MersenneTwister(14)

    h = 0.01
    dt = h^2
    V = make_V(h, 0, 2)
    A, B = make_AB(h, dt, V)
    n = size(A.diag)[1]

    timesteps = 20

    ws = collect(range(0.7, 1.1, step=0.1))
    is = similar(ws)

    for k in 1:length(ws)
        println("omega = ", ws[k])
        x = randn(rng, ComplexF64, n)
        i_ = zeros(timesteps)
        for t in 1:timesteps
            b = B * x
            x, i_[t] = SOR(A, b, omega=ws[k], initial_guess=x)
            print(i_[t], " ")
        end
        is[k] = sum(i_[2:end]) / length(i_[2:end])
        println("Average: ", is[k])
    end
    println("For M = ", sqrt(n), ", the omega giving fewest SOR-iterations over ", timesteps, " was ", ws[argmin(is)], " averaging ", minimum(is), " iterations per SOR")
end

function test_sparse()
    h = 0.005
    dt = h^2
    V = make_V(h, 0, 2)
    A, B = make_AB(h, dt, V)
    # A_ds = to_dense(A)
    A_sp = to_sparse(A)

    n = size(A.diag)[1]
    println(sqrt(n))
    b = randn(ComplexF64, n)

    @time x_or, i = SOR(A, b, omega=0.9)
    # @time x_ds = inv(A_ds) * b
    @time x_sp = A_sp \ b
    @time x_spor, j = SOR(A_sp, b, omega=0.9)

    println("iterations: ", i)
    println("iterations: ", j)
    # r = sum(x_or .- x_ds)
    # s = sum(x_ds .- x_sp)
    t = sum(x_sp .- x_or)
    r = sum(x_spor .- x_or)
    s = sum(x_spor .- x_sp)

    println(r)
    println(s)
    println(t)

end

# test_multiplication()
# test_SOR()
# determine_best_omega()
test_sparse()


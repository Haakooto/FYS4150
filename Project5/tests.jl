using Random
include("sparse_matrices.jl")

function denseMat(d, r)
    """
    Set up the equivalent dense matrix as the sparse one
    """
    n = length(d)
    m = Int(sqrt(n))
    A = zeros(ComplexF64, n, n)
    for i in 1:n
        A[i, i] = d[i]
        if i > 1
            if mod(i, m) != 1
                A[i, i - 1] = -r
                A[i - 1, i] = -r
            end
            if i > m
                A[i, i - m] = -r
                A[i - m, i] = -r
            end
        end
    end
    return A
end

function test_multiplication()
    """
    verify that multiplaction is done correctly with the sparse matrix
    """
    n = 16
    r = pi * 1im
    d = ones(n) .+ 4r
    A_sp = SparseMat(d, r)
    A_ds = denseMat(d, -r)

    x = randn(ComplexF64, n)

    show_SparseMat(A_sp)

    # for j in [-A_sp.m, -1, 0, 1, A_sp.m]
    #     println(get_at(A_sp, 1, 1 + j))
    # end
    Ax_sp = A_sp * x
    Ax_ds = A_ds * x

    r = sum(Ax_sp .- Ax_ds)
    @assert real(r) < 1e-15
    @assert imag(r) < 1e-15

    println("Multiplication works")

end

function test_SOR()
    """
    Verify that the successive over relaxation for
    solving matrix problems is implemented correctly,
    using the sparse matricies
    """
    rng = MersenneTwister(14)

    n = 9
    r = 1im / 2
    d = ones(n) .+ 4r
    # d[1] *= 10

    A_sp = SparseMat(d, r)
    A_ds = denseMat(d, r)
    # z = zeros(ComplexF64, n)
    z = 0

    b = randn(rng, ComplexF64, n)

    x_SOR, i = SOR(A_sp, b, initial_guess=z, omega=0.5, tol=1e-10)
    println(norm(b .- (A_sp * x_SOR)))
    # x_SOR, i = SOR(A_sp, b, omega=0.0001, tol=1e-7)
    # x_SOR, i = SOR_wiki(A_ds, b, omega=1.5, tol=1e-7)

    x_INV = inv(A_ds) * b
    println()
    println()
    println(i)
    display(x_SOR)
    println()
    display(x_INV)
    println()
    println(norm(b .- A_ds * x_INV))
    println(norm(x_INV .- x_SOR))
    # r = sum(A_sp .- A_ds)
    # @assert real(r) < 1e-15
    # @assert imag(r) < 1e-15
    # println("SOR works")
end


# test_multiplication()

function SOR_wiki(A, b; omega=1.5, tol=1e-4)
    n = length(b)
    step = 0
    x = similar(b)
    randn!(x)
    residual = 1
    # println(x)
    while residual > tol
        # for what in 1:10
        residual = 0
        for i in 1:n
            sigma = 0
            for j in 1:n
                if i != j
                    sigma += A[i, j] * x[j]
                    # println(sigma)
                end
            end
            update = (omega / A[i, i]) * (b[i] - sigma) - omega * x[i]
            x[i] = x[i] + update
            #println(update, " ", x)
            #println()
            if abs(update) > residual
                residual = abs(update)
            end
        end
        step += 1
        if step > 1e5
            break
        end
    end
    return x, step
end
test_SOR()



function SOR_o(A, b, w=1.5, tol=1e-13)
    m = length(b)
    x = similar(b)  #  make vector of same type
    randn!(x)  #  initiate random vector
    max_change = 1  #  something to get it going
    iters = 0  #  count iterations needed
    while max_change > tol  # run until convergence
        iters += 1
        max_change = 0
        for i in 1:m
            Sum = 0  # sum(Aij * xj)
            for j in 1:m  # This way also does Gauss-Seidel
                Sum += A[i, j] * x[j]
            end
            update = w * (b[i] - Sum) / A[i, i]
            x[i] = x[i] + update
            if abs(update) > max_change
                max_change = abs(update)
            end
        end
    end
    return x, iters
end

# panic()

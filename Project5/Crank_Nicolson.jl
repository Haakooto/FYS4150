using LinearAlgebra
using Random
using NPZ
# using DataFrames
# using Plotly

println("Done loading packages")

ij_to_k(i, j, m) = return (i - 1) * m + j

function make_AB_(r, a, b)
    n = length(a)
    m = Int(sqrt(n))
    A = zeros(ComplexF64, n, n)
    for i in 1:n
        A[i, i] = a[i]
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
    if b == 0
        return A
    else
        return A, make_AB_(-r, b, 0)
    end
end


function make_AB(M, h, dt, V)
    m = M - 2
    a = ones(ComplexF64, m*m)
    b = ones(ComplexF64, m*m)
    r = 1im * dt / 2
    for k in 1:m*m
        ch = 4 * r / h^2 + r * V[k]
        a[k] += ch
        b[k] -= ch
    end
    return make_AB_(r, a, b)
end

function advance(u, A, B)
    b = B * u
    u, i = SOR(A, b)
    # u = inv(A) * b
    return u
end

function SOR(A, b, w=1.5, tol=1e-13)
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
                Sum += A[j, i] * x[i]
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

function make_V(M)
    m = M - 2
    V = zeros(m, m)
    return V
end

function make_u(M)
    m = M - 2
    # u = zeros(m, m)
    # u = reshape(u, m * m)
    u = Array{ComplexF64}(1:m*m)

    randn!(u)  # for now, fill u with random
    u /= sqrt(sum(abs2.(u)))  # normalize u
    return u
end

function initialize_system(M, dt, h)
    V = make_V(M)
    A, B = make_AB(M, h, dt, V)
    u = make_u(M)
    return A, B, u
end

function main()
    M = 40
    h = 1 / (M - 3)
    T = 1
    dt = 0.01
    nT = Int(round(T / dt) + 1)
    A, B, u = initialize_system(M, dt, h)
    save_all = zeros(Float64, nT, (M - 2)^2 + 1)
    # display(A)
    # println(size(A), size(B), size(u))
    save_all[:, 1] = collect(range(0, T, nT))
    save_all[1, 2:end] = abs2.(u)
    # display(save_all[1:10, 1:10])
    for n in 2:nT
        u = advance(u, A, B)
        save_all[n, 2:end] = abs2.(u)
    end
    println(sum(abs2.(u)))
    npzwrite("test.npz", save_all)
    # plot(save_all, M-2)
end

# function plot(A, m)
#     # df = DataFrame(A, :auto)
#     x = collect(range(0, 1, m))
#     plot(contour(
#         z=reshape(A[1, 2:end], m, m),
#         x=x, y=x,
#     ))
# end


main()
# m = 4
# n = m ^ 2
# a = Array(1:n)
# r = 1.4
# A, B = make_AB(r, a, a)
# display(B)
# println()

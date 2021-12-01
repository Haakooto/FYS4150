using LinearAlgebra
import Base.:*

struct SparseMat
    diag::Array
    r::ComplexF64
    m::Int
end

SparseMat(d, r) = SparseMat(d, r, Int(sqrt(length(d))))

function get_at(A::SparseMat, i::Int, j::Int)
    if i <= 0 || j <= 0 || i > A.m^2 || j > A.m^2
        return 0
    else
        if i == j
            return A.diag[i]
        elseif abs(i - j) == 1
            t = min(i, j)
            return (mod(t, A.m) != 0) * A.r
        else
            return (abs(i - j) == A.m) * A.r
        end
    end
end

function over_index(a::Array, i::Int)
    if i <=0 || i > length(a)
        return 0
    else
        return a[i]
    end
end

function to_dense(SM::SparseMat)
    n = SM.m^2
    out = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n
            out[i, j] = get_at(SM, i, j)
        end
    end
    return out
end

function *(A::SparseMat, x::Array)
    b = zeros(ComplexF64, length(x))
    # b = zeros(typeof(x[1]), size(x))
    indices = [-A.m, -1, 0, 1, A.m]
    for i in 1:A.m^2
        for j in indices
            b[i] += get_at(A, i, i + j) * over_index(x, i + j)
        end
    end
    return b
end

function SOR(A::SparseMat, b::Array; initial_guess=0, omega=1.5, tol=1e-13)
    max_iter = 1e5
    n = length(b)
    if initial_guess == 0
        x = similar(b)  #  make vector of same type
    else
        x = initial_guess
    end
    iters = 0  #  count iterations needed
    indices = [-A.m, -1, 1, A.m]
    while iters < max_iter # run until convergence
        iters += 1
        for i in 1:n
            Sigma = 0  # sum(Aij * xj)
            for j in indices
                Sigma += get_at(A, i, j + i) * over_index(x, i + j)
            end
            x[i] = (1 - omega) * x[i] + (omega / A.diag[i]) * (b[i] - Sigma)
        end
        residual = maximum(abs.((A * x) .- b))
        if residual < tol
            break
        end
    end
    return x, iters
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

using LinearAlgebra
import Base.:*

struct SparseMat
    """
    Because the shape of A and B, 
    only need to store the diagonal and off-diag element r 
    """
    diag::Array
    r::ComplexF64
    m::Int
end

SparseMat(d, r) = SparseMat(d, r, Int(sqrt(length(d))))  # Functor ('constructor function')

function Base.getindex(A::SparseMat, i::Int, j::Int)
    """
    For indexing SparseMat
    Returns matrix element at index (i, j)
    Used as normal matrix indexing, A[i, j]
    """
    if i <= 0 || j <= 0 || i > A.m^2 || j > A.m^2
        return 0  # outside of matrix
    else
        if i == j  # diagonal
            return A.diag[i]
        elseif abs(i - j) == 1  # first off-diagonal is r, except every mth element
            t = min(i, j)  # to go down or left
            return (mod(t, A.m) != 0) * A.r
        else  # mth off-diagonal is r
            return (abs(i - j) == A.m) * A.r 
        end
    end
end

function over_index(a::Array, i::Int)
    """
    To not deal with pesky indexing errors, 
    let indexing outside of array return 0
    """
    if i <=0 || i > length(a)
        return 0
    else
        return a[i]
    end
end

function to_dense(SM::SparseMat)
    """
    Convert a SparseMat to dense
    Must be done to print sparse matrix
    Loops over all indexes and fills empty matrix
    """
    n = SM.m^2
    out = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n
            out[i, j] = SM[i, j]
        end
    end
    return out
end

function *(A::SparseMat, x::Array)
    """
    Performs matrix multiplication between SparseMat and Array
    Uses the fact that A is sparse to only perform nessisary calculations
    """
    b = zeros(ComplexF64, length(x))
    indices = [-A.m, -1, 0, 1, A.m]  # non-zero elements of A
    for i in 1:A.m^2
        for j in indices
            b[i] += A[i, i + j] * over_index(x, i + j)
        end
    end
    return b
end

function SOR(A::SparseMat, b::Array; initial_guess=0, omega=1.5, tol=1e-13, max_iter=1e4)
    """
    Performs Successive Over Relaxation to solve matrix problem Ax=b
    Starting with an initial guess of x, iteratively update it until convergence
    for each iteration calculate each new element through the following equation
        x_i^(k + 1) = (1 - w)x_i^k + (w / a_ii)*[b_i - sum from j=1 to j=N, skipping j = i, of a_ij * x_j]
        
    See https://en.wikipedia.org/wiki/Successive_over-relaxation for derivation and discussion of algorithm
    
    This particular implementation uses the fact that A is SparseMat to only loop over j where A is non-zero

    Arguments
        A: SparseMat
            Matrix in equation
        b: Array{ComplexF}
            RHS of equation
        initial_guess: Array{ComplexF}
            Initial guess to start from. By default 0, in which case a zeros is used.
            A good initial array can significantly speed up convergence
        omega: Float
            relaxation factor
        tol: Float
            Convergence criteria
        max_iter: Int
            Maximum number of iterations to try before giving up
    Return:
        x: Array{ComplexF}
            Final solution to Ax=b
        iters: Int
            Number of iterations required to reach convergence
    """
    n = length(b)
    if initial_guess == 0
        x = zeros(ComplexF64, n)  #  make vector of same type
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
                Sigma += A[i, j + i] * over_index(x, i + j)
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
    """
    return a meshgrid of x and y
    """
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

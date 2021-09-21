using LinearAlgebra

function create_A(N, h)
    A = zeros((N, N))
    for i in 1:N
        A[i, i] = 2 * h
        if i != 1
            A[i, i - 1] = -h
            A[i - 1, i] = -h
        end
    end
    return A
end

function create_A(N)
    h = float(N + 1) ^ 2
    a = zeros(N - 1) .- h
    b = zeros(N) .+ 2 .* h
    A = SymTridiagonal(b, a)
    return Array(A)
end

function analytic_solutions(N, M, d, a)
    vals = [d + 2 * a * cos(i * pi / (N + 1)) for i=1:M]
    vecs = [sin(i * j * pi / (N + 1)) for j=1:M, i in 1:N]
    MatrixNormalize!(vecs)
    return vals, vecs
end

function MatrixNormalize!(A)
    n = size(A, 1)
    for i = 1:n
        A[i, :] = normalize(A[i, :])
    end
end

function analytic_solutions(N, M)
    h = float(N + 1) ^ 2
    return analytic_solutions(N, M, 2 * h, -h)
end

function sort_mat_by_vec!(A, b)
    indices = sortperm(b)
    A[:] = A[:, indices]
    b[:] = b[indices]


    for i = 1:length(b)
        if A[1, i] < 0
            A[:, i] .*= -1
        end
    end
end
using LinearAlgebra

ij_to_k(i, j, m) = return (i - 1) * m + j

function make_AB_(r, a, b)
    n = length(a)
    m = Int(sqrt(n))
    A = zeros(n, n)
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
    a = ones(m)
    b = ones(m)
    r = 1im * dt / 2
    for k in 1:m
        a[k] += 4 * r / h^2 + r * v[k_to_ij(k)]
        b[k] -= 4 * r / h^2 + r * v[k_to_ij(k)]
    end
    return make_AB_(r, a, b)
end

m = 4
n = m ^ 2
a = Array(1:n)
r = 1.4
A, B = make_AB(r, a, a)
display(B)
println()

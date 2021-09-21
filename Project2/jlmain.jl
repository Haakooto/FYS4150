using Formatting
using Plotly
include("utils.jl")

#=
    This is a solution of the problem in Julia-lang.
    This was done for fun, and is not a part of the submitted project.
    The code is very similar to the C++ code, but adapted to Julia features.
    Plotting is done in-lang. The produced plots are incomplete, and will stay so.
=#

function max_offdiag_symmetric(A, N)
    max = 0
    k = 0
    l = 0
    for i = 1:N
        for j = i + 1:N
            if abs(A[i, j]) > max
                max = abs(A[i, j])
                k = i
                l = j
            end
        end
    end
    return max, k, l
end

function Rotation!(A, R, k, l, N)
    tau = (A[l, l] - A[k, k]) / 2 / A[k, l]
    if tau > 0
        t = 1 / (sqrt(1 + tau ^ 2) + tau)
    else
        t = - 1 / (sqrt(1 + tau ^ 2) - tau)
    end
    c = 1 / sqrt(1 + t ^ 2)
    s = c * t

    akk = A[k, k]
    all = A[l, l]
    A[k, k] = akk * c ^ 2 - 2 * A[k, l] * c * s + all * s ^ 2
    A[l, l] = all * c ^ 2 + 2 * A[k, l] * c * s + akk * s ^ 2
    A[k, l] = A[l, k] = 0

    for i = 1:N
        if (i != k) && (i != l)
            A[i, k] = A[i, k] * c - A[i, l] * s
            A[i, l] = A[i, l] * c + A[k, i] * s
            A[k, i] = A[i, k]
            A[l, i] = A[i, l]
        end
        rik = R[i, k]
        R[i, k] = rik * c - R[i, l] * s
        R[i, l] = R[i, l] * c + rik * s
    end
end


function Jacobi!(A, R, tol, stop)
    N = size(A, 1)
    max = tol + 1
    iter = 0
    while max > tol
        max, k, l = max_offdiag_symmetric(A, N)
        Rotation!(A, R, k, l, N)
        iter += 1
        if iter == stop
            return false, 0
        end
    end
    return true, iter
end

function test_analyticity()
    N = 6
    tol = 1e-10
    A = create_A(N)

    D = eigvals(A)
    S = eigvecs(A)

    a_vals, a_vecs = analytic_solutions(N, N)
    sort_mat_by_vec!(S, D)
    sort_mat_by_vec!(a_vecs, a_vals)

    @assert maximum(broadcast(abs, (S .- a_vecs))) < tol
    @assert maximum(broadcast(abs, (D .- a_vals))) < tol
end

function tests()
    test_analyticity()

    println("All tests passed")
end

function run_Jacobi(N, maxiter, plotthis=false)
    tol = 1e-10

    A = create_A(N)
    eig_vecs = Matrix{Float64}(I, N, N)


    succ, iters = Jacobi!(A, eig_vecs, tol, maxiter)
    if succ
        printfmtln("Algorithm converged after {:d} iterations", iters)

        eig_vals = diag(A)
        MatrixNormalize!(eig_vecs)
        sort_mat_by_vec!(eig_vecs, eig_vals)

        if plotthis
            MyPlot(eig_vals, eig_vecs, 3)
        end

    else
        printfmtln("Algorithm did not converge. Terminated after {:d} iterations", maxiter)
    end
end

function MyPlot(vals, vecs, K)
    N = length(vals)
    h = 1 / (N + 1)
    x = collect(range(1, N, step=1)) .* h
    avals, avecs = analytic_solutions(N, K)


    line1 = scatter(x=x, y=vecs[:, 1], mode="lines",
                    line=Dict(:color => "firebrick",
                             :width => 4,
                             :dash => "solid"
                             ),
                    markers=true, marker_size=10,
                    name="Computed 1"
                    )

    line2 = scatter(x=x, y=vecs[:, 2], mode="lines",
                    line=Dict(:color => "darkgreen",
                             :width => 4,
                             :dash => "solid"
                             ),
                    markers=true, marker_size=10,
                    name="Computed 2"
                    )

    line3 = scatter(x=x, y=vecs[:, 3], mode="lines",
                    line=Dict(:color => "marine",
                             :width => 4,
                             :dash => "solid"
                             ),
                    markers=true, marker_size=10,
                    name="Computed 3"
                    )

    line4 = scatter(x=x, y=avecs[1, :], mode="lines",
                    line=Dict(:color => "firebrick",
                             :width => 4,
                             :dash => "dash"
                             ),
                    markers=true, marker_size=10,
                    name="Analytical 1"
                    )

    line5 = scatter(x=x, y=avecs[2, :], mode="lines",
                    line=Dict(:color => "darkgreen",
                             :width => 4,
                             :dash => "dash"
                             ),
                    markers=true, marker_size=10,
                    name="Analytical 2"
                    )

    line6 = scatter(x=x, y=avecs[3, :], mode="lines",
                    line=Dict(:color => "marine",
                             :width => 4,
                             :dash => "dash"
                             ),
                    markers=true, marker_size=10,
                    name="Analytical 3"
                    )


    lines = [line1, line2, line3, line4, line5, line6]

    layout = Layout(title=format("First {:d} solutions the the differential equation for N = {:d}", K, N),
                    width=1800, height=1200,
                    legend=attr(font_size=10)
                    )
    p = plot(lines, layout)
    savefig(p, format("results/test_{:d}.svg", N))
end

function main()
    tests()

    maxiter = 100000
    # N = 100
    # for i = 5:N
    #     run_Jacobi(i, maxiter)
    # end

    run_Jacobi(10, maxiter, true)
    run_Jacobi(100, maxiter, true)
end

main()
println()
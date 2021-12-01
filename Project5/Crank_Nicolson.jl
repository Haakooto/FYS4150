using LinearAlgebra
using Random
using SparseArrays
# using Pkg
# Pkg.add("NPZ")
# Pkg.add("ProgressMeter")
using NPZ
using ProgressMeter

include("sparse_matrices.jl")

function make_AB(h, dt, V)
    """
    Sets up the sparse matrices A and B,
    with diagonals a and b, r as non-diagonal elements

    Arguments:
        h: Float
            Step length in spatial discretisation
        dt: Float
            Step length in temporal discretisation
        V: Matrix{Float}
            Matrix of potential values in each point
    Returns 
        A, B: SparseMat
            Sparse matrices, as defined in assignment
    """
    m = Int(1 / h) - 2
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

function make_V(h, v0, slits)
    """
    Sets up potential matrix V
    Must still be written, currently no potential is used

    Arguments:
        h: Float
            Step length in spatial discretisation
        v0: Float
            Value of potential where not zero
        slits: Int
            Number of slits in the slit experiment
    Returns:
        V: Matrix{Float}
            Potential matrix
    """
    m = Int(1 / h) - 2
    V = zeros(m, m)
    return V
end

function Gauss_wave_packet(h, xc, yc, px, py, sx, sy)
    """
    Initializes the system with a 2D gaussian wave packet

    Arguments:
        h: Float
            Step length in spatial discretisation 
        xc, yc: Float
          Centre of wave packet in x and y coordinate
        px, py: Float
            wave packet momentum in x and y coordinate
        sx, sy: Float
            Initial width of wave packet in x and y coordinate
    Returns
        u:  Array{ComplexF}
            Normalised, flattened array of initial wave packet
    """
    x = [i * h for i in 0:Int(1/h)-1]
    X, Y = meshgrid(x, x)

    t1x = -((X .- xc).^2) ./ (2 * sx^2)
    t1y = -((Y .- yc).^2) ./ (2 * sy^2)
    t2x = 1im * px .* (X .- xc)
    t2y = 1im * py .* (Y .- yc)

    u = exp.(t1x .+ t1y .+ t2x .+ t2y)

    u /= sqrt(sum(abs2.(u)))  # normalize u
    u = vec(u[2:end-1, 2:end-1])  # Flatten, chop off boundaries
    return u
end

function initialize_system(dt, h, v0, slits, P)
    """
    Wrapper, initializes the system

    Arguments
        dt: Float
            Step length in temporal discretisation
        h: Float
            Step length in spatial discretisation
        v0: Float
            Value of potential where not zero
        slits: Int
            number of slits in the slit experiment
        P: Array{Float}
            Container for parameters to gaussian wave packet
    Returns:
        A, B: SparseMat
            A and B as specified in assigmnet
        u: Array{ComplexF}
            Normalized, flattened complex array of initial wave packet
    """
    V = make_V(h, v0, slits)
    u = Gauss_wave_packet(h, P...)
    A, B = make_AB(h, dt, V)
    return A, B, u
end

function main()
    # Discretisation free parameters
    h = 0.005
    dt = 2.5e-5
    T = 0.002

    # System free parameters
    xc = 0.25
    sx = 0.005
    px = 200
    yc = 0.5
    sy = 0.05
    py = 0
    v0 = 0
    slits = 2

    # Dependent parameters
    gauss_params = [xc, yc, px, py, sx, sy]
    M = Int(1 / h)
    nT = Int(T / dt)

    A, B, u = initialize_system(dt, h, v0, slits, gauss_params)

    # Probs contains the probability at all points, at all times, as well as time array
    # Must be changed, depending on what we want show
    probs = zeros(Float64, nT, (M - 2)^2 + 1)
    probs[:, 1] = collect(range(0, T, nT))
    probs[1, 2:end] = abs2.(u)

    pbar = Progress(nT; showspeed=true, barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=10)
    for time in 2:nT
        b = B * u  # Problem 3 step 1
        u, i = SOR(A, b, initial_guess=u, omega=1, tol=1e-13)  # Problem 3 step 2
        P = abs2.(u)  # record
        probs[time, 2:end] = P 
        ProgressMeter.next!(pbar; showvalues=[(:time,time), (:i,i), (:prob,sum(P))])
    end
    println("\n\n\n\n")
    println(sum(abs2.(u)))  # final total probability
    npzwrite("evolve3.npz", probs)
end

main()

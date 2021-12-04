function make_AB(h, dt, V)
    #=
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
    =#
    m = Int(1 / h) - 2
    a = ones(ComplexF64, m*m)
    b = ones(ComplexF64, m*m)
    r = 1im * dt / 2 / h^2
    for k in 1:m*m
        ch = 4r + r * V[k] * h^2
        a[k] += ch
        b[k] -= ch
    end
    return SparseMat(a, -r), SparseMat(b, r)
end

function make_V(h, v0=0, slits=2; thickness=0.02, position=0.5, aperture=0.05, seperator=0.05 )
    #=
    Sets up potential matrix V
    Must still be written, currently no potential is used
    Use parameters to calculate

    It is possible the matrix has to be transposed,
    as Julia is column-major,
    but I have put no brain power into contemplating
    this further than the potential nessescity

    Arguments:
        h: Float
           Step length in spatial discretisation
        v0: Float
            Value of potential where not zero
        slits: Int
               Number of slits in the slit experiment
        thickness: Float
                   width of wall
        position: Float
                  position along x-axis of wall
        aperture: Float
                  width along y-axis of slit
        seperator: Float
                   width of wall between slits
    Returns:
        V: Matrix{Float}
           Potential matrix
    =#
    m = Int(1 / h) - 2
    V = zeros(m, m)
    return V
end

function Gauss_wave_packet(h, xc, yc, px, py, sx, sy)
    #=
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
    =#
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
    #=
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
    =#
    V = make_V(h, v0, slits)
    u = Gauss_wave_packet(h, P...)
    A, B = make_AB(h, dt, V)
    return A, B, u
end

function simulate(args, name, realimag=false)
    #=
    Main function in program, sets up system and evolves it through time
    Takes 11 arguments in a list, and a name to save file as

    Arguments:
        args: Array
            Contains the following parameters:
                1: h, step length. Never change from 0.005
                2: dt, time step. Never change from 2.5 * 10^-5
                3: T, total simulation time. Set at own perogative
                4: xc, wave packer parameter. Never change from 0.25
                5: yc, wave packer parameter. Never change from 0.5
                6: px, wave packer parameter. Never change from 200
                7: py, wave packer parameter. Never change from 0
                8: sx, wave packer parameter. Never change from 0.05
                9: sy, wave packer parameter. Set at own perogative
                10: v0, potential wall constant. Set at own perogative
                11: slits, number of slits in wall. Set at own perogative
        name: String
              name to save datafile as
        realimag: Bool
                  to save real and imag part seperately, in addition to popbability
    Returns
        saves datafile with name name
    =#
    # unpack parameters
    h = args[1]
    dt = args[2]
    T = args[3]
    gauss_params = args[4:9]
    v0 = args[10]
    slits = args[11]

    file = "npz/" * name * ".npz"

    A, B, u = initialize_system(dt, h, v0, slits, gauss_params)
    println("Finished initialized system")

    M = Int(1 / h)
    points = (M - 2) ^ 2
    total = points
    if realimag
        total *= 3
    t = range(0, T, step=dt)

    # set up storage mat, place in initial conditions
    storage = zeros(Float32, length(t), total + 2)
    storage[:, 1] = t
    storage[1, 2] = sum(abs2.(u))
    storage[1, 3:points] = abs2.(u)

    println("Starting simulation...")
    pbar = Progress(nT; showspeed=true, barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=10)
    for time in 2:length(t)
        b = B * u  # Problem 3 part 1
        u, i = SOR(A, b, initial_guess=u, omega=0.9)  # problem 3 part 2
        P = abs2.(u)  # Born rule

        storage[time, 2] = sum(P)  # total probability at time
        storage[time, 3:points+2] = P  # probability as each points
        if realimag
            storage[time, points+1:2*points] = real.(u)
            storage[time, 2*points+1:3*points] = imag.(u)
        ProgressMeter.next!(pbar; showvalues=[(:time, time * dt), (:SOR_iters, i), (:Total_probability, sum(P))])
    end
    println("\n\n\n\n")
    npzwrite(file, storage)
    println("Saved data to file: ", file)
end

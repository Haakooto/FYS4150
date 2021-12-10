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
    m = Int(1 / h) - 1
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


function position_to_index(position, h=0.005)
    #=
    Translates a spatial position to the corresponding index for the internal points
    Arguments:
        position: Float
            position of whatever is going on there
        h: Float
            Step length in spatial discretisation
    =#
    index = round(position/h)
    return Int(index)
end


function make_V(h, v0=0, slits=1; thickness=0.02, position=0.5, aperture=0.05, separator=0.05)
    #=
    Sets up potential matrix V

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
    m = Int(1 / h) - 1
    V = zeros(m, m)

    x_start = position_to_index(position - thickness/2, h)
    x_stop = position_to_index(position + thickness/2, h)

    if slits == 1
        y_start = position_to_index(0.5 - aperture/2, h)
        y_stop = position_to_index(0.5 + aperture/2, h)
        V[1:y_start-1, x_start:x_stop] .= v0
        V[y_stop:end, x_start:x_stop] .= v0

    elseif slits == 2
        y_start1 = position_to_index(0.5 - separator/2 - aperture, h)
        y_stop1 = position_to_index(0.5 - separator/2, h)

        y_start2 = y_start1 + position_to_index(separator + aperture, h)
        y_stop2 = y_stop1 + position_to_index(separator + aperture, h)

        V[1 : y_start1-1, x_start:x_stop] .= v0
        V[y_stop1 : y_start2-1, x_start:x_stop] .= v0
        V[y_stop2 : end, x_start:x_stop] .= v0

    else slits == 3
        y_start1 = position_to_index(0.5 - aperture*3/2 - separator, h)
        y_stop1 = position_to_index(0.5 - aperture*1/2 - separator, h)

        y_start2 = y_start1 + position_to_index(separator + aperture, h)
        y_stop2 = y_stop1 + position_to_index(separator + aperture, h)

        y_start3 = y_start2 + position_to_index(separator + aperture, h)
        y_stop3 = y_stop2 + position_to_index(separator + aperture, h)

        V[1:y_start1-1, x_start:x_stop] .= v0
        V[y_stop1:y_start2-1, x_start:x_stop] .= v0
        V[y_stop2:y_start3-1, x_start:x_stop] .= v0
        V[y_stop3:end, x_start:x_stop] .= v0

    end

    file = "npz/potential_" * string(Int(slits)) * "_slits.npz"   #write the potential to file to check what it looks like
    npzwrite(file, V)

    return transpose(V)

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
    x = [i * h for i in 0:Int(1/h)]
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
        File contains matrix of size (length(t), points + 2),
        where points = (M-2)^2.
        The first column is the time-values, the second column is the total probability in domain at that time
        The rest of each row is the (M-2)^2 probabilities corresponding to each point. Must be reshaped as (M-2 X M-2) for plotting

        For problem 8, two additional save files are created, one for real part, one for imag part
        These do NOT contain time nor probability, so size is (length(t), points).
        Can simply be reshaped and plotted for each time as is.
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
    A = to_sparse(A)  # to speed-up code, use SparseArrays, not our SparseMat
    println("Finished initialized system")

    M = Int(1 / h) + 1
    points = (M - 2) ^ 2
    t = range(0, T, step=dt)

    # set up storage mat, place in initial conditions
    storage = zeros(Float32, length(t), points + 2)
    storage[:, 1] = t
    storage[1, 2] = (1-sum(abs2.(u)))
    storage[1, 3:end] = abs2.(u)

    # @assert abs2.(u) == (real.(u) - 1im * imag.(u)) .* (real.(u) + 1im * imag.(u))

    # For problem 8, also save real and imag parts, seperatrely
    if realimag
        store_real = zeros(Float64, length(t), points)
        store_imag = zeros(Float64, length(t), points)
        store_real[1, :] = real.(u)
        store_imag[1, :] = imag.(u)
    end

    println("Starting simulation...")
    pbar = Progress(length(t); showspeed=true, barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=10)
    i = 30  # must be sufficiently high to give good results, but low enough to not take forever
    for time in 2:length(t)
        b = B * u  # Problem 3 part 1
        # u, i = SOR(A, b, initial_guess=u, omega=0.9)  # problem 3 part 2
        ssor!(u, A, b, 0.9, maxiter=i)  # this is faster than our sor-solver
        P = abs2.(u)  # Born rule

        storage[time, 2] = 1-sum(P)  # deviation from 1 of total probability at time
        storage[time, 3:end] = P  # probability as each points

        if realimag
            store_real[time, 1:end] = real.(u)
            store_imag[time, 1:end] = imag.(u)
        end
        ProgressMeter.next!(pbar; showvalues=[(:time, time * dt), (:SOR_iters, i), (:Total_probability, sum(P))])
    end
    println("\n\n\n\n")
    npzwrite(file, storage)
    println("Saved data to file: ", file)

    if realimag
        file_r = "npz/" * name * "_real" * ".npz"
        file_i = "npz/" * name * "_imag" * ".npz"
        npzwrite(file_r, store_real)
        npzwrite(file_i, store_imag)
        println("Saved real and imaginary parts to seperate files.")
    end
end

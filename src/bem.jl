using .Propeller

# using ForwardDiff
# noNaNs(x::Real) = true
# noNaNs(x::ForwardDiff.Dual) = !any(isnan, ForwardDiff.partials(x))

"""
    _hbem(
        rotor::Rotor,
        flight_condition::FlightCondition,
        rpm,
        vRel_R,
        T;
        tol_newton = 1.0e-9,
        maxiter_newton = 1000,
    )

# Returns 
- `error_CT`: error in thrust coefficient
- `T_BET`: thrust
- `Q_BET`: torque
- `Np`: normal loads
- `Tp`: tangential loads
- `α`: angle of attack
- `cl`: lift coefficient
- `cd`: drag coefficient
- `CT_BET`: thrust coefficient
- `U`: inflow velocity
- `Φ`: inflow angle
- `reynolds`: Reynolds number
"""
function _hbem(
    rotor::Rotor,
    flight_condition::FlightCondition,
    rpm,
    vRel_R,
    T;
    tol_newton=1.0e-9,
    maxiter_newton=1000,
)

    # ---------------------------------
    # Load Rotor Object
    # ---------------------------------
    r_tip = rotor.r_tip
    n_blades = rotor.n_blades
    n_radial = rotor.n_radial
    n_azimuthal = rotor.n_azimuthal

    Vx, Vy, Vz = vRel_R
    # Vinf = [Vx, Vy, Vz] # relative incoming velocity

    n = rpm / 60.0
    Ω = 2π * n # rad per second

    Vtip = Ω * r_tip # tip velocity

    # rotor disk area 
    A = π * r_tip^2


    # ---------------------------------
    # Calculate μ and λ_c
    # ---------------------------------

    # println("--------------------")
    # println(Vx, " ", Vy, " ", Vtip)
    # println("--------------------")

    # μ_z = -Vz / Vtip # advance ratio, perpendicular to rotor plane
    μ = (Vx^2 + Vy^2)^0.5 / Vtip # advance ratio, parallel to rotor plane
    # println("μ: ", μ)
    λ_c = -Vz / Vtip # rotor climb ratio, perpendicular to rotor plane

    # Copter angle of attack
    # AoA = atan(Vinf[3], norm_cs_safe(Vinf))

    # ---------------------------------
    # CT_momentum
    # ---------------------------------

    # Thrust coefficient
    CT_momentum = T / (flight_condition.ρ * A * Vtip^2)


    # ---------------------------------
    # Finding λ_0 (Newton-Raphson)
    # ---------------------------------

    # Initial guess for lambda (inflow ratio)
    λ_h = (0.5 * CT_momentum)^0.5 # inflow ratio in hover
    λ_0 = λ_h

    # println("CT_momentum: ", CT_momentum)


    # Iteratively find the inflow ratio
    f(λ_0) = λ_0 - λ_c - 0.5 * CT_momentum * (μ^2 + λ_0^2)^(-0.5)
    fp(λ_0) = 1.0 + 0.5 * CT_momentum * (μ^2 + λ_0^2)^(-1.5) * λ_0

    err = 1.0
    iter = 0
    while err > tol_newton && iter < maxiter_newton
        f = λ_0 - λ_c - 0.5 * CT_momentum * (μ^2 + λ_0^2)^(-0.5)
        fp = 1 + 0.5 * CT_momentum * (μ^2 + λ_0^2)^(-1.5) * λ_0
        λ_new = λ_0 - f / fp
        err = (λ_new - λ_0)^2 / λ_h^2
        λ_0 = λ_new
        iter += 1
    end





    # ---------------------------------
    # Use λ_0 to Calculate CT_BET 
    # ---------------------------------

    # Advance Ratio
    # μ_par = μ # parallel to the rotor
    # μ_perp = λ_c # prependicular to the rotor (lambda_c)

    # ---- Inflow Models ----
    # Empirically/vortex-theory-based assumptions about the distribution 
    # of the induced velocity

    # Wake skew angle
    # χ = atan(μ_par / (μ_perp + λ_0))
    χ = atan(μ, λ_0)

    # Linear inflow models
    if rotor.inflow_model == :pittpeters
        k_x = (15π / 23) * tan(χ / 2)
        k_y = 0
    elseif rotor.inflow_model == :drees
        k_x = 4 / 3 * (1 - cos(χ) - 1.8μ^2) / sin(χ)
        k_y = -2μ
    elseif rotor.inflow_model == :payne
        k_x = 4 / 3 * (μ / λ_0 / (1.2 + μ / λ_0))
        k_y = 0
    end


    # Fixing ψ given a Vy velocity, ψ is defined when x is aligned with the free stream
    # ψ = hcat(rotor.psi)' .- atan(-Vy, Vx)
    ψ = repeat(rotor.psi', rotor.n_radial, 1) # n_radial x n_azimuthal
    r_nondim = repeat(hcat(rotor.r_nondim), 1, n_azimuthal) # n_radial x n_azimuthal

    # Constructing inflow ratio based on a linear model
    λ = λ_0 .* (1 .+ k_x .* r_nondim .* cos.(ψ) + k_y .* r_nondim .* sin.(ψ))


    # Integrate to find roll and pitch for a rotor
    s = repeat(rotor.s, 1, n_azimuthal)

    β = 0 # flapping angle
    β_dot = 0

    # The incoming velocity seen by the blade (nondimensional!)
    u_p = λ + μ * β * cos.(ψ) + r_nondim * β_dot / Ω

    u_t = r_nondim + μ * sin.(ψ)

    u = (u_t .^ 2 + u_p .^ 2) .^ 0.5 # resultant velocity
    U = u .* Vtip


    Φ = atan.(u_p, u_t) # relative inflow angle (induced angle of attack)
    θ = repeat(hcat(rotor.theta), 1, n_azimuthal)
    chord = repeat(hcat(rotor.chord), 1, n_azimuthal)
    reynolds =
        flight_condition.ρ .* (u .* Vtip) .* chord ./ flight_condition.μ

    if rotor.use_polar_model
        slope_cl = repeat(hcat(rotor.cl_a_slope), 1, n_azimuthal)
        α_eff_0 = θ - Φ
        cl1 = slope_cl .* α_eff_0
        α = (θ .+ rotor.alpha_cl0 .- Φ)
        alp = α .- rotor.alpha_cl0
        sigma =
            (
                1 .+ exp.(-rotor.stall_transition_rate .* (alp .- rotor.alpha_cutoff)) .+
                exp.(rotor.stall_transition_rate .* (alp .+ rotor.alpha_cutoff))
            ) ./ (
                (1 .+ exp.(-rotor.stall_transition_rate .* (alp .- rotor.alpha_cutoff))) .*
                (1 .+ exp.(rotor.stall_transition_rate .* (alp .+ rotor.alpha_cutoff)))
            )
        cl =
            (1 .- sigma) .* cl1 .+ sigma .* (2 .* sign.(alp) .* sin.(alp) .^ 2 .* cos.(alp))
        cd = 0.0 .* similar(cl)
        _ = α[1, 1]

    else
        # Evaluate Lift Coefficient
        α = θ - Φ

        polar_interpolator = nothing

        if rotor.online_correction3D
            polar_data_3D = correction3D(rotor.polar_data, rotor.r, rotor.chord, mean(U, dims=2)[:], Ω)
            polar_interpolator = polar_interpolators(polar_data_3D)
        else
            polar_interpolator = rotor.polar_interpolator
        end

        cl = similar(u)
        cd = similar(u)
        for i = 1:n_radial
            for j = 1:n_azimuthal
                cl[i, j] = polar_interpolator[1](α[i, j], reynolds[i, j], s[i, j])
                cd[i, j] = polar_interpolator[2](α[i, j], reynolds[i, j], s[i, j])
            end
        end
    end

    # Velocity, dynamic pressure 
    q = 0.5 .* flight_condition.ρ .* U .^ 2

    # Mach correction (Glauert)
    if rotor.mach_correction
        mach = abs.(U) ./ flight_condition.c

        if any(mach .> 1.0)
            error(
                "Local Mach number is greater than 1.0. Check your rotor and
                operating conditions or disable Mach correction 
                (`rotor.mach_correction=false`)."
            )
        end

        cl = cl ./ (1 .- mach .^ 2) .^ 0.5
    end

    # Hub/tip loss correction
    Ftip = similar(u)
    Fhub = similar(u)
    if rotor.hubtip_correction == :prandtl
        for i = 1:n_radial
            for j = 1:n_azimuthal
                Ftip[i, j], Fhub[i, j] = hubtip_correction_prandtl(
                    r_tip,
                    rotor.r_hub,
                    rotor.r[i],
                    Φ[i, j],
                    n_blades,
                )
            end
        end
    elseif rotor.hubtip_correction == :none
        Ftip .= 1
        Fhub .= 1
    end

    # Normal and Tangential Force Coefficients
    cn = Ftip .* Fhub .* cl .* cos.(Φ) .- cd .* sin.(Φ) # normal force coefficient
    ct = Ftip .* Fhub .* cl .* sin.(Φ) .+ cd .* cos.(Φ) # tangential force coefficient

    # Normal and Tangential Forces / Blade
    Np = cn .* q .* chord
    Tp = ct .* q .* chord

    # Average loads over azimuthal direction to approx. the thrust and torque
    #   Note: In reality, the thrust/torque varies as the blades have different 
    #   azimuthal positions in time
    Np_avg_azimuthal = Trapz.trapz(rotor.psi, Np, Val(2)) ./ 2π
    Tp_avg_azimuthal = Trapz.trapz(rotor.psi, Tp .* r_nondim .* r_tip, Val(2)) ./ 2π
    T_BET = n_blades .* Trapz.trapz(rotor.r, Np_avg_azimuthal)
    Q_BET = n_blades .* Trapz.trapz(rotor.r, Tp_avg_azimuthal)

    # Thrust coefficient BET
    CT_BET = T_BET ./ (flight_condition.ρ * A * Vtip^2)

    # ---------------------------------
    # Calculate Error: 
    #   |CT_BET - CT_Momentum|
    # ---------------------------------
    error_CT = abs(CT_BET - CT_momentum)

    return error_CT, T_BET, Q_BET, Np, Tp, α, cl, cd, CT_BET, U, Φ, reynolds
end


"""
thrust2rpm_hbem(T, vRel_R, rotor; lb=100.0, ub=16000.0)

Computes the RPM based on the specified thrust using the hybrid BEM developed
by Davoudi [1].

# Arguments
- `T`: thrust
- `vRel_R`: relative velocity in {R}-coordinates
- `rotor::Rotor`: rotor struct
- `lb`: lower bound of thrust for optimizer
- `ub`: lower bound of thrust for optimizer

# References
[1] B. Davoudi, “A Hybrid Blade Element Momentum Model for Flight Simulation
    of Rotary Wing Unmanned Aerial Vehicles,” AIAA Paper, 2019.
"""
# function thrust2rpm(
#     rotor::Rotor,
#     flight_condition::FlightCondition,
#     T,
#     vRel_R;
#     lb_rpm=100.0,
#     ub_rpm=16000.0,
#     online_correction_3D=false,
#     mach_correction=false,
# )
#     obj = rpm -> _hbem(
#         rpm, flight_condition, vRel_R, T, rotor, online_correction_3D=online_correction_3D, mach_correction=mach_correction)[1]
#     res = Optim.optimize(obj, lb_rpm, ub_rpm, Optim.Brent()) # Brent apparently has guaranteed convergence
#     return res.minimizer
# end



"""
    rpm2thrust(rotor::Rotor, rpm, vRel_R; flight_condition::FlightCondition=FlightCondition(0.0), lb_T=0.0, ub_T=20.0, method_optim=:brent, trace_optim=false, tol_optim=1e-5, maxiter_optim=1000, tol_newton=1.0e-9, maxiter_newton=1000q

Computes the thrust based on the specified RPM using the hybrid BEM developed
by Davoudi [1].

# Arguments
- `rotor::Rotor`: rotor struct
- `rpm`: RPM
- `vRel_R`: relative velocity in {R}-coordinates
- `flight_condition::FlightCondition`: flight condition
- `lb_T`: lower bound of thrust for optimizer
- `ub_T`: lower bound of thrust for optimizer
- `method_optim`: optimization method (`:brent`, `:golden`)
- `trace_optim`: show trace of optimizer
- `tol_optim`: tolerance for optimizer
- `maxiter_optim`: maximum number of iterations for optimizer
- `tol_newton`: tolerance for Newton-Raphson
- `maxiter_newton`: maximum number of iterations for Newton-Raphson

# References
[1] B. Davoudi, “A Hybrid Blade Element Momentum Model for Flight Simulation
    of Rotary Wing Unmanned Aerial Vehicles,” AIAA Paper, 2019.
"""
function rpm2thrust(
    rotor::Rotor,
    rpm,
    vRel_R;
    flight_condition::FlightCondition=FlightCondition(0.0),
    lb_T=0.0,
    ub_T=20.0,
    method_optim=:brent,
    trace_optim=false,
    tol_optim=1e-5,
    maxiter_optim=1000,
    tol_newton=1.0e-9,
    maxiter_newton=1000,
)

    obj = T -> _hbem(
        rotor, flight_condition, rpm, vRel_R, T; tol_newton=tol_newton, maxiter_newton=maxiter_newton
    )[1]


    x = NaN
    # INFO: the optimization algorithm allocates most of the memory here...
    if method_optim == :brent
        res_optim = Optim.optimize(obj, lb_T, ub_T, Optim.Brent(); iterations=maxiter_optim, show_trace=trace_optim)
        # https://github.com/JuliaNLSolvers/Optim.jl/blob/master/src/univariate/solvers/brent.jl
        x = res_optim.minimizer
    elseif method_optim == :golden
        res_optim = Optim.optimize(obj, lb_T, ub_T, Optim.GoldenSection(); iterations=maxiter_optim, show_trace=trace_optim)
        x = res_optim.minimizer
    end

    # Calculate all values
    err, _T, Q, Np, Tp, α, cl, cd, CT, U, phi, re =
        _hbem(
            rotor, flight_condition, rpm, vRel_R, x; tol_newton=tol_newton, maxiter_newton=maxiter_newton
        )


    # Check convergence

    # println(err)

    if err > tol_optim
        converged = false
    else
        converged = true
    end

    # Nondimensional quantities
    n = rpm / 60.0
    D = 2.0 * rotor.r_tip
    CT = _T / (flight_condition.ρ * n^2 * D^4)
    CQ = Q / (flight_condition.ρ * n^2 * D^5)
    # Ω = 2π * n
    # Pin = Ω * Q

    return _T, Q, converged, Np, Tp, α, cl, cd, CT, CQ, U, phi, re

end

"""
    struct AeroResult

A struct to hold results of aerodynamic calculations.

# Fields

Vectors: n_τ
- `T::Vector{Float64}`: Thrust.
- `Q::Vector{Float64}`: Torque.

- `CT::Vector{Float64}`: Thrust coefficients.
- `CQ::Vector{Float64}`: Torque coefficients.
- `η::Vector{Float64}`: Propulsive efficiency.

- `converged::Vector{Bool}`: Convergence status.

Arrays: n_τ x n_radial x n_azimuthal
- `Np::Array{Float64}`: Normal loads, N'' [N/m^2].
- `Tp::Array{Float64}`: Tangential loads, T'' [N/m^2].
- `α::Array{Float64}`: Angle of attack.
- `cl::Array{Float64}`: Lift coefficient.
- `cd::Array{Float64}`: Drag coefficient.
- `re::Array{Float64}`: Reynolds numbers.
- `V::Array{Float64}`: Inflow velocity.	
"""
struct AeroResult
    T::Vector{Float64}
    Q::Vector{Float64}

    CT::Vector{Float64}
    CQ::Vector{Float64}
    # η::Vector{Float64}

    converged::Vector{Bool}

    Np::Array{Float64}
    Tp::Array{Float64}
    α::Array{Float64}
    cl::Array{Float64}
    cd::Array{Float64}
    re::Array{Float64}
    V::Array{Float64}
end



"""
    simulate_rpm2thrust(
        rotor::Rotor,
        rpm::Vector,
        vRel_R::Matrix;
        n_threads=nothing,
        lb_T=0.0,
        ub_T=10.0,
        parallel::Bool=true,
        flight_condition::FlightCondition=FlightCondition(0.0),
        tol_optim=1e-5,
        maxiter_optim::Int=500,
        trace_optim::Bool=false,
        method_optim::Symbol=:brent,
        tol_newton=1.0e-9,
        maxiter_newton::Int=1000,
    )

Computes the thrust based on the specified RPM using the hybrid BEM developed
by Davoudi [1].

# Arguments
- `rotor::Rotor`: rotor struct
- `rpm::Vector`: RPM
- `vRel_R::Matrix`: relative velocity in {R}-coordinates
- `lb`: lower bound of thrust for optimizer
- `ub`: lower bound of thrust for optimizer
- `flight_condition::FlightCondition`: flight condition
- `tol_optim`: tolerance for optimizer
- `maxiter_optim`: maximum number of iterations for optimizer
- `trace_optim::Bool`: show trace of optimizer
- `method_optim::Symbol`: optimization method (`:brent`, `:golden`)
- `tol_newton`: tolerance for Newton-Raphson
- `maxiter_newton`: maximum number of iterations for Newton-Raphson
"""
function simulate_rpm2thrust(
    rotor::Rotor,
    rpm::Vector,
    vRel_R::Matrix;
    n_threads=nothing,
    lb_T=0.0,
    ub_T=10.0,
    parallel::Bool=true,
    flight_condition::FlightCondition=FlightCondition(0.0),
    tol_optim=1e-5,
    maxiter_optim::Int=500,
    trace_optim::Bool=false,
    method_optim::Symbol=:brent,
    tol_newton=1.0e-9,
    maxiter_newton::Int=1000,
)::AeroResult


    n_τ = length(rpm)
    n_radial = rotor.n_radial
    n_azimuthal = rotor.n_azimuthal

    p = Progress(n_τ, dt=1, barlen=50, color=:yellow)

    T = zeros(n_τ)
    converged = zeros(n_τ)
    Np = zeros((n_τ, n_radial, n_azimuthal))
    Tp = zeros((n_τ, n_radial, n_azimuthal))
    α = zeros((n_τ, n_radial, n_azimuthal))
    cl = zeros((n_τ, n_radial, n_azimuthal))
    cd = zeros((n_τ, n_radial, n_azimuthal))
    phi = zeros((n_τ, n_radial, n_azimuthal))
    Q = zeros(n_τ)
    CT = zeros(n_τ)
    CQ = zeros(n_τ)
    u = zeros((n_τ, n_radial, n_azimuthal))
    re = zeros((n_τ, n_radial, n_azimuthal))


    # Check if all rpm and vRel_R are the same  
    same = false
    if all(rpm .≈ rpm[1]) && all(vRel_R[:, 1] .≈ vRel_R[1, 1]) && all(vRel_R[:, 2] .≈ vRel_R[1, 2]) && all(vRel_R[:, 3] .≈ vRel_R[1, 3])
        println("INFO: All RPM and vRel_R are approx. the same.")
        same = true
    end

    if same
        i = 1
        T[i],
        Q[i],
        converged[i],
        Np[i, :, :],
        Tp[i, :, :],
        α[i, :, :],
        cl[i, :, :],
        cd[i, :, :],
        CT[i],
        CQ[i],
        u[i, :, :],
        phi[i, :, :],
        re[i, :, :] = rpm2thrust(
            rotor,
            rpm[1],
            vRel_R[1, :];
            flight_condition=flight_condition,
            lb_T=lb_T,
            ub_T=ub_T,
            tol_optim=tol_optim,
            maxiter_optim=maxiter_optim,
            trace_optim=trace_optim,
            method_optim=method_optim,
            tol_newton=tol_newton,
            maxiter_newton=maxiter_newton,
        )
        T .= T[1]
        Q .= Q[1]
        CT .= CT[1]
        CQ .= CQ[1]
        converged .= converged[1]

        for i = 1:n_τ
            Np[i, :, :] = Np[1, :, :]
            Tp[i, :, :] = Tp[1, :, :]
            α[i, :, :] = α[1, :, :]
            cl[i, :, :] = cl[1, :, :]
            cd[i, :, :] = cd[1, :, :]

            u[i, :, :] = u[1, :, :]
            phi[i, :, :] = phi[1, :, :]
            re[i, :, :] = re[1, :, :]
        end
    else

        if parallel
            function process_chunk(start, stop)
                for i = start:stop
                    T[i],
                    Q[i],
                    converged[i],
                    Np[i, :, :],
                    Tp[i, :, :],
                    α[i, :, :],
                    cl[i, :, :],
                    cd[i, :, :],
                    CT[i],
                    CQ[i],
                    u[i, :, :],
                    phi[i, :, :],
                    re[i, :, :] = rpm2thrust(
                        rotor,
                        rpm[i],
                        vRel_R[i, :];
                        flight_condition=flight_condition,
                        lb_T=lb_T,
                        ub_T=ub_T,
                        tol_optim=tol_optim,
                        maxiter_optim=maxiter_optim,
                        trace_optim=trace_optim,
                        method_optim=method_optim,
                        tol_newton=tol_newton,
                        maxiter_newton=maxiter_newton,
                    )
                    next!(p)
                end
            end

            if n_threads === nothing
                num_threads = nthreads()
            else
                num_threads = 1
            end
            println(f"Threads: {num_threads}")

            # Calculate chunk size per thread
            chunk_size = div(n_τ, num_threads)
            rest = mod(n_τ, num_threads)
            threads = []

            # Split the work among threads
            for i = 1:num_threads
                start_idx = (i - 1) * chunk_size + 1
                stop_idx = min(i * chunk_size, n_τ)
                thread = @spawn process_chunk(start_idx, stop_idx)
                push!(threads, thread)
            end

            # Wait for all threads to finish
            for thread in threads
                wait(thread)
            end

            i_serial_start = n_τ - rest + 1
            if rest > 0
                println("Running the remaining $rest job in serial mode...")
            end
        else
            i_serial_start = 1
            println("Serial computing...")
        end

        # Run the rest of the jobs
        for i = i_serial_start:n_τ

            println("\n---- $i: $(rpm[i]), $(vRel_R[i, :]) -----")

            T[i],
            Q[i],
            converged[i],
            Np[i, :, :],
            Tp[i, :, :],
            α[i, :, :],
            cl[i, :, :],
            cd[i, :, :],
            CT[i],
            CQ[i],
            u[i, :, :],
            phi[i, :, :],
            re[i, :, :] = rpm2thrust(
                rotor,
                rpm[i],
                vRel_R[i, :];
                flight_condition=flight_condition,
                lb_T=lb_T,
                ub_T=ub_T,
                tol_optim=tol_optim,
                maxiter_optim=maxiter_optim,
                trace_optim=trace_optim,
                method_optim=method_optim,
                tol_newton=tol_newton,
                maxiter_newton=maxiter_newton,
            )
            if !parallel
                next!(p)
            end
        end

        finish!(p)


    end

    aero = AeroResult(
        T,
        Q,
        CT,
        CQ,
        converged,
        Np,
        Tp,
        α,
        cl,
        cd,
        re,
        u,
    )

    return aero
end

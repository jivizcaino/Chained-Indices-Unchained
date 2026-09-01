#------------------------------------------------------------------------------
#Replication code for: Chained Indices Unchained: Structural Transformation and the Welfare Foundations of Income Growth Measurement
#Benchmark HRV 2021 Model Under NHCES Preferences
#By:                   Omar Licandro and Juan I. Vizcaino
#This Version:         08/2026
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Set the Working Directory
currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Activate local project (Code/Project.toml)
using Pkg
#Pkg.activate(@__DIR__)
#Pkg.instantiate()

using XLSX, DataFrames, BlackBoxOptim ,  Statistics
using OrderedCollections, OrderedCollections, LsqFit,  Random, PrettyTables, LaTeXStrings, Plots; gr()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Configuration 
# Set to true to save figures, false to only display
save_figures = false

# Set to true to run SMM estimation, false to use pre-estimated parameters
run_smm      = false

# σ held fixed during SMM estimation (only ηs is estimated; ωc is calibrated
# internally to the 1980 goods share)
σ_fixed      = 0.20

# Choose σ when run_smm = false. Available values: 0.01, 0.10, 0.20, 0.30, 0.34
σ_choice     = 0.20
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Set the Parameter Values for Parameters Taken From HRV
ωx  = 0.650

εx  = 0.000
εc  = 0.000

θ   = 1/3
ρ   = 0.040
δ   = 0.080
ψ   = 1.500

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Import the HRV Data
file_name  = "HRV2021_Data.xlsx"
sheet_name = "Data"
data_range = "A1:DH78"

# Read the data from the Excel file
HRV_data   = XLSX.readtable(joinpath(datadir, file_name), sheet_name)
HRV_df     = DataFrame(HRV_data)
#HRV_df    = DataFrame(HRV_data...)

#Get the shares of Goods and Services in Investment and Consumption (X,C) respectively
VAX_GOOD_SHARE = HRV_df[HRV_df.year .>= 1980, "VAX_GOOD_S"]
VAX_SERV_SHARE = HRV_df[HRV_df.year .>= 1980, "VAX_SERV_S"]
VAC_GOOD_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_GOOD_S"]
VAC_SERV_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_SERV_S"]

#Get the Population
Nt             = HRV_df[HRV_df.year .>= 1980, "POP"]

#Get Total Efficiency Units of Labor
Lt             = HRV_df[HRV_df.year .>= 1980, "LAB_TOT_QI"]

#Get the Initial and Final Values for the TFP Indices
Ag_1980        = HRV_df[HRV_df.year .== 1980, "TFPVA_GOOD_I"][1]
Ag_2023        = HRV_df[HRV_df.year .== 2023, "TFPVA_GOOD_I"][1]
As_1980        = HRV_df[HRV_df.year .== 1980, "TFPVA_SERV_I"][1]
As_2023        = HRV_df[HRV_df.year .== 2023, "TFPVA_SERV_I"][1]
calAx_1980     = HRV_df[HRV_df.year .== 1980, "calA_X_I_TD"][1]
calAx_2023     = HRV_df[HRV_df.year .== 2023, "calA_X_I_TD"][1]

#Get the Initial and Final Values for Population and Efficiency Units of Labor
N_1980         = HRV_df[HRV_df.year .== 1980, "POP"][1]
N_2023         = HRV_df[HRV_df.year .== 2023, "POP"][1]

L_1980         = HRV_df[HRV_df.year .== 1980,"LAB_TOT_QI"][1]
L_2023         = HRV_df[HRV_df.year .== 2023,"LAB_TOT_QI"][1]

C_GOOD_P       = HRV_df[HRV_df.year .>= 1980, "C_GOOD_P"]./HRV_df[HRV_df.year .== 1980, "C_GOOD_P"]
C_SERV_P       = HRV_df[HRV_df.year .>= 1980, "C_SERV_P"]./HRV_df[HRV_df.year .== 1980, "C_SERV_P"]
X_TOT_P        = HRV_df[HRV_df.year .>= 1980, "X_TOT_P"] ./HRV_df[HRV_df.year .== 1980, "X_TOT_P"]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Define the Functions Characterizing the Equilibrium Conditions
# Recall that L = N*h, where N is the population and h is the avg efficiency units of labor per person
##Relative Prices
Ps(calAx,As,wedge_s)  = (calAx/As)*wedge_s
Pg(calAx,Ag,wedge_g)  = (calAx/Ag)*wedge_g

## Investment Productivity
calAx(Ax,Ag,As,ωx,εx) = Ax*(ωx * (Ag^(εx-1)) + (1-ωx) * (As^(εx-1)) )^(1/(εx-1))

# Aggregate Variables
## Output
Y(calAx,K,L,θ)        = calAx*(K^θ)*(L^(1-θ))

## Capital accumulation
dKdt(calAx,K,L,E,δ,θ) = calAx*(K^θ)*(L^(1-θ)) - E - δ*K

## Consumption Expenditure
E(Ps,Pg,Cs,Cg)        = Ps*Cs + Pg*Cg

# Per Capita Variables
## Output per Capita
y(calAx,k,h,θ)        = calAx*(k^θ)*(h^(1-θ))

## Capital accumulation
dkdt(calAx_hat,k,e,δ,n,θ)        = calAx_hat*(k^θ) - e - (δ+n)*k
dedt_e(calAx_hat,k,θ,ρ,δ,χ,g_ps) = (1 / (1-χ)) * (θ*calAx_hat*k^(θ-1) - ρ - δ - χ*g_ps )

# Cumulate a vector of interval growth rates into a level index, 1980 = 0.
# g[1] is the padding element for 1980 and must be zero. Each g[i] is the growth
# rate over the interval from i-1 to i, that is, already an integral over that
# interval, so the elements are summed. Passing them through a quadrature rule
# would average adjacent intervals and drop half of the first and last year.
function chain_index(g::AbstractVector)
    @assert iszero(g[1]) "chain_index expects g[1] = 0 (the 1980 padding element)."
    return cumsum(g)
end

# Tornqvist share weights for the interval [t-1, t]. Returns a vector of length T-1.
chain_share(s::AbstractVector) = 0.5 .* (s[1:end-1] .+ s[2:end])
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# NHCES Preference Block
function E_nhces(Pg, Ps, C, ωc, σ, ηg, ηs)
    """
    Consumption expenditure function:
    E_t = e(P_t, C_t)
    """
    inside =
        ωc       * Pg^(1 - σ) * C^((1 - σ) * ηg) +
        (1 - ωc) * Ps^(1 - σ) * C^((1 - σ) * ηs)

    return inside^(1 / (1 - σ))
end

function dEdC_nhces(Pg, Ps, C, ωc, σ, ηg, ηs)
    """
    Marginal expenditure cost:
    ∂E_t / ∂C_t
    """
    E = E_nhces(Pg, Ps, C, ωc, σ, ηg, ηs)

    bracket =
        ωc       * ηg * Pg^(1 - σ) * C^((1 - σ) * ηg - 1) +
        (1 - ωc) * ηs * Ps^(1 - σ) * C^((1 - σ) * ηs - 1)

    return E^σ * bracket
end

function sg_nhces(Pg, Ps, C, E, ωc, σ, ηg, ηs)
    """
    Goods share in consumption expenditure.
    """
    return ωc * Pg^(1 - σ) * C^((1 - σ) * ηg) * E^(σ - 1)
end

function ss_nhces(Pg, Ps, C, E, ωc, σ, ηg, ηs)
    """
    Services share in consumption expenditure.
    """
    return (1 - ωc) * Ps^(1 - σ) * C^((1 - σ) * ηs) * E^(σ - 1)
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# NHCES Equivalent Variation Functions (Section 1.5)
#
# EV income measure (eq. 50):
#   m̂_{t,z,τ} = m_τ + (Φ_t(P_{g,τ},P_{s,τ}) − Φ_t(P_{g,z},P_{s,z})) / ν_t
#
# where Φ_t(Pg,Ps) = max_C { U(C) − ν_t · E(Pg,Ps,C) }
# and the optimal Ĉ satisfies:  Ĉ^{−ψ} = ν_t · ∂E/∂C(Pg,Ps,Ĉ)
#------------------------------------------------------------------------------

# Solve for Bellman-optimal C̃: C̃^{−ψ} = ν_val · ∂E/∂C(Pg,Ps,C̃)
function solve_C_bellman_nhces(ν_val, Pg, Ps, ωc, σ, ηg, ηs, ψ;
                                C_low=1e-10, C_high=1e8,
                                tol=1e-12, maxiter=500)
    residual(C) = C^(-ψ) - ν_val * dEdC_nhces(Pg, Ps, C, ωc, σ, ηg, ηs)
    lo, hi = C_low, C_high
    while residual(lo) < 0
        lo /= 10
        lo < 1e-300 && return NaN
    end
    while residual(hi) > 0
        hi *= 10
        hi > 1e300 && return NaN
    end
    for _ in 1:maxiter
        mid = sqrt(lo * hi)          # geometric mean: fast on wide log-scale intervals
        val = residual(mid)
        abs(val) < tol && return mid
        val > 0 ? (lo = mid) : (hi = mid)   # f(lo)≥0, f(hi)≤0 convention
    end
    return sqrt(lo * hi)
end

# Φ_t(Pg, Ps) = U(C̃) − ν_t · E(Pg,Ps,C̃)
function Phi_t_nhces(ν_val, Pg, Ps, ωc, σ, ηg, ηs, ψ)
    C_star = solve_C_bellman_nhces(ν_val, Pg, Ps, ωc, σ, ηg, ηs, ψ)
    !isfinite(C_star) && return NaN
    E_star = E_nhces(Pg, Ps, C_star, ωc, σ, ηg, ηs)
    U_star = (C_star^(1 - ψ) - 1) / (1 - ψ)
    return U_star - ν_val * E_star
end

# EV income: m̂_{t,z,τ} = y_τ + (Φ_t(P_τ) − Φ_t(P_z)) / ν_t  (eq. 50)
# y_τ = gross per-capita income at τ (matches the PIGL convention)
function mhat_nhces(t::Int, z::Int, τ::Int; ν_t, Pgt, Pst, yt, ωc, σ, ηg, ηs, ψ)
    ν_val = ν_t[t]
    Phi_τ = Phi_t_nhces(ν_val, Pgt[τ], Pst[τ], ωc, σ, ηg, ηs, ψ)
    Phi_z = Phi_t_nhces(ν_val, Pgt[z], Pst[z], ωc, σ, ηg, ηs, ψ)
    return yt[τ] + (Phi_τ - Phi_z) / ν_val
end

# Hypothetical expenditure Ê_{t,z} = E(P_{g,z},P_{s,z},Ĉ_{t,z})  (eq. 45)
function Ehat_nhces(t::Int, z::Int;
                    ν_t, Pgt, Pst, ωc, σ, ηg, ηs, ψ)
    C_hat = solve_C_bellman_nhces(ν_t[t], Pgt[z], Pst[z], ωc, σ, ηg, ηs, ψ)
    return E_nhces(Pgt[z], Pst[z], C_hat, ωc, σ, ηg, ηs)
end

# EV goods expenditure share ŝ_{g,t,z}  (eq. 52)
function sg_hat_nhces(t::Int, z::Int;
                      ν_t, Pgt, Pst, ωc, σ, ηg, ηs, ψ)
    C_hat = solve_C_bellman_nhces(ν_t[t], Pgt[z], Pst[z], ωc, σ, ηg, ηs, ψ)
    E_hat = E_nhces(Pgt[z], Pst[z], C_hat, ωc, σ, ηg, ηs)
    return sg_nhces(Pgt[z], Pst[z], C_hat, E_hat, ωc, σ, ηg, ηs)
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Solve C from E = e(P,C)
function solve_C_from_E(E_target, Pg, Ps, ωc, σ, ηg, ηs;
                        C_low=1e-12,
                        C_high=1e12,
                        tol=1e-12,
                        maxiter=500)

    f(C) = E_nhces(Pg, Ps, C, ωc, σ, ηg, ηs) - E_target

    lo = C_low
    hi = C_high

    # Expand lower bracket if necessary
    while f(lo) > 0
        lo /= 10
        if lo < 1e-300
            return NaN
        end
    end

    # Expand upper bracket if necessary
    while f(hi) < 0
        hi *= 10
        if hi > 1e300
            return NaN
        end
    end

    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        val = f(mid)

        if abs(val) < tol
            return mid
        elseif val > 0
            hi = mid
        else
            lo = mid
        end
    end

    return 0.5 * (lo + hi)
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Solve C_{t+1} from the NHCES Euler Equation
function solve_C_next(C_now,Pg_now, Ps_now, Pg_next, Ps_next, R_now,
                      ρ, δ, ψ,ωc, σ, ηg, ηs;
                      C_low_factor=1e-4,C_high_factor=1e4,tol=1e-12,maxiter=500)

    M_now = dEdC_nhces(Pg_now, Ps_now, C_now, ωc, σ, ηg, ηs)

    function euler_residual(C_next)
        M_next = dEdC_nhces(Pg_next, Ps_next, C_next, ωc, σ, ηg, ηs)

        lhs =
            ψ * log(C_next / C_now) +
            log(M_next / M_now)

        rhs = R_now - δ - ρ

        return lhs - rhs
    end

    lo = C_now * C_low_factor
    hi = C_now * C_high_factor

    while euler_residual(lo) > 0
        lo /= 10
        if lo < 1e-300
            return NaN
        end
    end

    while euler_residual(hi) < 0
        hi *= 10
        if hi > 1e300
            return NaN
        end
    end

    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        val = euler_residual(mid)

        if abs(val) < tol
            return mid
        elseif val > 0
            hi = mid
        else
            lo = mid
        end
    end

    return 0.5 * (lo + hi)
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Calibrate C0 and ωc to match the initial goods consumption expenditure share
function calibrate_C0_ωc_from_initial_share(E0,Pg0,Ps0,sg0_data,σ,ηg,ηs;z_low=-50.0,z_high=50.0,tol=1e-12,maxiter=500)
                                             
    if !(E0 > 0 && Pg0 > 0 && Ps0 > 0 && 0 < sg0_data < 1)
        error("Invalid inputs in calibrate_C0_ωc_from_initial_share.")
    end

    # From the expenditure share conditions:
    #
    # A = s_g0 * E0^(1-σ)
    # B = s_s0 * E0^(1-σ)
    #
    # where
    # A = ωc Pg0^(1-σ) C0^((1-σ)ηg)
    # B = (1-ωc) Ps0^(1-σ) C0^((1-σ)ηs)

    A = sg0_data * E0^(1 - σ)
    B = (1 - sg0_data) * E0^(1 - σ)

    ag = (1 - σ) * ηg
    as = (1 - σ) * ηs

    # Solve in z = log(C0) for numerical stability
    function residual_logC(z)
        C = exp(z)

        term_g = A / (Pg0^(1 - σ) * C^ag)
        term_s = B / (Ps0^(1 - σ) * C^as)

        return term_g + term_s - 1
    end

    lo = z_low
    hi = z_high

    f_lo = residual_logC(lo)
    f_hi = residual_logC(hi)

    expand_iter = 0
    while sign(f_lo) == sign(f_hi) && expand_iter < 100
        lo -= 10.0
        hi += 10.0

        f_lo = residual_logC(lo)
        f_hi = residual_logC(hi)

        expand_iter += 1
    end

    if sign(f_lo) == sign(f_hi)
        return (NaN, NaN)
    end

    z_mid = NaN
    f_mid = Inf

    iter = 0
    while abs(f_mid) > tol && iter < maxiter
        z_mid = 0.5 * (lo + hi)
        f_mid = residual_logC(z_mid)

        if sign(f_mid) == sign(f_lo)
            lo = z_mid
            f_lo = f_mid
        else
            hi = z_mid
            f_hi = f_mid
        end

        iter += 1
    end

    C0 = exp(z_mid)

    ωc =
        A /
        (
            Pg0^(1 - σ) *
            C0^((1 - σ) * ηg)
        )

    ωc = clamp(ωc, eps(), 1.0 - eps())

    return C0, ωc
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Calibrate the growth rates of Agt,Ast,Axt (these values are from the HRV TFP series)
g_Ag        = log(Ag_2023/Ag_1980)/(2023-1980)
g_As        = log(As_2023/As_1980)/(2023-1980)
g_calAx     = log(calAx_2023/calAx_1980)/(2023-1980)
g_n         = (log(N_2023/N_1980))/(2023-1980)
g_l         = (log(L_2023/L_1980))/(2023-1980)
g_h         = g_l - g_n
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Build the series for Agt, Ast, Axt (Normalized to 1 in 1980)
Ag_0        = 1.000
As_0        = 1.000
calAx_0     = 1.000
Ag_t        = Ag_0.*exp.( g_Ag.*((1980:2023) .- 1980))
As_t        = As_0.*exp.( g_As.*((1980:2023) .- 1980))
calAx_t     = calAx_0.*exp.( g_calAx.*((1980:2023) .- 1980))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Build Wedges Such that the Prices in the Model Match the Prices in the Data
###Build Series for Prices in the Data First

#Express Data Prices as Relative to the Investment Numeraire
Pg_t_data   = C_GOOD_P ./ X_TOT_P
Ps_t_data   = C_SERV_P ./ X_TOT_P

###Undistorted Model-Implied Prices
Ps_t_undist = Ps.(calAx_t,As_t,1.000)
Pg_t_undist = Pg.(calAx_t,Ag_t,1.000)

##Wedge: Data vs Model
wedge_Pg_Ps = (Pg_t_data./Ps_t_data)./(Pg_t_undist./Ps_t_undist)

##Fit an exponential trend with a fixed intercept: wedge_Pg_Ps ≈ exp(b * (year - 1980))
years = 1980:2023
xdata = collect(years) .- 1980
ydata = wedge_Pg_Ps

##Only fit the exponent b, force a=1
exp_model_fixed_a(x, p) = exp.(p[1] .* x)
p0    = [0.0]  # initial guess for b
fit   = curve_fit(exp_model_fixed_a, xdata, ydata, p0)
b_fit = fit.param[1]

g_wedge_Pg_Ps     = b_fit  #Annual Growth Rate of the Wedge
wedge_Pg_Ps_trend = exp_model_fixed_a(xdata, fit.param)

# Compute Equilibrium Prices (relative to the investment numeraire)
Ps_t  = Ps.(calAx_t,As_t,1.000)
Pg_t  = Pg.(calAx_t,Ag_t,wedge_Pg_Ps_trend)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Simulate the Model Under NHCES Preferences
#   Initialize the economy on the aggregate BGP in 1980,
#   then simulate the true NHCES dynamics thereafter.
#------------------------------------------------------------------------------
function sim_model_full_nhces(params::Vector{Float64};Pg_t::Vector{Float64},Ps_t::Vector{Float64},
                                θ::Float64,ρ::Float64,δ::Float64,ψ::Float64,g_calAx::Float64,g_As::Float64,g_Ag::Float64,
                                g_h::Float64,g_l::Float64,g_n::Float64,sg0_data::Float64)

    σ  = params[1]
    ηs = params[2]
    ηg = 1.0

    T = length(1980:2023)

    # Basic admissibility checks
    if !(σ > 0 && ψ > 0 && ηg > 0 && ηs > 0)
        eq = Dict{String, Vector{Float64}}()
        eq["sgt"] = fill(1e6, T)
        return eq
    end

    # Exogenous series
    N_t         = 1.000 .* exp.(g_n .* ((1980:2023) .- 1980))
    h_t         = 1.000 .* exp.(g_h .* ((1980:2023) .- 1980))
    L_t         = 1.000 .* exp.(g_l .* ((1980:2023) .- 1980))

    Ag_t        = 1.000 .* exp.(g_Ag    .* ((1980:2023) .- 1980))
    As_t        = 1.000 .* exp.(g_As    .* ((1980:2023) .- 1980))
    calAx_t     = 1.000 .* exp.(g_calAx .* ((1980:2023) .- 1980))

    calAx_hat_t = calAx_t .* (h_t .^ (1 - θ))

    n           = g_n
    g_abgp      = g_calAx / (1 - θ) + g_h

    #--------------------------------------------------------------------------
    # Step 1: Find k0 such that initial expenditure growth equals the ABGP rate
    #--------------------------------------------------------------------------

    function E0_from_k0(k0)
        return calAx_hat_t[1] * k0^θ - (g_abgp + δ + n) * k0
    end

    function initial_expenditure_growth_error(k0)

        E0 = E0_from_k0(k0)

        # E0 <= 0 means k0 is too large (over-capitalised): growth is below target.
        if !isfinite(E0) || E0 <= 0
            return -1e6
        end

        C0, ωc_local = calibrate_C0_ωc_from_initial_share(E0,Pg_t[1],Ps_t[1],sg0_data,σ,ηg,ηs)

        if !(0 < ωc_local < 1)
            return -1e6
        end

        R0 = θ * calAx_hat_t[1] * k0^(θ - 1)

        C1 = solve_C_next(C0,
                          Pg_t[1], Ps_t[1],
                          Pg_t[2], Ps_t[2],
                          R0,
                          ρ, δ, ψ,
                          ωc_local, σ, ηg, ηs)

        if !isfinite(C1) || C1 <= 0
            return -1e6
        end

        E1 = E_nhces(Pg_t[2], Ps_t[2], C1, ωc_local, σ, ηg, ηs)

        if !isfinite(E1) || E1 <= 0
            return -1e6
        end

        g_E0 = log(E1 / E0)

        return g_E0 - g_abgp
    end

    # Bracket k0
    k0_min = 1e-4
    k0_max = 100.0

    err_min = initial_expenditure_growth_error(k0_min)
    err_max = initial_expenditure_growth_error(k0_max)

    iter_expand = 0
    while sign(err_min) == sign(err_max) && iter_expand < 100
        k0_min /= 2
        k0_max *= 2

        err_min = initial_expenditure_growth_error(k0_min)
        err_max = initial_expenditure_growth_error(k0_max)

        iter_expand += 1
    end

    if sign(err_min) == sign(err_max)
        eq = Dict{String, Vector{Float64}}()
        eq["sgt"] = fill(1e6, T)
        return eq
    end

    # Bisection for k0
    k0_mid = NaN
    err_mid = Inf

    iter = 0
    while abs(err_mid) > 1e-10 && iter < 500
        k0_mid = 0.5 * (k0_min + k0_max)
        err_mid = initial_expenditure_growth_error(k0_mid)

        if sign(err_mid) == sign(err_min)
            k0_min = k0_mid
            err_min = err_mid
        else
            k0_max = k0_mid
            err_max = err_mid
        end

        iter += 1
    end

    k0 = k0_mid

    #--------------------------------------------------------------------------
    # Step 2: Compute E0 from the ABGP capital accumulation condition
    #--------------------------------------------------------------------------

    e0 = E0_from_k0(k0)

    #--------------------------------------------------------------------------
    # Step 3: Calibrate ωc and C0 jointly from E0 and the initial goods share
    #--------------------------------------------------------------------------

    c0, ωc = calibrate_C0_ωc_from_initial_share(e0, Pg_t[1], Ps_t[1], sg0_data, σ, ηg, ηs)

    if !(0 < ωc < 1)
        eq = Dict{String, Vector{Float64}}()
        eq["sgt"] = fill(1e6, T)
        return eq
    end

    #--------------------------------------------------------------------------
    # Storage
    #--------------------------------------------------------------------------

    vars = ["Yt", "Xt", "Et", "Kt",
            "yt", "et", "ct", "kt", "ht",
            "sst", "sgt",
            "cst", "cgt", "Cst", "Cgt",
            "Xst", "Xgt",
            "Pst", "Pgt",
            "Kst", "Kgt",
            "Lt", "Nt", "Ht",
            "Lst", "Lgt",
            "Yst", "Ygt",
            "Wt", "Rt",
            "Agt", "Ast",
            "Axt", "calAxt", "calAxhat",
            "ωc"]

    eq = Dict{String, Vector{Any}}()
    for var in vars
        eq[var] = Vector{Any}()
    end

    push!(eq["ωc"], ωc)
    
    # Initial values
    push!(eq["et"], e0)
    push!(eq["ct"], c0)
    push!(eq["kt"], k0)

    push!(eq["ht"], h_t[1])
    push!(eq["Agt"], Ag_t[1])
    push!(eq["Ast"], As_t[1])

    push!(eq["calAxt"], calAx_t[1])
    push!(eq["calAxhat"], calAx_hat_t[1])

    push!(eq["Pst"], Ps_t[1])
    push!(eq["Pgt"], Pg_t[1])

    push!(eq["Lt"], L_t[1])
    push!(eq["Nt"], N_t[1])

    push!(eq["yt"], y(eq["calAxt"][1], eq["kt"][1], eq["ht"][1], θ))

    #--------------------------------------------------------------------------
    # Step 4: Simulate forward
    #--------------------------------------------------------------------------

    for t in 1:T

        # Consumption expenditure shares
        sg_t = sg_nhces(eq["Pgt"][t], eq["Pst"][t],
                        eq["ct"][t], eq["et"][t],
                        ωc, σ, ηg, ηs)

        ss_t = ss_nhces(eq["Pgt"][t], eq["Pst"][t],
                        eq["ct"][t], eq["et"][t],
                        ωc, σ, ηg, ηs)

        if !isfinite(sg_t) || !isfinite(ss_t) ||
           sg_t < 0 || sg_t > 1 ||
           ss_t < 0 || ss_t > 1
            return Dict("sgt" => fill(1e6, T))
        end

        push!(eq["sgt"], sg_t)
        push!(eq["sst"], ss_t)

        # Aggregate quantities
        push!(eq["Yt"], eq["yt"][t] * eq["Nt"][t])
        push!(eq["Kt"], eq["kt"][t] * eq["Nt"][t])
        push!(eq["Et"], eq["et"][t] * eq["Nt"][t])

        X_t = eq["Yt"][t] - eq["Et"][t]
        push!(eq["Xt"], X_t)

        # Investment allocation
        Xs_t =
            eq["Xt"][t] /
            (
                1 +
                (ωx / (1 - ωx)) *
                (eq["Agt"][t] / eq["Ast"][t])^εx
            )

        Xg_t = X_t - Xs_t

        push!(eq["Xst"], Xs_t)
        push!(eq["Xgt"], Xg_t)

        # Sectoral consumption quantities
        cs_t = (eq["sst"][t] * eq["et"][t]) / eq["Pst"][t]
        cg_t = (eq["sgt"][t] * eq["et"][t]) / eq["Pgt"][t]

        push!(eq["cst"], cs_t)
        push!(eq["cgt"], cg_t)

        push!(eq["Cst"], cs_t * eq["Nt"][t])
        push!(eq["Cgt"], cg_t * eq["Nt"][t])

        # Factor prices. R is the object that enters the Euler equation, so it is
        # stored in the same form used there: R = θ calAxhat k^(θ-1), which equals
        # θ calAx K^(θ-1) L^(1-θ) with L = N h.
        Rt_t = θ * eq["calAxhat"][t] * eq["kt"][t]^(θ - 1)
        Wt_t = (1 - θ) * eq["calAxhat"][t] * eq["kt"][t]^θ

        push!(eq["Rt"], Rt_t)
        push!(eq["Wt"], Wt_t)

        # Labor allocation.
        # Xst and Xgt are expenditure values, since Xgt = Xt - Xst, so the ratio
        # of investment expenditures is Xgt/Xst and must not be multiplied by
        # prices again. Cgt and Cst are quantities, so the consumption term does
        # carry prices.
        Ls =
            (eq["Xt"][t] / eq["Yt"][t]) *
            (
                1 /
                (
                    (eq["Xgt"][t] / eq["Xst"][t])
                    + 1
                )
            ) +
            (eq["Et"][t] / eq["Yt"][t]) *
            (
                1 /
                (
                    ((eq["Pgt"][t] * eq["Cgt"][t]) /
                     (eq["Pst"][t] * eq["Cst"][t]))
                    + 1
                )
            )

        if !isfinite(Ls) || Ls < 0 || Ls > 1
            return Dict("sgt" => fill(1e6, T))
        end

        Lg = 1 - Ls

        push!(eq["Lst"], Ls)
        push!(eq["Lgt"], Lg)

        Ks = eq["Kt"][t] * Ls
        Kg = eq["Kt"][t] * Lg

        push!(eq["Kst"], Ks)
        push!(eq["Kgt"], Kg)

        Ys = eq["Ast"][t] * (Ks^θ) * (Ls^(1 - θ))
        Yg = eq["Agt"][t] * (Kg^θ) * (Lg^(1 - θ))

        push!(eq["Yst"], Ys)
        push!(eq["Ygt"], Yg)

        # Stop after 2023
        if t == T
            break
        end

        # Push next-period exogenous objects
        push!(eq["Pst"], Ps_t[t+1])
        push!(eq["Pgt"], Pg_t[t+1])

        push!(eq["calAxt"], calAx_t[t+1])
        push!(eq["calAxhat"], calAx_hat_t[t+1])

        push!(eq["Agt"], Ag_t[t+1])
        push!(eq["Ast"], As_t[t+1])

        push!(eq["ht"], h_t[t+1])
        push!(eq["Lt"], L_t[t+1])
        push!(eq["Nt"], N_t[t+1])

        #--------------------------------------------------------------
        # ABGP step. Section 4.1 assumes the economy has been on an ABGP
        # since 1980, so per capita consumption expenditure grows at
        # g_abgp each period and C_{t+1} is recovered by inverting the
        # NHCES expenditure function. The k0 bisection above imposes the
        # Euler equation at 1980. With ηg ≠ ηs there is no exact BGP in C,
        # so the Euler equation holds at 1980 and the residual thereafter
        # measures the size of the ABGP restriction. It is reported below.
        #--------------------------------------------------------------

        e_next = eq["et"][t] * exp(g_abgp)

        c_next = solve_C_from_E(e_next,
                                Pg_t[t+1], Ps_t[t+1],
                                ωc, σ, ηg, ηs)

        if !isfinite(c_next) || !isfinite(e_next) ||
           c_next <= 0 || e_next <= 0
            return Dict("sgt" => fill(1e6, T))
        end

        push!(eq["ct"], c_next)
        push!(eq["et"], e_next)

        #--------------------------------------------------------------
        # Capital accumulation
        #--------------------------------------------------------------

        g_k =
            dkdt(eq["calAxhat"][t],
                 eq["kt"][t],
                 eq["et"][t],
                 δ, n, θ) / eq["kt"][t]

        if !isfinite(g_k) || abs(g_k) > 1.0
            return Dict("sgt" => fill(1e6, T))
        end

        k_next = eq["kt"][t] * exp(g_k)

        push!(eq["kt"], k_next)

        push!(
            eq["yt"],
            y(eq["calAxt"][t+1],
              eq["kt"][t+1],
              eq["ht"][t+1],
              θ)
        )
    end

    return eq
end
#------------------------------------------------------------------------------
if run_smm

println("Running SMM Estimation...")
println("Estimating ηs only; fixing σ = $σ_fixed, ηg = 1, and calibrating ωc to match the 1980 goods share.")

    ηg       = 1.000
    ### Initial goods expenditure share in the data
    sg0_data = VAC_GOOD_SHARE[1]

    ### Pass fixed parameters to the simulation function
    ### params[1] = ηs  (σ is fixed at σ_fixed)
    sim_model(params) =
        sim_model_full_nhces([σ_fixed , params[1]];
            Pg_t=Pg_t,Ps_t=Ps_t,θ=θ,ρ=ρ,δ=δ,ψ=ψ,g_calAx=g_calAx,g_As=g_As,g_Ag=g_Ag,g_h=g_h,g_l=g_l,g_n=g_n,sg0_data=sg0_data)

    ### Define the SSE function to minimize
    function SSE(params; cons_goods_share_data)

        sim_eq = sim_model(params)

        if !haskey(sim_eq, "sgt") ||
           length(sim_eq["sgt"]) != length(cons_goods_share_data) ||
           any(x -> !isfinite(x) || x > 1.0 || x < 0.0, sim_eq["sgt"])

            @show "Invalid simulation detected."
            return 1e10
        end

        goods_errors = (sim_eq["sgt"] .- cons_goods_share_data)

        SSE = sum(goods_errors .^ 2)

        return SSE
    end

    ### Parameter space
    ###
    ### params[1] = ηs
    ###
    ### σ is fixed at σ_fixed.
    ### ηg is fixed at 1.
    ParamSpace = [(1.00001, 20.000)] # ηs
            
    ### Seeds for repeated global optimization
    seeds = [1234]
    optimization_results = OrderedDict()

    ### Run optimization for different seeds
    for seed in seeds
        println("\n" * "="^80)
        println("Running optimization with seed: $seed")
        println("="^80)

        opt_problem =
            bbsetup( params -> SSE(params;cons_goods_share_data=VAC_GOOD_SHARE);                    
                SearchRange=ParamSpace,TraceMode=:compact,Method=:adaptive_de_rand_1_bin,                
                PopulationSize = 100, MaxFuncEvals = 50_000,FitnessTolerance = 1e-20,
                MaxStepsWithoutProgress = 20_000, rng=MersenneTwister(seed) )
               
        opt_results = bboptimize(opt_problem)

        best_params_seed = best_candidate(opt_results)

        optimization_results[seed] = OrderedDict(
            "params"  => best_params_seed,
            "fitness" => best_fitness(opt_results),
            "seed"    => seed,
            "ηs"      => best_params_seed[1]
        )

        println("Seed $seed: Fitness = $(best_fitness(opt_results))")
        println(
            "Parameters: " *
            "ηs=$(round(best_params_seed[1], digits=6))"
        )
    end

    ### Find best result across seeds
    best_seed_index =
        argmin([optimization_results[s]["fitness"] for s in seeds])

    best_result =
        optimization_results[seeds[best_seed_index]]

    ### Re-run the model at the best parameters to recover calibrated ωc
    best_sim_eq =
        sim_model(best_result["params"])

    ωc_best =
        best_sim_eq["ωc"][1]

    println("""
    $(repeat("=", 80))
    BEST RESULT ACROSS ALL SEEDS
    $(repeat("=", 80))
    Best seed:    $(best_result["seed"])
    Best fitness: $(best_result["fitness"])

    Estimated parameters:
    ηs = $(round(best_result["ηs"], digits=6))

    Fixed parameters:
    σ  = $(σ_fixed)
    ηg = 1.000000

    Calibrated parameter:
    ωc = $(round(ωc_best, digits=6))
    $(repeat("=", 80))
    """)

    ### Extract best parameters
    σ  = σ_fixed
    ηs = round(best_result["params"][1], digits=5)

    ηg = 1.0
    ωc = round(ωc_best, digits=5)

    @assert σ == σ_fixed "σ was not set from σ_fixed; check the run_smm branch."

    #------------------------------------------------------------------------------
    # TABLE 2 - Calibrated / Estimated Parameters
    println("""
    $(repeat("=", 80))
    NHCES Preference Parameters
    $(repeat("=", 80))

    Estimated via SMM:
    ηs = $(round(best_result["ηs"], digits=6))

    Fixed:
    σ  = $(σ_fixed)
    ηg = 1.000000

    Calibrated to match the 1980 goods consumption expenditure share:
    ωc = $(round(ωc_best, digits=6))
    $(repeat("=", 80))
    """)
#------------------------------------------------------------------------------

else
    println("Skipping SMM estimation, using pre-estimated NHCES parameters...")

    ηg = 1.0

    # Pre-estimated parameters for each σ value
    # σ_choice (set in Configuration block) selects the calibration.
    # Available σ values: 0.01, 0.10, 0.20, 0.30, 0.34
    nhces_params = Dict(
        0.01 => (σ=0.010000, ωc=0.243124, ηs=1.201051),
        0.10 => (σ=0.100000, ωc=0.245369, ηs=1.418465),
        0.20 => (σ=0.200000, ωc=0.24797741620073313, ηs=1.854580),   # benchmark
        0.30 => (σ=0.300000, ωc=0.250759, ηs=2.931970),
        0.34 => (σ=0.340000, ωc=0.251939, ηs=3.991059),
    )

    if !haskey(nhces_params, σ_choice)
        error("σ_choice = $σ_choice is not available. Choose from: $(sort(collect(keys(nhces_params))))")
    end

    p_sel = nhces_params[σ_choice]
    σ     = p_sel.σ
    ωc    = p_sel.ωc
    ηs    = p_sel.ηs

    println("Selected σ = $σ  →  ωc = $ωc,  ηs = $ηs")

    # ωc is also calibrated internally by sim_model_full_nhces to match VAC_GOOD_SHARE[1].
end
#------------------------------------------------------------------------------
# Simulate the Model with the Calibrated Parameters
sim_eq = sim_model_full_nhces([σ, ηs];
    Pg_t=Pg_t, Ps_t=Ps_t,
    θ=θ, ρ=ρ, δ=δ, ψ=ψ,
    g_calAx=g_calAx, g_As=g_As, g_Ag=g_Ag,
    g_h=g_h, g_l=g_l, g_n=g_n,
    sg0_data=VAC_GOOD_SHARE[1]
)

### Stop with a readable message if the simulation failed
if !haskey(sim_eq, "ωc") || any(x -> !isfinite(x) || x > 1e5, sim_eq["sgt"])
    error("""
    The NHCES simulation failed at σ = $σ, ηs = $ηs. The returned dictionary
    holds only the sentinel "sgt" series, so this parameter combination is not
    admissible.
    """)
end

ωc    = Float64(sim_eq["ωc"][1])

### ωc is calibrated inside the simulation to match VAC_GOOD_SHARE[1]; the value
### stored in nhces_params is a record of the estimation, not an input.

#Compute the Marginal Value of Wealth
ct_f  = Float64.(sim_eq["ct"])
et_f  = Float64.(sim_eq["et"])
sg_C  = Float64.(sim_eq["sgt"])   # goods expenditure share (= P_g c_g / E)
ss_C  = Float64.(sim_eq["sst"])   # services expenditure share (= P_s c_s / E)
ν_t   = ct_f.^(1 - ψ) ./ (et_f .* (ηg .* sg_C .+ ηs .* ss_C))

#-----------------------------------------------------------------------------------------------------

#---------------------------------------------------------------
#Check that the Model is in ABGP
kt_f    = Float64.(sim_eq["kt"])
et_f    = Float64.(sim_eq["et"])
Yt_f    = Float64.(sim_eq["Yt"])
Et_agg  = Float64.(sim_eq["Et"])    # aggregate consumption expenditure
Pgt_f   = Float64.(sim_eq["Pgt"])
Pst_f   = Float64.(sim_eq["Pst"])

gk      = log.(kt_f[2:end])   .- log.(kt_f[1:end-1])
ge      = log.(et_f[2:end])   .- log.(et_f[1:end-1])
gY      = log.(Yt_f[2:end])   .- log.(Yt_f[1:end-1])
gE      = log.(Et_agg[2:end]) .- log.(Et_agg[1:end-1])

g_pg    = log.(Pgt_f[2:end]) .- log.(Pgt_f[1:end-1])
g_ps    = log.(Pst_f[2:end]) .- log.(Pst_f[1:end-1])

Rt_f      = Float64.(sim_eq["Rt"])
dEdC_i(i) = dEdC_nhces(Pgt_f[i], Pst_f[i], ct_f[i], ωc, σ, ηg, ηs)
M_f       = [dEdC_i(i) for i in eachindex(ct_f)]
res_euler = ψ .* diff(log.(ct_f)) .+ diff(log.(M_f)) .- (Rt_f[1:end-1] .- δ .- ρ)

gY_f      = diff(log.(Yt_f))
ky_f      = kt_f ./ Float64.(sim_eq["yt"])
se_f      = et_f ./ Float64.(sim_eq["yt"])

println("""
$(repeat("=", 62))
ABGP check and Euler residual
$(repeat("=", 62))
  g_Y  1981 / 2023             : $(round(gY_f[1], digits=6)) / $(round(gY_f[end], digits=6))
  g_Y  spread over the sample  : $(round(maximum(gY_f) - minimum(gY_f), sigdigits=3))
  g_abgp + g_n                 : $(round(g_calAx/(1-θ) + g_h + g_n, digits=6))
  k/y  1980 / 2023             : $(round(ky_f[1], digits=5)) / $(round(ky_f[end], digits=5))
  s_e  1980 / 2023             : $(round(se_f[1], digits=5)) / $(round(se_f[end], digits=5))
  Euler residual, 1980 / 2023  : $(round(res_euler[1], sigdigits=3)) / $(round(res_euler[end], sigdigits=3))
  Euler residual, max |.|      : $(round(maximum(abs.(res_euler)), sigdigits=3))
  implied error in g_C, pp/yr  : $(round(100*maximum(abs.(res_euler))/ψ, digits=3))
$(repeat("=", 62))
""")

@assert maximum(gY_f) - minimum(gY_f) < 1e-8 "Aggregate output growth is not constant: the economy is not on the ABGP."
#----------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
#CHARTS
# Define font settings for all plots
tickfont   = font(12)
guidefont  = font(12)
legendfont = font(12)
#-----------------------------------------------------------------------------------------------------

#----------------------------------------------------------
### Figure 3 (b)
year = 1980:2023
plot(year, sim_eq["sgt"], label=L"s_{g,t} \textrm{- \ Model}", linestyle=:solid, color=:black, lw=2.00, 
    ylabel="Share of Goods in Consumption\nExpenditure",
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,legendfont=legendfont,
    legend=(0.775, 0.900),
    xticks=1980:5:2025, yticks=0.0:0.02:0.80,xrotation=45)

plot!(year, VAC_GOOD_SHARE, label=L"s_{g,t} \textrm{- \ Data}",  linestyle=:dot, lw=2.00,
    minorgrid=false, color=:black,left_margin=5Plots.mm,framestyle=:box)

if save_figures 
    savefig(joinpath(figuresdir, "sg_t_Model_Fit_NHCES.png"))
    println("Figure saved to: ", joinpath(figuresdir, "sg_t_Model_Fit_NHCES.png"))
end
#----------------------------------------------------------

#---------------------------------------------------------
### Figure 4 (a)
GDP_data = HRV_df[HRV_df.year .>= 1980, "GDP_QI"]./ HRV_df[HRV_df.year .>= 1980, "GDP_QI"][1]

#Compute the Divisia Index
x_t         = sim_eq["Xt"]./sim_eq["Nt"]
s_e         = sim_eq["Et"]./sim_eq["Yt"]
sg_t        = sim_eq["sgt"]
ss_t        = sim_eq["sst"]

g_cg        = sim_eq["cgt"][2:end]./sim_eq["cgt"][1:end-1] .- 1
g_cs        = sim_eq["cst"][2:end]./sim_eq["cst"][1:end-1] .- 1
g_x         = x_t[2:end]./x_t[1:end-1] .- 1

sg_w       = chain_share(Float64.(sim_eq["sgt"]))
se_w       = chain_share(Float64.(s_e))

gD         = se_w.*( sg_w.*g_cg .+ (1 .- sg_w).*g_cs ) .+ (1 .- se_w).*g_x
gD_agg     = [ 0.0 ; (gD.+ g_n) ]
FS         = chain_index(gD_agg)

plot(1980:2023, FS , ylabel="GDP Index\n(log scale; 1980=0)",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1980:5:2024, yticks=0.0:0.2:1.4,
    ylim=(0.00, 1.40),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    label = "Chained Index - Model",
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

plot!(1980:2023, log.(GDP_data),label = "Chained Index - Data",
    linestyle=:dot, lw=2.0, color=:black,
    legend=(0.15, 0.90))

if save_figures 
    savefig(joinpath(figuresdir,"GDP_Model_vs_Data_NHCES.png"))
    println("Figure saved to: ", joinpath(figuresdir, "GDP_Model_vs_Data_NHCES.png"))
end
#---------------------------------------------------------

#---------------------------------------------------------------
### Figure 4 (b)
g_FS          = gD .+ g_n 
GDP_data_long = HRV_df[!, "GDP_QI"]

# Compute the growth rate (log difference)
g_GDP_long = [NaN; GDP_data_long[2:end] ./ GDP_data_long[1:end-1] .- 1]

# Compute 10-period centered moving average (5 before, 5 after)
function moving_average(x, window)
    n = length(x)
    half = div(window, 2)
    ma = similar(x, Float64)
    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        ma[i] = mean(skipmissing(x[lo:hi]))
    end
    return ma
end

g_GDP_long_ma10 = moving_average(g_GDP_long, 11)

plot(1981:2023, g_GDP_long_ma10[1981-1947+1:end],
    label = "Data - 10-year Moving Average", 
    linestyle=:dot, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1980:5:2025,
    yticks=0.00:0.01:0.05,
    ylim=(0.00, 0.05),    
    xtickfont=tickfont, ytickfont=guidefont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

plot!(1981:2023, g_FS ,label = "Model - FS Chained Index", 
    ylabel="GDP Growth Rate",linestyle=:solid, 
    lw=2.0, color=:black, minorgrid=false,
    xticks=1980:5:2024, yticks=0.010:0.010:0.50,
    ylim=(0.010,0.050),legend=(0.15, 0.90),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

if save_figures 
    savefig(joinpath(figuresdir,"GDP_Growth_Model_vs_Data_NHCES.png"))
    println("Figure saved to: ", joinpath(figuresdir, "GDP_Growth_Model_vs_Data_NHCES.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 5 (a)
sc   = sim_eq["et"]./sim_eq["yt"]
sg_z = sim_eq["sgt"]
ss_z = sim_eq["sst"]
se_z = sc[1:end]

#Compute the Marginal Value of Capital (NHCES formula, same as above)
ν_t  = ct_f.^(1 - ψ) ./ (et_f .* (ηg .* sg_C .+ ηs .* ss_C))

#Gross Income Per Capita entering the EV measure (m_t = y_t), matching the
#gross concept behind the GDP Divisia index, which uses s_e = e/y.
yt_f = Float64.(sim_eq["yt"])
m_t  = yt_f
# Consumption expenditure per capita. et_f already holds E_t/N_t, so this is a
# consistency check rather than a redefinition. Note that Et_agg below keeps the
# aggregate series distinct from the per-capita one.
s_e_t = Float64.(sim_eq["Et"] ./ sim_eq["Yt"])
@assert maximum(abs.(s_e_t .* yt_f .- et_f)) < 1e-8

# Convenience closures capturing NHCES parameters
_mhat(t, z, τ) = mhat_nhces(t, z, τ;
    ν_t=ν_t, Pgt=Pgt_f, Pst=Pst_f, yt=yt_f,
    ωc=ωc, σ=σ, ηg=ηg, ηs=ηs, ψ=ψ)
_Ehat(t, z) = Ehat_nhces(t, z;
    ν_t=ν_t, Pgt=Pgt_f, Pst=Pst_f, ωc=ωc, σ=σ, ηg=ηg, ηs=ηs, ψ=ψ)
_sg_hat(t, z) = sg_hat_nhces(t, z;
    ν_t=ν_t, Pgt=Pgt_f, Pst=Pst_f, ωc=ωc, σ=σ, ηg=ηg, ηs=ηs, ψ=ψ)

τ   = (1980:1:2023) .- 1980 .+ 1
T   = length(τ)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# NHCES Fisher-Shell Indices — EV definition
#
# In m̂_{t,z,τ} the first index is the base preference order, the second the base
# prices and the third the date being evaluated. A fixed-base index sets t = z at
# the base year and lets the evaluated date run:
#
#   P_{t1,τ} = log m̂_{t1,t1,τ} − log m̂_{t1,t1,t0}     (paper, eq. 31)
#   L_{t0,τ} = log m̂_{t0,t0,τ} − log m_{t0}           (paper, eq. 32)
#
# so the calls below are _mhat(t_base, t_base, τ) with τ rolling. The paper writes
# the rolling index as z in equations (31) and (32); it sits in the third slot.
# LV_main.jl rolls the second slot instead, which fixes income at y_{t_base} and
# is not the same object.
#------------------------------------------------------------------------------

# Base year t0 = 2023 — m̂_{t,τ,z}, z rolling (t=τ=2023)
t_base           = 2023
t_prime          = t_base - 1980 + 1
mhat_2023_2023_τ = [_mhat(t_prime, t_prime,zi) for zi in τ]
@assert isapprox(mhat_2023_2023_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_2023_2023_τ   = log.(mhat_2023_2023_τ) .- log.(mhat_2023_2023_τ[1])
g_2023_2023_τ    = diff(FS_2023_2023_τ) 

# Base year t0 = 2010 — m̂_{t,τ,z}, z rolling (t=τ=2010)
t_base           = 2010
t_prime          = t_base - 1980 + 1
mhat_2010_2010_τ = [_mhat(t_prime, t_prime,zi) for zi in τ]
@assert isapprox(mhat_2010_2010_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_2010_2010_τ   = log.(mhat_2010_2010_τ) .- log.(mhat_2010_2010_τ[1]) 
g_2010_2010_τ    = diff(FS_2010_2010_τ) 

# Base year t0 = 2000 — m̂_{t,τ,z}, z rolling (t=τ=2000)
t_base           = 2000
t_prime          = t_base - 1980 + 1
mhat_2000_2000_τ = [_mhat(t_prime, t_prime, zi ) for zi in τ]
@assert isapprox(mhat_2000_2000_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_2000_2000_τ   = log.(mhat_2000_2000_τ) .- log.(mhat_2000_2000_τ[1]) 

# Base year t0 = 1990 — m̂_{t,τ,z}, z rolling (t=τ=1990)
t_base           = 1990
t_prime          = t_base - 1980 + 1
mhat_1990_1990_τ = [_mhat(t_prime, t_prime, zi) for zi in τ]
@assert isapprox(mhat_1990_1990_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_1990_1990_τ   = log.(mhat_1990_1990_τ) .- log.(mhat_1990_1990_τ[1]) 

# Base year t0 = 1980 — m̂_{t,τ,z}, z rolling (t=τ=1980)
t_base           = 1980
t_prime          = t_base - 1980 + 1
mhat_1980_1980_τ = [_mhat(t_prime, t_prime,zi) for zi in τ]
@assert isapprox(mhat_1980_1980_τ[1], yt_f[t_prime]; rtol=1e-10)    
FS_1980_1980_τ   = log.(mhat_1980_1980_τ) .- log.(mhat_1980_1980_τ[1]) 
g_1980_1980_τ    = diff(FS_1980_1980_τ)

#------------------------------------------------------------------------------
# Population index and aggregate EV indices
# FS_b_agg = FS_b_pc + log(N_z / N_1980);  FS (Divisia) is already aggregate
#------------------------------------------------------------------------------
pop_index = log.(sim_eq["Nt"] ./ sim_eq["Nt"][1])
g_N       = diff(log.(sim_eq["Nt"]))

FS_2023_2023_τ_agg = FS_2023_2023_τ .+ pop_index
FS_2010_2010_τ_agg = FS_2010_2010_τ .+ pop_index
FS_2000_2000_τ_agg = FS_2000_2000_τ .+ pop_index
FS_1990_1990_τ_agg = FS_1990_1990_τ .+ pop_index
FS_1980_1980_τ_agg = FS_1980_1980_τ .+ pop_index

g_2023_2023_τ_agg = diff(FS_2023_2023_τ_agg)
g_2010_2010_τ_agg = diff(FS_2010_2010_τ_agg)
g_2000_2000_τ_agg = diff(FS_2000_2000_τ_agg)
g_1990_1990_τ_agg = diff(FS_1990_1990_τ_agg)
g_1980_1980_τ_agg = diff(FS_1980_1980_τ_agg)

@assert maximum(abs.((FS_2023_2023_τ_agg .- FS_2023_2023_τ) .- pop_index)) < 1e-10
@assert maximum(abs.(g_2023_2023_τ_agg .- (g_2023_2023_τ .+ g_N))) < 1e-10

p = plot(1980:2023, FS_2023_2023_τ_agg, label=L"\mathcal{P}_{2023,z} - \textrm{2023‑base \ Fisher‑Shell \ Index}",
    ylabel="Cumulative Growth", linestyle=:dash, lw=2.0,
    xticks=1980:5:2023, yticks=0.00:0.20:1.40,
    ylims=(0.00, 1.40),
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.10, 0.920), legendfont=legendfont,
    xrotation=45, framestyle=:box)

plot!(p, 1980:2023, FS_1980_1980_τ_agg, label=L"\mathcal{L}_{1980,z} - \textrm{1980‑base \  Fisher‑Shell \ Index}",
    linestyle=:dot, lw=2, color=:black)

plot!(p, 1980:2023, FS, label=L"\mathcal{D}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    linestyle=:solid, lw=2, color=:black)

if save_figures
    savefig(joinpath(figuresdir,"FS_BBEV_NHCES.png"))
    println("Figure saved to: ", joinpath(figuresdir, "FS_BBEV_NHCES.png"))
end

#Differences Between Fixed-Base and Chained Indices
println("""
$(repeat("=", 60))
Differences Between Fixed-Base and Chained Indices (2023)
$(repeat("=", 60))
FS Index (base 1980) - Chained Index: $(round(FS_1980_1980_τ_agg[end] - FS[end], digits=3))
FS Index (base 2023) - Chained Index: $(round(FS_2023_2023_τ_agg[end] - FS[end], digits=3))
$(repeat("=", 60))
""")
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 6 (a)
# Equivalent variation measure of the consumption expenditure share.
#
# Notation. In m̂_{t,z,τ} the first index is the base preference order, the
# second the base prices and the third the date being evaluated. A fixed-base
# index sets t = z at the base year and lets τ run over 1980:2023. The running
# index in this block is therefore τ, not z.
#
# Appendix D, equation (40):
#
#     ŝ_{e,t,z,τ} = ê_{t,z} / m̂_{t,z,τ}
#
# The numerator is constant in τ, since quasi-linearity pins hypothetical
# consumption expenditure independently of the date evaluated. The denominator is
# the equivalent variation measure of income, not actual income. For a fixed base
# t = z the ratio therefore starts at ê_{t,t}/m̂_{t,t,t} = e_t/m_t = s_{e,t} and,
# for the current-base index, declines in τ converging from above to s_{e,t1}.
#
# mhat_1980_1980_τ and mhat_2023_2023_τ are the fixed-base EV income series
# already built for Figure 5(a), so they are reused here.
Ehat_1980_1980     = _Ehat(1, 1)     # ê_{1980,1980} = e_1980
Ehat_2023_2023     = _Ehat(T, T)     # ê_{2023,2023} = e_2023

se_1980_τ          = Ehat_1980_1980 ./ mhat_1980_1980_τ
se_2023_τ          = Ehat_2023_2023 ./ mhat_2023_2023_τ

se                 = et_f ./ yt_f

### At its own base the hypothetical objects collapse onto the actual ones,
### because ν_t is defined so that C_t solves the Bellman first order condition
### at prices P_t. Ê_{t,t} = e_t and m̂_{t,t,t} = m_t, so ŝ_{e,t,t,t} = s_{e,t}.
@assert isapprox(Ehat_1980_1980, et_f[1]; rtol=1e-8)
@assert isapprox(Ehat_2023_2023, et_f[T]; rtol=1e-8)
@assert isapprox(mhat_1980_1980_τ[1], yt_f[1]; rtol=1e-8)
@assert isapprox(mhat_2023_2023_τ[T], yt_f[T]; rtol=1e-8)
@assert isapprox(se_1980_τ[1], et_f[1]/yt_f[1]; rtol=1e-8)
@assert isapprox(se_2023_τ[T], et_f[T]/yt_f[T]; rtol=1e-8)

### Appendix D: the current-base share declines in τ, converging from above to s_e.
@assert all(diff(se_2023_τ) .< 0)

### Axis limits from the three series, so the panel is not compressed against a
### fixed range. Rounded outwards to the tick step.
se_step = 0.25
se_lo   = 0.25
se_hi   = 1.35

plot(1980:2023,se,
    label=L"s_{e,\tau}",
    ylabel="Share of Consumption Expenditure\nin Gross Income",
    linestyle=:solid, lw=2,
    xticks=1980:5:2025, yticks=se_lo:se_step:se_hi, ylims=(se_lo, se_hi),
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.75, 0.92),
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box) 

plot!(1980:2023,se_1980_τ,label=L"\hat{s}_{e,1980,1980,\tau}",
    linestyle=:dot, lw=2,color=:black)

plot!(1980:2023,se_2023_τ,label=L"\hat{s}_{e,2023,2023,\tau}",
    linestyle=:dash, lw=2,color=:black)

if save_figures
    savefig(joinpath(figuresdir,"se_t_NHCES.png"))
    println("Figure saved to: ", joinpath(figuresdir, "se_t_NHCES.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 6 (b)
#
# EV goods expenditure share ŝ_{g,t,z}. Appendix D, equation (41): applying Roy's
# identity to V(ê_{t,z}, P_{g,z}, P_{s,z}) gives an object that does not depend on
# the date τ being evaluated, because quasi-linearity pins ê_{t,z} independently
# of τ and base prices are fixed at z. For a fixed base t = z it further collapses
# to the actual share at the base year, ŝ_{g,t,t} = s_{g,t}. The two lines are
# therefore horizontal in τ by construction, which is the object plotted here.
sg_2023_τ  = _sg_hat(T, T)
sg_1980_τ  = _sg_hat(1, 1)

@assert isapprox(sg_1980_τ, sg_C[1]; rtol=1e-8)
@assert isapprox(sg_2023_τ, sg_C[T]; rtol=1e-8)

plot(1980:2023, sg_C, label=L"s_{g,\tau}",linestyle=:solid, lw=2,color=:black)
plot!(1980:2023, fill(sg_1980_τ, length(1980:2023)), label=L"\hat{s}_{g,1980,1980,\tau}",
    linestyle=:dot, lw=2,color=:black)

plot!(1980:2023, fill(sg_2023_τ, length(1980:2023)), label=L"\hat{s}_{g,2023,2023,\tau}",
    ylabel="Share of Goods in Consumption\nExpenditure",
    linestyle=:dash, lw=2,
    xticks=1980:5:2025, yticks=0.000:0.05:0.300,
    ylims=(0.000, 0.300),
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.75, 0.65),
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir, "sg_t_z_NHCES.png"))
    println("Figure saved to: ", joinpath(figuresdir, "sg_t_z_NHCES.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 7
g_FS              = gD .+ g_n
g_FS_2023_2023_τ  = g_2023_2023_τ_agg  
g_FS_2010_2010_τ  = g_2010_2010_τ_agg  
g_FS_2000_2000_τ  = g_2000_2000_τ_agg  
g_FS_1990_1990_τ  = g_1990_1990_τ_agg  
g_FS_1980_1980_τ  = g_1980_1980_τ_agg  

p = plot(1981:2023, g_FS_2023_2023_τ, label=L"g^{D}_{2023,z} - \textrm{2023‑base \  FS \ Index}",
    ylabel="Growth Rate",
    linestyle=:dash, lw=2.5,
    xticks=1980:5:2025, yticks=0.00:0.005:0.05,
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.475, 0.500),
    legendfont=legendfont,
    xrotation=45,
    framestyle=:box)

plot!(p, 1981:2023, g_FS_2010_2010_τ, label=L"g^{D}_{2010,z} - \textrm{\ 2010‑base \  FS \ Index}",
    linestyle=:dash, lw=2, color=:black)

plot!(p, 1981:2023, g_FS_2000_2000_τ, label=L"g^{D}_{2000,z}- \textrm{2000‑base \  FS \ Index}",
    linestyle=:dash, lw=1.5, minorgrid=false, color=:black)

plot!(p, 1981:2023, g_FS_1990_1990_τ, label=L"g^{D}_{1990,z} - \textrm{1990‑base \  FS \ Index}",
    linestyle=:dash, lw=1.0, minorgrid=false, color=:black)

plot!(p, 1981:2023, g_FS_1980_1980_τ, label=L"g^{D}_{1980,z} - \textrm{1980‑base \  FS \ Index}",
    linestyle=:dash, lw=0.5, minorgrid=false, color=:black)

plot!(p, 1981:2023, g_FS , label=L"g^{D}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    linestyle=:solid, lw=2, minorgrid=false, color=:black)

xaxis!(p, minor_ticks=true, minor_tick_step=1.00)
yaxis!(p, minor_ticks=true, minor_tick_step=0.01)

if save_figures
    savefig(joinpath(figuresdir, "FS_GrowthRates_1980_2023_NHCES.png"))
    println("Figure saved to: ", joinpath(figuresdir, "FS_GrowthRates_1980_2023_NHCES.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
#Decline in the Growth Rate Between 1981 and 2023 for the Chained Index
decline_growth_rate = (g_FS[end] - g_FS[1]) * 100

println("""
$(repeat("=", 60))
Decline in GDP Growth Rate (1981-2023)
$(repeat("=", 60))
Chained Index (2023 - 1981): $(round(decline_growth_rate, digits=2)) percentage points
$(repeat("=", 60))
""")
#---------------------------------------------------------------

#------------------------------------------------------------------------------
# Real Consumption Expenditure Indices (NHCES)
#
# Konus expenditure at base-t prices reaching the period-z consumption index:
#     Etilde_{t,z} = E(P_{g,t}, P_{s,t}, C_z)
# The marginal value of wealth does not enter: this is the purely intratemporal
# object, the counterpart of Appendix C in the PIGL code.
#
# Growth-rate correction factor:
#     g_e(t,z)/gD_e(z) = elasE(t,z) / elasE(z,z),
#     elasE(t,z) = dlog E / dlog C at prices P_t and index C_z = eta_g*sg + eta_s*ss
# It equals one when eta_g = eta_s, i.e. under homothetic preferences.
#------------------------------------------------------------------------------
function elasE(t::Int, z::Int)
    C = ct_f[z]
    E = E_nhces(Pgt_f[t], Pst_f[t], C, ωc, σ, ηg, ηs)
    s = sg_nhces(Pgt_f[t], Pst_f[t], C, E, ωc, σ, ηg, ηs)
    return ηg*s + ηs*(1 - s)
end
aux_e_nhces(t::Int, z::Int) = elasE(t, z) / elasE(z, z)

# Chained Divisia for real consumption expenditure. The weighting and cumulation
# follow the chain_weights / chain_cumulate flags set in the Configuration block,
# so this block and the GDP block above share one convention. Setting
# (:endpoint, :trap) reproduces the construction in LV_main.jl and LV_main_HCD.jl.
cgt_f    = Float64.(sim_eq["cgt"])
cst_f    = Float64.(sim_eq["cst"])
sg_w_e   = chain_share(sg_C)
lg_cg    = log.(cgt_f[2:end] ./ cgt_f[1:end-1])
lg_cs    = log.(cst_f[2:end] ./ cst_f[1:end-1])
gD_e_z   = sg_w_e .* lg_cg .+ (1 .- sg_w_e) .* lg_cs
D_e_z    = chain_index([0.0; gD_e_z .+ g_n])

# Fixed-base Fisher-Shell index with base t (aggregate, 1980 = 0).
# The correction factor is averaged over each interval to match the share weights.
aux_mid(t::Int, z::Int) = 0.5 * (aux_e_nhces(t, z-1) + aux_e_nhces(t, z))
FS_e(t::Int) = chain_index([0.0; gD_e_z .* [aux_mid(t, z) for z in 2:T] .+ g_n])

P_e_2023_z = FS_e(2023 - 1980 + 1)
P_e_2010_z = FS_e(2010 - 1980 + 1)
P_e_2000_z = FS_e(2000 - 1980 + 1)
P_e_1990_z = FS_e(1990 - 1980 + 1)
L_e_1980_z = FS_e(1980 - 1980 + 1)

g_e_2023_z = diff(P_e_2023_z)
g_e_1980_z = diff(L_e_1980_z)

### Checks
@assert all(isapprox.([aux_e_nhces(b, b) for b in (1, 11, 21, 31, T)], 1.0; atol=1e-12))
@assert isapprox(elasE(1, 1), ηg*sg_C[1] + ηs*ss_C[1]; rtol=1e-10)
@assert isapprox(E_nhces(Pgt_f[1], Pst_f[1], ct_f[1], ωc, σ, ηg, ηs), et_f[1]; rtol=1e-8)

println("""
$(repeat("=", 60))
Real Consumption Expenditure Indices (NHCES), σ = $σ
$(repeat("=", 60))
  P_e(2023) - D_e : $(round(P_e_2023_z[end] - D_e_z[end], digits=4))
  L_e(1980) - D_e : $(round(L_e_1980_z[end] - D_e_z[end], digits=4))
  spread          : $(round(P_e_2023_z[end] - L_e_1980_z[end], digits=4))
$(repeat("=", 60))
""")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Export the model-implied series for the HCD run
#------------------------------------------------------------------------------
let path = joinpath(datadir, "NHCES_model_series.csv")
    open(path, "w") do io
        println(io, "year,Pg,Ps,sg,e,C")
        for i in 1:T
            println(io, join([1979 + i, Pgt_f[i], Pst_f[i], sg_C[i], et_f[i], ct_f[i]], ","))
        end
    end
    println("Model series written to: $path  (σ = $σ, ηs = $ηs, ωc = $ωc)")
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compensated goods share underlying the fixed-base consumption expenditure indices
#
#   sg_kon(t,z) = goods share at base-t prices and the period-z utility level C_z,
#                 i.e. the share implicit in Etilde_{t,z} = E(P_t, C_z).
#
# Contrast with the intertemporal measure. There the hypothetical share is constant
# in tau, because quasilinearity pins the consumption bundle and the whole adjustment
# to the utility target runs through investment. Here there is no investment margin,
# so attaining the period-z utility level requires moving along the Engel curve and
# the compensated share varies with z even though prices are fixed at t.
#------------------------------------------------------------------------------
function sg_kon_nhces(t::Int, z::Int)
    C = ct_f[z]
    E = E_nhces(Pgt_f[t], Pst_f[t], C, ωc, σ, ηg, ηs)
    return sg_nhces(Pgt_f[t], Pst_f[t], C, E, ωc, σ, ηg, ηs)
end

sg_kon_2023_z = [sg_kon_nhces(T, z)                 for z in 1:T]
sg_kon_2010_z = [sg_kon_nhces(2010 - 1980 + 1, z)   for z in 1:T]
sg_kon_1980_z = [sg_kon_nhces(1, z)                 for z in 1:T]

### Check: at its own base the compensated share equals the actual one
@assert isapprox(sg_kon_nhces(1, 1), sg_C[1]; rtol=1e-10)
@assert isapprox(sg_kon_nhces(T, T), sg_C[T]; rtol=1e-10)

println(repeat("=", 62))
println("Goods share in consumption expenditure: actual vs compensated")
println(repeat("=", 62))
println("  actual       s_g      1980 / 2023 : ",
        round(sg_C[1], digits=4), " / ", round(sg_C[end], digits=4))
println("  2023-base    sg_kon   1980 / 2023 : ",
        round(sg_kon_2023_z[1], digits=4), " / ", round(sg_kon_2023_z[end], digits=4))
println("  1980-base    sg_kon   1980 / 2023 : ",
        round(sg_kon_1980_z[1], digits=4), " / ", round(sg_kon_1980_z[end], digits=4))
println(repeat("=", 62))

#------------------------------------------------------------------------------
# Figure - Goods share in consumption expenditure: actual vs compensated
plot(1980:2023, sg_C,
    label  = L"s_{g,z} \ \ \ \ \ \ \ - \textrm{Actual}",
    ylabel = "Share of Goods in\nConsumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.475, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(1980:2023, sg_kon_2023_z,
    label = L"\hat{s}_{g,2023,z} - \textrm{2023-base \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)

plot!(1980:2023, sg_kon_1980_z,
    label = L"\hat{s}_{g,1980,z} - \textrm{1980-base \ Index}",
    linestyle = :dot, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "NHCES_goods_share_consumption_expenditure.png"))
    println("Figure saved to: ", joinpath(figuresdir, "NHCES_goods_share_consumption_expenditure.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure - Cumulative growth: Chained Divisia vs 1980-base vs 2023-base
plot(1980:2023, D_e_z,
    label  = L"\mathcal{D_{e}}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel = "Cumulative Growth in\nReal Consumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.100, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(1980:2023, P_e_2023_z,
    label = L"\mathcal{P_{e}}_{2023,z} - \textrm{2023-base  \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)

plot!(1980:2023, L_e_1980_z,
    label = L"\mathcal{L_{e}}_{1980,z} - \textrm{1980-base  \ Index}",
    linestyle = :dot, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "NHCES_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "NHCES_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure - Growth rates for several base years vs the chained Divisia
gP_e_2023_z = diff(P_e_2023_z)[2:end]
gQ_e_2010_z = diff(P_e_2010_z)[2:end]
gQ_e_2000_z = diff(P_e_2000_z)[2:end]
gQ_e_1990_z = diff(P_e_1990_z)[2:end]
gL_e_1980_z = diff(L_e_1980_z)[2:end]
gD_e_plot   = diff(D_e_z)[2:end]

plot(1982:2023, gP_e_2023_z,
    label = L"g^{D}_{e}_{2023,z} - \textrm{2023-base \ Index}",
    linestyle = :dash, lw = 2.5, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.400, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(1982:2023, gQ_e_2010_z, label = L"g^{D}_{e}_{2010,z} - \textrm{2010-base  \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)
plot!(1982:2023, gQ_e_2000_z, label = L"g^{D}_{e}_{2000,z} - \textrm{2000-base  \ Index}",
    linestyle = :dash, lw = 1.5, color = :black)
plot!(1982:2023, gQ_e_1990_z, label = L"g^{D}_{e}_{1990,z} - \textrm{1990-base  \ Index}",
    linestyle = :dash, lw = 1.0, color = :black)
plot!(1982:2023, gL_e_1980_z, label = L"g^{D}_{e}_{1980,z} - \textrm{1980-base  \ Index}",
    linestyle = :dash, lw = 0.5, color = :black)
plot!(1982:2023, gD_e_plot,
    label  = L"g^{D}_{e}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel = "Growth Rate of Real \n Consumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "NHCES_growth_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "NHCES_growth_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# σ SWEEP — Real consumption expenditure indices across NHCES calibrations
# Requires nhces_params, i.e. the branch with run_smm_nhces = false.
#------------------------------------------------------------------------------
# Runs only when run_smm = false, since nhces_params is defined in that branch.
if !run_smm

function cons_indices_nhces(σ_i, ηs_i)
    eq = sim_model_full_nhces([σ_i, ηs_i];
             Pg_t=Pg_t, Ps_t=Ps_t, θ=θ, ρ=ρ, δ=δ, ψ=ψ,
             g_calAx=g_calAx, g_As=g_As, g_Ag=g_Ag,
             g_h=g_h, g_l=g_l, g_n=g_n, sg0_data=VAC_GOOD_SHARE[1])
    (!haskey(eq, "ct") || any(x -> x > 1e5, eq["sgt"])) && return nothing

    ωc_i = Float64(eq["ωc"][1])
    Pg   = Float64.(eq["Pgt"]);  Ps = Float64.(eq["Pst"])
    ct   = Float64.(eq["ct"])
    sg   = Float64.(eq["sgt"]);  ss = Float64.(eq["sst"])
    Tl   = length(ct)

    # ∂log E/∂log C = ηg·sg + ηs·ss, evaluated at prices P_t and consumption index C_z
    function ela(t, z)
        E = E_nhces(Pg[t], Ps[t], ct[z], ωc_i, σ_i, ηg, ηs_i)
        s = sg_nhces(Pg[t], Ps[t], ct[z], E, ωc_i, σ_i, ηg, ηs_i)
        return ηg*s + ηs_i*(1 - s)
    end

    cg  = Float64.(eq["cgt"]);  cs = Float64.(eq["cst"])
    smid = chain_share(sg)
    gD  = smid .* log.(cg[2:end]./cg[1:end-1]) .+ (1 .- smid) .* log.(cs[2:end]./cs[1:end-1])
    D   = chain_index([0.0; gD .+ g_n])

    aux(t, z) = 0.5 * (ela(t,z-1)/ela(z-1,z-1) + ela(t,z)/ela(z,z))
    FS(t) = chain_index([0.0; gD .* [aux(t,z) for z in 2:Tl] .+ g_n])

    # HCD limiting case: the base-period bias when the entire share change is
    # attributed to income effects equals Δs_g · Δln(P_g/P_s), exactly.
    Δsg  = sg[end] - sg[1]
    Δlnp = log(Pg[end]/Ps[end]) - log(Pg[1]/Ps[1])

    return (σ = σ_i, ηs = ηs_i, ωc = ωc_i,
            gapP = FS(Tl)[end] - D[end],
            gapL = FS(1)[end]  - D[end],
            Δsg = Δsg, Δlnp = Δlnp, hcd = Δsg*Δlnp)
end

sweep_nhces = [cons_indices_nhces(nhces_params[k].σ, nhces_params[k].ηs)
               for k in sort(collect(keys(nhces_params)))]

println("\n" * "="^112)
println(rpad("σ",7), rpad("ηs",9), rpad("(1-σ)(ηs-ηg)",15),
        rpad("Δs_g",10), rpad("Δln p",10),
        rpad("P_e-D_e",10), rpad("L_e-D_e",10), rpad("spread",10),
        rpad("HCD bound",11), "% of HCD")
println("="^112)
for r in sweep_nhces
    r === nothing && continue
    spread = r.gapP - r.gapL
    println(rpad(round(r.σ,  digits=3), 7),
            rpad(round(r.ηs, digits=4), 9),
            rpad(round((1-r.σ)*(r.ηs-ηg), digits=4), 15),
            rpad(round(r.Δsg,  digits=4), 10),
            rpad(round(r.Δlnp, digits=4), 10),
            rpad(round(r.gapP, digits=4), 10),
            rpad(round(r.gapL, digits=4), 10),
            rpad(round(spread, digits=4), 10),
            rpad(round(r.hcd, digits=4), 11),
            round(100*spread/r.hcd, digits=1))
end
println("="^112)
println("HCD bound = Δs_g · Δln(P_g/P_s), the base-period bias when all structural")
println("transformation is income-driven (σ = 1). Exact, parameter-free.")

# CSV for the paper table
open(joinpath(datadir, "NHCES_sigma_sweep.csv"), "w") do io
    println(io, "sigma,eta_s,income_effect,dsg,dlnp,gapP,gapL,spread,hcd_bound,pct_of_hcd")
    for r in sweep_nhces
        r === nothing && continue
        sp = r.gapP - r.gapL
        println(io, join([r.σ, r.ηs, (1-r.σ)*(r.ηs-ηg), r.Δsg, r.Δlnp,
                          r.gapP, r.gapL, sp, r.hcd, 100*sp/r.hcd], ","))
    end
end
println("Sweep written to: ", joinpath(datadir, "NHCES_sigma_sweep.csv"))

end  # if !run_smm  (σ sweep)
#------------------------------------------------------------------------------
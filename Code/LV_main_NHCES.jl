#------------------------------------------------------------------------------
#Replication code for: Chained Indices Unchained: Structural Transformation and the Welfare Foundations of Income Growth Measurement
#By:                   Omar Licandro and Juan I. Vizcaino
#This Version:         14/04/2026
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
using OrderedCollections, OrderedCollections, LsqFit,  Random, PrettyTables
#using MathJaxRenderer
using LaTeXStrings
using Plots; gr()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Configuration 
# Set to true to save figures, false to only display
save_figures = false

# Set to true to run SMM estimation, false to use pre-estimated parameters
run_smm      = false

# Choose which moment/specification to use for SMM estimation:
#   :baseline  — matches goods consumption expenditure share (s_g) in levels;
#                σ fixed, ωc calibrated to initial share, ηs estimated.
#   :log_ratio — matches log(s_sl / s_gl) using Eq. (7);
#                σ fixed at σ_fixed_alt, ωc and ηs estimated simultaneously.
#   :both      — runs both specifications sequentially.
smm_spec     = :baseline
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Set the Parameter Values for Parameters Taken From HRV
ωx  = 0.650
ωx  = 0.290

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
HRV_data = XLSX.readtable(joinpath(datadir, file_name), sheet_name)
HRV_df   = DataFrame(HRV_data)
#HRV_df  = DataFrame(HRV_data...)

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

function integrate_trap(g::AbstractVector)
    T = length(g)
    out = zeros(eltype(g), T)
    for i in 2:T
        out[i] = out[i-1] + 0.5*(g[i-1] + g[i])
    end
    return out
end
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
# Option A:
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

        # Factor prices
        Rt_t = θ * eq["calAxt"][t] * eq["Kt"][t]^(θ - 1)
        Wt_t = (1 - θ) * eq["calAxt"][t] * eq["Kt"][t]^θ

        push!(eq["Rt"], Rt_t)
        push!(eq["Wt"], Wt_t)

        # Labor allocation
        Ls =
            (eq["Xt"][t] / eq["Yt"][t]) *
            (
                1 /
                (
                    ((eq["Pgt"][t] * eq["Xgt"][t]) /
                     (eq["Pst"][t] * eq["Xst"][t]))
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
        # ABGP step: on the Asymptotic Balanced Growth Path, per-capita
        # consumption expenditure grows at exactly g_abgp each period
        # (budget constraint + constant K/Y ratio).  We impose this
        # directly and recover the composite consumption index C_{t+1}
        # by inverting the NHCES expenditure function.
        # The k0 bisection above (which uses the Euler equation) already
        # pins down the unique initial capital stock consistent with this
        # BGP, so no information is lost.
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
# Alternative simulation: ωc is a FREE parameter (not calibrated to match
# the initial goods share).  C0 is recovered via solve_C_from_E directly.
# Used for the alternative calibration based on equation (7).
#
# params[1] = ωc   (share weight on goods in consumption aggregator)
# params[2] = σ    (elasticity parameter)
# params[3] = ηs   (income effect for services)
function sim_model_full_nhces_free_ωc(params::Vector{Float64};
                                       Pg_t::Vector{Float64},
                                       Ps_t::Vector{Float64},
                                       θ::Float64, ρ::Float64,
                                       δ::Float64, ψ::Float64,
                                       g_calAx::Float64, g_As::Float64,
                                       g_Ag::Float64, g_h::Float64,
                                       g_l::Float64, g_n::Float64)

    ωc = params[1]
    σ  = params[2]
    ηs = params[3]
    ηg = 1.0

    T = length(1980:2023)

    if !(0 < ωc < 1 && σ > 0 && ψ > 0 && ηg > 0 && ηs > 0)
        return Dict("sgt" => fill(1e6, T), "sst" => fill(1e6, T))
    end

    # Exogenous series
    N_t         = 1.000 .* exp.(g_n     .* ((1980:2023) .- 1980))
    h_t         = 1.000 .* exp.(g_h     .* ((1980:2023) .- 1980))
    L_t         = 1.000 .* exp.(g_l     .* ((1980:2023) .- 1980))
    Ag_t        = 1.000 .* exp.(g_Ag    .* ((1980:2023) .- 1980))
    As_t        = 1.000 .* exp.(g_As    .* ((1980:2023) .- 1980))
    calAx_t     = 1.000 .* exp.(g_calAx .* ((1980:2023) .- 1980))
    calAx_hat_t = calAx_t .* (h_t .^ (1 - θ))

    n      = g_n
    g_abgp = g_calAx / (1 - θ) + g_h

    #--------------------------------------------------------------------------
    # Step 1: Find k0 — same bisection as before, but use solve_C_from_E
    #         (ωc is now known, so no need to calibrate it)
    #--------------------------------------------------------------------------

    function E0_from_k0(k0)
        return calAx_hat_t[1] * k0^θ - (g_abgp + δ + n) * k0
    end

    function initial_expenditure_growth_error(k0)
        E0 = E0_from_k0(k0)
        (!isfinite(E0) || E0 <= 0) && return -1e6

        C0 = solve_C_from_E(E0, Pg_t[1], Ps_t[1], ωc, σ, ηg, ηs)
        (!isfinite(C0) || C0 <= 0) && return -1e6

        R0 = θ * calAx_hat_t[1] * k0^(θ - 1)
        C1 = solve_C_next(C0, Pg_t[1], Ps_t[1], Pg_t[2], Ps_t[2], R0,
                          ρ, δ, ψ, ωc, σ, ηg, ηs)
        (!isfinite(C1) || C1 <= 0) && return -1e6

        E1 = E_nhces(Pg_t[2], Ps_t[2], C1, ωc, σ, ηg, ηs)
        (!isfinite(E1) || E1 <= 0) && return -1e6

        return log(E1 / E0) - g_abgp
    end

    k0_min, k0_max = 1e-4, 100.0
    err_min = initial_expenditure_growth_error(k0_min)
    err_max = initial_expenditure_growth_error(k0_max)

    iter_expand = 0
    while sign(err_min) == sign(err_max) && iter_expand < 100
        k0_min /= 2; k0_max *= 2
        err_min = initial_expenditure_growth_error(k0_min)
        err_max = initial_expenditure_growth_error(k0_max)
        iter_expand += 1
    end

    if sign(err_min) == sign(err_max)
        return Dict("sgt" => fill(1e6, T), "sst" => fill(1e6, T))
    end

    k0_mid, err_mid = NaN, Inf
    iter = 0
    while abs(err_mid) > 1e-10 && iter < 500
        k0_mid  = 0.5 * (k0_min + k0_max)
        err_mid = initial_expenditure_growth_error(k0_mid)
        if sign(err_mid) == sign(err_min)
            k0_min = k0_mid; err_min = err_mid
        else
            k0_max = k0_mid; err_max = err_mid
        end
        iter += 1
    end

    k0 = k0_mid

    #--------------------------------------------------------------------------
    # Step 2: E0 and C0
    #--------------------------------------------------------------------------

    e0 = E0_from_k0(k0)
    c0 = solve_C_from_E(e0, Pg_t[1], Ps_t[1], ωc, σ, ηg, ηs)

    if !isfinite(c0) || c0 <= 0
        return Dict("sgt" => fill(1e6, T), "sst" => fill(1e6, T))
    end

    #--------------------------------------------------------------------------
    # Storage (same layout as sim_model_full_nhces)
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
    for var in vars; eq[var] = Vector{Any}(); end

    push!(eq["ωc"],      ωc)
    push!(eq["et"],      e0)
    push!(eq["ct"],      c0)
    push!(eq["kt"],      k0)
    push!(eq["ht"],      h_t[1])
    push!(eq["Agt"],     Ag_t[1])
    push!(eq["Ast"],     As_t[1])
    push!(eq["calAxt"],  calAx_t[1])
    push!(eq["calAxhat"],calAx_hat_t[1])
    push!(eq["Pst"],     Ps_t[1])
    push!(eq["Pgt"],     Pg_t[1])
    push!(eq["Lt"],      L_t[1])
    push!(eq["Nt"],      N_t[1])
    push!(eq["yt"],      y(calAx_t[1], k0, h_t[1], θ))

    #--------------------------------------------------------------------------
    # Step 3: Forward simulation (identical to sim_model_full_nhces)
    #--------------------------------------------------------------------------

    for t in 1:T
        sg_t = sg_nhces(eq["Pgt"][t], eq["Pst"][t], eq["ct"][t], eq["et"][t], ωc, σ, ηg, ηs)
        ss_t = ss_nhces(eq["Pgt"][t], eq["Pst"][t], eq["ct"][t], eq["et"][t], ωc, σ, ηg, ηs)

        if !isfinite(sg_t) || !isfinite(ss_t) ||
           sg_t < 0 || sg_t > 1 || ss_t < 0 || ss_t > 1
            return Dict("sgt" => fill(1e6, T), "sst" => fill(1e6, T))
        end

        push!(eq["sgt"], sg_t)
        push!(eq["sst"], ss_t)

        push!(eq["Yt"], eq["yt"][t] * eq["Nt"][t])
        push!(eq["Kt"], eq["kt"][t] * eq["Nt"][t])
        push!(eq["Et"], eq["et"][t] * eq["Nt"][t])

        X_t  = eq["Yt"][t] - eq["Et"][t]
        push!(eq["Xt"], X_t)

        Xs_t = eq["Xt"][t] / (1 + (ωx / (1 - ωx)) * (eq["Agt"][t] / eq["Ast"][t])^εx)
        Xg_t = X_t - Xs_t
        push!(eq["Xst"], Xs_t)
        push!(eq["Xgt"], Xg_t)

        cs_t = (eq["sst"][t] * eq["et"][t]) / eq["Pst"][t]
        cg_t = (eq["sgt"][t] * eq["et"][t]) / eq["Pgt"][t]
        push!(eq["cst"], cs_t)
        push!(eq["cgt"], cg_t)
        push!(eq["Cst"], cs_t * eq["Nt"][t])
        push!(eq["Cgt"], cg_t * eq["Nt"][t])

        Rt_t = θ * eq["calAxt"][t] * eq["Kt"][t]^(θ - 1)
        Wt_t = (1 - θ) * eq["calAxt"][t] * eq["Kt"][t]^θ
        push!(eq["Rt"], Rt_t)
        push!(eq["Wt"], Wt_t)

        Ls = (eq["Xt"][t] / eq["Yt"][t]) *
             (1 / (((eq["Pgt"][t] * eq["Xgt"][t]) / (eq["Pst"][t] * eq["Xst"][t])) + 1)) +
             (eq["Et"][t] / eq["Yt"][t]) *
             (1 / (((eq["Pgt"][t] * eq["Cgt"][t]) / (eq["Pst"][t] * eq["Cst"][t])) + 1))

        if !isfinite(Ls) || Ls < 0 || Ls > 1
            return Dict("sgt" => fill(1e6, T), "sst" => fill(1e6, T))
        end

        Lg = 1 - Ls
        push!(eq["Lst"], Ls)
        push!(eq["Lgt"], Lg)

        Ks = eq["Kt"][t] * Ls;  Kg = eq["Kt"][t] * Lg
        push!(eq["Kst"], Ks);   push!(eq["Kgt"], Kg)

        Ys = eq["Ast"][t] * (Ks^θ) * (Ls^(1 - θ))
        Yg = eq["Agt"][t] * (Kg^θ) * (Lg^(1 - θ))
        push!(eq["Yst"], Ys);   push!(eq["Ygt"], Yg)

        t == T && break

        push!(eq["Pst"],      Ps_t[t+1])
        push!(eq["Pgt"],      Pg_t[t+1])
        push!(eq["calAxt"],   calAx_t[t+1])
        push!(eq["calAxhat"], calAx_hat_t[t+1])
        push!(eq["Agt"],      Ag_t[t+1])
        push!(eq["Ast"],      As_t[t+1])
        push!(eq["ht"],       h_t[t+1])
        push!(eq["Lt"],       L_t[t+1])
        push!(eq["Nt"],       N_t[t+1])

        e_next = eq["et"][t] * exp(g_abgp)
        c_next = solve_C_from_E(e_next, Pg_t[t+1], Ps_t[t+1], ωc, σ, ηg, ηs)

        if !isfinite(c_next) || !isfinite(e_next) || c_next <= 0 || e_next <= 0
            return Dict("sgt" => fill(1e6, T), "sst" => fill(1e6, T))
        end

        push!(eq["ct"], c_next)
        push!(eq["et"], e_next)

        g_k = dkdt(eq["calAxhat"][t], eq["kt"][t], eq["et"][t], δ, n, θ) / eq["kt"][t]
        if !isfinite(g_k) || abs(g_k) > 1.0
            return Dict("sgt" => fill(1e6, T), "sst" => fill(1e6, T))
        end

        k_next = eq["kt"][t] * exp(g_k)
        push!(eq["kt"], k_next)
        push!(eq["yt"], y(calAx_t[t+1], k_next, h_t[t+1], θ))
    end

    return eq
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Diagnostic helper: run the simulation step-by-step with verbose output.
# Call from the REPL as:  debug_sim_nhces([σ_val, ηs_val])
function debug_sim_nhces(params::Vector{Float64};
                         Pg_t=Pg_t, Ps_t=Ps_t, θ=θ, ρ=ρ, δ=δ, ψ=ψ,
                         g_calAx=g_calAx, g_As=g_As, g_Ag=g_Ag,
                         g_h=g_h, g_l=g_l, g_n=g_n,
                         sg0_data=VAC_GOOD_SHARE[1])

    σ  = params[1]
    ηs = params[2]
    ηg = 1.0

    println("─── debug_sim_nhces ───────────────────────────────────")
    println("  σ=$(round(σ,digits=6))  ηs=$(round(ηs,digits=6))  ηg=$(ηg)")

    T         = length(1980:2023)
    h_t       = 1.0 .* exp.(g_h     .* ((1980:2023) .- 1980))
    calAx_t   = 1.0 .* exp.(g_calAx .* ((1980:2023) .- 1980))
    calAx_hat = calAx_t .* (h_t .^ (1 - θ))
    n         = g_n
    g_abgp    = g_calAx / (1 - θ) + g_h

    E0_from_k0(k0) = calAx_hat[1] * k0^θ - (g_abgp + δ + n) * k0

    function init_growth_err(k0)
        E0 = E0_from_k0(k0)
        (!isfinite(E0) || E0 <= 0) && return -1e6
        C0, ωc_l = calibrate_C0_ωc_from_initial_share(E0, Pg_t[1], Ps_t[1], sg0_data, σ, ηg, ηs)
        !(0 < ωc_l < 1) && return -1e6
        R0  = θ * calAx_hat[1] * k0^(θ - 1)
        C1  = solve_C_next(C0, Pg_t[1], Ps_t[1], Pg_t[2], Ps_t[2], R0, ρ, δ, ψ, ωc_l, σ, ηg, ηs)
        (!isfinite(C1) || C1 <= 0) && return -1e6
        E1  = E_nhces(Pg_t[2], Ps_t[2], C1, ωc_l, σ, ηg, ηs)
        (!isfinite(E1) || E1 <= 0) && return -1e6
        return log(E1 / E0) - g_abgp
    end

    # --- Step 1: bracket k0 ---
    k0_min, k0_max = 1e-4, 100.0
    err_lo, err_hi = init_growth_err(k0_min), init_growth_err(k0_max)
    iexp = 0
    while sign(err_lo) == sign(err_hi) && iexp < 100
        k0_min /= 2; k0_max *= 2
        err_lo = init_growth_err(k0_min); err_hi = init_growth_err(k0_max)
        iexp += 1
    end
    if sign(err_lo) == sign(err_hi)
        println("  FAIL  k0 bracketing failed after $iexp expansions.")
        println("        err_lo=$(round(err_lo,digits=6))  err_hi=$(round(err_hi,digits=6))")
        return nothing
    end
    println("  OK    k0 bracketed in [$k0_min, $k0_max] after $iexp expansions")

    # --- Step 2: bisect k0 ---
    k0_mid, err_mid = NaN, Inf
    for _ in 1:500
        k0_mid = 0.5*(k0_min + k0_max)
        err_mid = init_growth_err(k0_mid)
        abs(err_mid) < 1e-10 && break
        sign(err_mid) == sign(err_lo) ? (k0_min = k0_mid; err_lo = err_mid) :
                                        (k0_max = k0_mid; err_hi = err_mid)
    end
    println("  OK    k0=$(round(k0_mid,digits=8))  residual=$(round(err_mid,sigdigits=4))")

    # --- Step 3: calibrate ωc ---
    e0 = E0_from_k0(k0_mid)
    c0, ωc = calibrate_C0_ωc_from_initial_share(e0, Pg_t[1], Ps_t[1], sg0_data, σ, ηg, ηs)
    if !(0 < ωc < 1)
        println("  FAIL  ωc=$(round(ωc,sigdigits=6)) is outside (0,1).")
        return nothing
    end
    println("  OK    e0=$(round(e0,digits=6))  c0=$(round(c0,digits=6))  ωc=$(round(ωc,digits=6))")

    # --- Step 4: first few periods ---
    println("  Simulating $(T) periods...")
    Ag_t    = 1.0 .* exp.(g_Ag .* ((1980:2023) .- 1980))
    As_t    = 1.0 .* exp.(g_As .* ((1980:2023) .- 1980))
    ct, et, kt  = c0, e0, k0_mid
    ok = true
    for t in 1:T
        sg_t = sg_nhces(Pg_t[t], Ps_t[t], ct, et, ωc, σ, ηg, ηs)
        ss_t = ss_nhces(Pg_t[t], Ps_t[t], ct, et, ωc, σ, ηg, ηs)
        if !isfinite(sg_t) || !isfinite(ss_t) || sg_t < 0 || sg_t > 1 || ss_t < 0 || ss_t > 1
            println("  FAIL  t=$t: shares out of range  sg=$(round(sg_t,sigdigits=4))  ss=$(round(ss_t,sigdigits=4))")
            ok = false; break
        end
        t == T && break
        Rt = θ * calAx_hat[t] * kt^(θ - 1)
        c_next = solve_C_next(ct, Pg_t[t], Ps_t[t], Pg_t[t+1], Ps_t[t+1], Rt, ρ, δ, ψ, ωc, σ, ηg, ηs)
        e_next = E_nhces(Pg_t[t+1], Ps_t[t+1], c_next, ωc, σ, ηg, ηs)
        if !isfinite(c_next) || !isfinite(e_next) || c_next <= 0 || e_next <= 0
            println("  FAIL  t=$t: c_next=$(c_next)  e_next=$(e_next)")
            ok = false; break
        end
        g_k = dkdt(calAx_hat[t], kt, et, δ, n, θ) / kt
        kt  = kt * exp(g_k)
        ct, et = c_next, e_next
    end
    if ok
        sg1   = sg_nhces(Pg_t[1],   Ps_t[1],   c0, e0, ωc, σ, ηg, ηs)
        sgend = sg_nhces(Pg_t[end], Ps_t[end], ct, et, ωc, σ, ηg, ηs)
        ss1   = ss_nhces(Pg_t[1],   Ps_t[1],   c0, e0, ωc, σ, ηg, ηs)
        ssend = ss_nhces(Pg_t[end], Ps_t[end], ct, et, ωc, σ, ηg, ηs)
        println("  OK    all $T periods completed")
        println("        sg[1]=$(round(sg1,digits=4))  sg[end]=$(round(sgend,digits=4))  data[1]=$(round(VAC_GOOD_SHARE[1],digits=4))  data[end]=$(round(VAC_GOOD_SHARE[end],digits=4))")
        println("        ss[1]=$(round(ss1,digits=4))  ss[end]=$(round(ssend,digits=4))")
        println("  NOTE: both optimizer bounds hit (σ at lower, ηs at upper) → consider expanding ParamSpace.")
    end
    println("───────────────────────────────────────────────────────")
    return nothing
end
#------------------------------------------------------------------------------

run_smm = false
#------------------------------------------------------------------------------

if run_smm

    if smm_spec ∈ (:baseline, :both)
    println("Running SMM Estimation...")
    println("Estimating ηs only; fixing σ = 0.05, ηg = 1, and calibrating ωc to match the 1980 goods share.")

    σ_fixed  = 0.00001
    ηg       = 1.00000

    ### Initial goods expenditure share in the data
    sg0_data = VAC_GOOD_SHARE[1]

    ### Pass fixed parameters to the simulation function
    ### params[1] = ηs  (σ is fixed at σ_fixed)
    sim_model(params) =
        sim_model_full_nhces(
            [σ_fixed , params[1]];
            Pg_t=Pg_t,
            Ps_t=Ps_t,
            θ=θ,
            ρ=ρ,
            δ=δ,
            ψ=ψ,
            g_calAx=g_calAx,
            g_As=g_As,
            g_Ag=g_Ag,
            g_h=g_h,
            g_l=g_l,
            g_n=g_n,
            sg0_data=sg0_data
        )

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
    ParamSpace = [
        (1.00001, 3.000)    # ηs
    ]

    ### Seeds for repeated global optimization
    seeds = [1234]
    optimization_results = OrderedDict()

    ### Quick diagnostic at the centre of the parameter space before optimizing
    println("\nDiagnostic run at parameter-space midpoint:")
    debug_sim_nhces([0.5*(0.050+5.000), 0.5*(1.001+5.000)])

    ### Run optimization for different seeds
    for seed in seeds
        println("\n" * "="^80)
        println("Running optimization with seed: $seed")
        println("="^80)

        opt_problem =
            bbsetup(
                params ->
                    SSE(
                        params;
                        cons_goods_share_data=VAC_GOOD_SHARE
                    );
                SearchRange=ParamSpace,
                TraceMode=:compact,
                Method=:adaptive_de_rand_1_bin,
                PopulationSize=100,
                MaxFuncEvals=50_000,
                FitnessTolerance=1e-20,
                MaxStepsWithoutProgress=500_000,
                rng=MersenneTwister(seed)
            )

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

    end # smm_spec ∈ (:baseline, :both)

    #------------------------------------------------------------------------------

    if smm_spec ∈ (:log_ratio, :both)
    # ALTERNATIVE CALIBRATION  (Equation 7 — log share-ratio moments)
    #
    # Matches the log ratio of services-to-goods consumption expenditure shares
    # across all periods, using the analytical NHCES relationship:
    #
    #   log(s_sl^C / s_gl^C) = log((1-ωc)/ωc)
    #                        + (1-σ) log(P_sl / P_gl)
    #                        + (1-σ)(ηs - ηg) log(C_t)
    #
    # Like the baseline calibration, σ is fixed at σ_fixed.  The approach
    # optimises over ωc and ηs simultaneously; ωc is a free parameter (not
    # calibrated to match the initial share) and C_t is taken from the model.
    #------------------------------------------------------------------------------

    σ_fixed_alt = 0.0001

    println("\n" * "="^80)
    println("Running Alternative SMM Estimation (Eq. 7 – Log Share-Ratio Moments)...")
    println("Estimating ωc and ηs simultaneously; σ fixed at $(σ_fixed_alt).")
    println("="^80)

    ### Simulation wrapper for the alternative calibration
    ### params[1] = ωc, params[2] = ηs  (σ fixed at σ_fixed_alt)
    sim_model_alt(params) =
        sim_model_full_nhces_free_ωc(
            [params[1], σ_fixed_alt, params[2]];  # [ωc, σ, ηs]
            Pg_t=Pg_t, Ps_t=Ps_t,
            θ=θ, ρ=ρ, δ=δ, ψ=ψ,
            g_calAx=g_calAx, g_As=g_As, g_Ag=g_Ag,
            g_h=g_h, g_l=g_l, g_n=g_n
        )

    ### SSE based on equation (7): log(s_sl / s_gl) in model vs. data
    data_log_ratio_alt = log.(VAC_SERV_SHARE ./ VAC_GOOD_SHARE)

    function SSE_log_ratio(params)
        sim_eq_alt = sim_model_alt(params)

        T_alt = length(VAC_SERV_SHARE)

        if !haskey(sim_eq_alt, "sgt") ||
           length(sim_eq_alt["sgt"]) != T_alt ||
           any(x -> !isfinite(x) || x <= 0.0 || x >= 1.0, sim_eq_alt["sgt"]) ||
           any(x -> !isfinite(x) || x <= 0.0 || x >= 1.0, sim_eq_alt["sst"])
            return 1e10
        end

        model_log_ratio = log.(Float64.(sim_eq_alt["sst"]) ./ Float64.(sim_eq_alt["sgt"]))

        if any(!isfinite, model_log_ratio)
            return 1e10
        end

        return sum((model_log_ratio .- data_log_ratio_alt) .^ 2)
    end

    ### Parameter space for alternative calibration
    ### params[1] = ωc ∈ (0,1), params[2] = ηs > 1  (σ fixed at σ_fixed)
    ParamSpace_alt = [
        (0.05,    0.95),     # ωc
        (1.00001, 5.000)     # ηs
    ]

    seeds_alt = [1234]
    optimization_results_alt = OrderedDict()

    for seed in seeds_alt
        println("\n" * "="^80)
        println("Alt. calibration — seed: $seed")
        println("="^80)

        opt_problem_alt =
            bbsetup(
                params -> SSE_log_ratio(params);
                SearchRange=ParamSpace_alt,
                TraceMode=:compact,
                Method=:adaptive_de_rand_1_bin,
                PopulationSize=100,
                MaxFuncEvals=50_000,
                FitnessTolerance=1e-20,
                MaxStepsWithoutProgress=500_000,
                rng=MersenneTwister(seed)
            )

        opt_results_alt = bboptimize(opt_problem_alt)

        best_params_alt = best_candidate(opt_results_alt)

        optimization_results_alt[seed] = OrderedDict(
            "params"  => best_params_alt,
            "fitness" => best_fitness(opt_results_alt),
            "seed"    => seed,
            "ωc"      => best_params_alt[1],
            "ηs"      => best_params_alt[2]
        )

        println("Seed $seed: Fitness = $(best_fitness(opt_results_alt))")
        println(
            "Parameters: " *
            "ωc=$(round(best_params_alt[1], digits=6))  " *
            "ηs=$(round(best_params_alt[2], digits=6))"
        )
    end

    ### Best result across seeds
    best_seed_index_alt =
        argmin([optimization_results_alt[s]["fitness"] for s in seeds_alt])

    best_result_alt = optimization_results_alt[seeds_alt[best_seed_index_alt]]

    ωc_alt = round(best_result_alt["ωc"], digits=5)
    σ_alt  = σ_fixed_alt
    ηs_alt = round(best_result_alt["ηs"], digits=5)
    ηg_alt = 1.0

    println("""
    $(repeat("=", 80))
    ALTERNATIVE CALIBRATION — BEST RESULT (Eq. 7 Log Share-Ratio)
    $(repeat("=", 80))
    Best seed:    $(best_result_alt["seed"])
    Best fitness: $(best_result_alt["fitness"])

    Estimated parameters (ωc and ηs simultaneously):
    ωc = $(ωc_alt)
    ηs = $(ηs_alt)

    Fixed:
    σ  = $(σ_alt)
    ηg = 1.000000
    $(repeat("=", 80))
    """)
    #--------------------------------------------------------------------------

    end # smm_spec ∈ (:log_ratio, :both)

else
    println("Skipping SMM estimation, using pre-estimated NHCES parameters...")

    ηg = 1.0

    # TODO: update σ and ηs after running the optimizer (run_smm = true)
    ωc = 0.242602      # pre-estimated; 
    σ  = 0.0001        # pre-estimated; 
    ηs = 1.182618      # pre-estimated; 

    # ωc is calibrated internally by sim_model_full_nhces to match VAC_GOOD_SHARE[1].
end
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Simulate the Model with the Calibrated Parameters
sim_eq = sim_model_full_nhces([σ, ηs];
    Pg_t=Pg_t, Ps_t=Ps_t,
    θ=θ, ρ=ρ, δ=δ, ψ=ψ,
    g_calAx=g_calAx, g_As=g_As, g_Ag=g_Ag,
    g_h=g_h, g_l=g_l, g_n=g_n,
    sg0_data=VAC_GOOD_SHARE[1]
)

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
Et_f    = Float64.(sim_eq["Et"])
Pgt_f   = Float64.(sim_eq["Pgt"])
Pst_f   = Float64.(sim_eq["Pst"])

gk      = log.(kt_f[2:end])  .- log.(kt_f[1:end-1])
ge      = log.(et_f[2:end])  .- log.(et_f[1:end-1])
gY      = log.(Yt_f[2:end])  .- log.(Yt_f[1:end-1])
gE      = log.(Et_f[2:end])  .- log.(Et_f[1:end-1])

g_pg    = log.(Pgt_f[2:end]) .- log.(Pgt_f[1:end-1])
g_ps    = log.(Pst_f[2:end]) .- log.(Pst_f[1:end-1])


println("""
$(repeat("=", 60))
NHCES ABGP (ABGP = $(round(g_abgp, digits=5)))
$(repeat("=", 60))
         t=1980   t=2023
  gk :   $(round(gk[1],  digits=5))   $(round(gk[end],  digits=5))
  ge :   $(round(ge[1],  digits=5))   $(round(ge[end],  digits=5))
  gY :   $(round(gY[1],  digits=5))   $(round(gY[end],  digits=5))
  gE :   $(round(gE[1],  digits=5))   $(round(gE[end],  digits=5))
$(repeat("=", 60))
""")
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
    savefig(joinpath(figuresdir, "sg_t_Model_Fit.png"))
    println("Figure saved to: ", joinpath(figuresdir, "sg_t_Model_Fit.png"))
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

gD         = s_e[2:end].*( sim_eq["sgt"][2:end].*g_cg .+ sim_eq["sst"][2:end].*g_cs  ) .+ (1 .- s_e[2:end]).*g_x
gD_agg     = [ 0 ; (gD.+ g_n) ]
FS         = integrate_trap(gD_agg)

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
    savefig(joinpath(figuresdir,"GDP_Model_vs_Data.png"))
    println("Figure saved to: ", joinpath(figuresdir, "GDP_Model_vs_Data.png"))
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
    savefig(joinpath(figuresdir,"GDP_Growth_Model_vs_Data.png"))
    println("Figure saved to: ", joinpath(figuresdir, "GDP_Growth_Model_vs_Data.png"))
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

#Gross Income Per Capita  (m_t = y_t)
m_t  = Float64.(sim_eq["yt"])

# Float64 income array for EV computations
yt_f  = Float64.(sim_eq["yt"])
# Consumption expenditure per capita: E_t = s_e · y_t
# Only E_t (not savings) enters utility, so the EV formula uses E_τ as the reference.
s_e_t = Float64.(sim_eq["Et"] ./ sim_eq["Yt"])
Et_f  = s_e_t .* yt_f

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
# NHCES Fisher-Shell Indices — EV definition (eq. 50), Object A
# Fix t = τ at the base year; let z roll from 1980 to 2023:
#   m̂_{t,z,t} = y_t + (Φ_t(P_t) − Φ_t(P_z)) / ν_t,  t=τ=base year, z∈1980:2023
# (Object A: rolling reference-price year z, fixed income and preferences at base year t)
#
# At z = base year the price correction vanishes: m̂_{t,t,t} = y_t exactly.
# m̂ is always positive in this calibration.
#
# FS index: P_{t,z} = log(m̂_{t,z,t}) − log(m̂_{t,z_0,t})  where z_0 = 1980
#------------------------------------------------------------------------------

# Base year t0 = 2023 — Object A: m̂_{t,z,t}, z rolling (t=τ=2023)
t_base           = 2023
t_prime          = t_base - 1980 + 1
mhat_2023_2023_τ = [_mhat(t_prime, t_prime,zi) for zi in τ]
@assert isapprox(mhat_2023_2023_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_2023_2023_τ   = log.(mhat_2023_2023_τ) .- log.(mhat_2023_2023_τ[1])
g_2023_2023_τ    = diff(FS_2023_2023_τ) 

# Base year t0 = 2010 — Object A: m̂_{t,z,t}, z rolling (t=τ=2010)
t_base           = 2010
t_prime          = t_base - 1980 + 1
mhat_2010_2010_τ = [_mhat(t_prime, t_prime,zi) for zi in τ]
@assert isapprox(mhat_2010_2010_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_2010_2010_τ   = log.(mhat_2010_2010_τ) .- log.(mhat_2010_2010_τ[1]) 
g_2010_2010_τ    = diff(FS_2010_2010_τ) 

# Base year t0 = 2000 — Object A: m̂_{t,z,t}, z rolling (t=τ=2000)
t_base           = 2000
t_prime          = t_base - 1980 + 1
mhat_2000_2000_τ = [_mhat(t_prime, t_prime, zi ) for zi in τ]
@assert isapprox(mhat_2000_2000_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_2000_2000_τ   = log.(mhat_2000_2000_τ) .- log.(mhat_2000_2000_τ[1]) 

# Base year t0 = 1990 — Object A: m̂_{t,z,t}, z rolling (t=τ=1990)
t_base           = 1990
t_prime          = t_base - 1980 + 1
mhat_1990_1990_τ = [_mhat(t_prime, t_prime, zi) for zi in τ]
@assert isapprox(mhat_1990_1990_τ[t_prime], yt_f[t_prime]; rtol=1e-10)
FS_1990_1990_τ   = log.(mhat_1990_1990_τ) .- log.(mhat_1990_1990_τ[1]) 

# Base year t0 = 1980 — Object A: m̂_{t,z,t}, z rolling (t=τ=1980)
t_base           = 1980
t_prime          = t_base - 1980 + 1
mhat_1980_1980_τ = [_mhat(t_prime, t_prime,zi) for zi in τ]
@assert isapprox(mhat_1980_1980_τ[1], yt_f[t_prime]; rtol=1e-10)    # at z=t=1980: m̂=y_1980
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
    savefig(joinpath(figuresdir,"FS_BBEV.png"))
    println("Figure saved to: ", joinpath(figuresdir, "FS_BBEV.png"))
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
# Alternative Expenditure Shares (NHCES EV)
# Ê_{t,z}: hypothetical expenditure at preferences t, prices of period z
z_prime            = (1980:1:2023) .- 1980 .+ 1

# At Preferences of 1980 (t=1, z=1): Ê_{1980,1980}
Ehat_1980_1980     = _Ehat(1, 1)
se_1980_z          = Ehat_1980_1980 ./ yt_f[z_prime]

# At Preferences of 2023 (t=T, z=T): Ê_{2023,2023}
Ehat_2023_2023     = _Ehat(T, T)
se_2023_z          = Ehat_2023_2023 ./ yt_f[z_prime]

se                 = sim_eq["et"]./sim_eq["yt"]
   
plot(1980:2023,se,
    label=L"s_{e,z}",
    ylabel="Share of Consumption Expenditure\nin Gross Income",
    linestyle=:solid, lw=2,
    xticks=1980:5:2025, yticks=0.0:0.50:5.00,
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.85, 0.92),
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box) 

plot!(1980:2023,se_1980_z,label=L"s_{e,1980,z}",
    linestyle=:dot, lw=2,color=:black)

plot!(1980:2023,se_2023_z,label=L"s_{e,2023,z}",
    linestyle=:dash, lw=2,color=:black)

if save_figures
    savefig(joinpath(figuresdir,"se_t.png"))
    println("Figure saved to: ", joinpath(figuresdir, "se_t.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 6 (b)

#Goods Expenditure Share: NHCES EV goods share ŝ_{g,t,z}
# At preferences=prices=2023
sg_2023_z  = _sg_hat(T, T)
# At preferences=prices=1980
sg_1980_z  = _sg_hat(1, 1)

plot(1980:2023, sim_eq["sgt"], label=L"s_{g,z}",linestyle=:solid, lw=2,color=:black)
plot!(1980:2023, fill(sg_1980_z, length(1980:2023)), label=L"s_{g,1980,z}",
    linestyle=:dot, lw=2,color=:black)

plot!(1980:2023, fill(sg_2023_z, length(1980:2023)), label=L"s_{g,2023,z}",
    ylabel="Share of Goods in Consumption\nExpenditure",
    linestyle=:dash, lw=2,
    xticks=1980:5:2025, yticks=0.000:0.05:0.300,
    ylims=(0.000, 0.300),
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.85, 0.65),
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir, "sg_t_z.png"))
    println("Figure saved to: ", joinpath(figuresdir, "sg_t_z.png"))
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
    savefig(joinpath(figuresdir, "FS_GrowthRates_1980_2023.png"))
    println("Figure saved to: ", joinpath(figuresdir, "FS_GrowthRates_1980_2023.png"))
end

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



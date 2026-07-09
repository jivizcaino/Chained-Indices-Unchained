# ======================================================================
# LV_CIU_HCD.jl
# Heterothetic Cobb-Douglas Preferences
#
#
# Purpose:
#   Compute welfare-based real income / quantity indices under
#   Heterothetic Cobb-Douglas (HCD) preferences.
#
# Core HCD objects:
#   log E(p_g, p_s, Ω) = log Ω + m_g(Ω) log p_g + [1 - m_g(Ω)] log p_s
#
#   log Ω_t = log y_t - m_{g,t} log p_{g,t} - (1 - m_{g,t}) log p_{s,t}
#
#   where m_g(Ω) is recovered non-parametrically from observed
#   pairs (log Ω_t, m_{g,t}).
#
# Notes:
#   - This first version follows the structure of LV_main_NHCES.jl.
#   - HCD has no SMM preference-parameter estimation block: the Engel
#     curve m_g(Ω) is recovered non-parametrically from the observed
#     consumption goods-share path.
# ======================================================================

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Uncomment only if running outside the project environment for the first time.
Pkg.add(["CSV", "DataFrames", "Statistics", "Interpolations", "Plots", "XLSX",
          "BlackBoxOptim", "OrderedCollections", "LsqFit", "Random",
          "PrettyTables", "LaTeXStrings"])

using CSV
using Statistics
using Interpolations
using Plots

using XLSX, DataFrames, BlackBoxOptim, Statistics
using OrderedCollections, LsqFit, Random, PrettyTables, LaTeXStrings, Plots; gr()

# ======================================================================
# 1. User settings
# ======================================================================

currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))

#------------------------------------------------------------------------------
#Import the HRV Data
file_name  = "HRV2021_Data.xlsx"
sheet_name = "Data"
data_range = "A1:DH78"

# Read the data from the Excel file
HRV_data   = XLSX.readtable(joinpath(datadir, file_name), sheet_name)
HRV_df     = DataFrame(HRV_data)

# Output paths
currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Configuration
# Set to true to save figures, false to only display
save_figures = false
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
#Import the HRV Data Objects
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
# HCD Preference Block
#------------------------------------------------------------------------------

function clamp_share(x; eps_share=1e-8)
    return min(max(x, eps_share), 1.0 - eps_share)
end

function natural_cubic_spline_functions(x_raw::AbstractVector, y_raw::AbstractVector)
    df = DataFrame(x = Float64.(collect(x_raw)), y = Float64.(collect(y_raw)))
    sort!(df, :x)

    # If duplicate x values arise numerically, average their y values.
    grouped = combine(groupby(df, :x), :y => mean => :y)

    x = Float64.(grouped.x)
    y = Float64.(grouped.y)

    n = length(x)

    if n < 4
        error("Need at least four distinct observations for the HCD Engel-curve spline.")
    end

    if any(diff(x) .<= 0)
        error("Spline x-grid must be strictly increasing.")
    end

    h = diff(x)

    α = zeros(Float64, n)
    for i in 2:n-1
        α[i] =
            (3.0 / h[i])     * (y[i+1] - y[i]) -
            (3.0 / h[i-1])   * (y[i]   - y[i-1])
    end

    l  = ones(Float64, n)
    μ  = zeros(Float64, n)
    z  = zeros(Float64, n)

    for i in 2:n-1
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1] * μ[i-1]
        μ[i] = h[i] / l[i]
        z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
    end

    c = zeros(Float64, n)
    b = zeros(Float64, n-1)
    d = zeros(Float64, n-1)
    a = y[1:n-1]

    for j in (n-1):-1:1
        c[j] = z[j] - μ[j] * c[j+1]
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2.0*c[j]) / 3.0
        d[j] = (c[j+1] - c[j]) / (3.0*h[j])
    end

    function interval_index(x0)
        if x0 <= x[1]
            return 1
        elseif x0 >= x[end]
            return n - 1
        else
            j = searchsortedlast(x, x0)
            return min(max(j, 1), n - 1)
        end
    end

    function spline_eval(x0)
        # Flat extrapolation outside the observed support.
        if x0 <= x[1]
            return y[1]
        elseif x0 >= x[end]
            return y[end]
        end

        j  = interval_index(x0)
        dx = x0 - x[j]

        return a[j] + b[j]*dx + c[j]*dx^2 + d[j]*dx^3
    end

    function spline_deriv(x0)
        # Flat extrapolation outside the observed support.
        if x0 <= x[1] || x0 >= x[end]
            return 0.0
        end

        j  = interval_index(x0)
        dx = x0 - x[j]

        return b[j] + 2.0*c[j]*dx + 3.0*d[j]*dx^2
    end

    return spline_eval, spline_deriv
end

function make_mg_functions(logΩ::AbstractVector, mg::AbstractVector)
    # Smooth non-parametric Engel curve: natural cubic spline through
    # observed pairs (log Ω_t, m_{g,t}), with flat extrapolation.
    spline_eval, spline_deriv = natural_cubic_spline_functions(logΩ, mg)

    mg_hat(z)  = clamp_share(spline_eval(z))
    dmg_hat(z) = spline_deriv(z)

    return mg_hat, dmg_hat
end

function logE_hcd(Pg, Ps, Ω, mg_hat)
    z  = log(Ω)
    mg = mg_hat(z)
    return z + mg*log(Pg) + (1.0 - mg)*log(Ps)
end

function E_hcd(Pg, Ps, Ω, mg_hat)
    return exp(logE_hcd(Pg, Ps, Ω, mg_hat))
end

function sg_hcd(Pg, Ps, Ω, mg_hat)
    return mg_hat(log(Ω))
end

function ss_hcd(Pg, Ps, Ω, mg_hat)
    return 1.0 - sg_hcd(Pg, Ps, Ω, mg_hat)
end

function dEdΩ_hcd(Pg, Ps, Ω, mg_hat, dmg_hat)
    # If z = log Ω, then
    #
    #   log E = z + m_g(z) log Pg + [1-m_g(z)] log Ps
    #
    # so
    #
    #   ∂E/∂Ω = (E/Ω) * [1 + m_g'(z) * (log Pg - log Ps)].
    z          = log(Ω)
    E_val      = E_hcd(Pg, Ps, Ω, mg_hat)
    mg_prime_z = dmg_hat(z)
    elasticity = 1.0 + mg_prime_z*(log(Pg) - log(Ps))

    if !isfinite(elasticity) || elasticity <= 0
        return NaN
    end

    return (E_val / Ω) * elasticity
end

function U_crra(Ω, ψ)
    if isapprox(ψ, 1.0)
        return log(Ω)
    else
        return (Ω^(1.0 - ψ) - 1.0) / (1.0 - ψ)
    end
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Solve Ω from E = e(P, Ω)
#------------------------------------------------------------------------------
function solve_Ω_from_E_hcd(E_target, Pg, Ps, mg_hat;
                            z_low=-50.0, z_high=50.0,
                            tol=1e-12, maxiter=500)

    if !(E_target > 0 && Pg > 0 && Ps > 0)
        return NaN
    end

    f(z) = logE_hcd(Pg, Ps, exp(z), mg_hat) - log(E_target)

    lo = z_low
    hi = z_high
    f_lo = f(lo)
    f_hi = f(hi)

    expand_iter = 0
    while sign(f_lo) == sign(f_hi) && expand_iter < 100
        lo -= 10.0
        hi += 10.0
        f_lo = f(lo)
        f_hi = f(hi)
        expand_iter += 1
    end

    if sign(f_lo) == sign(f_hi)
        return NaN
    end

    mid = NaN
    f_mid = Inf

    iter = 0
    while abs(f_mid) > tol && iter < maxiter
        mid = 0.5*(lo + hi)
        f_mid = f(mid)

        if sign(f_mid) == sign(f_lo)
            lo = mid
            f_lo = f_mid
        else
            hi = mid
            f_hi = f_mid
        end

        iter += 1
    end

    return exp(mid)
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Solve Ω_{t+1} from the HCD Euler Equation
#------------------------------------------------------------------------------
function solve_Ω_next_hcd(Ω_now, Pg_now, Ps_now, Pg_next, Ps_next, R_now,
                          ρ, δ, ψ, mg_hat, dmg_hat;
                          z_width=20.0, tol=1e-12, maxiter=500)

    M_now = dEdΩ_hcd(Pg_now, Ps_now, Ω_now, mg_hat, dmg_hat)

    if !isfinite(M_now) || M_now <= 0
        return NaN
    end

    z_now = log(Ω_now)

    function euler_residual_z(z_next)
        Ω_next = exp(z_next)
        M_next = dEdΩ_hcd(Pg_next, Ps_next, Ω_next, mg_hat, dmg_hat)

        if !isfinite(M_next) || M_next <= 0
            return NaN
        end

        lhs = ψ*(z_next - z_now) + log(M_next / M_now)
        rhs = R_now - δ - ρ

        return lhs - rhs
    end

    lo = z_now - z_width
    hi = z_now + z_width

    f_lo = euler_residual_z(lo)
    f_hi = euler_residual_z(hi)

    expand_iter = 0
    while (!isfinite(f_lo) || !isfinite(f_hi) || sign(f_lo) == sign(f_hi)) && expand_iter < 100
        lo -= 5.0
        hi += 5.0
        f_lo = euler_residual_z(lo)
        f_hi = euler_residual_z(hi)
        expand_iter += 1
    end

    if !isfinite(f_lo) || !isfinite(f_hi) || sign(f_lo) == sign(f_hi)
        return NaN
    end

    mid = NaN
    f_mid = Inf

    iter = 0
    while abs(f_mid) > tol && iter < maxiter
        mid = 0.5*(lo + hi)
        f_mid = euler_residual_z(mid)

        if !isfinite(f_mid)
            return NaN
        end

        if sign(f_mid) == sign(f_lo)
            lo = mid
            f_lo = f_mid
        else
            hi = mid
            f_hi = f_mid
        end

        iter += 1
    end

    return exp(mid)
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# HCD Equivalent Variation Functions
#------------------------------------------------------------------------------
function solve_Ω_bellman_hcd(ν_val, Pg, Ps, ψ, mg_hat, dmg_hat;
                             z_low=-50.0, z_high=50.0,
                             tol=1e-12, maxiter=500)

    residual_z(z) =
        exp(-ψ*z) - ν_val*dEdΩ_hcd(Pg, Ps, exp(z), mg_hat, dmg_hat)

    lo = z_low
    hi = z_high
    f_lo = residual_z(lo)
    f_hi = residual_z(hi)

    expand_iter = 0
    while (!isfinite(f_lo) || !isfinite(f_hi) || sign(f_lo) == sign(f_hi)) && expand_iter < 100
        lo -= 10.0
        hi += 10.0
        f_lo = residual_z(lo)
        f_hi = residual_z(hi)
        expand_iter += 1
    end

    if !isfinite(f_lo) || !isfinite(f_hi) || sign(f_lo) == sign(f_hi)
        return NaN
    end

    mid = NaN
    f_mid = Inf

    iter = 0
    while abs(f_mid) > tol && iter < maxiter
        mid = 0.5*(lo + hi)
        f_mid = residual_z(mid)

        if !isfinite(f_mid)
            return NaN
        end

        if sign(f_mid) == sign(f_lo)
            lo = mid
            f_lo = f_mid
        else
            hi = mid
            f_hi = f_mid
        end

        iter += 1
    end

    return exp(mid)
end

function Phi_t_hcd(ν_val, Pg, Ps, ψ, mg_hat, dmg_hat)
    Ω_star = solve_Ω_bellman_hcd(ν_val, Pg, Ps, ψ, mg_hat, dmg_hat)
    !isfinite(Ω_star) && return NaN

    E_star = E_hcd(Pg, Ps, Ω_star, mg_hat)
    U_star = U_crra(Ω_star, ψ)

    return U_star - ν_val*E_star
end

# EV income: m̂_{t,z,τ} = y_τ + [Φ_t(P_τ) − Φ_t(P_z)] / ν_t
function mhat_hcd(t::Int, z::Int, τ::Int; ν_t, Pgt, Pst, yt, ψ, mg_hat, dmg_hat)
    ν_val = ν_t[t]
    Phi_τ = Phi_t_hcd(ν_val, Pgt[τ], Pst[τ], ψ, mg_hat, dmg_hat)
    Phi_z = Phi_t_hcd(ν_val, Pgt[z], Pst[z], ψ, mg_hat, dmg_hat)

    return yt[τ] + (Phi_τ - Phi_z) / ν_val
end

# Hypothetical expenditure Ê_{t,z}
function Ehat_hcd(t::Int, z::Int; ν_t, Pgt, Pst, ψ, mg_hat, dmg_hat)
    Ω_hat = solve_Ω_bellman_hcd(ν_t[t], Pgt[z], Pst[z], ψ, mg_hat, dmg_hat)
    return E_hcd(Pgt[z], Pst[z], Ω_hat, mg_hat)
end

# EV goods expenditure share ŝ_{g,t,z}
function sg_hat_hcd(t::Int, z::Int; ν_t, Pgt, Pst, ψ, mg_hat, dmg_hat)
    Ω_hat = solve_Ω_bellman_hcd(ν_t[t], Pgt[z], Pst[z], ψ, mg_hat, dmg_hat)
    return sg_hcd(Pgt[z], Pst[z], Ω_hat, mg_hat)
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
p0    = [0.0]
fit   = curve_fit(exp_model_fixed_a, xdata, ydata, p0)
b_fit = fit.param[1]

g_wedge_Pg_Ps     = b_fit
wedge_Pg_Ps_trend = exp_model_fixed_a(xdata, fit.param)

# Compute Equilibrium Prices (relative to the investment numeraire)
Ps_t  = Ps.(calAx_t,As_t,1.000)
Pg_t  = Pg.(calAx_t,Ag_t,wedge_Pg_Ps_trend)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# HCD Engel-Curve Recovery Helper
#------------------------------------------------------------------------------
function build_hcd_engel_curve_from_e0(e0::Float64,
                                       Pg_t::Vector{Float64},
                                       Ps_t::Vector{Float64},
                                       mg_data::AbstractVector,
                                       g_abgp::Float64)

    T      = length(mg_data)
    xgrid  = collect(0:T-1)
    e_path = e0 .* exp.(g_abgp .* xgrid)

    mg_vec = Float64.(mg_data)

    logΩ_fit =
        log.(e_path) .-
        mg_vec .* log.(Pg_t) .-
        (1.0 .- mg_vec) .* log.(Ps_t)

    mg_hat, dmg_hat = make_mg_functions(logΩ_fit, mg_vec)
    mg_fit          = [mg_hat(z) for z in logΩ_fit]

    return mg_hat, dmg_hat, logΩ_fit, exp.(logΩ_fit), mg_fit
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Simulate the Model Under HCD Preferences
#------------------------------------------------------------------------------
function sim_model_full_hcd(;Pg_t::Vector{Float64},Ps_t::Vector{Float64},
                              θ::Float64,ρ::Float64,δ::Float64,ψ::Float64,
                              g_calAx::Float64,g_As::Float64,g_Ag::Float64,
                              g_h::Float64,g_l::Float64,g_n::Float64,
                              sg_data::AbstractVector)

    T = length(1980:2023)

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

        if !isfinite(E0) || E0 <= 0
            return -1e6
        end

        mg_hat_local, dmg_hat_local, _, _, _ =
            build_hcd_engel_curve_from_e0(Float64(E0), Pg_t, Ps_t, sg_data, g_abgp)

        Ω0 = solve_Ω_from_E_hcd(E0, Pg_t[1], Ps_t[1], mg_hat_local)

        if !isfinite(Ω0) || Ω0 <= 0
            return -1e6
        end

        R0 = θ * calAx_hat_t[1] * k0^(θ - 1)

        Ω1 = solve_Ω_next_hcd(Ω0,
                              Pg_t[1], Ps_t[1],
                              Pg_t[2], Ps_t[2],
                              R0,
                              ρ, δ, ψ,
                              mg_hat_local, dmg_hat_local)

        if !isfinite(Ω1) || Ω1 <= 0
            return -1e6
        end

        E1 = E_hcd(Pg_t[2], Ps_t[2], Ω1, mg_hat_local)

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
        return (
            eq       = eq,
            mg_hat   = nothing,
            dmg_hat  = nothing,
            logΩ_fit = fill(NaN, T),
            Ω_fit    = fill(NaN, T),
            mg_fit   = fill(NaN, T)
        )
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
    # Step 2: Compute E0 and recover the final HCD Engel curve
    #--------------------------------------------------------------------------

    e0 = E0_from_k0(k0)

    mg_hat, dmg_hat, logΩ_fit, Ω_fit, mg_fit =
        build_hcd_engel_curve_from_e0(Float64(e0), Pg_t, Ps_t, sg_data, g_abgp)

    Ω0 = solve_Ω_from_E_hcd(e0, Pg_t[1], Ps_t[1], mg_hat)

    if !isfinite(Ω0) || Ω0 <= 0
        eq = Dict{String, Vector{Float64}}()
        eq["sgt"] = fill(1e6, T)
        return (
            eq       = eq,
            mg_hat   = mg_hat,
            dmg_hat  = dmg_hat,
            logΩ_fit = logΩ_fit,
            Ω_fit    = Ω_fit,
            mg_fit   = mg_fit
        )
    end

    #--------------------------------------------------------------------------
    # Storage
    #--------------------------------------------------------------------------

    vars = ["Yt", "Xt", "Et", "Kt",
            "yt", "et", "ct", "Ωt", "kt", "ht",
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
            "Axt", "calAxt", "calAxhat"]

    eq = Dict{String, Vector{Any}}()
    for var in vars
        eq[var] = Vector{Any}()
    end

    # Initial values
    push!(eq["et"], e0)
    push!(eq["ct"], Ω0)       # keep ct name for downstream comparability with NHCES
    push!(eq["Ωt"], Ω0)
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
    # Step 3: Simulate forward
    #--------------------------------------------------------------------------

    for t in 1:T

        # Consumption expenditure shares
        sg_t = sg_hcd(eq["Pgt"][t], eq["Pst"][t], eq["Ωt"][t], mg_hat)
        ss_t = 1.0 - sg_t

        if !isfinite(sg_t) || !isfinite(ss_t) ||
           sg_t < 0 || sg_t > 1 ||
           ss_t < 0 || ss_t > 1
            return (
                eq       = Dict("sgt" => fill(1e6, T)),
                mg_hat   = mg_hat,
                dmg_hat  = dmg_hat,
                logΩ_fit = logΩ_fit,
                Ω_fit    = Ω_fit,
                mg_fit   = mg_fit
            )
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
            return (
                eq       = Dict("sgt" => fill(1e6, T)),
                mg_hat   = mg_hat,
                dmg_hat  = dmg_hat,
                logΩ_fit = logΩ_fit,
                Ω_fit    = Ω_fit,
                mg_fit   = mg_fit
            )
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
        # ABGP step: as in the NHCES file, per-capita consumption
        # expenditure grows at g_abgp. In HCD, we recover Ω_{t+1}
        # by inverting the HCD expenditure function.
        #--------------------------------------------------------------

        e_next = eq["et"][t] * exp(g_abgp)

        Ω_next = solve_Ω_from_E_hcd(e_next,
                                    Pg_t[t+1], Ps_t[t+1],
                                    mg_hat)

        if !isfinite(Ω_next) || !isfinite(e_next) ||
           Ω_next <= 0 || e_next <= 0
            return (
                eq       = Dict("sgt" => fill(1e6, T)),
                mg_hat   = mg_hat,
                dmg_hat  = dmg_hat,
                logΩ_fit = logΩ_fit,
                Ω_fit    = Ω_fit,
                mg_fit   = mg_fit
            )
        end

        push!(eq["ct"], Ω_next)
        push!(eq["Ωt"], Ω_next)
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
            return (
                eq       = Dict("sgt" => fill(1e6, T)),
                mg_hat   = mg_hat,
                dmg_hat  = dmg_hat,
                logΩ_fit = logΩ_fit,
                Ω_fit    = Ω_fit,
                mg_fit   = mg_fit
            )
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

    return (
        eq       = eq,
        mg_hat   = mg_hat,
        dmg_hat  = dmg_hat,
        logΩ_fit = logΩ_fit,
        Ω_fit    = Ω_fit,
        mg_fit   = mg_fit
    )
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Simulate the Model
#------------------------------------------------------------------------------
sim_result = sim_model_full_hcd(
    Pg_t=Float64.(Pg_t), Ps_t=Float64.(Ps_t),
    θ=θ, ρ=ρ, δ=δ, ψ=ψ,
    g_calAx=g_calAx, g_As=g_As, g_Ag=g_Ag,
    g_h=g_h, g_l=g_l, g_n=g_n,
    sg_data=VAC_GOOD_SHARE
)

sim_eq   = sim_result.eq
mg_hat   = sim_result.mg_hat
dmg_hat  = sim_result.dmg_hat
logΩ_fit = sim_result.logΩ_fit
Ω_fit    = sim_result.Ω_fit
mg_fit   = sim_result.mg_fit

if !haskey(sim_eq, "sgt") || any(x -> !isfinite(x) || x > 1.0 || x < 0.0, sim_eq["sgt"])
    error("Invalid HCD simulation detected.")
end

#Compute the Marginal Value of Wealth
Ωt_f  = Float64.(sim_eq["Ωt"])
ct_f  = Float64.(sim_eq["ct"])
et_f  = Float64.(sim_eq["et"])
Pgt_f = Float64.(sim_eq["Pgt"])
Pst_f = Float64.(sim_eq["Pst"])

ν_t   =
    Ωt_f.^(-ψ) ./
    [dEdΩ_hcd(Pgt_f[t], Pst_f[t], Ωt_f[t], mg_hat, dmg_hat) for t in eachindex(Ωt_f)]
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
    savefig(joinpath(figuresdir, "sg_t_Model_Fit_HCD.png"))
    println("Figure saved to: ", joinpath(figuresdir, "sg_t_Model_Fit_HCD.png"))
end
#----------------------------------------------------------

#---------------------------------------------------------
### Figure 4 (a)
GDP_data = HRV_df[HRV_df.year .>= 1980, "GDP_QI"]./ HRV_df[HRV_df.year .>= 1980, "GDP_QI"][1]

#Compute the Divisia Index
x_t         = Float64.(sim_eq["Xt"]) ./ Float64.(sim_eq["Nt"])
s_e         = Float64.(sim_eq["Et"]) ./ Float64.(sim_eq["Yt"])
sg_t        = Float64.(sim_eq["sgt"])
ss_t        = Float64.(sim_eq["sst"])

cgt_f       = Float64.(sim_eq["cgt"])
cst_f       = Float64.(sim_eq["cst"])

g_cg        = cgt_f[2:end]./cgt_f[1:end-1] .- 1
g_cs        = cst_f[2:end]./cst_f[1:end-1] .- 1
g_x         = x_t[2:end]./x_t[1:end-1] .- 1

gD          = s_e[2:end].*( sg_t[2:end].*g_cg .+ ss_t[2:end].*g_cs ) .+ (1 .- s_e[2:end]).*g_x
gD_agg      = [ 0 ; (gD .+ g_n) ]
FS          = integrate_trap(gD_agg)

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
    savefig(joinpath(figuresdir,"GDP_Model_vs_Data_HCD.png"))
    println("Figure saved to: ", joinpath(figuresdir, "GDP_Model_vs_Data_HCD.png"))
end
#---------------------------------------------------------

#---------------------------------------------------------------
### Figure 4 (b)
g_FS          = gD .+ g_n
GDP_data_long = HRV_df[!, "GDP_QI"]

# Compute the growth rate
g_GDP_long = [NaN; GDP_data_long[2:end] ./ GDP_data_long[1:end-1] .- 1]

# Compute 10-period centered moving average
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
    savefig(joinpath(figuresdir,"GDP_Growth_Model_vs_Data_HCD.png"))
    println("Figure saved to: ", joinpath(figuresdir, "GDP_Growth_Model_vs_Data_HCD.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 5 (a)
sc   = Float64.(sim_eq["et"]) ./ Float64.(sim_eq["yt"])
sg_z = Float64.(sim_eq["sgt"])
ss_z = Float64.(sim_eq["sst"])
se_z = sc[1:end]

#Gross Income Per Capita
m_t  = Float64.(sim_eq["yt"])

# Float64 income array for EV computations
yt_f  = Float64.(sim_eq["yt"])

# Convenience closures capturing HCD objects
_mhat(t, z, τ) = mhat_hcd(t, z, τ;
    ν_t=ν_t, Pgt=Pgt_f, Pst=Pst_f, yt=yt_f,
    ψ=ψ, mg_hat=mg_hat, dmg_hat=dmg_hat)

_Ehat(t, z) = Ehat_hcd(t, z;
    ν_t=ν_t, Pgt=Pgt_f, Pst=Pst_f,
    ψ=ψ, mg_hat=mg_hat, dmg_hat=dmg_hat)

_sg_hat(t, z) = sg_hat_hcd(t, z;
    ν_t=ν_t, Pgt=Pgt_f, Pst=Pst_f,
    ψ=ψ, mg_hat=mg_hat, dmg_hat=dmg_hat)

τ   = (1980:1:2023) .- 1980 .+ 1
T   = length(τ)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# HCD Fisher-Shell Indices — EV definition
#------------------------------------------------------------------------------
function FS_base_index_hcd(t_base; τ, yt_f, _mhat)
    t_prime = t_base - 1980 + 1
    mhat_b_b_τ = [_mhat(t_prime, t_prime, zi) for zi in τ]

    @assert isapprox(mhat_b_b_τ[t_prime], yt_f[t_prime]; rtol=1e-8)

    FS_b_b_τ = log.(mhat_b_b_τ) .- log.(mhat_b_b_τ[1])

    return FS_b_b_τ
end

FS_2023_2023_τ = FS_base_index_hcd(2023; τ=τ, yt_f=yt_f, _mhat=_mhat)
FS_2010_2010_τ = FS_base_index_hcd(2010; τ=τ, yt_f=yt_f, _mhat=_mhat)
FS_2000_2000_τ = FS_base_index_hcd(2000; τ=τ, yt_f=yt_f, _mhat=_mhat)
FS_1990_1990_τ = FS_base_index_hcd(1990; τ=τ, yt_f=yt_f, _mhat=_mhat)
FS_1980_1980_τ = FS_base_index_hcd(1980; τ=τ, yt_f=yt_f, _mhat=_mhat)

g_2023_2023_τ = diff(FS_2023_2023_τ)
g_2010_2010_τ = diff(FS_2010_2010_τ)
g_2000_2000_τ = diff(FS_2000_2000_τ)
g_1990_1990_τ = diff(FS_1990_1990_τ)
g_1980_1980_τ = diff(FS_1980_1980_τ)

#------------------------------------------------------------------------------
# Population index and aggregate EV indices
#------------------------------------------------------------------------------
pop_index = log.(Float64.(sim_eq["Nt"]) ./ Float64.(sim_eq["Nt"])[1])
g_N       = diff(log.(Float64.(sim_eq["Nt"])))

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

p = plot(1980:2023, FS_2023_2023_τ_agg, label=L"\mathcal{P}_{2023,z} - \textrm{2023-base \ Fisher-Shell \ Index}",
    ylabel="Cumulative Growth", linestyle=:dash, lw=2.0,
    xticks=1980:5:2023, yticks=0.00:0.20:1.40,
    ylims=(0.00, 1.40),
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.10, 0.920), legendfont=legendfont,
    xrotation=45, framestyle=:box)

plot!(p, 1980:2023, FS_1980_1980_τ_agg, label=L"\mathcal{L}_{1980,z} - \textrm{1980-base \  Fisher-Shell \ Index}",
    linestyle=:dot, lw=2, color=:black)

plot!(p, 1980:2023, FS, label=L"\mathcal{D}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    linestyle=:solid, lw=2, color=:black)

if save_figures
    savefig(joinpath(figuresdir,"FS_BBEV_HCD.png"))
    println("Figure saved to: ", joinpath(figuresdir, "FS_BBEV_HCD.png"))
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
z_prime = (1980:1:2023) .- 1980 .+ 1

# At Preferences of 1980 (t=1, z=1): Ê_{1980,1980}
Ehat_1980_1980 = _Ehat(1, 1)
se_1980_z      = Ehat_1980_1980 ./ yt_f[z_prime]

# At Preferences of 2023 (t=T, z=T): Ê_{2023,2023}
Ehat_2023_2023 = _Ehat(T, T)
se_2023_z      = Ehat_2023_2023 ./ yt_f[z_prime]

se = Float64.(sim_eq["et"]) ./ Float64.(sim_eq["yt"])

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
    savefig(joinpath(figuresdir,"se_t_HCD.png"))
    println("Figure saved to: ", joinpath(figuresdir, "se_t_HCD.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 6 (b)
sg_2023_z = _sg_hat(T, T)
sg_1980_z = _sg_hat(1, 1)

plot(1980:2023, Float64.(sim_eq["sgt"]), label=L"s_{g,z}",linestyle=:solid, lw=2,color=:black)

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
    savefig(joinpath(figuresdir, "sg_t_z_HCD.png"))
    println("Figure saved to: ", joinpath(figuresdir, "sg_t_z_HCD.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 7
g_FS              = gD .+ g_n
g_FS_2023_2023_τ = g_2023_2023_τ_agg
g_FS_2010_2010_τ = g_2010_2010_τ_agg
g_FS_2000_2000_τ = g_2000_2000_τ_agg
g_FS_1990_1990_τ = g_1990_1990_τ_agg
g_FS_1980_1980_τ = g_1980_1980_τ_agg

p = plot(1981:2023, g_FS_2023_2023_τ, label=L"g^{D}_{2023,z} - \textrm{2023-base \  FS \ Index}",
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

plot!(p, 1981:2023, g_FS_2010_2010_τ, label=L"g^{D}_{2010,z} - \textrm{\ 2010-base \  FS \ Index}",
    linestyle=:dash, lw=2, color=:black)

plot!(p, 1981:2023, g_FS_2000_2000_τ, label=L"g^{D}_{2000,z}- \textrm{2000-base \  FS \ Index}",
    linestyle=:dash, lw=1.5, minorgrid=false, color=:black)

plot!(p, 1981:2023, g_FS_1990_1990_τ, label=L"g^{D}_{1990,z} - \textrm{1990-base \  FS \ Index}",
    linestyle=:dash, lw=1.0, minorgrid=false, color=:black)

plot!(p, 1981:2023, g_FS_1980_1980_τ, label=L"g^{D}_{1980,z} - \textrm{1980-base \  FS \ Index}",
    linestyle=:dash, lw=0.5, minorgrid=false, color=:black)

plot!(p, 1981:2023, g_FS , label=L"g^{D}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    linestyle=:solid, lw=2, minorgrid=false, color=:black)

xaxis!(p, minor_ticks=true, minor_tick_step=1.00)
yaxis!(p, minor_ticks=true, minor_tick_step=0.01)

if save_figures
    savefig(joinpath(figuresdir, "FS_GrowthRates_1980_2023_HCD.png"))
    println("Figure saved to: ", joinpath(figuresdir, "FS_GrowthRates_1980_2023_HCD.png"))
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

#---------------------------------------------------------------
# Export HCD index results
results_HCD = DataFrame(
    year                  = collect(1980:2023),
    sgt_model             = Float64.(sim_eq["sgt"]),
    sgt_data              = Float64.(VAC_GOOD_SHARE),
    logΩ_fit              = Float64.(logΩ_fit),
    Ω_fit                 = Float64.(Ω_fit),
    Ωt_model              = Float64.(sim_eq["Ωt"]),
    FS_Divisia            = Float64.(FS),
    FS_1980_base          = Float64.(FS_1980_1980_τ_agg),
    FS_1990_base          = Float64.(FS_1990_1990_τ_agg),
    FS_2000_base          = Float64.(FS_2000_2000_τ_agg),
    FS_2010_base          = Float64.(FS_2010_2010_τ_agg),
    FS_2023_base          = Float64.(FS_2023_2023_τ_agg)
)

CSV.write(joinpath(figuresdir, "HCD_indices_results.csv"), results_HCD)
println("Results saved to: ", joinpath(figuresdir, "HCD_indices_results.csv"))
#---------------------------------------------------------------

#------------------------------------------------------------------------------
# End of file
#------------------------------------------------------------------------------
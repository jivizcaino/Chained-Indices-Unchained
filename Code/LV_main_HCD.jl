#------------------------------------------------------------------------------
# Replication code for: Chained Indices Unchained: Structural Transformation and
#                       the Welfare Foundations of Income Growth Measurement
# HCD (Heterothetic Cobb–Douglas) — CONSUMPTION EXPENDITURE INDICES
# By:           Omar Licandro and Juan I. Vizcaino
# Adapted from: LV_main_NHCES.jl (NHCES) and LV_main.jl (PIGL)
# Preferences:  Bohr–Mestieri–Robert-Nicoud (2026); measurement per Licandro's note
# This Version: 14/07/2026
#------------------------------------------------------------------------------
#
# SMOOTHING (two layers, parallel to LV_main.jl / LV_main_NHCES.jl)
# -----------------------------------------------------------------
# 1. Prices  — fit exponential trends to Pg_t_data and Ps_t_data through the
#    1980 origin (Pg_1980 = Ps_1980 = 1 by normalisation).  This mirrors the
#    wedge-fitted model prices in the PIGL/NHCES pipelines.
# 2. Shares — fit a polynomial Engel curve sg_t = m_g(Ω_t) in ln Ω_t, iterating
#    because Ω_t depends on sg_t through eq. (5).  With smooth prices the
#    iteration converges to a smooth m_g.
# Toggles: `smooth_prices`, `smooth_shares`.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Working directories
currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Packages
using XLSX, DataFrames, Statistics
using LaTeXStrings
using Plots; gr()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# CONFIGURATION
save_figures        = false
aggregate_index     = true
share_source        = :data          # :data or :model
expenditure_source  = :abgp          # :abgp or :model

smooth_prices       = true
smooth_shares       = true
poly_degree         = 3
smooth_maxiter      = 200
smooth_tol          = 1e-12
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Parameters taken from HRV
θ = 1/3
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Helper: cumulative trapezoidal integration
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
# Import the HRV data
file_name  = "HRV2021_Data.xlsx"
sheet_name = "Data"

HRV_data = XLSX.readtable(joinpath(datadir, file_name), sheet_name)
HRV_df   = DataFrame(HRV_data)

VAC_GOOD_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_GOOD_S"]
VAC_SERV_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_SERV_S"]

Nt = HRV_df[HRV_df.year .>= 1980, "POP"]
Lt = HRV_df[HRV_df.year .>= 1980, "LAB_TOT_QI"]

calAx_1980 = HRV_df[HRV_df.year .== 1980, "calA_X_I_TD"][1]
calAx_2023 = HRV_df[HRV_df.year .== 2023, "calA_X_I_TD"][1]

N_1980 = HRV_df[HRV_df.year .== 1980, "POP"][1]
N_2023 = HRV_df[HRV_df.year .== 2023, "POP"][1]
L_1980 = HRV_df[HRV_df.year .== 1980, "LAB_TOT_QI"][1]
L_2023 = HRV_df[HRV_df.year .== 2023, "LAB_TOT_QI"][1]

C_GOOD_P = HRV_df[HRV_df.year .>= 1980, "C_GOOD_P"] ./ HRV_df[HRV_df.year .== 1980, "C_GOOD_P"]
C_SERV_P = HRV_df[HRV_df.year .>= 1980, "C_SERV_P"] ./ HRV_df[HRV_df.year .== 1980, "C_SERV_P"]
X_TOT_P  = HRV_df[HRV_df.year .>= 1980, "X_TOT_P"]  ./ HRV_df[HRV_df.year .== 1980, "X_TOT_P"]

Pg_t_data = Float64.(C_GOOD_P ./ X_TOT_P)
Ps_t_data = Float64.(C_SERV_P ./ X_TOT_P)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# ABGP growth rates
g_calAx = log(calAx_2023 / calAx_1980) / (2023 - 1980)
g_n     = log(N_2023 / N_1980)         / (2023 - 1980)
g_l     = log(L_2023 / L_1980)         / (2023 - 1980)
g_h     = g_l - g_n
g_abgp  = g_calAx / (1 - θ) + g_h
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Time grid
years = 1980:2023
τidx  = collect(years) .- 1980
Tn    = length(years)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Smoothed prices — exponential trend through 1980 origin
# ln P_{g,t} = g_Pg · τ,  ln P_{s,t} = g_Ps · τ, OLS without intercept so that
# P_{g,1980} = P_{s,1980} = 1 by construction.
#------------------------------------------------------------------------------
if smooth_prices
    lnPg_data = log.(Pg_t_data)
    lnPs_data = log.(Ps_t_data)

    g_Pg = (τidx' * lnPg_data) / (τidx' * τidx)
    g_Ps = (τidx' * lnPs_data) / (τidx' * τidx)

    Pg_t = exp.(g_Pg .* τidx)
    Ps_t = exp.(g_Ps .* τidx)
else
    g_Pg = log(Pg_t_data[end] / Pg_t_data[1]) / (Tn - 1)
    g_Ps = log(Ps_t_data[end] / Ps_t_data[1]) / (Tn - 1)
    Pg_t = Pg_t_data
    Ps_t = Ps_t_data
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Goods share (raw / initial) and expenditure path
if share_source == :data
    sg_t = Float64.(VAC_GOOD_SHARE)
elseif share_source == :model
    @assert (@isdefined sim_eq) "share_source = :model requires a PIGL/NHCES `sim_eq` in scope."
    sg_t = Float64.(sim_eq["sgt"])
else
    error("share_source must be :data or :model")
end
ss_t = 1 .- sg_t

if expenditure_source == :abgp
    lnE_t = g_abgp .* τidx
elseif expenditure_source == :model
    @assert (@isdefined sim_eq) "expenditure_source = :model requires a PIGL/NHCES `sim_eq` in scope."
    lnE_t = log.(Float64.(sim_eq["et"]))
else
    error("expenditure_source must be :abgp or :model")
end
E_t = exp.(lnE_t)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Smoothed Engel curve  (HCD theory: sg_t = m_g(Ω_t))
#
# Polynomial-in-log Ω fitted to the RAW share sg_raw, iterating on ln Ω_t
# because Ω_t depends on sg_t through eq. (5).  With smooth prices, ln Ω_t
# is smooth and the fitted sg is smooth.
#------------------------------------------------------------------------------
sg_raw = share_source == :data ? Float64.(VAC_GOOD_SHARE) : copy(sg_t)

if share_source == :data && smooth_shares
    sg_iter = copy(sg_raw)

    local β_engel
    for it in 1:smooth_maxiter
        lnΩ_iter = lnE_t .- sg_iter .* log.(Pg_t) .- (1 .- sg_iter) .* log.(Ps_t)
        X        = hcat([lnΩ_iter .^ p for p in 0:poly_degree]...)
        β_engel  = X \ sg_raw
        sg_new   = X * β_engel

        Δ       = maximum(abs.(sg_new .- sg_iter))
        sg_iter = sg_new

        if Δ < smooth_tol
            @info "Engel-curve smoothing converged in $it iterations (Δ = $(round(Δ, sigdigits=3)))."
            break
        end
        if it == smooth_maxiter
            @warn "Engel-curve smoothing hit maxiter = $smooth_maxiter (Δ = $(round(Δ, sigdigits=3)))."
        end
    end

    sg_t = sg_iter
    ss_t = 1 .- sg_t
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Sectoral per-capita consumption quantities (using smoothed share and prices)
cg_t = sg_t .* E_t ./ Pg_t
cs_t = ss_t .* E_t ./ Ps_t
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# HCD real consumption  (note eq. 5)
lnΩ_t = lnE_t .- sg_t .* log.(Pg_t) .- ss_t .* log.(Ps_t)
lnΩ0  = lnΩ_t[1]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# INDEX 1 — Chained Divisia
g_cg = cg_t[2:end] ./ cg_t[1:end-1] .- 1
g_cs = cs_t[2:end] ./ cs_t[1:end-1] .- 1
gD_e = sg_t[2:end] .* g_cg .+ ss_t[2:end] .* g_cs

g_cg

addn = aggregate_index ? g_n : 0.0
D_e  = integrate_trap([0;gD_e .+ addn])
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# INDEX 2/3 — Fixed-base true quantity indices (note eq. 9)
popcum = aggregate_index ? g_n .* τidx : zeros(Tn)

function Q_base(zstar_year::Int)
    zi = zstar_year - 1980 + 1
    return (lnΩ_t .- lnΩ0) .+
           (sg_t .- sg_t[1]) .* log(Pg_t[zi]) .+
           (ss_t .- ss_t[1]) .* log(Ps_t[zi]) .+
           popcum
end

L_e_1980 = Q_base(1980)
Q_e_1990 = Q_base(1990)
Q_e_2000 = Q_base(2000)
Q_e_2010 = Q_base(2010)
P_e_2023 = Q_base(2023)

gL_e_1980 = diff(L_e_1980)[2:end]
gQ_e_1990 = diff(Q_e_1990)[2:end]
gQ_e_2000 = diff(Q_e_2000)[2:end]
gQ_e_2010 = diff(Q_e_2010)[2:end]
gP_e_2023 = diff(P_e_2023)[2:end]
gD_e_plot = diff(D_e)[2:end]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Diagnostics
println("""
$(repeat("=", 64))
HCD — Consumption Expenditure Indices (cumulative growth, 1980→2023)
share_source = $(share_source) | expenditure_source = $(expenditure_source) | aggregate = $(aggregate_index)
smooth_prices = $(smooth_prices) | smooth_shares = $(smooth_shares) | poly_degree = $(poly_degree)
g_Pg (fit) = $(round(g_Pg, digits=5)) | g_Ps (fit) = $(round(g_Ps, digits=5))
$(repeat("=", 64))
  Chained Divisia    D_e                : $(round(D_e[end],       digits=4))
  1980-base (Laspeyres) L_e             : $(round(L_e_1980[end],  digits=4))
  2023-base (Paasche)   P_e             : $(round(P_e_2023[end],  digits=4))
$(repeat("-", 64))
  L_e(1980-base) − Chained Divisia      : $(round(L_e_1980[end] - D_e[end], digits=4))
  P_e(2023-base) − Chained Divisia      : $(round(P_e_2023[end] - D_e[end], digits=4))
  Base-period bias  P_e − L_e           : $(round(P_e_2023[end] - L_e_1980[end], digits=4))
$(repeat("=", 64))
""")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# CHARTS
tickfont   = font(12)
guidefont  = font(12)
legendfont = font(12)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure 0a — Price smoothing diagnostic
if smooth_prices
    plot(years, Pg_t_data,
        label     = L"P_{g,t} \ \textrm{- \ Data}",
        ylabel    = "Consumption Prices\n(relative to investment)",
        linestyle = :dot, lw = 2.0, color = :black,
        xticks    = 1980:5:2025, minorgrid = false,
        xtickfont = tickfont, ytickfont = tickfont,
        xguidefont = guidefont, yguidefont = guidefont,
        legendfont = legendfont, legend = (0.100, 0.900),
        xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)
    plot!(years, Pg_t,      label = L"P_{g,t} \ \textrm{- \ Smoothed}",
        linestyle = :solid, lw = 2.0, color = :black)
    plot!(years, Ps_t_data, label = L"P_{s,t} \ \textrm{- \ Data}",
        linestyle = :dashdot, lw = 2.0, color = :black)
    plot!(years, Ps_t,      label = L"P_{s,t} \ \textrm{- \ Smoothed}",
        linestyle = :dash, lw = 2.0, color = :black)

    if save_figures
        savefig(joinpath(figuresdir, "HCD_price_smoothing.png"))
        println("Figure saved to: ", joinpath(figuresdir, "HCD_price_smoothing.png"))
    end
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure 0b — Engel-curve smoothing diagnostic
if share_source == :data && smooth_shares
    plot(years, sg_raw,
        label     = L"s_{g,t} \ \textrm{- \ Data}",
        ylabel    = "Goods Consumption\nExpenditure Share",
        linestyle = :dot, lw = 2.0, color = :black,
        xticks    = 1980:5:2025, minorgrid = false,
        xtickfont = tickfont, ytickfont = tickfont,
        xguidefont = guidefont, yguidefont = guidefont,
        legendfont = legendfont, legend = (0.775, 0.900),
        xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

    plot!(years, sg_t,
        label     = L"s_{g,t} \ \textrm{- \ Smoothed \ } m_g(\Omega_t)",
        linestyle = :solid, lw = 2.0, color = :black)

    if save_figures
        savefig(joinpath(figuresdir, "HCD_engel_curve_smoothing.png"))
        println("Figure saved to: ", joinpath(figuresdir, "HCD_engel_curve_smoothing.png"))
    end
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure — Cumulative growth: Chained Divisia vs 1980-base vs 2023-base
plot(years, D_e,
    label  = L"\mathcal{D_{e}}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel = "Cumulative Growth in\nReal Consumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.100, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(years, P_e_2023,
    label = L"\mathcal{P_{e}}_{2023,z} - \textrm{2023‑base  \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)

plot!(years, L_e_1980,
    label = L"\mathcal{L_{e}}_{1980,z} - \textrm{1980‑base  \ Index}",
    linestyle = :dot, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure — Growth rates for several base years vs the chained Divisia
plot(1982:2023, gP_e_2023,
    label = L"g^{D}_{e}_{2023,z} - \textrm{2023‑base \ Index}",
    linestyle = :dash, lw = 2.5, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.400, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(1982:2023, gQ_e_2010, label = L"g^{D}_{e}_{2010,z} - \textrm{2010‑base  \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)
plot!(1982:2023, gQ_e_2000, label = L"g^{D}_{e}_{2000,z} - \textrm{2000‑base  \ Index}",
    linestyle = :dash, lw = 1.5, color = :black)
plot!(1982:2023, gQ_e_1990, label = L"g^{D}_{e}_{1990,z} - \textrm{1990‑base  \ Index}",
    linestyle = :dash, lw = 1.0, color = :black)
plot!(1982:2023, gL_e_1980, label = L"g^{D}_{e}_{1980,z} - \textrm{1980‑base  \ Index}",
    linestyle = :dash, lw = 0.5, color = :black)
plot!(1982:2023, gD_e_plot,
    label  = L"g^{D}_{e}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel = "Growth Rate of Real \n Consumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_growth_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_growth_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------
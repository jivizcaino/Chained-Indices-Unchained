#------------------------------------------------------------------------------
# Replication code for: Chained Indices Unchained: Structural Transformation and
#                       the Welfare Foundations of Income Growth Measurement
# HCD (Heterothetic Cobb–Douglas) — CONSUMPTION EXPENDITURE INDICES
# By:           Omar Licandro and Juan I. Vizcaino
# Adapted from: LV_main_NHCES.jl (NHCES) and LV_main.jl (PIGL)
# Preferences:  Bohr–Mestieri–Robert-Nicoud (2026); measurement per Licandro's note
# This Version: 09/07/2026
#------------------------------------------------------------------------------
#
# WHAT CHANGES RELATIVE TO PIGL / NHCES
# -------------------------------------
# Under HCD, consumption-expenditure shares depend ONLY on real consumption Ω,
# not on relative prices (note, eq. 2).  Consequently the consumption welfare
# indices require NO preference parameters and NO estimation — they are pinned
# down nonparametrically by the observed series {E_t, Pg_t, Ps_t, sg_t}:
#
#   real consumption      ln Ω_t = ln E_t − sg_t ln Pg_t − ss_t ln Ps_t   (note eq. 5)
#
#   chained Divisia       g^D_t  = sg_t Δln c_{g,t} + ss_t Δln c_{s,t}
#                         (local growth rate of the EV measure; exact Konus COLI)
#
#   fixed-base (base z*)  ln Q^{z*}_t = (ln Ω_t − ln Ω_0)
#                                       + (sg_t − sg_0) ln Pg_{z*}
#                                       + (ss_t − ss_0) ln Ps_{z*}          (note eq. 9)
#
#   • z* = 1980  → Laspeyres / 1980-base Fisher–Shell index  (L_e)
#   • z* = 2023  → Paasche  / 2023-base Fisher–Shell index    (P_e)
#   • cumulated g^D → Chained Divisia index                    (D_e)
#
# Because prices are normalised to 1980 = 1, ln Pg_{1980} = ln Ps_{1980} = 0, so
# the 1980-base index reduces to the cumulative change in real consumption Ω.
# The 2023-base index adds the base-period-bias wedge
#   (sg_t − sg_0)(ln Pg_{2023} − ln Pg_{1980}) + (ss_t − ss_0)(ln Ps_{2023} − …),
# which is exact and is the object the paper isolates.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Working directories
currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Packages  (no BlackBoxOptim needed: HCD requires no preference estimation.
#            LsqFit is used only by the :logistic Engel-curve smoother.)
using XLSX, DataFrames, Statistics
using LsqFit
using LaTeXStrings
using Plots; gr()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# CONFIGURATION
save_figures        = false

# Add population growth g_n so the indices measure AGGREGATE real consumption
# (matches the PIGL/NHCES convention).  Set false for per-capita indices.
aggregate_index     = true

# Source of the goods consumption-expenditure share sg_t:
#   :data   →  VAC_GOOD_SHARE  (HCD nonparametric m_g(Ω) recovery — recommended)
#   :model  →  Float64.(sim_eq["sgt"])  (requires a PIGL/NHCES `sim_eq` in scope)
share_source        = :data

# Source of consumption expenditure E_t (its LEVEL cancels in every index):
#   :abgp   →  ln E_t = g_abgp·(t−1980)  (the model's ABGP consumption path)
#   :model  →  log.(Float64.(sim_eq["et"]))  (requires `sim_eq` in scope)
expenditure_source  = :abgp

# SMOOTHED ENGEL CURVE  m̂_g(Ω).  As in Licandro's note, the goods share fed into
# the indices is a smooth monotone trend fitted through the (Ω_t, s_g,t) pairs,
# NOT the raw annual share — this removes the business-cycle / measurement wiggles
# that would otherwise contaminate the Divisia growth rate and the fixed-base
# share-change terms.
#   :none      →  use the raw share (no smoothing)
#   :poly      →  least-squares polynomial of degree `poly_deg` in lnΩ  (base Julia)
#   :logistic  →  decreasing 4-parameter logistic in lnΩ  (monotone; needs LsqFit)
smooth_share_method = :logistic
poly_deg            = 3        # only used when smooth_share_method = :poly
n_share_iter        = 3        # Ω↔m̂_g(Ω) fixed-point passes (2–3 is plenty; Ω barely moves)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Parameters taken from HRV (only θ and the growth rates are needed here)
θ = 1/3
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Helper: cumulative trapezoidal integration (identical to LV_main_NHCES.jl)
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
# Engel-curve smoothing:  m̂_g(Ω) fitted through the (lnΩ_t, s_g,t) scatter.
#
# _fit_share returns a closure  φ(lnΩ) → share  for the chosen method.
# Polynomial: centred/scaled design for conditioning.
# Logistic:   decreasing 4-parameter S-curve, monotone by construction.
#------------------------------------------------------------------------------
function _fit_share(x::AbstractVector, y::AbstractVector, method::Symbol, deg::Int)
    if method == :poly
        μ, σ = mean(x), std(x)
        z    = (x .- μ) ./ σ
        V    = reduce(hcat, [z .^ p for p in 0:deg])     # Vandermonde
        β    = V \ collect(y)                            # OLS
        return lnΩ -> begin
            zz = (lnΩ - μ) / σ
            sum(β[p+1] * zz^p for p in 0:deg)
        end

    elseif method == :logistic
        # share(x) = c + L / (1 + exp(k (x - x0))),  L>0, k>0  ⇒  decreasing
        model(x, p) = p[4] .+ p[1] ./ (1 .+ exp.(p[2] .* (x .- p[3])))
        p0  = [max(y[1] - y[end], 1e-3), 1.0, mean(x), min(y[end], y[1])]
        fit = curve_fit(model, collect(x), collect(y), p0)
        p   = fit.param
        return lnΩ -> p[4] + p[1] / (1 + exp(p[2] * (lnΩ - p[3])))

    else
        error("smooth_share_method must be :none, :poly, or :logistic")
    end
end

# Fit m̂_g(Ω) and evaluate it, iterating the Ω ↔ share fixed point a few times
# (Ω depends on the share through eq. 5, so re-fit as the smoothed share updates).
function smooth_engel_share(sg_raw, lnE, lnPg, lnPs;
                            method::Symbol, deg::Int, n_iter::Int)
    method == :none && return Float64.(sg_raw)
    sg = Float64.(sg_raw)
    local φ
    for _ in 1:max(n_iter, 1)
        lnΩ = lnE .- sg .* lnPg .- (1 .- sg) .* lnPs          # eq. 5, current share
        φ   = _fit_share(lnΩ, sg_raw, method, deg)           # fit trend through DATA
        sg  = clamp.(φ.(lnΩ), 1e-6, 1 - 1e-6)                # evaluate smoothed trend
    end
    return sg
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Import the HRV data  (same block as LV_main_NHCES.jl)
file_name  = "HRV2021_Data.xlsx"
sheet_name = "Data"

HRV_data = XLSX.readtable(joinpath(datadir, file_name), sheet_name)
HRV_df   = DataFrame(HRV_data)

# Goods / services shares in consumption expenditure
VAC_GOOD_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_GOOD_S"]
VAC_SERV_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_SERV_S"]

# Population and efficiency units of labour
Nt = HRV_df[HRV_df.year .>= 1980, "POP"]
Lt = HRV_df[HRV_df.year .>= 1980, "LAB_TOT_QI"]

# TFP endpoints (for g_abgp)
calAx_1980 = HRV_df[HRV_df.year .== 1980, "calA_X_I_TD"][1]
calAx_2023 = HRV_df[HRV_df.year .== 2023, "calA_X_I_TD"][1]

# Population / labour endpoints
N_1980 = HRV_df[HRV_df.year .== 1980, "POP"][1]
N_2023 = HRV_df[HRV_df.year .== 2023, "POP"][1]
L_1980 = HRV_df[HRV_df.year .== 1980, "LAB_TOT_QI"][1]
L_2023 = HRV_df[HRV_df.year .== 2023, "LAB_TOT_QI"][1]

# Consumption prices relative to the investment numeraire (1980 = 1)
C_GOOD_P = HRV_df[HRV_df.year .>= 1980, "C_GOOD_P"] ./ HRV_df[HRV_df.year .== 1980, "C_GOOD_P"]
C_SERV_P = HRV_df[HRV_df.year .>= 1980, "C_SERV_P"] ./ HRV_df[HRV_df.year .== 1980, "C_SERV_P"]
X_TOT_P  = HRV_df[HRV_df.year .>= 1980, "X_TOT_P"]  ./ HRV_df[HRV_df.year .== 1980, "X_TOT_P"]

Pg_t_data = Float64.(C_GOOD_P ./ X_TOT_P)   # goods consumption price  (numeraire = investment)
Ps_t_data = Float64.(C_SERV_P ./ X_TOT_P)   # services consumption price
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# ABGP growth rates (identical construction to LV_main_NHCES.jl)
g_calAx = log(calAx_2023 / calAx_1980) / (2023 - 1980)
g_n     = log(N_2023 / N_1980)         / (2023 - 1980)
g_l     = log(L_2023 / L_1980)         / (2023 - 1980)
g_h     = g_l - g_n
g_abgp  = g_calAx / (1 - θ) + g_h        # per-capita consumption-expenditure growth on the ABGP
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Build the common economy series {Pg_t, Ps_t, sg_t, E_t}
years = 1980:2023
τidx  = collect(years) .- 1980           # 0,1,…,43
Tn    = length(years)

Pg_t  = Pg_t_data
Ps_t  = Ps_t_data
lnPg  = log.(Pg_t)
lnPs  = log.(Ps_t)

# Consumption expenditure per capita (level is irrelevant — it cancels in every index)
if expenditure_source == :abgp
    lnE_t = g_abgp .* τidx
elseif expenditure_source == :model
    @assert (@isdefined sim_eq) "expenditure_source = :model requires a PIGL/NHCES `sim_eq` in scope."
    lnE_t = log.(Float64.(sim_eq["et"]))
else
    error("expenditure_source must be :abgp or :model")
end
E_t = exp.(lnE_t)

# RAW goods consumption-expenditure share (the noisy data)
if share_source == :data
    sg_raw = Float64.(VAC_GOOD_SHARE)
elseif share_source == :model
    @assert (@isdefined sim_eq) "share_source = :model requires a PIGL/NHCES `sim_eq` in scope."
    sg_raw = Float64.(sim_eq["sgt"])
else
    error("share_source must be :data or :model")
end

# SMOOTHED share = fitted Engel curve m̂_g(Ω) evaluated along the sample.
# All indices below use sg_t (the trend), so the wiggles never enter.
sg_t = smooth_engel_share(sg_raw, lnE_t, lnPg, lnPs;
                          method = smooth_share_method, deg = poly_deg, n_iter = n_share_iter)
ss_t = 1 .- sg_t

# Fit quality (smoothed trend vs raw share)
share_rmse = sqrt(mean((sg_t .- sg_raw).^2))
share_mono = all(diff(sg_t) .<= 1e-10)   # true if the trend is (weakly) monotone decreasing



# Sectoral per-capita consumption quantities (built from the SMOOTHED share)
cg_t = sg_t .* E_t ./ Pg_t
cs_t = ss_t .* E_t ./ Ps_t
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# HCD real consumption  (note eq. 5):  ln Ω_t = ln E_t − sg_t ln Pg_t − ss_t ln Ps_t
# (built from the SMOOTHED share, consistent with the last smoothing pass)
lnΩ_t = lnE_t .- sg_t .* lnPg .- ss_t .* lnPs
lnΩ0  = lnΩ_t[1]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# INDEX 1 — Chained Divisia (real consumption expenditure)
# g^D_t = sg_t Δln c_{g,t} + ss_t Δln c_{s,t}
# (built exactly as in LV_main.jl / LV_main_NHCES.jl: net rates + integrate_trap;
#  +g_n makes it aggregate.)
g_cg = cg_t[2:end] ./ cg_t[1:end-1] .- 1
g_cs = cs_t[2:end] ./ cs_t[1:end-1] .- 1
gD_e = sg_t[2:end] .* g_cg .+ ss_t[2:end] .* g_cs     # Divisia quantity growth (per capita)

addn = aggregate_index ? g_n : 0.0
D_e  = integrate_trap([0.0; gD_e .+ addn])            # cumulative Chained Divisia index
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# INDEX 2/3 — Fixed-base true quantity indices  (note eq. 9)
# ln Q^{z*}_t = (ln Ω_t − ln Ω_0) + (sg_t − sg_0) ln Pg_{z*} + (ss_t − ss_0) ln Ps_{z*}
popcum = aggregate_index ? g_n .* τidx : zeros(Tn)    # cumulative population growth

function Q_base(zstar_year::Int)
    zi = zstar_year - 1980 + 1
    return (lnΩ_t .- lnΩ0) .+
           (sg_t .- sg_t[1]) .* log(Pg_t[zi]) .+
           (ss_t .- ss_t[1]) .* log(Ps_t[zi]) .+
           popcum
end

L_e_1980 = Q_base(1980)     # Laspeyres / 1980-base Fisher–Shell
Q_e_1990 = Q_base(1990)
Q_e_2000 = Q_base(2000)
Q_e_2010 = Q_base(2010)
P_e_2023 = Q_base(2023)     # Paasche / 2023-base Fisher–Shell

# Growth rates (first differences of the cumulative levels)
gL_e_1980 = diff(L_e_1980)
gQ_e_1990 = diff(Q_e_1990)
gQ_e_2000 = diff(Q_e_2000)
gQ_e_2010 = diff(Q_e_2010)
gP_e_2023 = diff(P_e_2023)
gD_e_plot = diff(D_e)                                  # chained growth (aggregate if addn≠0)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Diagnostics: base-period bias in cumulative real consumption growth by 2023
println("""
$(repeat("=", 64))
HCD — Consumption Expenditure Indices (cumulative growth, 1980→2023)
share_source = $(share_source) | expenditure_source = $(expenditure_source) | aggregate = $(aggregate_index)
Engel-curve smoothing = $(smooth_share_method)$(smooth_share_method == :poly ? " (deg $(poly_deg))" : "")  |  trend RMSE vs data = $(round(share_rmse, digits=4))  |  monotone↓ = $(share_mono)
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
# CHARTS  (styling parallel to the PIGL "Real Consumption Expenditure Indices" block)
tickfont   = font(12)
guidefont  = font(12)
legendfont = font(12)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure — Engel curve: raw goods share vs the smoothed trend m̂_g(Ω) used downstream
plot(years, sg_raw,
    label = L"s_{g,t} \textrm{ - \ Data}",
    ylabel = "Share of Goods in Consumption\nExpenditure",
    seriestype = :scatter, markersize = 3, markercolor = :gray, markerstrokewidth = 0,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.775, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(years, sg_t,
    label = L"\hat{m}_{g}(\Omega_t) \textrm{ - \ Smoothed \ trend}",
    linestyle = :solid, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_goods_share_smoothed.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_goods_share_smoothed.png"))
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
    label = L"\mathcal{P_{e}}_{2023,z} - \textrm{2023‑base \  FS \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)

plot!(years, L_e_1980,
    label = L"\mathcal{L_{e}}_{1980,z} - \textrm{1980‑base \  FS \ Index}",
    linestyle = :dot, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure — Growth rates for several base years vs the chained Divisia
plot(1981:2023, gP_e_2023,
    label = L"g^{D}_{e}_{2023,z} - \textrm{2023‑base \  FS \ Index}",
    linestyle = :dash, lw = 2.5, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.100, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(1981:2023, gQ_e_2010, label = L"g^{D}_{e}_{2010,z} - \textrm{2010‑base \  FS \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)
plot!(1981:2023, gQ_e_2000, label = L"g^{D}_{e}_{2000,z} - \textrm{2000‑base \  FS \ Index}",
    linestyle = :dash, lw = 1.5, color = :black)
plot!(1981:2023, gQ_e_1990, label = L"g^{D}_{e}_{1990,z} - \textrm{1990‑base \  FS \ Index}",
    linestyle = :dash, lw = 1.0, color = :black)
plot!(1981:2023, gL_e_1980, label = L"g^{D}_{e}_{1980,z} - \textrm{1980‑base \  FS \ Index}",
    linestyle = :dash, lw = 0.5, color = :black)
plot!(1981:2023, gD_e_plot,
    label  = L"g^{D}_{e}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel = "Growth Rate of Real \n Consumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_growth_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_growth_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Replication code for: Chained Indices Unchained: Structural Transformation and
#                       the Welfare Foundations of Income Growth Measurement
# HCD (Heterothetic Cobb–Douglas) — CONSUMPTION EXPENDITURE INDICES
# By:           Omar Licandro and Juan I. Vizcaino
# Adapted from: LV_main_NHCES.jl (NHCES) and LV_main.jl (PIGL)
# Preferences:  Bohr–Mestieri–Robert-Nicoud (2026); measurement per Licandro's note
# This Version: 09/07/2026
#------------------------------------------------------------------------------
#
# WHAT CHANGES RELATIVE TO PIGL / NHCES
# -------------------------------------
# Under HCD, consumption-expenditure shares depend ONLY on real consumption Ω,
# not on relative prices (note, eq. 2).  Consequently the consumption welfare
# indices require NO preference parameters and NO estimation — they are pinned
# down nonparametrically by the observed series {E_t, Pg_t, Ps_t, sg_t}:
#
#   real consumption      ln Ω_t = ln E_t − sg_t ln Pg_t − ss_t ln Ps_t   (note eq. 5)
#
#   chained Divisia       g^D_t  = s̄g_t Δln c_{g,t} + s̄s_t Δln c_{s,t}
#                         (local growth rate of the EV measure; exact Konus COLI)
#
#   current-base (t1)     ln Q^CB_t = (ln Ω_t − ln Ω_0)
#                                     + Σ_i [m_i(Ω_t) − m_i(Ω_0)] ln p_{i,t1}   (note eq. 8–9)
#                         with reference prices fixed at the terminal date t1 = 2023.
#
# These are the ONLY two welfare objects in Licandro's note: the chained Divisia
# quantity index and the current-base true quantity index.  (No separate 1980-base
# index: with prices normalised to 1980 = 1 that would collapse to ln Ω growth.)
# The gap  Q^CB − D_e  is the chained-vs-fixed-base wedge the paper isolates.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Working directories
currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Packages  (no BlackBoxOptim needed: HCD requires no preference estimation.
#            LsqFit is used only by the :logistic Engel-curve smoother.)
using XLSX, DataFrames, Statistics
using LsqFit
using LaTeXStrings
using Plots; gr()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# CONFIGURATION
save_figures        = false

# Add population growth g_n so the indices measure AGGREGATE real consumption
# (matches the PIGL/NHCES convention).  Set false for per-capita indices.
aggregate_index     = true

# Source of the goods consumption-expenditure share sg_t:
#   :data   →  VAC_GOOD_SHARE  (HCD nonparametric m_g(Ω) recovery — recommended)
#   :model  →  Float64.(sim_eq["sgt"])  (requires a PIGL/NHCES `sim_eq` in scope)
share_source        = :data

# Source of consumption expenditure E_t (its LEVEL cancels in every index):
#   :abgp   →  ln E_t = g_abgp·(t−1980)  (the model's ABGP consumption path)
#   :model  →  log.(Float64.(sim_eq["et"]))  (requires `sim_eq` in scope)
expenditure_source  = :abgp

# SMOOTHED ENGEL CURVE  m̂_g(Ω).  As in Licandro's note, the goods share fed into
# the indices is a smooth monotone trend fitted through the (Ω_t, s_g,t) pairs,
# NOT the raw annual share — this removes the business-cycle / measurement wiggles
# that would otherwise contaminate the Divisia growth rate and the fixed-base
# share-change terms.
#   :none      →  use the raw share (no smoothing)
#   :poly      →  least-squares polynomial of degree `poly_deg` in lnΩ  (base Julia)
#   :logistic  →  decreasing 4-parameter logistic in lnΩ  (monotone; needs LsqFit)
smooth_share_method = :poly
poly_deg            = 3        # only used when smooth_share_method = :poly
n_share_iter        = 3        # Ω↔m̂_g(Ω) fixed-point passes (2–3 is plenty; Ω barely moves)

# Price series entering the indices:
#   :model  →  smooth wedge-adjusted prices (exp. TFP trends × fitted exp. wedge) — RECOMMENDED
#   :prod   →  smooth productivity-implied prices Pg=calAx/Ag, Ps=calAx/As (no wedge)
#   :data   →  raw annual data prices (wiggly; makes ln Ω and the Engel fit wiggle)
price_source        = :model
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Parameters taken from HRV (only θ and the growth rates are needed here)
θ = 1/3
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Engel-curve smoothing:  m̂_g(Ω) fitted through the (lnΩ_t, s_g,t) scatter.
#
# _fit_share returns a closure  φ(lnΩ) → share  for the chosen method.
# Polynomial: centred/scaled design for conditioning.
# Logistic:   decreasing 4-parameter S-curve, monotone by construction.
#------------------------------------------------------------------------------
function _fit_share(x::AbstractVector, y::AbstractVector, method::Symbol, deg::Int)
    if method == :poly
        μ, σ = mean(x), std(x)
        z    = (x .- μ) ./ σ
        V    = reduce(hcat, [z .^ p for p in 0:deg])     # Vandermonde
        β    = V \ collect(y)                            # OLS
        return lnΩ -> begin
            zz = (lnΩ - μ) / σ
            sum(β[p+1] * zz^p for p in 0:deg)
        end

    elseif method == :logistic
        # share(x) = c + L / (1 + exp(k (x - x0))),  L>0, k>0  ⇒  decreasing
        model(x, p) = p[4] .+ p[1] ./ (1 .+ exp.(p[2] .* (x .- p[3])))
        p0  = [max(y[1] - y[end], 1e-3), 1.0, mean(x), min(y[end], y[1])]
        fit = curve_fit(model, collect(x), collect(y), p0)
        p   = fit.param
        return lnΩ -> p[4] + p[1] / (1 + exp(p[2] * (lnΩ - p[3])))

    else
        error("smooth_share_method must be :none, :poly, or :logistic")
    end
end

# Fit m̂_g(Ω) and evaluate it, iterating the Ω ↔ share fixed point a few times
# (Ω depends on the share through eq. 5, so re-fit as the smoothed share updates).
function smooth_engel_share(sg_raw, lnE, lnPg, lnPs;
                            method::Symbol, deg::Int, n_iter::Int)
    method == :none && return Float64.(sg_raw)
    sg = Float64.(sg_raw)
    local φ
    for _ in 1:max(n_iter, 1)
        lnΩ = lnE .- sg .* lnPg .- (1 .- sg) .* lnPs          # eq. 5, current share
        φ   = _fit_share(lnΩ, sg_raw, method, deg)           # fit trend through DATA
        sg  = clamp.(φ.(lnΩ), 1e-6, 1 - 1e-6)                # evaluate smoothed trend
    end
    return sg
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Import the HRV data  (same block as LV_main_NHCES.jl)
file_name  = "HRV2021_Data.xlsx"
sheet_name = "Data"

HRV_data = XLSX.readtable(joinpath(datadir, file_name), sheet_name)
HRV_df   = DataFrame(HRV_data)

# Goods / services shares in consumption expenditure
VAC_GOOD_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_GOOD_S"]
VAC_SERV_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_SERV_S"]

# Population and efficiency units of labour
Nt = HRV_df[HRV_df.year .>= 1980, "POP"]
Lt = HRV_df[HRV_df.year .>= 1980, "LAB_TOT_QI"]

# TFP endpoints
calAx_1980 = HRV_df[HRV_df.year .== 1980, "calA_X_I_TD"][1]
calAx_2023 = HRV_df[HRV_df.year .== 2023, "calA_X_I_TD"][1]
Ag_1980    = HRV_df[HRV_df.year .== 1980, "TFPVA_GOOD_I"][1]
Ag_2023    = HRV_df[HRV_df.year .== 2023, "TFPVA_GOOD_I"][1]
As_1980    = HRV_df[HRV_df.year .== 1980, "TFPVA_SERV_I"][1]
As_2023    = HRV_df[HRV_df.year .== 2023, "TFPVA_SERV_I"][1]

# Population / labour endpoints
N_1980 = HRV_df[HRV_df.year .== 1980, "POP"][1]
N_2023 = HRV_df[HRV_df.year .== 2023, "POP"][1]
L_1980 = HRV_df[HRV_df.year .== 1980, "LAB_TOT_QI"][1]
L_2023 = HRV_df[HRV_df.year .== 2023, "LAB_TOT_QI"][1]

# Consumption prices relative to the investment numeraire (1980 = 1)
C_GOOD_P = HRV_df[HRV_df.year .>= 1980, "C_GOOD_P"] ./ HRV_df[HRV_df.year .== 1980, "C_GOOD_P"]
C_SERV_P = HRV_df[HRV_df.year .>= 1980, "C_SERV_P"] ./ HRV_df[HRV_df.year .== 1980, "C_SERV_P"]
X_TOT_P  = HRV_df[HRV_df.year .>= 1980, "X_TOT_P"]  ./ HRV_df[HRV_df.year .== 1980, "X_TOT_P"]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# ABGP growth rates (identical construction to LV_main_NHCES.jl)
g_calAx = log(calAx_2023 / calAx_1980) / (2023 - 1980)
g_Ag    = log(Ag_2023 / Ag_1980)       / (2023 - 1980)
g_As    = log(As_2023 / As_1980)       / (2023 - 1980)
g_n     = log(N_2023 / N_1980)         / (2023 - 1980)
g_l     = log(L_2023 / L_1980)         / (2023 - 1980)
g_h     = g_l - g_n
g_abgp  = g_calAx / (1 - θ) + g_h        # per-capita consumption-expenditure growth on the ABGP
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# SMOOTH consumption prices (relative to the investment numeraire), built as in
# LV_main_NHCES.jl:
#   Ps(calAx,As,wedge_s) = (calAx/As)·wedge_s
#   Pg(calAx,Ag,wedge_g) = (calAx/Ag)·wedge_g
# where the goods–services wedge is fitted to the data as a smooth exponential
# trend exp(b·(t−1980)).  Two smooth variants are produced:
#   *_model  : with the fitted wedge (matches the data relative price)  → price_source = :model
#   *_prod   : no wedge, pure relative productivities                   → price_source = :prod
Ps_price(calAx, As, wedge_s) = (calAx / As) * wedge_s
Pg_price(calAx, Ag, wedge_g) = (calAx / Ag) * wedge_g

τidx0   = collect((1980:2023) .- 1980)
Ag_t    = exp.(g_Ag    .* τidx0)
As_t    = exp.(g_As    .* τidx0)
calAx_t = exp.(g_calAx .* τidx0)

# Raw annual data prices (numeraire = investment) — used to fit the wedge / :data switch
Pg_t_data = Float64.(C_GOOD_P ./ X_TOT_P)
Ps_t_data = Float64.(C_SERV_P ./ X_TOT_P)

# Undistorted (no-wedge) productivity-implied prices
Ps_t_undist = Ps_price.(calAx_t, As_t, 1.000)
Pg_t_undist = Pg_price.(calAx_t, Ag_t, 1.000)
Pg_t_prod   = Pg_t_undist
Ps_t_prod   = Ps_t_undist

# Fit the relative-price wedge as exp(b·(t−1980)) and apply it to the goods price
wedge_Pg_Ps       = (Pg_t_data ./ Ps_t_data) ./ (Pg_t_undist ./ Ps_t_undist)
exp_model_fixed_a(x, p) = exp.(p[1] .* x)
wedge_fit         = curve_fit(exp_model_fixed_a, τidx0, Float64.(wedge_Pg_Ps), [0.0])
wedge_Pg_Ps_trend = exp_model_fixed_a(τidx0, wedge_fit.param)

Pg_t_model = Pg_price.(calAx_t, Ag_t, wedge_Pg_Ps_trend)   # smooth, wedge-adjusted goods price
Ps_t_model = Ps_price.(calAx_t, As_t, 1.000)               # smooth services price
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Build the common economy series {Pg_t, Ps_t, sg_t, E_t}
years = 1980:2023
τidx  = collect(years) .- 1980           # 0,1,…,43
Tn    = length(years)

if price_source == :model
    Pg_t = Pg_t_model            # wedge-adjusted smooth prices (match data relative price)
    Ps_t = Ps_t_model
elseif price_source == :prod
    Pg_t = Pg_t_prod             # productivity-implied smooth prices (no wedge)
    Ps_t = Ps_t_prod
elseif price_source == :data
    Pg_t = Pg_t_data             # raw annual data prices (wiggly)
    Ps_t = Ps_t_data
else
    error("price_source must be :model, :prod, or :data")
end
lnPg  = log.(Pg_t)
lnPs  = log.(Ps_t)

# Consumption expenditure per capita (level is irrelevant — it cancels in every index)
if expenditure_source == :abgp
    lnE_t = g_abgp .* τidx
elseif expenditure_source == :model
    @assert (@isdefined sim_eq) "expenditure_source = :model requires a PIGL/NHCES `sim_eq` in scope."
    lnE_t = log.(Float64.(sim_eq["et"]))
else
    error("expenditure_source must be :abgp or :model")
end
E_t = exp.(lnE_t)

# RAW goods consumption-expenditure share (the noisy data)
if share_source == :data
    sg_raw = Float64.(VAC_GOOD_SHARE)
elseif share_source == :model
    @assert (@isdefined sim_eq) "share_source = :model requires a PIGL/NHCES `sim_eq` in scope."
    sg_raw = Float64.(sim_eq["sgt"])
else
    error("share_source must be :data or :model")
end

# SMOOTHED share = fitted Engel curve m̂_g(Ω) evaluated along the sample.
# All indices below use sg_t (the trend), so the wiggles never enter.
sg_t = smooth_engel_share(sg_raw, lnE_t, lnPg, lnPs;
                          method = smooth_share_method, deg = poly_deg, n_iter = n_share_iter)
ss_t = 1 .- sg_t

# Fit quality (smoothed trend vs raw share)
share_rmse = sqrt(mean((sg_t .- sg_raw).^2))
share_mono = all(diff(sg_t) .<= 1e-10)   # true if the trend is (weakly) monotone decreasing

# Sectoral per-capita consumption quantities (built from the SMOOTHED share)
cg_t = sg_t .* E_t ./ Pg_t
cs_t = ss_t .* E_t ./ Ps_t
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# HCD real consumption  (note eq. 5):  ln Ω_t = ln E_t − sg_t ln Pg_t − ss_t ln Ps_t
# (built from the SMOOTHED share, consistent with the last smoothing pass)
lnΩ_t = lnE_t .- sg_t .* lnPg .- ss_t .* lnPs
lnΩ0  = lnΩ_t[1]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Population term added identically to EVERY index (per-capita → aggregate).
# Using one common additive g_n·(t−1980) keeps all indices on the same footing.
popcum = aggregate_index ? g_n .* τidx : zeros(Tn)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# INDEX 1 — Chained Divisia (real consumption expenditure), per capita then +pop.
# Törnqvist (superlative) discretisation of  g^D_t = Σ_i s̄_{i,t} Δln c_{i,t},
# with average shares s̄ = ½(s_t + s_{t−1}).  Log-based and cumulated by simple
# summation, so it is directly comparable to the fixed-base log levels below.
Δln_cg = log.(cg_t[2:end]) .- log.(cg_t[1:end-1])
Δln_cs = log.(cs_t[2:end]) .- log.(cs_t[1:end-1])
s̄g     = 0.5 .* (sg_t[2:end] .+ sg_t[1:end-1])
s̄s     = 1 .- s̄g
gD_e   = s̄g .* Δln_cg .+ s̄s .* Δln_cs                    # Divisia quantity growth (per capita)

D_e = [0.0; cumsum(gD_e)] .+ popcum                       # cumulative Chained Divisia index
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# INDEX 2 — Current-base true quantity index  Q^CB  (note eq. 8–9)
# ln Q^CB_t = (ln Ω_t − ln Ω_0) + Σ_i [m_i(Ω_t) − m_i(Ω_0)] ln p_{i,t1},
# with reference prices fixed at the TERMINAL date t1 = 2023.  This is the only
# fixed-base welfare measure in the note (there is no separate 1980-base index;
# with 1980-normalised prices that would collapse to ln Ω growth itself).
i_t1 = length(years)                                     # terminal period index (2023)
Q_CB = (lnΩ_t .- lnΩ0) .+
       (sg_t .- sg_t[1]) .* log(Pg_t[i_t1]) .+
       (ss_t .- ss_t[1]) .* log(Ps_t[i_t1]) .+
       popcum

# Growth rates (first differences of the cumulative levels)
gQ_CB     = diff(Q_CB)
gD_e_plot = diff(D_e)                                  # chained per-period growth (incl. population)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Diagnostics: chained Divisia vs current-base true quantity index (note)
println("""
$(repeat("=", 64))
HCD — Consumption Expenditure Indices (cumulative growth, 1980→2023)
share_source = $(share_source) | expenditure_source = $(expenditure_source) | aggregate = $(aggregate_index)
Engel-curve smoothing = $(smooth_share_method)$(smooth_share_method == :poly ? " (deg $(poly_deg))" : "")  |  trend RMSE vs data = $(round(share_rmse, digits=4))  |  monotone↓ = $(share_mono)
$(repeat("=", 64))
  Chained Divisia            D_e         : $(round(D_e[end], digits=4))
  Current-base (2023) Q^CB               : $(round(Q_CB[end], digits=4))
$(repeat("-", 64))
  Chaining wedge  Q^CB − D_e             : $(round(Q_CB[end] - D_e[end], digits=4))
  (real income   ln Ω_t − ln Ω_0 + pop   : $(round((lnΩ_t[end]-lnΩ0)+popcum[end], digits=4)))
$(repeat("=", 64))
""")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# CHARTS  (styling parallel to the PIGL "Real Consumption Expenditure Indices" block)
tickfont   = font(12)
guidefont  = font(12)
legendfont = font(12)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure — Engel curve: raw goods share vs the smoothed trend m̂_g(Ω) used downstream
plot(years, sg_raw,
    label = L"s_{g,t} \textrm{ - \ Data}",
    ylabel = "Share of Goods in Consumption\nExpenditure",
    seriestype = :scatter, markersize = 3, markercolor = :gray, markerstrokewidth = 0,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.775, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(years, sg_t,
    label = L"\hat{m}_{g}(\Omega_t) \textrm{ - \ Smoothed \ trend}",
    linestyle = :solid, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_goods_share_smoothed.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_goods_share_smoothed.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure — Cumulative growth: Chained Divisia vs Current-base (2023) true index
plot(years, D_e,
    label  = L"\mathcal{D}_{e} \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel = "Cumulative Growth in\nReal Consumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.100, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(years, Q_CB,
    label = L"\mathcal{Q}^{CB}_{e} - \textrm{Current-base \ (2023) \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Figure — Growth rates: Chained Divisia vs Current-base (2023)
plot(1981:2023, gD_e_plot,
    label  = L"g^{D}_{e} \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel = "Growth Rate of Real \n Consumption Expenditure",
    linestyle = :solid, lw = 2.0, color = :black,
    xticks = 1980:5:2025, minorgrid = false,
    xtickfont = tickfont, ytickfont = tickfont,
    xguidefont = guidefont, yguidefont = guidefont,
    legendfont = legendfont, legend = (0.100, 0.900),
    xrotation = 45, left_margin = 5Plots.mm, framestyle = :box)

plot!(1981:2023, gQ_CB,
    label = L"g^{Q^{CB}}_{e} - \textrm{Current-base \ (2023) \ Index}",
    linestyle = :dash, lw = 2.0, color = :black)

if save_figures
    savefig(joinpath(figuresdir, "HCD_growth_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "HCD_growth_consumption_expenditure_indices.png"))
end
#------------------------------------------------------------------------------

run_konus_decomp = false
#-==============================================================================
# Konus INCOME-vs-SUBSTITUTION decomposition (consumption, 1980→2023)
# Same object/split as the NHCES and PIGL append-blocks, writing one row to the
# shared CSV.  Under HCD the compensated share m̂_g(Ω) is PRICE-INDEPENDENT, so
# G_income = G_total and G_subst ≡ 0 by construction — a numerical check.
#-==============================================================================
if run_konus_decomp
  const DECOMP_CSV = abspath(joinpath(@__DIR__, "consumption_index_decomposition.csv"))
  function upsert_decomp_row(path, model, G_total, G_inc_1980, G_inc_2023, G_inc_avg, G_sub_avg)
      header = "model,G_total,G_income_1980,G_income_2023,G_income_avg,G_subst_avg"
      newline = join(string.((model, G_total, G_inc_1980, G_inc_2023, G_inc_avg, G_sub_avg)), ",")
      keep = String[]
      if isfile(path)
          for l in readlines(path)
              (isempty(strip(l)) || l == header || startswith(l, string(model, ","))) && continue
              push!(keep, l)
          end
      end
      open(path, "w") do io
          println(io, header)
          foreach(l -> println(io, l), keep)
          println(io, newline)
      end
      return nothing
  end

  i0 = 1
  iT = length(sg_t)                                     # 1980 → 2023
  Δln_PgPs = (lnPg[iT] - lnPg[i0]) - (lnPs[iT] - lnPs[i0])   # Δ ln(P_g/P_s)
  Δsg      = sg_t[iT] - sg_t[i0]                        # smoothed Engel-share change

  # HCD compensated share depends only on Ω → income term is ref-price-invariant
  G_income_hcd = Δsg * Δln_PgPs
  # total = current-base − base-base;  base-base ≡ ln Ω growth (prices 1980-normalised)
  G_total_hcd  = Q_CB[iT] - ((lnΩ_t[iT] - lnΩ0) + popcum[iT])
  G_subst_hcd  = G_total_hcd - G_income_hcd             # ≡ 0 up to rounding

  upsert_decomp_row(DECOMP_CSV, "HCD", G_total_hcd, G_income_hcd, G_income_hcd, G_income_hcd, G_subst_hcd)

  println("""
  $(repeat("=", 68))
  HCD — Konus income/substitution decomposition (consumption, 1980→2023)
  $(repeat("=", 68))
    Δ ln(P_g/P_s)                        = $(round(Δln_PgPs, digits=4))
    G_total  (current-base − base-base)  = $(round(G_total_hcd, digits=4))
    G_income (= Δs_g · Δln(P_g/P_s))     = $(round(G_income_hcd, digits=4))
    G_subst  (residual, should be ≈ 0)   = $(round(G_subst_hcd, digits=6))
    → wrote row to $(DECOMP_CSV)
  $(repeat("=", 68))
  """)
end
#------------------------------------------------------------------------------
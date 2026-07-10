#------------------------------------------------------------------------------
#Replication code for: Chained Indices Unchained: Structural Transformation and the Welfare Foundations of Income Growth Measurement
#By:                   Omar Licandro and Juan I. Vizcaino
#This Version:         04/2026
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Set the Working Directory
currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Configuration 
# Set to true to save figures, false to only display
save_figures = false

# Set to true to run SMM estimation, false to use pre-estimated parameters
run_smm      = false

# Set to true to re-estimate η and γ for a given fixed value of χ
run_smm_fixed_chi = false
chi_fixed         = 0.36

# Set to true to run the Konus decomposition of the fixed-base bias in the consumption expenditure indices
run_konus_decomp = false  
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Activate local project (Code/Project.toml)
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using XLSX, DataFrames, BlackBoxOptim ,  Statistics
#using MathJaxRenderer, 
using LaTeXStrings, LsqFit, Random, OrderedCollections, PrettyTables
using Plots; gr()
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

## Goods' Share in Consumption Expenditure
sg(Pg,e,Ps,η,χ,γ) = η*((e/Ps)^(-χ))*((Pg/Ps)^γ)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Functions Used to Compute the Equivalent Variation Measures
e_t_x(t,x; χ, ν_t, Pst)    = (ν_t[t]*(Pst[x]^χ))^(1/(χ-1))
sg_t_x(t,x;χ,η,γ,Pst,Pgt)  = η*( ( e_t_x(t,x; χ, ν_t, Pst)/Pst[x] )^(-χ) )*((Pgt[x]/Pst[x])^γ)
cal_v_t_x(t,x;χ,Pst)       = ( Pst[t]/Pst[x] )^(  (χ/(1-χ)) )

function B_t_z_τ(t::Int, z::Int, τ::Int; χ, η, γ, Pst::AbstractVector, Pgt::AbstractVector)
    """
    B_t_z_τ(t, z, τ; χ, η, γ, Pst, Pgt)
    Compute the bias correction coefficient used in the equivalent variation measure `mhat_t_z_τ`.

    Given preferences evaluated at period `t`, a base price period `z`, and a welfare comparison 
    period `τ`, this function returns the scalar `B` such that:

        mhat(t, z, τ) = y(τ) + B * ehat(t, z)

    where `ehat(t, z)` is the compensating expenditure at preferences `t` and prices `z`.

    The coefficient is defined as:

        B = (1/χ - 1 - ŝg(t,τ)/γ) * v(z,τ) - (1/χ - 1 - ŝg(t,z)/γ)

    where `ŝg(t,x)` is the goods share in consumption expenditure evaluated at preferences `t`
    and prices at period `x`, and `v(z,τ) = (Ps[z]/Ps[τ])^(χ/(1-χ))` is the service-price
    scaling factor.

    # Arguments
    - `t::Int`: base preference period index
    - `z::Int`: base price period index
    - `τ::Int`: welfare comparison period index

    # Keyword Arguments
    - `χ`: curvature parameter of the expenditure function (0 < χ < γ < 1)
    - `η`: taste parameter for goods in consumption
    - `γ`: substitution parameter (χ < γ < 1)
    - `Pst::AbstractVector`: time series of service prices
    - `Pgt::AbstractVector`: time series of goods prices
    """              
    sghat_t_τ = sg_t_x(t, τ; χ=χ, η=η, γ=γ, Pst=Pst, Pgt=Pgt)
    sghat_t_z = sg_t_x(t, z; χ=χ, η=η, γ=γ, Pst=Pst, Pgt=Pgt)
    v_z_τ     = cal_v_t_x(z, τ; χ=χ, Pst=Pst)

    return (1/χ - 1 - sghat_t_τ/γ) * v_z_τ - (1/χ - 1 - sghat_t_z/γ)
end

function mhat_t_z_τ(t::Int, z::Int, τ::Int;  χ, η, γ, ν_t::AbstractVector, Pst::AbstractVector,  Pgt::AbstractVector,  yt::AbstractVector)
    """
    mhat_t_z_τ(t, z, τ; χ, η, γ, ν_t, Pst, Pgt, yt)

    Compute the equivalent variation measure of welfare, expressed in units of net income at
    period `τ`.

    This is the generalized Fisher-Shell welfare index that adjusts net income `y(τ)` for the
    difference between actual preferences (period `t`) and base prices (period `z`). The
    measure is defined as:

    mhat(t, z, τ) = y(τ) + B(t, z, τ) * ehat(t, z)

    where `ehat(t, z) = e(t, z)` is the compensating expenditure evaluated at preferences `t`
    and prices at period `z`, and `B(t, z, τ)` is the bias correction coefficient from
    `B_t_z_τ`.

    When `t = z = τ`, this reduces to net income `y(τ)`.

    # Arguments
    - `t::Int`: base preference period index
    - `z::Int`: base price period index
    - `τ::Int`: welfare comparison period index

    # Keyword Arguments
    - `χ`: curvature parameter of the expenditure function (0 < χ < γ < 1)
    - `η`: taste parameter for goods in consumption
    - `γ`: substitution parameter (χ < γ < 1)
    - `ν_t::AbstractVector`: time series of the marginal utility of wealth
    - `Pst::AbstractVector`: time series of service prices
    - `Pgt::AbstractVector`: time series of goods prices
    - `yt::AbstractVector`: time series of net income per capita
    """

    ehat_t_z = e_t_x(t, z; χ=χ, ν_t=ν_t, Pst=Pst)
    B        = B_t_z_τ(t, z, τ; χ=χ, η=η, γ=γ, Pst=Pst, Pgt=Pgt)

    return yt[τ] + B * ehat_t_z
end

#Functions Used to Compute Indices for Real Consumption Expenditure
g_D_e(sg,g_cg,g_cs)                 = ( sg*g_cg + (1 - sg)*g_cs  ) 
g_e_t_x(t,x;χ,ν_t,Pst,sg,g_cg,g_cs) = g_D_e(sg,g_cg,g_cs)*((e_t_x(x,x; χ, ν_t, Pst)/e_t_x(t,x; χ, ν_t, Pst))*(Pst[t]/Pst[x]))^(χ) 
aux_e_t_x(t,x;χ,ν_t,Pst)            = ( (e_t_x(x,x; χ, ν_t, Pst) / e_t_x(t,x; χ,ν_t, Pst) )*(Pst[t]/Pst[x]))^(χ)

function dev_t_z_τ(t::Int, z::Int, τ::Int; χ, η, γ, ν_t::AbstractVector, Pst::AbstractVector, Pgt::AbstractVector, yt::AbstractVector, g_pg, g_ps)
        """
    dev_t_z_τ(t, z, τ; χ, η, γ, ν_t, Pst, Pgt, yt, g_pg, g_ps)
    Compute the deviation term for the generalized Fisher-Shell growth rate.

    The deviation is defined as:

    dev = (se_τ * sg_τ - v(z,τ) * se(t,z) * sg(t,τ)) * g_pg + (se_τ * (1 - sg_τ) - v(z,τ) * se(t,z) * (1 - sg(t,τ))) * g_ps

    where:
    - `se_τ     = e(τ,τ) / y(τ)`       – expenditure share at base (τ,τ)
    - `sg_τ     = sg(τ,τ)`             – goods share at base
    - `se(t,z)  = e(t,z) / y(z)`       – scaled expenditure at preferences `t`, prices `z`
    - `sg(t,τ)  = sg(t,τ)`             – goods share at preferences `t`, prices `τ`
    - `v(z,τ)   = (Ps[z]/Ps[τ])^(χ/(1-χ))` – service-price scaling factor

    # Arguments
    - `t::Int`: base preference period index
    - `z::Int`: base price period index
    - `τ::Int`: welfare comparison period index

    # Keyword Arguments
    - `χ`: parameter governing income effects  (0 < χ < γ < 1)
    - `η`: taste parameter for goods in consumption
    - `γ`: substitution parameter (χ < γ < 1)
    - `ν_t::AbstractVector`: time series of the marginal utility of wealth
    - `Pst::AbstractVector`: time series of service prices
    - `Pgt::AbstractVector`: time series of goods prices
    - `yt::AbstractVector`: time series of net income per capita
    - `g_pg`: growth rate of the goods price
    - `g_ps`: growth rate of the services price
    """
    se_τ      = e_t_x(τ, τ; χ=χ, ν_t=ν_t, Pst=Pst) / yt[τ]
    sg_τ      = sg_t_x(τ, τ; χ=χ, η=η, γ=γ, Pst=Pst, Pgt=Pgt)
    se_t_z_τ  = e_t_x(t, z; χ=χ, ν_t=ν_t, Pst=Pst) / yt[τ]
    #se_t_z_τ  = e_t_x(t, z; χ=χ, ν_t=ν_t, Pst=Pst) / yt[z]
    cal_v_z_τ = cal_v_t_x(z, τ; χ=χ, Pst=Pst)
    sg_t_τ    = sg_t_x(t, τ; χ=χ, η=η, γ=γ, Pst=Pst, Pgt=Pgt)
    return ( se_τ*sg_τ - cal_v_z_τ*se_t_z_τ*sg_t_τ ) * g_pg  + ( se_τ*(1 - sg_τ) - cal_v_z_τ*se_t_z_τ*(1 - sg_t_τ) ) * g_ps
end
#------------------------------------------------------------------------------

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
#TABLE 1
growth_rates = DataFrame(
Parameter    = ["θ (Capital Share)","ρ (Discount Rate)","δ (Depreciation Rate)", "g_n (Population)", "g_h (Human Capital)", 
                "g_calAx (Investment TFP)", "g_As (Services TFP)", "g_Ag (Goods TFP)", "ζ (Price Wedge)"],

Value       = [round(θ, digits=4), round(ρ, digits=4), round(δ, digits=4), round(g_n, digits=4), round(g_h, digits=4), 
               round(g_calAx, digits=4), round(g_As, digits=4), round(g_Ag, digits=4), round(g_wedge_Pg_Ps, digits=4)])     

io = IOBuffer()
println(io, """
$(repeat("=", 60))
Calibrated Parameters and Growth Rates
$(repeat("=", 60))
""")
pretty_table(io, growth_rates, alignment=:l)
println(io, "$(repeat("=", 60))\n")
table_string = String(take!(io))
println(table_string)
#------------------------------------------------------------------------------

### Define the Function that Simulates the Model
### Note 1: We set the model to be in ABGP from 1980 onwards
#------------------------------------------------------------------------------
function sim_model_full(params::Vector{Float64};Pg_t::Vector{Float64},Ps_t::Vector{Float64},
                                        θ::Float64,ρ::Float64,δ::Float64,g_calAx::Float64,g_As::Float64,g_Ag::Float64,
                                        g_h::Float64,g_l::Float64,g_n::Float64)
        
    χ = params[1]
    η = params[2]
    γ = params[3]

    # Impose the Parameter Restriction From HRV
    if !(1 > γ > χ > 0)
        eq = Dict{String, Vector{Float64}}()
        for var in ["Yt" , "Xt" , "Et" , "Kt", 
                    "yt",  "et" , "kt" , "ht", "sst", "sgt", 
                    "cst", "cgt", "Cst", "Cgt","Xst", "Xgt", "Pst", "Pgt", "Kst", "Kgt", 
                    "Lt" , "Nt" , "Ht" ,
                    "Lst", "Lgt", "Yst", "Ygt", "Wt", "Rt", 
                    "Agt", "Ast", "Axt", "calAxt", "calAxhat"]
            eq[var] = fill(1e6, length(1980:2023))
        end
        return eq
    end


    N_t         = 1.000.*exp.(g_n.*((1980:2023) .- 1980))
    h_t         = 1.000.*exp.(g_h.*((1980:2023) .- 1980))
    L_t         = 1.000.*exp.(g_l.*((1980:2023) .- 1980))
    g_abgp      = g_calAx/(1-θ) + g_h
    Ag_t        = 1.000.*exp.( g_Ag.*((1980:2023)    .- 1980))
    As_t        = 1.000.*exp.( g_As.*((1980:2023)    .- 1980))
    calAx_t     = 1.000.*exp.( g_calAx.*((1980:2023) .- 1980))

    calAx_hat_t = calAx_t.*(h_t.^(1-θ))

    g_ps        = g_calAx - g_As 
    n           = g_n

    #Compute the k0,e0 that puts the Dynamic System in ABGP from t0
    k0_max     = 100.0 
    k0_min     = 0.001 
    k0_mid     = []
    err_ge     = Inf
    iter       = 0 

    while abs(err_ge) .> 1e-10 && iter < 100
        k0_mid = (k0_max + k0_min)/2
        err_ge = dedt_e(calAx_hat_t[1],k0_mid,θ,ρ,δ,χ,g_ps) - g_abgp

        if err_ge > 0
            k0_min = k0_mid
        else
            k0_max = k0_mid
        end
        k0_mid = (k0_max + k0_min)/2
        iter += 1
    end
    k0 = k0_mid

    e0_max = k0*10 
    e0_min = k0/10 
    e0_mid = []

    err_gk = Inf
    iter   = 0 

    #Find e0
    while abs(err_gk) .> 1e-10 && iter < 100
        e0_mid = (e0_max + e0_min)/2
        err_gk = dkdt(calAx_hat_t[1],k0,e0_mid,δ,n,θ)/k0  - g_abgp

        if err_gk > 0
            e0_min = e0_mid
        else
            e0_max = e0_mid
        end

        iter += 1
    end

    e0   = e0_mid

    vars = ["Yt" , "Xt" , "Et" , "Kt" , "yt",  "et" , "kt" , "ht", "sst", "sgt", 
            "cst", "cgt", "Cst", "Cgt","Xst", "Xgt", "Pst", "Pgt", "Kst", "Kgt", 
            "Lt" , "Nt" , "Ht" , "Lst", "Lgt", "Yst", "Ygt", "Wt", "Rt", 
            "Agt", "Ast", "Axt", "calAxt", "calAxhat"]

    eq = Dict{String, Vector{Any}}()
    for var in vars
        eq[var] = Vector{Any}()
    end

    push!(eq["et"] , e0)
    push!(eq["kt"] , k0)
    push!(eq["ht"] , h_t[1])
    push!(eq["Agt"], Ag_t[1])
    push!(eq["Ast"], As_t[1])
    push!(eq["calAxt"],calAx_t[1])
    push!(eq["calAxhat"],calAx_hat_t[1])
    push!(eq["Pst"], Ps_t[1])
    push!(eq["Pgt"], Pg_t[1])
    push!(eq["Lt"], L_t[1])
    push!(eq["Nt"], N_t[1])
    push!(eq["yt"], y(eq["calAxt"][1],eq["kt"][1],eq["ht"][1],θ))
    
    t = 1
    for t in 1:(2024-1980) 

        sg_t = sg(eq["Pgt"][t],eq["et"][t],eq["Pst"][t],η,χ,γ)
        ss_t = 1 - sg_t
        push!(eq["sst"],ss_t)
        push!(eq["sgt"],sg_t)

        push!(eq["Yt"], eq["yt"][t]*eq["Nt"][t])
        push!(eq["Kt"], eq["kt"][t]*eq["Nt"][t])
        push!(eq["Et"], eq["et"][t]*eq["Nt"][t])
        X_t   = eq["Yt"][t] - eq["Et"][t]
        push!(eq["Xt"],X_t)

        Xs_t   = eq["Xt"][t]/((1 + (ωx/(1-ωx))*(eq["Agt"][t]/eq["Ast"][t])^(εx)))
        push!(eq["Xst"],Xs_t)
        Xg_t   = X_t - Xs_t
        push!(eq["Xgt"],Xg_t)

        cs_t   = (eq["sst"][t]*eq["et"][t])/eq["Pst"][t] 
        cg_t   = (eq["sgt"][t]*eq["et"][t])/eq["Pgt"][t]  

        push!(eq["cst"],cs_t)
        push!(eq["cgt"],cg_t)
        push!(eq["Cst"],eq["cst"][t]*eq["Nt"][t])
        push!(eq["Cgt"],eq["cgt"][t]*eq["Nt"][t])

        push!(eq["Rt"],θ*eq["calAxt"][t]*eq["Kt"][t]^(θ-1))
        push!(eq["Wt"],(1-θ)*eq["calAxt"][t]*eq["Kt"][t]^θ)

        Ls = (eq["Xt"][t]/eq["Yt"][t])*(1/(((eq["Pgt"][t]*eq["Xgt"][t])/(eq["Pst"][t]*eq["Xst"][t]))+1))+(eq["Et"][t]/eq["Yt"][t])*(1/(((eq["Pgt"][t]*eq["Cgt"][t])  /(eq["Pst"][t]*eq["Cst"][t]))    +1))   
        if !isfinite(Ls) || Ls < 0 || Ls > 1
            return Dict("sgt" => fill(1e6, 44))
        end
        Lg = 1-Ls

        push!(eq["Lst"],Ls)
        push!(eq["Lgt"],Lg)

        Ks = eq["Kt"][t]*Ls
        Kg = eq["Kt"][t]*Lg

        push!(eq["Kst"],Ks)
        push!(eq["Kgt"],Kg)

        Ys = eq["Ast"][t]*(Ks^θ)*(Ls^(1-θ))
        Yg = eq["Agt"][t]*(Kg^θ)*(Lg^(1-θ))

        push!(eq["Yst"],Ys)
        push!(eq["Ygt"],Yg)

        # Break loop here
        if t == 44
            break
        end

        push!(eq["Pst"], Ps_t[t+1])
        push!(eq["Pgt"], Pg_t[t+1])

        push!(eq["calAxt"],calAx_t[t+1])
        push!(eq["calAxhat"],calAx_hat_t[t+1])

        push!(eq["Agt"], Ag_t[t+1])
        push!(eq["Ast"], As_t[t+1])

        push!(eq["ht"], h_t[t+1])
        push!(eq["Lt"], L_t[t+1])
        push!(eq["Nt"], N_t[t+1])

        g_e  = dedt_e(eq["calAxhat"][t],eq["kt"][t],θ,ρ,δ,χ,g_ps)  
        if !isfinite(g_e) || abs(g_e) > 1.0
            return Dict("sgt" => fill(1e6, 44))
        end
        e    = eq["et"][t]*exp(g_e)
        push!(eq["et"],e)

        g_k  = dkdt(eq["calAxhat"][t],eq["kt"][t],eq["et"][t],δ,n,θ)/eq["kt"][t]
        if !isfinite(g_k) || abs(g_k) > 1.0
            return Dict("sgt" => fill(1e6, 44))
        end
        k_t  = eq["kt"][t]*exp(g_k)

        push!(eq["kt"], k_t)
        push!(eq["yt"], y(eq["calAxt"][t+1],eq["kt"][t+1],eq["ht"][t+1],θ))
    end

    σ_t   = (1-γ) .- ((η.*((eq["Pgt"]./eq["Pst"]).^γ))./((eq["et"]./eq["Pst"]).^χ .- η.*((eq["Pgt"]./eq["Pst"]).^γ))).*(γ-χ)

    if any(x -> x < 0, σ_t)
        eq["sgt"] = fill(1e6, length(eq["sgt"]))
    end

    return eq
end
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
if run_smm
    println("Running SMM Estimation...")
    # Calibrate the Parameters of the Model Using a SMM
    ##Note 2: Recall that HRV calibration is: χ=0.550,η=0.440,γ=0.690


    


    ### Pass Some Parameters to the Simulation Function
    sim_model(params) = sim_model_full(params;Pg_t=Pg_t,Ps_t=Ps_t,θ=θ,ρ=ρ,δ=δ,g_calAx=g_calAx,g_As=g_As,g_Ag=g_Ag,g_h=g_h,g_l=g_l,g_n=g_n)

    ### Define the SSE Function to Minimize
    function SSE(params;cons_exp_share_data)
        sim_eq     = sim_model(params)
        # Check for NaN or Inf in sim_eq["sgt"]
        if any(x -> isnan(x) || isinf(x), sim_eq["sgt"])
            @show "NaN or Inf detected in sim_eq['sgt']"
            return 9e12
        else SSE =  1e20*( ( sum( 100.0.*(sim_eq["sgt"] .- cons_exp_share_data) ).^2 ) ) #Sum of Squared Errors 
            #@show SSE
            return SSE
        end
    end

    ### Define the Parameter Space for Optimization
    ParamSpace     = [(0.001,1.000),  #χ  
                    (0.001,1.000),  #η
                    (0.001,1.000)]  #γ

    ### Define different sees to start the optimization process
    ### Note 3: We will pick the parameter set with the lowest SSE across different seeds 
    seeds = [1234, 5678]
    optimization_results = OrderedDict()

    ### Run the Optimization for  Different Seeds
    for seed in seeds
        println("\n" * "="^80)
        println("Running optimization with seed: $seed")
        println("="^80)
        
        opt_problem = bbsetup(params -> SSE(params;cons_exp_share_data=VAC_GOOD_SHARE);
                            SearchRange=ParamSpace,
                            TraceMode=:compact,  # Use compact to reduce output clutter
                            Method                  =:adaptive_de_rand_1_bin,
                            PopulationSize          = 100,
                            MaxFuncEvals            = 2_000_000,
                            FitnessTolerance        = 1e-20,
                            MaxStepsWithoutProgress = 500_000,
                            rng  = MersenneTwister(seed))
        
        opt_results = bboptimize(opt_problem)
        
        # Store results
        optimization_results[seed] = OrderedDict(
            "params" => best_candidate(opt_results),
            "fitness" => best_fitness(opt_results),
            "seed"    => seed,
            "χ"       => best_candidate(opt_results)[1],
            "η"       => best_candidate(opt_results)[2],
            "γ"       => best_candidate(opt_results)[3]
        )
        
        println("Seed $seed: Fitness = $(best_fitness(opt_results))")
        println("Parameters: χ=$(round(best_candidate(opt_results)[1], digits=6)), "*
                "η=$(round(best_candidate(opt_results)[2], digits=6)), "*
                "γ=$(round(best_candidate(opt_results)[3], digits=6))")
    end

    ### Find the best result across all seeds
    best_seed   = argmin([optimization_results[s]["fitness"] for s in seeds])
    best_result = optimization_results[seeds[best_seed]]

    println("""
    $(repeat("=", 80))
    BEST RESULT ACROSS ALL SEEDS
    $(repeat("=", 80))
    Best seed:    $(best_result["seed"])
    Best fitness: $(best_result["fitness"])
    Best parameters:
    χ = $(round(best_result["χ"], digits=6))
    η = $(round(best_result["η"], digits=6))
    γ = $(round(best_result["γ"], digits=6))
    $(repeat("=", 80))
 """)

    ### Extract best parameters
    χ, η, γ = best_result["params"]
    χ, η, γ = round(χ,digits=3), round(η,digits=3), round(γ,digits=3)
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------
    #TABLE 2 - Calibrated Parameters
    println("""
    $(repeat("=", 80))
    Preference Parameters Estimated via SMM
    $(repeat("=", 80))

    Parameters:
    χ = $(round(best_result["χ"], digits=6))
    η = $(round(best_result["η"], digits=6))
    γ = $(round(best_result["γ"], digits=6))
    $(repeat("=", 80))
    """)
    #------------------------------------------------------------------------------
elseif run_smm_fixed_chi
    println("Running SMM Estimation with fixed χ = $chi_fixed ...")
    ##Note: χ is held at chi_fixed; only η and γ are optimised

    ### Pass Some Parameters to the Simulation Function (χ is fixed)
    sim_model_fixed_chi(params) = sim_model_full([chi_fixed; params]; Pg_t=Pg_t, Ps_t=Ps_t,
                                    θ=θ, ρ=ρ, δ=δ, g_calAx=g_calAx, g_As=g_As, g_Ag=g_Ag,
                                    g_h=g_h, g_l=g_l, g_n=g_n)

    ### Define the SSE Function to Minimize
    function SSE_fixed_chi(params; cons_exp_share_data)
        sim_eq = sim_model_fixed_chi(params)
        if any(x -> isnan(x) || isinf(x), sim_eq["sgt"])
            return 9e12
        else
            return 1e20 * (sum(100.0 .* (sim_eq["sgt"] .- cons_exp_share_data)).^2)
        end
    end

    ### Search over η and γ only; enforce γ > chi_fixed via the lower bound
    ParamSpace_fixed_chi = [(0.001, 1.000),                 # η
                            (chi_fixed + 0.001, 1.000)]     # γ  (must exceed χ)

    seeds = [1234, 5678, 9101, 1121, 3141, 5926, 5358, 9793, 2384, 6264]
    optimization_results_fc = OrderedDict()

    for seed in seeds
        println("\n" * "="^80)
        println("Running optimization with seed: $seed  (χ fixed at $chi_fixed)")
        println("="^80)

        opt_problem = bbsetup(params -> SSE_fixed_chi(params; cons_exp_share_data=VAC_GOOD_SHARE);
                              SearchRange             = ParamSpace_fixed_chi,
                              TraceMode               = :compact,
                              Method                  = :adaptive_de_rand_1_bin,
                              PopulationSize          = 100,
                              MaxFuncEvals            = 500_000,
                              FitnessTolerance        = 1e-20,
                              MaxStepsWithoutProgress = 100_000,
                              rng = MersenneTwister(seed))

        opt_results = bboptimize(opt_problem)

        optimization_results_fc[seed] = OrderedDict(
            "params"  => best_candidate(opt_results),
            "fitness" => best_fitness(opt_results),
            "seed"    => seed,
            "η"       => best_candidate(opt_results)[1],
            "γ"       => best_candidate(opt_results)[2]
        )

        println("Seed $seed: Fitness = $(best_fitness(opt_results))")
        println("Parameters: η=$(round(best_candidate(opt_results)[1], digits=6)), " *
                "γ=$(round(best_candidate(opt_results)[2], digits=6))")
    end

    best_seed_fc   = argmin([optimization_results_fc[s]["fitness"] for s in seeds])
    best_result_fc = optimization_results_fc[seeds[best_seed_fc]]

    println("""
    $(repeat("=", 80))
    BEST RESULT ACROSS ALL SEEDS (χ fixed at $chi_fixed)
    $(repeat("=", 80))
    Best seed:    $(best_result_fc["seed"])
    Best fitness: $(best_result_fc["fitness"])
    Best parameters:
    χ = $chi_fixed  (fixed)
    η = $(round(best_result_fc["η"], digits=6))
    γ = $(round(best_result_fc["γ"], digits=6))
    $(repeat("=", 80))
    """)

    χ        = chi_fixed
    η, γ     = best_result_fc["params"]
    χ, η, γ  = round(χ, digits=3), round(η, digits=3), round(γ, digits=3)
    #------------------------------------------------------------------------------
else
    println("Skipping SMM estimation, using pre-estimated parameters...")
    #=
    ================================================================================
    BEST RESULT ACROSS ALL SEEDS
    ================================================================================
    Best seed:    1121
    Best fitness: 8.697191480061655e-9
    Best parameters:
    χ = 0.360129
    η = 0.251054
    γ = 0.714119
    ================================================================================
    =#
    χ, η, γ = 0.360129,0.251054,0.714119
    #χ, η, γ  = 0.010000,0.251054,0.714119
    χ, η, γ  = round(χ,digits=3), round(η,digits=3), round(γ,digits=3)
end
#------------------------------------------------------------------------------


#-------------------------------------------------------------------
### Check that the Parameter Restrictions are Indeed Satisfied 
1 > γ > χ > 0
((1-γ)/(1-χ))
#-------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
# Simulate the Model with the Calibrated Parameters
sim_eq = sim_model_full([χ,η,γ];Pg_t=Pg_t,Ps_t=Ps_t,θ=θ,ρ=ρ,δ=δ,g_calAx=g_calAx,g_As=g_As,g_Ag=g_Ag,g_h=g_h,g_l=g_l,g_n=g_n)

#Compute the Marginal Value of Capital 
ν_t    = (sim_eq["et"].^(χ-1))./(sim_eq["Pst"].^χ)

#Check that the Elasticity of Substitution is Always Positive
σ_t    = (1-γ) .- ((η.*((sim_eq["Pgt"]./sim_eq["Pst"]).^γ))./((sim_eq["et"]./sim_eq["Pst"]).^χ .- η.*((sim_eq["Pgt"]./sim_eq["Pst"]).^γ))).*(γ-χ)

if any(x -> x < 0, σ_t)
  println("Warning: The elasticity of substitution is negative at some points.")
else
  println("The elasticity of substitution is positive at all points.")
end

if any(x -> x > ((1-γ)/(1-χ)) , sim_eq["sgt"]) 
    println("Warning: The Indirect Utility Function is not valid for this parameter set")
else
    println("The Indirect Utility Function IS VALID for this parameter set")
end

### Check that Assumption 3 is satisfied
if all(((1-γ)/(1-χ)) .> sim_eq["sgt"])
    println("""
    $(repeat("=", 60))
    ✓ Assumption 3 is SATISFIED
    $(repeat("=", 60))
    The condition (1-γ)/(1-χ) > sg(t) holds for all t
    
    Upper bound: (1-γ)/(1-χ) = $(round((1-γ)/(1-χ), digits=4))
    Max sg(t):   $(round(maximum(sim_eq["sgt"]), digits=4))
    Min sg(t):   $(round(minimum(sim_eq["sgt"]), digits=4))
    $(repeat("=", 60))
    """)
else
    println("""
    $(repeat("=", 60))
    ✗ WARNING: Assumption 3 is VIOLATED
    $(repeat("=", 60))
    The condition (1-γ)/(1-χ) > sg(t) does NOT hold for all t
    
    Upper bound: (1-γ)/(1-χ) = $(round((1-γ)/(1-χ), digits=4))
    Max sg(t):   $(round(maximum(sim_eq["sgt"]), digits=4))
    Min sg(t):   $(round(minimum(sim_eq["sgt"]), digits=4))
    Years where violated: $(findall(((1-γ)/(1-χ)) .<= sim_eq["sgt"]) .+ 1979)
    $(repeat("=", 60))
    """)
end
#-----------------------------------------------------------------------------------------------------

#---------------------------------------------------------------
#Check that the Economy is indeed in ABGP
gy  = log.(sim_eq["yt"][2:end]) .- log.(sim_eq["yt"][1:end-1])
gk  = log.(sim_eq["kt"][2:end]) .- log.(sim_eq["kt"][1:end-1])
ge  = log.(sim_eq["et"][2:end]) .- log.(sim_eq["et"][1:end-1])
gY  = log.(sim_eq["Yt"][2:end]) .- log.(sim_eq["Yt"][1:end-1])
gE  = log.(sim_eq["Et"][2:end]) .- log.(sim_eq["Et"][1:end-1])

g_pg = log.(sim_eq["Pgt"][2:end]) .- log.(sim_eq["Pgt"][1:end-1]) 
g_ps = log.(sim_eq["Pst"][2:end]) .- log.(sim_eq["Pst"][1:end-1]) 
#----------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
#CHARTS
# Define font settings for all plots
tickfont   = font(12)
guidefont  = font(12)
legendfont = font(12)
#-----------------------------------------------------------------------------------------------------

#---------------------------------------------------------------
### Figure 1 (a)
#Plot the Relative Price of Investment
X_P     = HRV_df[!,"X_TOT_P"]
C_P     = HRV_df[!,"C_TOT_P"]
rel_P_I = X_P ./ C_P

plot(1947:2023, rel_P_I ,
    ylabel="Relative Price of Investment\n(1947=1)",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1947:5:2025, 
    yticks=0.000:0.05:1.50,     
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir,"rel_Price_Investment.png"))
    println("Figure saved to: ", joinpath(figuresdir, "rel_Price_Investment.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 1 (b)
calA_X_I   = HRV_df[!,"calA_X_I_TD"]

plot(1947:2023, calA_X_I, ylabel="Effective Investment-Specific TFP\n(1947=1)",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1947:5:2025,ylim=(0.95, 2.35),   
    yticks=1.00:0.25:2.25, 
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir,"investment_TFP.png"))
    println("Figure saved to: ", joinpath(figuresdir, "investment_TFP.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 2 (a)
#Plot the Consumption Share
X_TOT      = HRV_df[!,"X_TOT"]
C_TOT      = HRV_df[!,"C_TOT"]
cons_share = C_TOT ./ (C_TOT+X_TOT)

plot(1947:2023, cons_share ,
    ylabel="Consumption Share in\nTotal Expenditure",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks = 1947:5:2025, 
    ylim   = (0.50, 1.00), 
    yticks = 0.40:0.10:1.00,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir,"consump_expend_share.png"))
    println("Figure saved to: ", joinpath(figuresdir, "consump_expend_share.png"))
end
#---------------------------------------------------------

#---------------------------------------------------------
### Figure 2 (b)
#Plot the Goods Consumption Value Added Share
VAC_GOOD_SHARE_long = HRV_df[HRV_df.year .>= 1947, "VAC_GOOD_S"]

plot(1947:2023, VAC_GOOD_SHARE_long,
    ylabel="Share of Goods in \n Total Consumption Value Added",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1947:5:2025,  
    yticks=0.000:0.05:0.400, 
    ylim=(0.08, 0.42), 
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    xrotation=45,
    legend=false,
    left_margin=5Plots.mm,  
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir,"goods_consumption_share.png"))
    println("Figure saved to: ", joinpath(figuresdir, "goods_consumption_share.png"))
end
#----------------------------------------------------------

#----------------------------------------------------------
### Figure 3 (a)
#Plot the Relative Price of Goods vs Services
plot(1980:2023,Pg_t_data./Ps_t_data, label=L"P_{g,t}/P_{s,t} \ \textrm{-} \ Data",
    ylabel="Relative Price of Goods vs Services \n (1980=1)",
    linestyle=:dash, lw=2,
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,legendfont=legendfont,
    legend=(0.55, 0.90),
    xticks=1980:5:2025, yticks=0.000:0.100:1.000,
    left_margin=6Plots.mm,xrotation=45,framestyle=:box)

#model WITHOUT wedge
plot!(1980:2023,Pg_t_undist./Ps_t_undist, label=L"P_{g,t}/P_{s,t} = A_{s,t}/A_{g,t}", linestyle=:dot, lw=2.00,color=:black)

#model WITH wedge
plot!(1980:2023,Pg_t./Ps_t, label=L"P_{g,t}/P_{s,t} = A_{s,t}/(A_{g,t}e^{\zeta t})", linestyle=:solid, lw=2.00,color=:black)

if save_figures 
    savefig(joinpath(figuresdir, "PgPs_Data_vs_HRV_vs_LV.png"))
    println("Figure saved to: ", joinpath(figuresdir, "PgPs_Data_vs_HRV_vs_LV.png"))
end

println("""
$(repeat("=", 60))
The exponential annual growth rate of the distortion is:
$(repeat("=", 60))
ζ = $(round(g_wedge_Pg_Ps, digits=4))
$(repeat("=", 60))
""")
#----------------------------------------------------------

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

#Compute the Marginal Value of Capital 
ν_t  = (sim_eq["et"].^(χ-1))./(sim_eq["Pst"].^χ)

#Compute Net Income Per Capita
m_t  = sim_eq["yt"] - δ*sim_eq["kt"]

#At Preferences of 2023
t_base             = 2023
t_prime            = t_base - 1980 + 1
τ                  = (1980:1:2023) .- 1980 .+ 1
m_τ                = sim_eq["yt"][τ]

mhat_2023_2023_τ = mhat_t_z_τ.(t_prime, t_prime, τ; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"],yt=sim_eq["yt"])
dev_2023_2023_τ  = dev_t_z_τ.( t_prime, t_prime, τ; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"],yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

#Sanity Checks
@assert isapprox(dev_2023_2023_τ[end], 0.0; atol=1e-8)
@assert isapprox(dev_2023_2023_τ[end], 0.0; atol=1e-8)

gD_2023_2023_τ   = ( m_τ[2:end] ./ mhat_2023_2023_τ[2:end] ) .* ( gD .+ g_n .+ dev_2023_2023_τ[2:end] )
FS_2023_2023_τ   = integrate_trap([0;gD_2023_2023_τ])

#Based on 2010
t_base             = 2010
t_prime            = t_base - 1980 + 1
τ_prime            = (1980:1:2023) .- 1980 .+ 1

mhat_2010_2010_τ = mhat_t_z_τ.(t_prime, t_prime, τ_prime;  χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"],  Pgt=sim_eq["Pgt"],  yt=sim_eq["yt"])
dev_2010_2010_τ  = dev_t_z_τ.(t_prime, t_prime, τ_prime; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

gD_2010_2010_τ   = ( m_τ[2:end] ./ mhat_2010_2010_τ[2:end] ) .* ( gD .+ g_n .+ dev_2010_2010_τ[2:end] )
FS_2010_2010_τ   = integrate_trap([0;gD_2010_2010_τ])

#Based on 2000
t_base           = 2000
t_prime          = t_base - 1980 + 1
τ_prime          = (1980:1:2023) .- 1980 .+ 1

mhat_2000_2000_τ = mhat_t_z_τ.(t_prime, t_prime, τ_prime;  χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"],  Pgt=sim_eq["Pgt"],  yt=sim_eq["yt"])
dev_2000_2000_τ  = dev_t_z_τ.(t_prime, t_prime, τ_prime; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

gD_2000_2000_τ   = ( m_τ[2:end] ./ mhat_2000_2000_τ[2:end] ) .* ( gD .+ g_n .+ dev_2000_2000_τ[2:end] )
FS_2000_2000_τ   = integrate_trap([0;gD_2000_2000_τ])

#Based on 1990
t_base             = 1990
t_prime            = t_base - 1980 + 1
τ_prime            = (1980:1:2023) .- 1980 .+ 1

mhat_1990_1990_τ = mhat_t_z_τ.(t_prime, t_prime, τ_prime;  χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"],  yt=sim_eq["yt"])
mhat_τ           = mhat_t_z_τ.(τ_prime, τ_prime, τ_prime;  χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"],  yt=sim_eq["yt"])
dev_1990_1990_τ  =   dev_t_z_τ.(t_prime, t_prime, τ_prime; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

gD_1990_1990_τ   = ( m_τ[2:end] ./ mhat_1990_1990_τ[2:end]) .* ( gD .+ g_n .+ dev_1990_1990_τ[2:end] )
FS_1990_1990_τ   = integrate_trap([0;gD_1990_1990_τ])

#At Preferences of 1980
t_base             = 1980
t_prime            = t_base - 1980 + 1
τ_prime            = (1980:1:2023) .- 1980 .+ 1

mhat_1980_1980_τ = mhat_t_z_τ.(t_prime,t_prime,τ; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"])
dev_1980_1980_τ  = dev_t_z_τ.( t_prime,t_prime,τ; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

gD_1980_1980_τ   = ( m_τ[2:end] ./ mhat_1980_1980_τ[2:end]) .* ( gD .+ g_n .+ dev_1980_1980_τ[2:end] )
FS_1980_1980_τ   = integrate_trap([0;gD_1980_1980_τ])

p = plot(1980:2023, FS_2023_2023_τ, label=L"\mathcal{P}_{2023,z} - \textrm{2023‑base \ Fisher‑Shell \ Index}",
    ylabel="Cumulative Growth", linestyle=:dash, lw=2.0,
    xticks=1980:5:2023, yticks=0.00:0.20:1.40,
    ylims=(0.00, 1.40),
    minorgrid=false, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.10, 0.920), legendfont=legendfont,
    xrotation=45, framestyle=:box)

plot!(p, 1980:2023, FS_1980_1980_τ, label=L"\mathcal{L}_{1980,z} - \textrm{1980‑base \  Fisher‑Shell \ Index}",
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
FS Index (base 1980) - Chained Index: $(round(FS_1980_1980_τ[end] - FS[end], digits=3))
FS Index (base 2023) - Chained Index: $(round(FS_2023_2023_τ[end] - FS[end], digits=3))
$(repeat("=", 60))
""")
#---------------------------------------------------------------
show(IOContext(stdout, :limit=>false), "text/plain",
     hcat(collect(1980:2023), FS, FS_1980_1980_τ, FS_2023_2023_τ))
#---------------------------------------------------------
### Figure 5 (b)
#Alternative Fisher-Shell Indices

###Price Chained Fisher-Shell Index
#At Preferences of 1980
t_base             = 1980
t_prime            = t_base - 1980 + 1
h                  = (1980:1:2023) .- 1980 .+ 1
m_h                = sim_eq["yt"][h]

mhat_1980_h_h      = mhat_t_z_τ.(t_prime, h, h; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"])
mhat_h             = mhat_t_z_τ.(τ_prime, h, h; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"])
dev_1980_h_h       =  dev_t_z_τ.(t_prime, h, h; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

#Sanity Checks
@assert isapprox(dev_1980_h_h[1], 0.0; atol=1e-8)

gD_1980_h_h        = ( m_h[2:end] ./ mhat_1980_h_h[2:end] ) .* ( gD .+ g_n .+ dev_1980_h_h[2:end] )
FS_1980_h_h        = integrate_trap([0;gD_1980_h_h])

#At Preferences of 2023
t_base             = 2023
t_prime            = t_base - 1980 + 1
h                  = (1980:1:2023) .- 1980 .+ 1
m_h                = sim_eq["yt"][h]

mhat_2023_h_h      = mhat_t_z_τ.(t_prime, h, h; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"])
mhat_h             = mhat_t_z_τ.(τ_prime, h, h; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"])
dev_2023_h_h       =  dev_t_z_τ.(t_prime, h, h; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

#Sanity Checks
@assert isapprox(dev_2023_h_h[end], 0.0; atol=1e-8)
gD_2023_h_h      = ( m_h[2:end] ./ mhat_2023_h_h[2:end]) .* (gD .+ g_n .+ dev_2023_h_h[2:end] )
FS_2023_h_h      = integrate_trap([0;gD_2023_h_h])

plot(1980:2023, FS_2023_h_h ,
    linestyle=:dash, lw=2.0, color=:black,
    xticks = 1980:5:2025,  
    label=L"\mathcal{\hat{P}}_{2023,z} - \textrm{2023‑base \  Price‑Chained \ FS \ Index}")

plot!(1980:2023, FS_1980_h_h ,
    linestyle=:dot, lw=2.0, color=:black,
    label=L"\mathcal{\hat{L}}_{1980,z} - \textrm{1980‑base \  Price‑Chained \ FS \ Index}",
    legend=(0.100, 0.920))

plot!(1980:2023, FS ,ylabel="Cumulative Growth",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1980:5:2024, yticks=0.0:0.2:1.4,
    ylim=(0.00, 1.40),
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    label=L"\mathcal{D}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    xrotation=45,
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir,"price_chained_FS.png"))
    println("Figure saved to: ", joinpath(figuresdir, "price_chained_FS.png"))
end
#---------------------------------------------------------

#---------------------------------------------------------
# Alternative Fisher-Shell Indices
#### Preferences Chained Fisher-Shell Index
## At Prices of 1980
z_base          = 1980
z_prime         = z_base - 1980 + 1

t               = (1980:1:2023) .- 1980 .+ 1
τ               = (1980:1:2023) .- 1980 .+ 1
m_τ             = sim_eq["yt"][τ]

mhat_t_1980_τ   = mhat_t_z_τ.(t, z_prime, τ ;  χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"],  Pgt=sim_eq["Pgt"],  yt=sim_eq["yt"])
dev_t_1980_τ    =  dev_t_z_τ.(t, z_prime, τ; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

gD_t_1980_τ     = ( m_τ[2:end] ./ mhat_t_1980_τ[2:end] ) .* ( gD .+ g_n .+ dev_t_1980_τ[2:end] )
FS_t_1980_τ     = integrate_trap([0;gD_t_1980_τ])

#At Prices of 2023
z_base          = 2023
z_prime         = z_base - 1980 + 1

t               = (1980:1:2023) .- 1980 .+ 1
τ               = (1980:1:2023) .- 1980 .+ 1
m_τ             = sim_eq["yt"][τ]

mhat_t_2023_τ   = mhat_t_z_τ.(t, z_prime, τ ;  χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"],  Pgt=sim_eq["Pgt"],  yt=sim_eq["yt"])
dev_t_2023_τ    =  dev_t_z_τ.(t, z_prime, τ; χ=χ, η=η, γ=γ, ν_t=ν_t, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"], yt=sim_eq["yt"], g_pg=g_pg[1], g_ps=g_ps[1])

gD_t_2023_τ      = ( m_τ[2:end] ./ mhat_t_2023_τ[2:end]) .* ( gD .+ g_n .+ dev_t_2023_τ[2:end] )
FS_t_2023_τ      = integrate_trap([0;gD_t_2023_τ])

plot(1980:2023, FS_t_2023_τ ,
    linestyle=:dash, lw=2.0,color=:black,
    label=L"\mathcal{\hat{P}}_{2023,z} - \textrm{2023‑base \  Pref‑Chained \ FS \ Index}")
plot!(1980:2023, FS_t_1980_τ ,
    linestyle=:dot, lw=2.0, color=:black,
    label=L"\mathcal{\hat{L}}_{1980,z} - \textrm{1980‑base \  Pref‑Chained \ FS \ Index}",
    legend=(0.100, 0.920))

plot!(1980:2023, FS , ylabel="Cumulative Growth",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1980:5:2024, yticks=0.0:0.2:1.4,
    ylim=(0.00, 1.40),
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    label=L"\mathcal{D}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    xrotation=45,
    framestyle=:box)

if save_figures
    savefig(joinpath(figuresdir,"prefs_chained_FS.png"))
    println("Figure saved to: ", joinpath(figuresdir, "prefs_chained_FS.png"))
end
#---------------------------------------------------------

#---------------------------------------------------------------
### Figure 6 (a)
# Alternative Expenditure Shares
# At Preferences of 1980
t_base             = 1980
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_1980_1980        = e_t_x.(t_prime,t_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"]) 
se_1980_z          = e_1980_1980./(sim_eq["yt"][ z_prime ])

#At Preferences of 2023
t_base             = 2023
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_2023_2023        = e_t_x.(t_prime,t_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"]) 
se_2023_z          = e_2023_2023./(sim_eq["yt"][ z_prime ])

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

#Goods Expenditure Share in ABGP vs a Fixed-Base Counterfactual
t_base     = 2023
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1
sg_2023_z  = sg_t_x.(t_prime,t_prime;χ=χ,η=η,γ=γ, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"])  

t_base     = 1980
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1
sg_1980_z  = sg_t_x.(t_prime,t_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"])    

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
g_FS_2023_2023_τ  = gD_2023_2023_τ
g_FS_2010_2010_τ  = gD_2010_2010_τ
g_FS_2000_2000_τ  = gD_2000_2000_τ
g_FS_1990_1990_τ  = gD_1990_1990_τ
g_FS_1980_1980_τ  = gD_1980_1980_τ

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

#------------------------------------------------------------------------
#APPENDIX A.1
#Decomposition of the Decline of the Growth Rate
πg  = sim_eq["Pgt"][2:end]./sim_eq["Pgt"][1:end-1] .- 1
πs  = sim_eq["Pst"][2:end]./sim_eq["Pst"][1:end-1] .- 1
πD  = s_e[2:end].*( sim_eq["sgt"][2:end].*πg .+ sim_eq["sst"][2:end].*πs ) 
(πg - πs)[end] 
s_e[end]
100*(πD[1] - πD[end])

#Alternative Decomposition of the Decline of the Growth Rate
x_t         = sim_eq["Xt"]./sim_eq["Nt"]
s_e         = sim_eq["Et"]./sim_eq["Yt"]
sg_t        = sim_eq["sgt"]
ss_t        = sim_eq["sst"]

g_cg        = sim_eq["cgt"][2:end]./sim_eq["cgt"][1:end-1] .- 1
g_cs        = sim_eq["cst"][2:end]./sim_eq["cst"][1:end-1] .- 1
g_x         = x_t[2:end]./x_t[1:end-1] .- 1

g_D         = s_e[2:end].*( sim_eq["sgt"][2:end].*g_cg .+ sim_eq["sst"][2:end].*g_cs  ) .+ (1 .- s_e[2:end]).*g_x
g_D[end]  - g_D[1]
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#APPENDIX A.5
#Cuantitative Evaluation of Condition (S_g,t0) in Step 1
#   s_g,t0 < min_{z∈(t0,t1)} [ π_s/(π_s-π_g) · exp(coeff·(t0-z)) ]
#   where coeff = χ/(1-χ)·π_s - γ·(π_s - π_g)

let
    z_range        = 1980:1:2023
    t0, t1         = 1980, 2023

    lhs            = sg_1980_z
    exponent_coeff = (χ/(1-χ))*πs[1] - γ*(πs[1] - πg[1])
    prefactor      = πs[1] / (πs[1] - πg[1])
    rhs_vec        = prefactor .* exp.(exponent_coeff .* (t0 .- z_range))
    rhs_min        = minimum(rhs_vec)

    println(""" 
        $(repeat("=", 60))
        CONDITION (S_g,t0)
        s_g,t0 < min_z [ π_s/(π_s-π_g) · exp(coeff·(t0-z)) ]
        
        $(repeat("=", 60))
        LHS
        s_g,t0                                    = $(round(lhs, digits=4))

        Exponent coefficient  [χ/(1-χ)·π_s - γ·(π_s-π_g)]
        = $(round(exponent_coeff, sigdigits=4))
        $(exponent_coeff < 0 ?
        "→ coeff < 0: min approached as z → t0⁺; RHS simplifies to π_s/(π_s-π_g)" :
        "→ coeff ≥ 0: min attained at z = t1 = $t1")

        RHS components
        Prefactor  π_s/(π_s-π_g)                  = $(round(prefactor, digits=4))
        min_z [ π_s/(π_s-π_g)·exp(coeff·(t0-z)) ] = $(round(rhs_min, digits=4))

        Conclusion
          Condition (S_g,t0) holds (LHS < RHS)?      $(lhs < rhs_min ? "YES ✓" : "NO ✗")
        $(repeat("=", 60))
        """)
end

#------------------------------------------------------------------------
#Cuantitative Evaluation of Condition (S_g,t1) in Step 1
#   s_g,t1 < min_{z∈(t0,t1)} [ π_s/(π_s-π_g) · exp(coeff·(t1-z)) ]
#   where coeff = χ/(1-χ)·π_s - γ·(π_s - π_g)

let
    z_range        = 1980:1:2023
    t0, t1         = 1980, 2023

    lhs            = sg_2023_z
    exponent_coeff = (χ/(1-χ))*πs[end] - γ*(πs[end] - πg[end])
    prefactor      = πs[end] / (πs[end] - πg[end])
    rhs_vec        = prefactor .* exp.(exponent_coeff .* (t1 .- z_range))
    rhs_min        = minimum(rhs_vec)

    println("""
        $(repeat("=", 60))
        CONDITION (S_g,t1)
          s_g,t1 < min_z [ π_s/(π_s-π_g) · exp(coeff·(t1-z)) ]
        $(repeat("=", 60))
        LHS
          s_g,t1                                    = $(round(lhs, digits=4))

        Exponent coefficient  [χ/(1-χ)·π_s - γ·(π_s-π_g)]  (at t1)
          = $(round(exponent_coeff, sigdigits=4))
          $(exponent_coeff < 0 ?
            "→ coeff < 0: min attained at z → t1⁻; RHS → π_s/(π_s-π_g)" :
            "→ coeff ≥ 0: min attained at z = t0 = $t0")

        RHS components
          Prefactor  π_s/(π_s-π_g)  (at t1)         = $(round(prefactor, digits=4))
          min_z [ π_s/(π_s-π_g)·exp(coeff·(t1-z)) ] = $(round(rhs_min, digits=4))

        Conclusion
          Condition (S_g,t1) holds (LHS < RHS)?      $(lhs < rhs_min ? "YES ✓" : "NO ✗")
        $(repeat("=", 60))
        """)
end
#------------------------------------------------------------------------
#APPENDIX C
#Real Consumption Expenditure Indices
gD_e_z     = g_D_e.(sim_eq["sgt"][2:end],g_cg,g_cs)
D_e_z      = integrate_trap([0;gD_e_z .+ g_n])

#Base 2023
t_base     = 2023
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1

g_e_2023_z = gD_e_z.*aux_e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])[2:end] 
P_e_2023_z = integrate_trap([0;g_e_2023_z .+ g_n])

#Base 2010
t_base     = 2010
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1

g_e_2010_z = gD_e_z.*aux_e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])[2:end] 
P_e_2010_z = integrate_trap([0;g_e_2010_z .+ g_n])

#Base 2000
t_base     = 2000
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1

g_e_2000_z = gD_e_z.*aux_e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])[2:end] 
P_e_2000_z = integrate_trap([0;g_e_2000_z .+ g_n])

#Base 1990
t_base     = 1990
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1

g_e_1990_z = gD_e_z.*aux_e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])[2:end] 
P_e_1990_z = integrate_trap([0;g_e_1990_z .+ g_n])

#Base 1980
t_base     = 1980
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1

g_e_1980_z = gD_e_z.*aux_e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])[2:end] 
L_e_1980_z = integrate_trap([0;g_e_1980_z .+ g_n])
#---------------------------------------------------------------

#---------------------------------------------------------------
### Real Consumption Expenditure Indices - D_e_z, P_e_2023_z,P_e_1980_z
 plot(1981:2023, g_e_2023_z, 
    label =L"g^{D}_{e}_{2023,z}  - \textrm{2023‑base \  FS \ Index}",
    linestyle=:dash, lw=2.5, 
    color=:black)

plot!(1981:2023, g_e_2010_z, 
    label =L"g^{D}_{e}_{2010,z}  - \textrm{2010‑base \  FS \ Index}",
    linestyle=:dash, lw=2.0, 
    color=:black)

plot!(1981:2023, g_e_2000_z, 
    label =L"g^{D}_{e}_{2000,z}  - \textrm{2000‑base \  FS \ Index}",
    linestyle=:dash, lw=1.5, 
    color=:black)

plot!(1981:2023, g_e_1990_z, 
    label =L"g^{D}_{e}_{1990,z}  - \textrm{1990‑base \  FS \ Index}",
    linestyle=:dash, lw=1.0, 
    color=:black)
    
plot!(1981:2023, g_e_1980_z, 
    label=L"g^{D}_{e}_{1980,z} - \textrm{1980‑base \  FS \ Index}",
    linestyle=:dash, lw=0.5, 
    color=:black)
plot!(1981:2023, gD_e_z, 
    label =L"g^{D}_{e}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel="Growth Rate of Real \n Consumption Expenditure",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1980:5:2025, 
    #yticks=0.0:0.2:1.4,
    ylim=(0.010, 0.025),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    legend=(0.100, 0.900),
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)
if save_figures
    savefig(joinpath(figuresdir, "growth_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "growth_consumption_expenditure_indices.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Real Consumption Expenditure Indices - D_e_z, P_e_2023_z,P_e_1980_z
plot(1980:2023, D_e_z, 
    label =L"\mathcal{D_{e}}_{z} \ \ \ \ \ \ - \textrm{Chained \ Divisia \ Index}",
    ylabel="Cumulative Growth in\nReal Consumption Expenditure",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1980:5:2025, 
    yticks=0.0:0.2:1.4,
    ylim=(0.00, 1.40),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    legend=(0.100, 0.900),
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

plot!(1980:2023, P_e_2023_z, 
    label=L"\mathcal{P_{e}}_{2023,z} - \textrm{2023‑base \  FS \ Index}",
    linestyle=:dash, lw=2.0, 
    color=:black)


plot!(1980:2023, L_e_1980_z, 
    label=L"\mathcal{L_{e}}_{1980,z} - \textrm{1980‑base \  FS \ Index}",
    linestyle=:dot, lw=2.0, 
    color=:black)
if save_figures
    savefig(joinpath(figuresdir, "consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "consumption_expenditure_indices.png"))
end
#---------------------------------------------------------------

#---------------------------------------------------------------
### Growth Rates of Real Consumption Expenditure Indices for Different Base Periods
 plot(1981:2023, g_e_2023_z, 
    label =L"g^{D}_{e}_{2023,z}",
    linestyle=:dash, lw=2.5, 
    color=:black)

plot!(1981:2023, g_e_2010_z, 
    label =L"g^{D}_{e}_{2010,z}",
    linestyle=:dash, lw=2.0, 
    color=:black)

plot!(1981:2023, g_e_2000_z, 
    label =L"g^{D}_{e}_{2000,z}",
    linestyle=:dash, lw=1.5, 
    color=:black)

plot!(1981:2023, g_e_1990_z, 
    label =L"g^{D}_{e}_{1990,z}",
    linestyle=:dash, lw=1.0, 
    color=:black)
    
plot!(1981:2023, g_e_1980_z, 
    label=L"g^{D}_{e}_{1980,z}",
    linestyle=:dash, lw=0.5, 
    color=:black)
plot!(1981:2023, gD_e_z, 
    label =L"g^{D}_{e}_{z}",
    ylabel="Growth Rate of Real \n Consumption Expenditure",
    linestyle=:solid, lw=2.0,
    minorgrid=false, color=:black,
    xticks=1980:5:2025, 
    #yticks=0.0:0.2:1.4,
    #ylim=(0.010, 0.025),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    legend=(0.800, 0.900),
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)
if save_figures
    savefig(joinpath(figuresdir, "growth_consumption_expenditure_indices.png"))
    println("Figure saved to: ", joinpath(figuresdir, "growth_consumption_expenditure_indices.png"))
end
#---------------------------------------------------------------


#-=============================================================================
# APPEND TO THE END OF  LV_main.jl   (PIGL)
#
# Konus INCOME-vs-SUBSTITUTION decomposition of the fixed-base bias in the
# consumption-expenditure quantity index.  Same object and split as the NHCES
# block (see that file for the derivation):
#
#   G_total       = current-base (2023) − base-base (1980) consumption qty index
#                 = P_e_2023_z[end] − L_e_1980_z[end]      (already built above)
#   G_income(ref) = [ s_g^H(p_ref,u_T) − s_g^H(p_ref,u_0) ] · Δln(P_g/P_s)
#                   with the compensated share's PRICE argument FROZEN at p_ref
#   G_subst       = G_total − G_income                      (share price response)
#
# PIGL compensated goods share at (utility period t, price period z):
#   sg_t_x(t, z; χ, η, γ, Pst, Pgt)      (already defined in this file; uses ν_t)
#
# Uses objects already defined in this file:
#   P_e_2023_z, L_e_1980_z, sg_t_x, χ, η, γ, ν_t, sim_eq
#-=============================================================================

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

  Pgt      = sim_eq["Pgt"]; Pst = sim_eq["Pst"]
  i0       = 1
  iT       = length(Pst)                                     # 1980 → 2023
  Δln_PgPs = (log(Pgt[iT]) - log(Pgt[i0])) - (log(Pst[iT]) - log(Pst[i0]))   # Δ ln(P_g/P_s)

  # PIGL compensated goods share at (utility period t, price period z)
  sgh(t, z)  = sg_t_x(t, z; χ=χ, η=η, γ=γ, Pst=Pst, Pgt=Pgt)

  # ---- TOTAL gap : current-base (2023) minus base-base (1980) quantity index ----
  # read off the exact Fisher–Shell consumption indices already constructed above
  G_total    = P_e_2023_z[end] - L_e_1980_z[end]

  # ---- INCOME part : freeze the compensated share's price argument at p_ref ------
  G_income_1980 = (sgh(iT, i0) - sgh(i0, i0)) * Δln_PgPs
  G_income_2023 = (sgh(iT, iT) - sgh(i0, iT)) * Δln_PgPs
  G_income_avg  = 0.5 * (G_income_1980 + G_income_2023)

  # ---- SUBSTITUTION part : residual ---------------------------------------------
  G_subst_avg = G_total - G_income_avg

  upsert_decomp_row(DECOMP_CSV, "PIGL", G_total, G_income_1980, G_income_2023, G_income_avg, G_subst_avg)

  println("""
  $(repeat("=", 68))
  PIGL — Konus income/substitution decomposition (consumption, 1980→2023)
  $(repeat("=", 68))
    Δ ln(P_g/P_s)                       = $(round(Δln_PgPs, digits=4))
  $(repeat("-", 68))
    G_total  (current-base − base-base)  = $(round(G_total, digits=4))
    G_income (freeze @1980 prices)       = $(round(G_income_1980, digits=4))
    G_income (freeze @2023 prices)       = $(round(G_income_2023, digits=4))
    G_income (average)                   = $(round(G_income_avg, digits=4))
    G_subst  (residual)                  = $(round(G_subst_avg, digits=4))
  $(repeat("-", 68))
    income sign  = $(sign(G_income_avg) > 0 ? "+ (as in HCD)" : "−")
    subst  sign  = $(sign(G_subst_avg)  > 0 ? "+" : "−")   |  flips total? $(sign(G_total) != sign(G_income_avg))
    → wrote row to $(DECOMP_CSV)
  $(repeat("=", 68))
  """)
end
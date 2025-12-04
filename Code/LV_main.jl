#------------------------------------------------------------------------------
#Replication code for: Chained Indices Unchained: Structural Transformation and the Welfare Foundations of Income Growth Measurement
#By:                   Omar Licandro and Juan I. Vizcaino
#This Version:         25/11/2025
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
Pkg.activate(@__DIR__)
Pkg.instantiate()

using XLSX, DataFrames, BlackBoxOptim , Plots; plotlyjs()
using Statistics,MathJaxRenderer, LaTeXStrings, LsqFit, Random, OrderedCollections, PrettyTables
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
data_range = "A1:DF78"

# Read the data from the Excel file
HRV_df = DataFrame(XLSX.readtable(joinpath(datadir, file_name), sheet_name, infer_eltypes=true))

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
L_2017         = HRV_df[HRV_df.year .== 2023,"LAB_TOT_QI"][1]

C_GOOD_P       = HRV_df[HRV_df.year .>= 1980, "C_GOOD_P"]./HRV_df[HRV_df.year .== 1980, "C_GOOD_P"]
C_SERV_P       = HRV_df[HRV_df.year .>= 1980, "C_SERV_P"]./HRV_df[HRV_df.year .== 1980, "C_SERV_P"]
X_TOT_P        = HRV_df[HRV_df.year .>= 1980, "X_TOT_P"] ./HRV_df[HRV_df.year .== 1980, "X_TOT_P"]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Define the Functions Characterizing the Equilibrium Conditions
# Recall that L = N*h, where N is the population and h is the avg efficiency units of labor per person

## Relative Prices
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
e_t_x(t,x; χ, ν_t, Pst)    = (ν_t[t] * (Pst[x]^χ))^(1/(χ-1))
sg_t_x(t,x;χ,η,γ,Pst,Pgt)  = η*( ( e_t_x(t,x; χ, ν_t, Pst)/Pst[x] )^(-χ) )*((Pgt[x]/Pst[x])^γ)
cal_v_t_x(t,x;χ,Pst)       = ( Pst[t]/Pst[x] )^(  (χ/(1-χ)) )

#Consumption Expenditure Indices
e_tilde_tz(t,z;χ,η,γ,e,Pgt,Pst) = (( (e[z]/Pst[z])^χ + (η*χ/γ)*( (Pgt[t]/Pst[t])^γ - (Pgt[z]/Pst[z])^γ ) )^(1/χ))* Pst[t]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Calibrate the growth rates of Agt,Ast,Axt (these values are from the HRV TFP series)
g_Ag        = log(Ag_2023/Ag_1980)/(2023-1980)
g_As        = log(As_2023/As_1980)/(2023-1980)
g_calAx     = log(calAx_2023/calAx_1980)/(2023-1980)
g_n         = (log(N_2023/N_1980))/(2023-1980)
g_l         = (log(L_2017/L_1980))/(2023-1980)

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
Parameter    = ["θ (Capital Share)",
                "ρ (Discount Rate)",
                "δ (Depreciation Rate)",
                "g_n (Population)", 
                "g_h (Human Capital)", 
                "g_calAx (Investment TFP)", 
                "g_As (Services TFP)", 
                "g_Ag (Goods TFP)",
                "ζ (Price Wedge)"],

Value       = [round(θ, digits=4),
               round(ρ, digits=4),
               round(δ, digits=4),   
               round(g_n, digits=4), 
               round(g_h, digits=4), 
               round(g_calAx, digits=4),
               round(g_As, digits=4), 
               round(g_Ag, digits=4),
               round(g_wedge_Pg_Ps, digits=4)])     

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

#------------------------------------------------------------------------------
# Calibrate the Parameters of the Model Using a SMM
##Note 1: We set the model to be in ABGP from 1980 onwards
##Note 2: Recall that HRV calibration is: χ=0.550,η=0.440,γ=0.690

### Define the Function that Simulates the Model
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
seeds = [1234, 5678, 9101, 1121, 3141, 5926, 5358, 9793, 2384, 6264]
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
χ, η, γ = round(χ,digits=3), round(η,digits=3), round(γ,digits=3)
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

#-------------------------------------------------------------------
### Check that the Parameter Restrictions are Indeed Satisfied 
1 > γ > χ > 0
((1-γ)/(1-χ))
#-------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
# Simulate the Model with the Calibrated Parameters
sim_eq = sim_model_full([χ,η,γ];Pg_t=Pg_t,Ps_t=Ps_t,θ=θ,ρ=ρ,δ=δ,g_calAx=g_calAx,g_As=g_As,g_Ag=g_Ag,g_h=g_h,g_l=g_l,g_n=g_n)

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
    
    Years where violated: $(findall(((1-γ)/(1-χ)) .<= sim_eq["sgt"]) .+ 1979)
    $(repeat("=", 60))
    """)
end

### Sanity Check: Do we get the same results as in the Analytical Solution?
g_abgp      = g_calAx/(1-θ) + g_h 
g_k         = g_abgp
g_pg        = g_calAx - g_Ag + g_wedge_Pg_Ps
g_ps        = g_calAx - g_As 
κ           = (ρ + δ + χ*g_ps + (1-χ)*g_k)/θ 

Ag_0        = 1.000
As_0        = 1.000
calAx_0     = 1.000

k0          = κ^(1/(θ-1))*(calAx_0^(1/(θ-1)))
sim_eq["kt"][1]

e0          = (κ - δ - g_n - g_k)*k0
sim_eq["et"][1]

y0         = κ*k0
sim_eq["yt"][1]

cg_t_an = η*(( e0/Ps_t[1] )^(1-χ))*(( Pg_t[1]/Ps_t[1])^γ ) 
sim_eq["cgt"][1]

cs_t_an = e0/Ps_t[1] - η*(( e0/Ps_t[1] )^(1-χ))*(( Pg_t[1]/Ps_t[1])^γ ) 
sim_eq["cst"][1]
#-----------------------------------------------------------------------------------------------------

#---------------------------------------------------------------
#Check that the Economy is indeed in ABGP
gy  = log.(sim_eq["yt"][2:end]) .- log.(sim_eq["yt"][1:end-1])
gk  = log.(sim_eq["kt"][2:end]) .- log.(sim_eq["kt"][1:end-1])
ge  = log.(sim_eq["et"][2:end]) .- log.(sim_eq["et"][1:end-1])
gY  = log.(sim_eq["Yt"][2:end]) .- log.(sim_eq["Yt"][1:end-1])
gE  = log.(sim_eq["Et"][2:end]) .- log.(sim_eq["Et"][1:end-1])

#Compute the Value of the Problem (in ABGP) -> Check this is correct (g_calAx=g_ABGP? sim_eq["et"][1] or sim_eq["Et"][1])
V_K0 = (1/χ)*(sim_eq["et"][1]/sim_eq["Pst"][1])*(1/(ρ - χ*(θ*g_calAx + g_As)))
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
    ylabel="Relative Price of Investment <br> (1947=1)", 
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.5, color=:black,
    xticks=1947:5:2025, 
    yticks=0.000:0.05:1.50,     
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

#savefig(joinpath(figuresdir,"rel_Price_Investment.png"))
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 1 (b)
calA_X_I   = HRV_df[!,"calA_X_I_TD"]

plot(1947:2023, calA_X_I, ylabel="Effective Investment-Specific TFP <br> (1947=1)", 
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.5, color=:black,
    xticks=1947:5:2025,ylim=(0.95, 2.35),   
    yticks=1.00:0.25:2.25, 
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

#savefig(joinpath(figuresdir,"investment_TFP.png"))
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 2 (a)
#Plot the Consumption Share
X_TOT      = HRV_df[!,"X_TOT"]
C_TOT      = HRV_df[!,"C_TOT"]
cons_share = C_TOT ./ (C_TOT+X_TOT)

plot(1947:2023, cons_share ,
    ylabel="Consumption Share in <br> Total Expenditure", 
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.50, color=:black,
    xticks = 1947:5:2025, 
    ylim   = (0.50, 1.00), 
    yticks = 0.40:0.10:1.00,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

#savefig(joinpath(figuresdir,"consump_expend_share.png"))
#---------------------------------------------------------

#---------------------------------------------------------
### Figure 2 (b)
#Plot the Goods Consumption Value Added Share
VAC_GOOD_SHARE_long = HRV_df[HRV_df.year .>= 1947, "VAC_GOOD_S"]

plot(1947:2023, VAC_GOOD_SHARE_long,
    ylabel="Share of Goods in Total Consumption <br> Value Added",    
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.5, color=:black,
    xticks=1947:5:2025,  
    yticks=0.000:0.05:0.400, 
    ylim=(0.08, 0.42), 
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    xrotation=45,
    legend=false,
    left_margin=5Plots.mm,  
    framestyle=:box)

#savefig(joinpath(figuresdir,"goods_consumption_share.png"))
#----------------------------------------------------------

#----------------------------------------------------------
### Figure 3 (a)
#Plot the Relative Price of Goods vs Services
plot(1980:2023,Pg_t_data./Ps_t_data, label="(Pg/Ps) - Data",
    ylabel="Relative Price of Goods vs Services <br> (Pg/Ps ; 1980=1)",
    linestyle=:dash, lw=2,
    minorgridalpha=0.5, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,legendfont=legendfont,
    legend=(0.60, 0.95),
    xticks=1980:5:2025, yticks=0.000:0.100:1.000,
    left_margin=6Plots.mm,xrotation=45,framestyle=:box)

#model without distortions
plot!(1980:2023,Pg_t_undist./Ps_t_undist, label="(Pg/Ps) = (As/Ag)", linestyle=:dot, lw=2.00,color=:black)

#model without distortions
plot!(1980:2023,Pg_t./Ps_t, label="(Pg/Ps) = (As/Ag*exp(ζt))", linestyle=:solid, lw=2.00,color=:black)

plot(1980:2023, wedge_Pg_Ps, label="Wedge ( (Pg/Ps)-Data / (Pg/Ps)-Model",
     ylabel="Relative Price of Goods vs Services <br> (Pg/Ps ; 1980=1)",
    linestyle=:dash, lw=2,xticks=1980:5:2025, 
    minorgridalpha=0.5, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.40, 0.95),legendfont=legendfont,
    left_margin=6Plots.mm,xrotation=45,framestyle=:box)

plot!(1980:2023, wedge_Pg_Ps_trend, label="Exponential Trend", lw=2, color=:black, linestyle=:solid)

#savefig(joinpath(figuresdir, "PgPs_Wedge_Fit.png"))

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
plot(year, sim_eq["sgt"], label="Model", linestyle=:solid, color=:black, lw=2.00, 
    ylabel="Share of Goods in Consumption <br> Expenditure", 
    #title="γ=$(round(γ, digits=3)), η=$(round(η, digits=3)), χ=$(round(χ, digits=3))",
    #title= L"\gamma".s,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,legendfont=legendfont,
    legend=(0.15, 0.95),
    xticks=1980:5:2025, yticks=0.0:0.02:0.80,xrotation=45)

plot!(year, VAC_GOOD_SHARE, label="Data", legend=(0.800, 0.950), linestyle=:dot, lw=2.00, 
    minorgrid=true, minorgridalpha=0.9, color=:black,left_margin=5Plots.mm,framestyle=:box)
#savefig(joinpath(figuresdir, "sg_t_Model_Fit.png"))
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

g_D         = s_e[2:end].*( sim_eq["sgt"][2:end].*g_cg .+ sim_eq["sst"][2:end].*g_cs  ) .+ (1 .- s_e[2:end]).*g_x
FS          = [0;cumsum( g_D .+ g_n )]

plot(1980:2023, FS , ylabel="GDP Index <br> (log scale;1980=0)", 
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.5, color=:black,
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
    legend=(0.15, 0.95))

#savefig(joinpath(figuresdir,"GDP_Model_vs_Data.png"))
#---------------------------------------------------------

#---------------------------------------------------------------
### Figure 4 (b)
g_FS          = g_D .+ g_n 
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

plot(1980:2023, g_GDP_long_ma10[1980-1947+1:end],
    label = "Data - 10-year Moving Average", 
    linestyle=:dot, lw=2.0,
    minorgridalpha=0.5, color=:black,
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
    lw=2.0, color=:black, minorgridalpha=0.5,
    xticks=1980:5:2024, yticks=0.010:0.010:0.50,
    ylim=(0.010,0.050),legend=(0.50, 0.95),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)

#savefig(joinpath(figuresdir,"GDP_Growth_Model_vs_Data_v2.png"))
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

#At Preferences of 1980
t_base             = 1980
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_1980_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"]) 
e_2023_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"]) 
e_1980_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])
se_1980_z          = e_1980_z./(sim_eq["yt"][ z_prime ])
sg_1980_z          = sg_t_x.(t_prime,z_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

#At Preferences of 2023
t_base             = 2023
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_2023_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])
se_2023_z          = e_2023_z./(sim_eq["yt"][ z_prime ])
sg_2023_z          = sg_t_x.(t_prime,z_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

#Compute the Equivalent Variation Measures
e_t_x(t,x; χ, ν_t, Pst)    = (ν_t[t] * (Pst[x]^χ))^(1/(χ-1))
sg_t_x(t,x;χ,η,γ,Pst,Pgt)  = η*( ( e_t_x(t,x; χ, ν_t, Pst)/Pst[x] )^(-χ) )*((Pgt[x]/Pst[x])^γ)
cal_v_t_x(t,x;χ,Pst)       = ( Pst[t]/Pst[x] )^(  (χ/(1-χ)) )

sg_z = sim_eq["sgt"]
ss_z = sim_eq["sst"]
se_z = s_e[1:end]

#Based on 2023
t_base      = 2023
t_prime     = t_base - 1980 + 1
z_prime     = (1980:1:2023) .- 1980 .+ 1
τ_base      = 2023
τ_prime     = τ_base - 1980 + 1

cal_v_z_2023  = cal_v_t_x.(z_prime,τ_prime;χ,Pst=sim_eq["Pst"])
se_2023_z     = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./sim_eq["yt"][z_prime]
sg_2023_z     = sg_t_x.(t_prime,τ_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_2023_z_2023 = (se_z.*sg_z .- cal_v_z_2023.*se_2023_z.*sg_2023_z).*g_pg .+ 
                  (se_z.*(1 .- sg_z) .- cal_v_z_2023.*se_2023_z.*(1 .- sg_2023_z)).*g_ps

#Based on 2010
t_base      = 2010
t_prime     = t_base - 1980 + 1
z_prime     = (1980:1:2023) .- 1980 .+ 1
τ_base      = 2010
τ_prime     = τ_base - 1980 + 1

cal_v_z_2010  = cal_v_t_x.(z_prime,τ_prime;χ,Pst=sim_eq["Pst"])
se_2010_z     = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./sim_eq["yt"][z_prime]
sg_2010_z     = sg_t_x.(t_prime,τ_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_2010_z_2010 = (se_z.*sg_z        .- cal_v_z_2010.*se_2010_z.*sg_2010_z       ).*g_pg .+ 
                  (se_z.*(1 .- sg_z) .- cal_v_z_2010.*se_2010_z.*(1 .- sg_2010_z)).*g_ps


#Based on 2000
t_base      = 2000
t_prime     = t_base - 1980 + 1
z_prime     = (1980:1:2023) .- 1980 .+ 1
τ_base      = 2000
τ_prime     = τ_base - 1980 + 1

cal_v_z_2000  = cal_v_t_x.(z_prime,τ_prime;χ,Pst=sim_eq["Pst"])
se_2000_z     = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./sim_eq["yt"][z_prime]
sg_2000_z     = sg_t_x.(t_prime,τ_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_2000_z_2000 = (se_z.*sg_z        .- cal_v_z_2000.*se_2000_z.*sg_2000_z       ).*g_pg .+ 
                  (se_z.*(1 .- sg_z) .- cal_v_z_2000.*se_2000_z.*(1 .- sg_2000_z)).*g_ps

#Based on 1990
t_base      = 1990
t_prime     = t_base - 1980 + 1
z_prime     = (1980:1:2023) .- 1980 .+ 1
τ_base      = 1990
τ_prime     = τ_base - 1980 + 1

cal_v_z_1990  = cal_v_t_x.(z_prime,τ_prime;χ,Pst=sim_eq["Pst"])
se_1990_z     = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./sim_eq["yt"][z_prime]
sg_1990_z     = sg_t_x.(t_prime,τ_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_1990_z_1990 = (se_z.*sg_z        .- cal_v_z_1990.*se_1990_z.*sg_1990_z       ).*g_pg .+ 
                  (se_z.*(1 .- sg_z) .- cal_v_z_1990.*se_1990_z.*(1 .- sg_1990_z)).*g_ps

#Based on 1980
t_base      = 1980
t_prime     = t_base - 1980 + 1
z_prime     = (1980:1:2023) .- 1980 .+ 1
τ_base      = 1980
τ_prime     = τ_base - 1980 + 1

cal_v_z_1980  = cal_v_t_x.(z_prime,τ_prime;χ,Pst=sim_eq["Pst"])
se_1980_z     = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./sim_eq["yt"][z_prime]
sg_1980_z     = sg_t_x.(t_prime,τ_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_1980_z_1980 = (se_z.*sg_z        .- cal_v_z_1980.*se_1980_z.*sg_1980_z       ).*g_pg .+ 
                  (se_z.*(1 .- sg_z) .- cal_v_z_1980.*se_1980_z.*(1 .- sg_1980_z)).*g_ps

FS              = [0;cumsum( g_D .+ g_n )]
FS_2023_z_2023  = [0;cumsum( g_D .+ g_n .+ dev_2023_z_2023[2:end] )]
FS_2010_z_2010  = [0;cumsum( g_D .+ g_n .+ dev_2010_z_2010[2:end] )]
FS_2000_z_2000  = [0;cumsum( g_D .+ g_n .+ dev_2000_z_2000[2:end] )]
FS_1990_z_1990  = [0;cumsum( g_D .+ g_n .+ dev_1990_z_1990[2:end] )]
FS_1980_z_1980  = [0;cumsum( g_D .+ g_n .+ dev_1980_z_1980[2:end] )]

p = plot(1980:2023, FS_2023_z_2023, label="FS Index (base 2023)",
    ylabel="Cummulative Growth",linestyle=:dash, lw=2.0,
    xticks=1980:5:2023, yticks=0.00:0.20:1.40,
    ylims = (0.00,1.40),
    minorgridalpha=0.2, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.15, 0.90),legendfont=legendfont,
    xrotation=45,framestyle=:box)

plot!(p, 1980:2023, FS_1980_z_1980 , label="FS Index (base 1980)",linestyle=:dot, lw=2, color=:black)

plot!(p, 1980:2023, FS, label="Chained Index",
    linestyle=:solid, lw=2, minorgrid=true, minorgridalpha=0.2, color=:black)

xaxis!(p, minor_ticks=true, minor_tick_step=1.00)
yaxis!(p, minor_ticks=true, minor_tick_step=0.01)

#savefig(joinpath(figuresdir, "FS_BBEV.png"))

#Differences Between Fixed-Base and Chained Indices
println("""
$(repeat("=", 60))
Differences Between Fixed-Base and Chained Indices (2023)
$(repeat("=", 60))
FS Index (base 1980) - Chained Index: $(round(FS_1980_z_1980[end] - FS[end], digits=3))
FS Index (base 2023) - Chained Index: $(round(FS_2023_z_2023[end] - FS[end], digits=3))
$(repeat("=", 60))
""")
#---------------------------------------------------------------

#---------------------------------------------------------
### Figure 5 (b)
#Alternative Fisher-Shell Indices
###Price Chained Fisher-Shell Index
sg_z = sim_eq["sgt"]
ss_z = sim_eq["sst"]
se_z = sc[1:end]

#At Preferences of 1980
t_base             = 1980
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_1980_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])
se_1980_z          = e_1980_z./(sim_eq["yt"][ z_prime ])
sg_1980_z          = sg_t_x.(t_prime,z_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

#At Preferences of 2023
t_base             = 2023
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_2023_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])
se_2023_z          = e_2023_z./(sim_eq["yt"][ z_prime ])
sg_2023_z          = sg_t_x.(t_prime,z_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_2023_z = ( se_z.*sg_z .-  se_2023_z.*sg_2023_z ).*g_pg  .+ ( se_z.*(1 .- sg_z) .-  se_2023_z.*(1 .- sg_2023_z) ).*g_ps
dev_1980_z = ( se_z.*sg_z .-  se_1980_z.*sg_1980_z ).*g_pg  .+ ( se_z.*(1 .- sg_z) .-  se_1980_z.*(1 .- sg_1980_z) ).*g_ps

FS          = [0;cumsum( g_D .+ g_n )]
FS_2023_z   = [0;cumsum( g_D .+ g_n .+ dev_2023_z[1:end-1] )]
FS_1980_z   = [0;cumsum( g_D .+ g_n .+ dev_1980_z[1:end-1] )]

plot(1980:2023, FS_2023_z ,
    linestyle=:dash, lw=2.0, color=:black,
    xticks = 1980:5:2025,  
    label = "FS Index (2023 Preferences)")

plot!(1980:2023, FS_1980_z ,
    linestyle=:dot, lw=2.0, color=:black,
    label = "FS Index (1980 Preferences)",
    legend=(0.15, 0.92))

plot!(1980:2023, FS ,ylabel="Cummulative Growth", 
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.5, color=:black,
    xticks=1980:5:2024, yticks=0.0:0.2:1.4,
    ylim=(0.00, 1.40),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    label = "Chained Index",
    xrotation=45,
    framestyle=:box)

#savefig(joinpath(figuresdir,"price_chained_FS.png"))
#---------------------------------------------------------

#---------------------------------------------------------
#Alternative Fisher-Shell Indices
####Preferences Chained Fisher-Shell Index

##At Prices of 1980
t_base          = 1980
t_prime         = t_base - 1980 + 1
z_prime         = (1980:1:2023) .- 1980 .+ 1
cal_v_1980_z    = cal_v_t_x.(t_prime,z_prime;χ,Pst=sim_eq["Pst"])

se_1980         = e_t_x.(t_prime,t_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./(sim_eq["yt"][ t_prime ])
sg_1980         = sg_t_x.(t_prime,t_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

se_1980_z       = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./(sim_eq["yt"][ z_prime ])
sg_z            = sg_t_x.(z_prime,z_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_z_1980      = (se_1980.*sg_1980 .- cal_v_1980_z.*se_1980_z.*sg_z).*g_pg .+ (se_1980.*(1- sg_1980) .- cal_v_1980_z.*se_1980_z.*(1 .- sg_z)).*g_ps    

#At Prices of 2023
t_base          = 2023
t_prime         = t_base - 1980 + 1   
z_prime         = (1980:1:2023) .- 1980 .+ 1

cal_v_2023_z    = cal_v_t_x.(t_prime,z_prime;χ,Pst=sim_eq["Pst"])

se_2023         = e_t_x.(t_prime,t_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./(sim_eq["yt"][ t_prime ])
sg_2023         = sg_t_x.(t_prime,t_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

se_2023_z       = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"])./(sim_eq["yt"][ z_prime ])
sg_z            = sg_t_x.(z_prime,t_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"],Pgt=sim_eq["Pgt"])

dev_z_2023      = (se_2023.*sg_2023 .- cal_v_2023_z.*se_2023_z.*sg_z).*g_pg .+ (se_2023.*(1- sg_2023) .- cal_v_2023_z.*se_2023_z.*(1 .- sg_z)).*g_ps    

FS              = [0;cumsum( g_D .+ g_n )]
FS_z_2023       = [0;cumsum( g_D .+ g_n .+ dev_z_2023[1:end-1] )]
FS_z_1980       = [0;cumsum( g_D .+ g_n .+ dev_z_1980[1:end-1] )]

plot(1980:2023, FS_z_2023 ,
    linestyle=:dash, lw=2.0,color=:black,
    label = "FS Index (2023 Prices)")
plot!(1980:2023, FS_z_1980 ,
    linestyle=:dot, lw=2.0, color=:black,
    label = "FS Index (1980 Prices)",
    legend=(0.15, 0.90))
    #title="Alternative Fisher-Shell Indices <br> (fixing prices)")

plot!(1980:2023, FS , ylabel="Cummulative Growth", 
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.5, color=:black,
    xticks=1980:5:2024, yticks=0.0:0.2:1.4,
    ylim=(0.00, 1.40),   
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legendfont=legendfont,
    label = "Chained Index",
    xrotation=45,
    framestyle=:box)

#savefig(joinpath(figuresdir,"prefs_chained_FS.png"))
#---------------------------------------------------------

#---------------------------------------------------------------
### Figure 6 (a)
#Alternative Expenditure Shares
#At Preferences of 1980
t_base             = 1980
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_1980_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"]) 
se_1980_z          = e_1980_z./(sim_eq["yt"][ z_prime ])

#At Preferences of 2023
t_base             = 2023
t_prime            = t_base - 1980 + 1
z_prime            = (1980:1:2023) .- 1980 .+ 1
e_2023_z           = e_t_x.(t_prime,z_prime;χ=χ,ν_t=ν_t,Pst=sim_eq["Pst"]) 
se_2023_z          = e_2023_z./(sim_eq["yt"][ z_prime ])

sc                 = sim_eq["et"]./sim_eq["yt"]
   
plot(1980:2023,sc, 
    label="se(t)",
    ylabel="Share of ConsumptionExpenditure <br> in Gross Income",
    linestyle=:solid, lw=2,
    xticks=1980:5:2025, yticks=0.0:0.50:5.00,
    minorgridalpha=0.5, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.75, 0.92),
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box) 

plot!(1980:2023,se_1980_z,label="se(1980,z)",
    linestyle=:dot, lw=2,color=:black)

plot!(1980:2023,se_2023_z,label="se(2023,z)",
    linestyle=:dash, lw=2,color=:black)

#savefig(joinpath(figuresdir,"se_t.png"))
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 6 (b)

#Goods Expenditure Share in ABGP vs a Fixed-Base Counterfactual
t_base     = 2023
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1
sg_2023_z  = sg_t_x.(t_prime,z_prime;χ=χ,η=η,γ=γ, Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"])  

t_base     = 1980
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1
sg_1980_z  = sg_t_x.(t_prime,z_prime;χ=χ,η=η,γ=γ,Pst=sim_eq["Pst"], Pgt=sim_eq["Pgt"])    

plot(1980:2023, sim_eq["sgt"], label="sg(t)",
    linestyle=:solid, lw=2,color=:black)

plot!(1980:2023, sg_1980_z, label="sg(1980,z)",
    linestyle=:dot, lw=2,color=:black)

plot!(1980:2023, sg_2023_z, label="sg(2023,z)",
    ylabel="Share of Goods in Consumption <br> Expenditure",
    linestyle=:dash, lw=2,
    xticks=1980:5:2025, yticks=0.000:0.05:0.300,
    ylims=(0.000, 0.300),
    minorgridalpha=0.5, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.75, 0.95),
    legendfont=legendfont,
    xrotation=45,
    left_margin=5Plots.mm,  
    framestyle=:box)
#savefig(joinpath(figuresdir, "sg_t_z.png"))
#---------------------------------------------------------------

#---------------------------------------------------------------
### Figure 7
g_FS              = g_D .+ g_n 
g_FS_2023_z_2023  = g_D .+ g_n .+ dev_2023_z_2023[2:end]
g_FS_2010_z_2010  = g_D .+ g_n .+ dev_2010_z_2010[2:end]
g_FS_2000_z_2000  = g_D .+ g_n .+ dev_2000_z_2000[2:end]
g_FS_1990_z_1990  = g_D .+ g_n .+ dev_1990_z_1990[2:end]
g_FS_1980_z_1980  = g_D .+ g_n .+ dev_1980_z_1980[2:end]

p = plot(1981:2023, g_FS_2023_z_2023, label="FS Index (base 2023)",
    ylabel="Growth Rate",
    linestyle=:dash, lw=2.5,
    xticks=1980:5:2025, yticks=0.00:0.005:0.05,
    minorgridalpha=0.2, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.620, 0.600),
    legendfont=legendfont,
    xrotation=45,
    framestyle=:box)

plot!(p, 1981:2023, g_FS_2010_z_2010, label="FS Index (base 2010)",
    linestyle=:dash, lw=2, color=:black)

plot!(p, 1981:2023, g_FS_2000_z_2000, label="FS Index (base 2000)",
    linestyle=:dash, lw=1.5, minorgrid=true, minorgridalpha=0.2, color=:black)

plot!(p, 1981:2023, g_FS_1990_z_1990, label="FS Index (base 1990)",
    linestyle=:dash, lw=1.0, minorgrid=true, minorgridalpha=0.2, color=:black)

plot!(p, 1981:2023, g_FS_1980_z_1980, label="FS Index (base 1980)",
    linestyle=:dash, lw=0.5, minorgrid=true, minorgridalpha=0.2, color=:black)

plot!(p, 1981:2023, g_FS , label="Chained Index",
    linestyle=:solid, lw=2, minorgrid=true, minorgridalpha=0.2, color=:black)

xaxis!(p, minor_ticks=true, minor_tick_step=1.00)
yaxis!(p, minor_ticks=true, minor_tick_step=0.01)
#savefig(joinpath(figuresdir, "FS_GrowthRates_1980_2023.png"))

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
#Decomposition of the Growth Rate Decline 
πg = sim_eq["Pgt"][2:end]./sim_eq["Pgt"][1:end-1] .- 1
πs = sim_eq["Pst"][2:end]./sim_eq["Pst"][1:end-1] .- 1
100*(s_e[2:end].*( (sim_eq["sgt"][end] .- sim_eq["sgt"][1]).*πg[1] .+ (sim_eq["sst"][end] .- sim_eq["sst"][1]).*πs[1] )) 

#Alternative Decomposition of the Growth Rate Decline
x_t         = sim_eq["Xt"]./sim_eq["Nt"]
s_e         = sim_eq["Et"]./sim_eq["Yt"]
sg_t        = sim_eq["sgt"]
ss_t        = sim_eq["sst"]

g_cg        = sim_eq["cgt"][2:end]./sim_eq["cgt"][1:end-1] .- 1
g_cs        = sim_eq["cst"][2:end]./sim_eq["cst"][1:end-1] .- 1
g_x         = x_t[2:end]./x_t[1:end-1] .- 1

g_D         = s_e[2:end].*( sim_eq["sgt"][2:end].*g_cg .+ sim_eq["sst"][2:end].*g_cs  ) .+ (1 .- s_e[2:end]).*g_x
g_D_agg     = round.(g_D .+ g_n, digits=5)
g_D_agg     = (g_D .+ g_n).*100

d_gD_add    = g_D_agg[2:end] .- g_D_agg[1:end-1]
d_sg_t      = sim_eq["sgt"][2:end] .- sim_eq["sgt"][1:end-1] 
ST_Effect   = (s_e[3:end].*(d_sg_t[2:end].*(g_cg[2:end] .- g_cs[2:end]))).*100
d_g_cs      = g_cs[2:end] .- g_cs[1:end-1] 
SGS_Effect  = (s_e[3:end].*(1 .- sim_eq["sgt"][3:end]).*d_g_cs).*100

final_SGS_Effect = (g_D_agg[1] .+ cumsum([0;SGS_Effect]))[end] .- (g_D_agg[1] .+ cumsum([0;SGS_Effect]))[1]
final_ST_Effect =  (g_D_agg[1] .+ cumsum([0;ST_Effect]))[end] .- (g_D_agg[1] .+ cumsum([0;ST_Effect]))[1]
#------------------------------------------------------------------------

#------------------------------------------------------------------------
round(final_SGS_Effect, digits=3) .+ round(final_ST_Effect, digits=3)

g_D_agg[end] - g_D_agg[1]


abs.(ST_Effect)./(abs.(ST_Effect) .+ abs.(SGS_Effect))

g_D_no_SGS       = s_e[2:end].*( sim_eq["sgt"][2:end].*g_cg .+ sim_eq["sst"][2:end].*g_cs[1]  ) .+ (1 .- s_e[2:end]).*g_x
g_D_agg_no_SGS    = (g_D_no_SGS .+ g_n).*100
round(g_D_agg_no_SGS[end] - g_D_agg_no_SGS[1],digits=3)

g_D_no_SCE       = s_e[2:end].*( sim_eq["sgt"][2].*g_cg .+ sim_eq["sst"][2].*g_cs  ) .+ (1 .- s_e[2]).*g_x
g_D_agg_no_SCE   = (g_D_no_SCE .+ g_n).*100
round(g_D_agg_no_SCE[end] - g_D_agg_no_SCE[1],digits=3)

g_D_agg .- g_D_agg_no_SCE 
SC_Effect  = round.( (g_D_agg .- g_D_agg_no_SCE), digits=5) 
SGS_Effect = round.((g_D_agg .- g_D_agg_no_SGS),digits=3)

SGS_Effect .+ SC_Effect

g_D_agg[end] - g_D_agg[1]
#------------------------------------------------------------------------

#------------------------------------------------------------------------
t_base     = 2023
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1

e_tilde_2023_z = e_tilde_tz.(t_prime,z_prime;χ=χ,η=η,γ=γ,e=sim_eq["et"],Pgt=sim_eq["Pgt"],Pst=sim_eq["Pst"])

#Alternative Decomposition of the Growth Rate Decline
gD_e         = sg_t[2:end].*( g_cg )  .+ (1 .- sg_t[2:end]).*( g_cs ) 
De           = [0;cumsum( gD_e  .+ g_n )]

[0;gD_e].*((sim_eq["et"]./e_tilde_2023_z).^(χ))
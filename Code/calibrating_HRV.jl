currentdir = @__DIR__
datadir    = abspath(joinpath(currentdir, "..", "Data"))
figuresdir = abspath(joinpath(currentdir, "..", "Figures"))

using Pkg
Pkg.activate(@__DIR__)  
Pkg.instantiate()
#Pkg.resolve()
#Pkg.precompile()

using XLSX, DataFrames
using BlackBoxOptim
using Plots; plotlyjs()

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Import the HRV Data
file_name  = "HRV2021_Data.xlsx"
sheet_name = "Data"
data_range = "A1:AU78"

# Read the data from the Excel file
HRV_data = XLSX.readtable(joinpath(datadir, file_name), sheet_name)
#HRV_df = DataFrame(HRV_data)
HRV_df = DataFrame(HRV_data...)

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
Ag_2017        = HRV_df[HRV_df.year .== 2017, "TFPVA_GOOD_I"][1]
As_1980        = HRV_df[HRV_df.year .== 1980, "TFPVA_SERV_I"][1]
As_2017        = HRV_df[HRV_df.year .== 2017, "TFPVA_SERV_I"][1]
calAx_1980     = HRV_df[HRV_df.year .== 1980, "calA_X_I_TD"][1]
calAx_2017     = HRV_df[HRV_df.year .== 2017, "calA_X_I_TD"][1]

#Get the Initial and Final Values for Population and Efficiency Units of Labor
N_1980         = HRV_df[HRV_df.year .== 1980, "POP"][1]
N_2023         = HRV_df[HRV_df.year .== 2023, "POP"][1]

L_1980         = HRV_df[HRV_df.year .== 1980,"LAB_TOT_QI"][1]
L_2017         = HRV_df[HRV_df.year .== 2017,"LAB_TOT_QI"][1]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Define the Functions Characterizing the Equilibrium Conditions
# Recall that L = N*h, where N is the population and h is the avg efficiency units of labor per person

## Relative Prices
Ps(calAx,As) = calAx/As
Pg(calAx,Ag) = calAx/Ag

## Investment Productivity
calAx(Ax,Ag,As,ωx,εx) = Ax*(ωx * (Ag^(εx-1)) + (1-ωx) * (As^(εx-1)) )^(1/(εx-1))

# Aggregate Variables
## Output
Y(calAx,K,L,θ)        = calAx*(K^θ)*(L^(1-θ))

## Capital accumulation
dKdt(calAx,K,L,E,δ,θ) = calAx*(K^θ)*(L^(1-θ)) - E - δ*K

## Consumption Expentiture
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
#Calibrate the growth rates of Agt,Ast,Axt (these values are from the HRV TFP series)
g_Ag        = log(Ag_2017/Ag_1980)/(2017-1980)
g_As        = log(As_2017/As_1980)/(2017-1980)
g_calAx     = log(calAx_2017/calAx_1980)/(2017-1980)
g_n         = (log(N_2023/N_1980))/(2023-1980)
g_l         = (log(L_2017/L_1980))/(2017-1980)

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
# Compute Equilibrium Prices (relative to the investment numeraire)
Ps_t        = Ps.(calAx_t,As_t)
Pg_t        = Pg.(calAx_t,Ag_t)
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
# Calibrate the Parameters of the Model Using a SMM
## Here we want the model to be in ABGP from 1980 onwards
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
                  "Agt", "Ast", "Axt", "calAxt", "calAxhatt"]
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

  while abs(err_ge) .> 1e-16 && iter < 100
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

  #Find the e0: here we use 
  while abs(err_gk) .> 1e-15 && iter < 100
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

  vars = ["Yt" , "Xt" , "Et" , "Kt", 
          "yt",  "et" , "kt" , "ht", "sst", "sgt", 
          "cst", "cgt", "Cst", "Cgt","Xst", "Xgt", "Pst", "Pgt", "Kst", "Kgt", 
          "Lt" , "Nt" , "Ht" ,
          "Lst", "Lgt", "Yst", "Ygt", "Wt", "Rt", 
          "Agt", "Ast", "Axt", "calAxt", "calAxhatt"]

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
  push!(eq["calAxhatt"],calAx_hat_t[1])
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
    push!(eq["calAxhatt"],calAx_hat_t[t+1])

    push!(eq["Agt"], Ag_t[t+1])
    push!(eq["Ast"], As_t[t+1])

    push!(eq["ht"], h_t[t+1])
    push!(eq["Lt"], L_t[t+1])
    push!(eq["Nt"], N_t[t+1])

    g_e  = dedt_e(eq["calAxhatt"][t],eq["kt"][t],θ,ρ,δ,χ,g_ps)  
    e    = eq["et"][t]*exp(g_e)
    push!(eq["et"],e)

    g_k  = dkdt(eq["calAxhatt"][t],eq["kt"][t],eq["et"][t],δ,n,θ)/eq["kt"][t]
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
sim_model(params) = sim_model_full(params;Pg_t=Pg_t,Ps_t=Ps_t,θ=θ,ρ=ρ,δ=δ,g_calAx=g_calAx,g_As=g_As,g_Ag=g_Ag,g_h=g_h,g_l=g_l,g_n=g_n)

function SSE(params;cons_exp_share_data)
    sim_eq = sim_model(params)
    return 100*sum( (sim_eq["sgt"] .- cons_exp_share_data).^2 )
end

#HRV cali: χ=0.550,η=0.440,γ=0.690
ParamSpace  = [(0.550,0.550),  #χ  
               (0.000,1.000),  #η
               (0.000,1.000)]  #γ


opt_problem    =  bbsetup(params -> SSE(params;cons_exp_share_data=VAC_GOOD_SHARE);
                                    SearchRange=ParamSpace,TraceMode=:verbose, TargetFitness=0.001,
                                    Method=:adaptive_de_rand_1_bin_radiuslimited)
                                
opt_results    = bboptimize(opt_problem, MaxFuncEvals= 20_000)

χ,η,γ          = best_candidate(opt_results) 

round(γ,digits=3), round(η,digits=3), round(χ,digits=3)

# Check that the Paramter Restrictions are Satisfied (1 > γ > χ > 0)
1 > γ > χ > 0

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

#Check: Do we get the same results as in the Analytical Solution?
g_abgp      = g_calAx/(1-θ) + g_h 
g_k         = g_abgp
g_ps        = g_calAx - g_As 
κ           = (ρ + δ + χ*g_ps + (1-χ)*g_k)/θ 

Ag_0        = 1.000
As_0        = 1.000
calAx_0     = 1.000

k0         = κ^(1/(θ-1))*(calAx_0^(1/(θ-1)))
sim_eq["kt"][1]

e0         = (κ - δ - g_n - g_k)*k0
sim_eq["et"][1]

y0         = κ*k0
sim_eq["yt"][1]

cg_t_an = η*(( e0/Ps_t[1] )^(1-χ))*(( Pg_t[1]/Ps_t[1])^γ ) 
sim_eq["cgt"][1]

cs_t_an = e0/Ps_t[1] - η*(( e0/Ps_t[1] )^(1-χ))*(( Pg_t[1]/Ps_t[1])^γ ) 
sim_eq["cst"][1]
#---------------------------------------------------------------

#-------------------------------------------------------
#Plot the Resulting Series
year = 1980:2023
plot(year, sim_eq["sgt"], label="Model", xlabel="Year", linestyle=:solid, color=:black, lw=2.00, 
    ylabel="Share of Goods in Consumption Expenditure", 
    title="γ=$(round(γ, digits=4)), η=$(round(η, digits=4)), χ=$(round(χ, digits=4))")
plot!(year, VAC_GOOD_SHARE, label="Data", legend=(0.800, 0.850), linestyle=:dot, lw=2.00, 
    minorgrid=true, minorgridalpha=0.9, color=:black)
savefig(joinpath(figuresdir, "Model_Fit.png"))

#plot(year, eq["egt"], label="Model", xlabel="Year", ylabel="Share of Goods in Consumption Expenditure", title="γ=0.990,η=0.272,χ=0.990")
#plot!(year, VAC_GOOD_SHARE, label="Data")
#savefig("HRV_chart_cali_alt.png")
#----------------------------------------------------------------------------------

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

#---------------------------------------------------------------
#Compute the Divisia Index
x_t   = sim_eq["Xt"]./sim_eq["Nt"]
sc    = sim_eq["Et"]./sim_eq["Yt"]
sg_t  = sim_eq["sgt"]
ss_t  = sim_eq["sst"]
g_cg  = log.(sim_eq["cgt"][2:end]) .- log.(sim_eq["cgt"][1:end-1])
g_cs  = log.(sim_eq["cst"][2:end]) .- log.(sim_eq["cst"][1:end-1])
g_x   = log.(x_t[2:end])  .- log.(x_t[1:end-1])

g_D   = sc[2:end].*(sim_eq["sgt"][2:end].*g_cg .+ sim_eq["sst"][2:end].*g_cs  ) .+ (1 .- sc[2:end]).*g_x
#---------------------------------------------------------------

#---------------------------------------------------------------
#Compute the Marginal Value of Capital
ν_t  = (sim_eq["et"].^(χ-1))./(sim_eq["Pst"].^χ)

ν_t./ν_t[1]

#Compute Net Income Per Capita
m_t  = sim_eq["yt"] - δ*sim_eq["kt"]
#---------------------------------------------------------------

#---------------------------------------------------------------
#Compute the Equivalent Variation Measures
V(e,Pg,Ps;χ,η,γ)          = (1/χ)*( (e/Ps)^χ ) - (η/γ)*((Pg/Ps)^(γ)) - (1/χ) + (η/γ) 
indU(m,Pg,Ps;χ,η,γ,ν)     = V( (( ν *(Ps^χ) )^(1/(χ-1))),Pg,Ps; χ,η,γ ) + ν*( m - (( ν *(Ps^χ) )^(1/(χ-1))) )  

exp_func(m,Pg,Ps;χ,η,γ,ν) = ( (ν *(Ps^χ) )^(1/(χ-1)) ) + (m/ν) - (1/ν)*V( ((ν *(Ps^χ) )^(1/(χ-1))) ,Pg,Ps ; χ,η,γ )

Mhat(t,z;χ,η,γ,ν)         = exp_func(indU(m_t[z],sim_eq["Pgt"][z],sim_eq["Pst"][z];χ,η,γ,ν),sim_eq["Pgt"][t],sim_eq["Pst"][t];χ,η,γ,ν) 
#---------------------------------------------------------------

#---------------------------------------------------------------
#Compute the 2023 based EV
t_base     = 2023
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:2023) .- 1980 .+ 1

EV_2023_pc = (Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])./
             (Mhat.( t_prime ,       1 ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ 1 ])


exp_share_2023_base = sim_eq["et"][end]./(Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])

#Compute the 1980 based EV
t_base     = 1980
t_prime    = t_base - 1980 + 1
z_prime    = (1980:1:2023) .- 1980 .+ 1
EV_1980_pc = (Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])./(Mhat.( t_prime , 1 ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ 1 ])

exp_share_1980_base = sim_eq["et"][1]./(Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])

FS         = cumsum( g_D .+ g_n ) 
FS         = [0;FS]

log_EV_2023 = log.(EV_2023_pc) .+ g_n.*(z_prime .- 1 )
log_EV_1980 = log.(EV_1980_pc) .+ g_n.*(z_prime .- 1 )

#---------------------------------------------------------------
tickfont   = font(10)
guidefont  = font(12)
legendfont = font(10)

p = plot(1980:2023, log_EV_2023, label="2023-based Fisher-Shell Index",
    ylabel="Cummulative Welfare Growth",
    linestyle=:dashdot, lw=2,
    xticks=1980:5:2023, yticks=0.0:0.2:1.4,
    minorgridalpha=0.2, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.15, 0.95),
    legendfont=legendfont
)

plot!(p, 1980:2023, log_EV_1980, label="1980-based Fisher-Shell Index",
    linestyle=:dash, lw=2, color=:black)

plot!(p, 1980:2023, FS, label="Chained Fisher-Shell Index",
    linestyle=:solid, lw=2, minorgrid=true, minorgridalpha=0.2, color=:black)

savefig(joinpath(figuresdir, "FS_BBEV.png"))
#---------------------------------------------------------------

#---------------------------------------------------------------
tickfont   = font(10)
guidefont  = font(12)
legendfont = font(10)

plot(1980:2023,exp_share_1980_base, 
    label="1980-based Fisher-Shell Index",
    ylabel="Expenditure Share of Income",
    linestyle=:solid, lw=2,
    xticks=1980:5:2025, yticks=0.0:0.05:0.80,
    minorgridalpha=0.2, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.50, 0.95),
    legendfont=legendfont)
savefig(joinpath(figuresdir,"Exp_Share_1980.png"))

plot(1980:2023,exp_share_2023_base, 
    label="2023-based Fisher-Shell Index",
    ylabel="Expenditure Share of Income",
    linestyle=:solid, lw=2,
    xticks=1980:5:2025, yticks=0.0:0.05:0.90,
    minorgridalpha=0.2, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.58, 0.95),
    legendfont=legendfont)
savefig(joinpath(figuresdir,"Exp_Share_2017.png"))
#---------------------------------------------------------------


#---------------------------------------------------------------
#Now Compute Different Base EV Measures
#Compute the 2023-based EV
t_base     = 2023
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:1990) .- 1980 .+ 1
EV_2023_pc = (Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])./
             (Mhat.( t_prime ,       1 ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ 1 ])

log_EV_2023 = log.(EV_2023_pc) .+ g_n.*(z_prime .- 1 )


#Compute the 2010-based EV
t_base     = 2010
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:1990) .- 1980 .+ 1
EV_2010_pc = (Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])./
             (Mhat.( t_prime ,       1 ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ 1 ])

log_EV_2010 = log.(EV_2010_pc) .+ g_n.*(z_prime .- 1 )

#Compute the 2000-based EV
t_base     = 2000
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:1990) .- 1980 .+ 1
EV_2000_pc = (Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])./
             (Mhat.( t_prime ,       1 ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ 1 ])

log_EV_2000 = log.(EV_2000_pc) .+ g_n.*(z_prime .- 1 )

#Compute the 1990-based EV
t_base     = 1990
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:1990) .- 1980 .+ 1
EV_1990_pc = (Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])./
             (Mhat.( t_prime ,       1 ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ 1 ])

log_EV_1990 = log.(EV_1990_pc) .+ g_n.*(z_prime .- 1 )

#Compute the 1980-based EV
t_base     = 1980
t_prime    = t_base - 1980 + 1   
z_prime    = (1980:1:1990) .- 1980 .+ 1
EV_1990_pc = (Mhat.( t_prime , z_prime ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ z_prime ])./
             (Mhat.( t_prime ,       1 ;χ=χ,η=η,γ=γ,ν=ν_t[t_prime]) .+ δ*sim_eq["kt"][ 1 ])

log_EV_1980 = log.(EV_1990_pc) .+ g_n.*(z_prime .- 1 )
#---------------------------------------------------------------

#---------------------------------------------------------------
tickfont   = font(10)
guidefont  = font(12)
legendfont = font(10)

p = plot(1980:1990, log_EV_2023, label="2023-based Index",
    ylabel="Accumulated Welfare Growth",
    linestyle=:dash, lw=2.5,
    xticks=1980:2:1992, yticks=0.00:0.05:0.30,
    minorgridalpha=0.2, color=:black,
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.15, 0.95),
    legendfont=legendfont
)

plot!(p, 1980:1990, log_EV_2010, label="2010-based Index",
    linestyle=:dash, lw=2, color=:black)

plot!(p, 1980:1990, log_EV_2000, label="2000-based Index",
    linestyle=:dash, lw=1.5, minorgrid=true, minorgridalpha=0.2, color=:black)

plot!(p, 1980:1990, log_EV_1990, label="1990-based Index",
    linestyle=:dash, lw=1.0, minorgrid=true, minorgridalpha=0.2, color=:black)

plot!(p, 1980:1990, log_EV_1980, label="1980-based Index",
    linestyle=:dash, lw=0.5, minorgrid=true, minorgridalpha=0.2, color=:black)

plot!(p, 1980:1990, FS[1:11], label="Chained Index",
    linestyle=:solid, lw=2, minorgrid=true, minorgridalpha=0.2, color=:black)

xaxis!(p, minor_ticks=true, minor_tick_step=1.00)
yaxis!(p, minor_ticks=true, minor_tick_step=0.01)
savefig(joinpath(figuresdir,"Different_Base_Indices.png"))

#---------------------------------------------------------------

#---------------------------------------------------------------
#Now Compute the Laspayres and Paasche GDP Indices
tickfont   = font(10)
guidefont  = font(12)
legendfont = font(10)

t_base      = 2023
T           = t_base - 1980 + 1   #-> 2023

t_base_2    = 1980
τ           = t_base_2 - 1980 + 1  #-> 1980 

e_Paasche = ( sim_eq["Pst"][T].*sim_eq["cst"] .+ sim_eq["Pgt"][T].*sim_eq["cgt"])./
( sim_eq["Pst"][T].*sim_eq["cst"][τ] .+ sim_eq["Pgt"][T].*sim_eq["cgt"][τ])

P_Paasche = sim_eq["et"]./e_Paasche

e_Laspeyres = ( sim_eq["Pst"][τ].*sim_eq["cst"] .+ sim_eq["Pgt"][τ].*sim_eq["cgt"])./( sim_eq["Pst"][τ].*sim_eq["cst"][τ] .+ sim_eq["Pgt"][τ].*sim_eq["cgt"][τ])
P_Laspeyres = sim_eq["et"]./e_Laspeyres


y_Laspeyres = (P_Paasche[τ].*e_Laspeyres .+ (sim_eq["yt"] - sim_eq["et"]) )./(P_Paasche[τ].*e_Laspeyres[τ] .+ (sim_eq["yt"][τ] - sim_eq["et"][τ]) )

y_Paasche   = (P_Laspeyres[T].*e_Laspeyres .+ (sim_eq["yt"] - sim_eq["et"]) )./(P_Laspeyres[T].*e_Laspeyres[τ] .+ (sim_eq["yt"][τ] - sim_eq["et"][τ]) )

plot(1980:2023, log.(y_Laspeyres) .+ g_n.*((1980 .- 1980 .+1):(2023 .- 1980 .+1)), 
    label="Laspeyres", ylabel="GDP Index (Log scale;1980=0)",    
    linestyle=:dash, lw=2.0,
    minorgridalpha=0.2, color=:black,
    xticks=1980:5:2025,  
    yticks=0.00:0.20:1.00,  
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=(0.15, 0.95),
    legendfont=legendfont)

plot!(1980:2023, log.(y_Paasche) .+ g_n.*(1980 .- 1980 .+1 : 2023 .- 1980 .+1 ), label="Paasche",
        linestyle=:dot, color=:black, lw=2.00)

plot!(1980:2023, 0.5.*log.(y_Paasche) .+ 0.5.*log.(y_Laspeyres) .+ g_n.*(1980 .- 1980 .+1 : 2023 .- 1980 .+1 ), label="Fisher-Ideal",
        linestyle=:solid, color=:black, lw=2.00)

savefig(joinpath(figuresdir,"GDP_sc.png"))
#---------------------------------------------------------------


#---------------------------------------------------------------
#Plot the Goods Consumption Value Added Share
tickfont   = font(10)
guidefont  = font(12)
legendfont = font(10)

plot(1980:2023, VAC_GOOD_SHARE,
    ylabel="% of Total Consumption Value Added",    
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.2, color=:black,
    xticks=1980:2:2025,  
    yticks=0.00:0.025:0.30,  
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    xrotation=45)

savefig(joinpath(figuresdir,"goods_consumption_share.png"))

#---------------------------------------------------------------
#Plot Investment TFP
tickfont   = font(10)
guidefont  = font(12)
#legendfont = font(10)

calA_X_I   = HRV_df[!,"calA_X_I_TD"]

plot(1947:2023, calA_X_I ,
    ylabel="TFP Index (1947=1)", 
    linestyle=:solid, lw=2.0,
    minorgridalpha=0.2, color=:black,
    xticks=1947:3:2025,    
    xtickfont=tickfont, ytickfont=tickfont,
    xguidefont=guidefont, yguidefont=guidefont,
    legend=:none,
    xrotation=45)

savefig(joinpath(figuresdir,"investment_TFP.png"))


stop here
#---------------------------------------------------------------
#=
#---------------------------------------------------------------
# Export the data to an Excel file

# Create a dictionary with the data to export
time_index = 1980:2017

data_to_export = Dict(
    "Year"    => time_index,
    "Agt"     => eq["Agt"],
    "Ast"     => eq["Ast"],
    "calAxt"  => eq["calAxt"],
    "Yt"      => eq["Yt"],
    "Et"      => eq["Et"],
    "Xt"      => eq["Xt"],
    "Xst"     => eq["Xst"],
    "Xgt"     => eq["Xgt"],
    "Kt"      => eq["Kt"],
    "Kst"     => eq["Kst"],
    "Kgt"     => eq["Kgt"],
    "Cgt"     => eq["Cgt"],    
    "Cst"     => eq["Cst"],    
    "Pst"     => eq["Pst"],
    "Pgt"     => eq["Pgt"],
    "egt"     => eq["egt"],
    "est"     => eq["est"],
    "Lgt"     => eq["Lgt"],
    "Lst"     => eq["Lst"],
    "Wt"      => eq["Wt"],
    "Rt"      => eq["Rt"],
    "nut"     => nu,
    "Divisia" => g_D,
    #"Value of the Problem" => V_K0,
    "σ_t"    => σ_t 
)

# Create a new Excel file and write the data
XLSX.openxlsx("HRV_sol_HRV_chi_Boppart_PopAadj.xlsx", mode="w") do xf
    sheet = xf[1]  # Get the first sheet
    row   = 1
    for (key, values) in data_to_export
        sheet[row, 1] = key  # Write the variable name
        for (i, value) in enumerate(values)
            sheet[row, i + 1] = value  # Write the variable values
        end
        row += 1
    end
end


using XLSX

# Create a new Excel file and write the data
XLSX.openxlsx("HRV_sol_HRV_chi_Boppart_PopAadj.xlsx", mode="w") do xf
    sheet = xf[1] 
    XLSX.rename!(sheet, "Variables")

    # Preserve the order of variables by iterating over the keys in the desired order
    ordered_keys = ["Year", "Agt", "Ast", "calAxt", "Yt", "Et", "Xt", "Xst", "Xgt", "Kt","Kst", "Kgt", 
                    "Cgt", "Cst", "Pst", "Pgt", "egt", "est", "Lgt", "Lst", "Wt", "Rt", "nut", "Divisia", "σ_t"]


    # Validate that all vectors are non-empty and have the same length
    num_rows = length(data_to_export["Year"])  # Use "Year" as the reference length
    for key in ordered_keys
        if !haskey(data_to_export, key)
            error("Key '$key' is missing in the data_to_export dictionary.")
        end
        if length(data_to_export[key]) != num_rows
            error("Key '$key' has a different length or is empty.")
        end
    end

    # Write the header row
    for (col, key) in enumerate(ordered_keys)
        sheet[1, col] = key  # Write the variable name in the first row
    end

    # Write the data column-wise
    for (col, key) in enumerate(ordered_keys)
        values = data_to_export[key]
        for row in 1:num_rows
            sheet[row + 1, col] = values[row]  # Write the data starting from the second row
        end
    end
    # Second sheet: Export parameter values
    sheet2 = XLSX.addsheet!(xf,"Parameters")

    # Define the parameters and their values
    parameters = Dict(
        "ωx" => ωx,
        "εx" => εx,
        "εc" => εc,
        "θ"  => θ,
        "ρ"  => ρ,
        "δ"  => δ,
        "χ"  => χ,
        "η"  => η,
        "γ"  => γ
    )

    # Write the parameters to the second sheet
    row = 1
    for (key, value) in parameters
        sheet2[row, 1] = key    
        sheet2[row, 2] = value  
        row += 1
    end
end

=#

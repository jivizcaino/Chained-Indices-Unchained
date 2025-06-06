cd("C:\\Users\\lezjv\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\Code")
#cd("C:\\Users\\Nacho\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\Code")
pwd()

using XLSX, DataFrames
using BlackBoxOptim
using Plots
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Import the HRV Data
file_path  = "C:\\Users\\lezjv\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\HRV2021_Data.xlsx"
#file_path  = "C:\\Users\\Nacho\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\HRV2021_Data.xlsx"
sheet_name = "Data"
data_range = "A1:CU72"

# Read the data from the Excel file
HRV_data = XLSX.readtable(file_path, sheet_name)
HRV_df = DataFrame(HRV_data)
#HRV_df = DataFrame(HRV_data...)

#Get the shares of Goods and Services in Investment and Consumption (X,C) respectively
VAX_GOOD_SHARE = HRV_df[HRV_df.year .>= 1980, "VAX_GOOD_S"]
VAX_SERV_SHARE = HRV_df[HRV_df.year .>= 1980, "VAX_SERV_S"]

VAC_GOOD_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_GOOD_S"]
VAC_SERV_SHARE = HRV_df[HRV_df.year .>= 1980, "VAC_SERV_S"]

#Get the Initial and Final Values for the TFP Indices
Ag_1980        = HRV_df[HRV_df.year .== 1980, "TFPVA_GOOD_I"][1]
Ag_2017        = HRV_df[HRV_df.year .== 2017, "TFPVA_GOOD_I"][1]
As_1980        = HRV_df[HRV_df.year .== 1980, "TFPVA_SERV_I"][1]
As_2017        = HRV_df[HRV_df.year .== 2017, "TFPVA_SERV_I"][1]

calAx_1980     = HRV_df[HRV_df.year .== 1980, "calA_X_I_BU"][1]
calAx_2017     = HRV_df[HRV_df.year .== 2017, "calA_X_I_BU"][1]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Define the Functions Characterizing the Equilibrium Conditions
## Relative Prices
P_st(calAxt,Ast) = calAxt/Ast
P_gt(calAxt,Agt) = calAxt/Agt

## Investment Productivity
calAx_t(Axt, Agt, Ast, ωx, εx) = Axt * (ωx * (Agt^(εx-1)) + (1-ωx) * (Ast^(εx-1)) )^(1/(εx-1))

## Aggregate Output 
Y_t(calAx_t, K_t, θ)           = calAx_t*K_t^θ

## Capital accumulation
dK_dt(calAx_t, Kt, Et, δ, θ)   = calAx_t*Kt^θ - Et - δ*Kt

## Euler Equation for Consumption Expentiture
dE_E(calAx_t,Kt,θ,ρ,δ,χ,dPst_Pst) = (1 / (1-χ)) * (θ*calAx_t*Kt^(θ-1) - ρ - δ - χ * dPst_Pst)

## Goods Share of Consumption Expenditure
Egt_share(Pgt, Et, Pst,ν, χ, γ) = ν*((Et / Pst)^(-χ)) * ((Pgt / Pst)^γ)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Calibrate the growth rates of Agt,Ast,Axt (these values are from the HRV TFP series)
g_Ag        = log(Ag_2017/Ag_1980)/(2017-1980)
g_As        = log(As_2017/As_1980)/(2017-1980)
g_calAx     = log(calAx_2017/calAx_1980)/(2017-1980)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Build the series for Agt, Ast, Axt (Normalized to 1 in 1980)
Ag         = Ag_1980.*exp.( g_Ag.*((1980:2017) .- 1980))
As         = As_1980.*exp.( g_As.*((1980:2017) .- 1980))
calAx      = calAx_1980.*exp.( g_calAx.*((1980:2017) .- 1980))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Compute Equilibrium Prices (relative to the investment numeraire)
Ps         = P_st.(calAx,As)
Pg         = P_gt.(calAx,Ag)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
ωx  = 0.650
ωc  = 0.290

εx  = 0.000
εc  = 0.000

θ   = 1/3
ρ   = 0.040
δ   = 0.080

ν   = 0.440
χ   = 0.550
γ   = 0.690

#χ = 0.42074919861636473
#ν = 0.1821488971899642
#γ = 1.4516543816138152 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Define a Dictionary to Store the Values of the Variables in Equilibrium
vars = ["Yt" , "Xt", "Et", "Kt", "est", "egt", "Cst", "Cgt", "Xst", "Xgt", "Pst", "Pgt", "Kst", "Kgt", 
        "Lst", "Lgt", "Wt", "Rt", "Agt", "Ast", "Axt", "calAxt"]

eq = Dict{String, Vector{Any}}()
for var in vars
    eq[var] = Vector{Any}()
end

#Compute the growth rate of the Aggregate Variables in ABGP
g_abgp = g_calAx/(1-θ) 


calAx0 = calAx_1980
#calAx0 = 1.000

#Compute the growth rate of the Aggregate Variables in ABGP
g_abgp     = g_calAx/(1-θ) 

#K0  = ((((1-χ)/(1-θ))*g_calAx + ρ + δ + χ*(g_calAx-g_As))/(θ*calAx0))^(1/(θ-1))
#E0  = calAx0*K0^θ - ( δ + g_abgp )*K0

#------------------------------------------------------------------------------
#Compute the K0,E0 that puts the Dynamic System in ABGP from t0
K0_max     = 100.0 
K0_min     = 0.001 

K0_mid = (K0_max + K0_min)/2

err_ge = Inf
iter   = 0 

while abs(err_ge) .> 1e-16 && iter < 100
    K0_mid = (K0_max + K0_min)/2
    err_ge = dE_E(calAx[1],K0_mid,θ,ρ,δ,χ,(Ps[2]-Ps[1])/Ps[1]) - g_abgp

    if err_ge > 0
        K0_min = K0_mid
    else
        K0_max = K0_mid
    end
    K0_mid = (K0_max + K0_min)/2
    iter += 1
end
K0 = K0_mid

E0_max = K0*10 
E0_min = K0/10 

err_gk = Inf
iter   = 0 

E0_mid = (E0_max + E0_min)/2

while abs(err_gk) .> 1e-15 && iter < 100

    err_gk = dK_dt(calAx[1],K0,E0_mid, δ, θ)/K0  - g_abgp

    if err_gk > 0
        E0_min = E0_mid
    else
        E0_max = E0_mid
    end

    E0_mid = (E0_max + E0_min)/2
    iter += 1
end

E0 = E0_mid

push!(eq["Et"] , E0)
push!(eq["Kt"] , K0)

push!(eq["Agt"], Ag[1])
push!(eq["Ast"], As[1])
push!(eq["calAxt"], calAx[1])
push!(eq["Pst"], Ps[1])
push!(eq["Pgt"], Pg[1])

t=0
for t in 1:(2018-1980) 
    eg_t = Egt_share(eq["Pgt"][t],eq["Et"][t],eq["Pst"][t],ν, χ, γ)
    es_t = 1 - eg_t
    push!(eq["est"],es_t)
    push!(eq["egt"],eg_t)

    Y    = Y_t(eq["calAxt"][t],eq["Kt"][t], θ)
    X    = Y - eq["Et"][t]
    push!(eq["Yt"],Y)
    push!(eq["Xt"],X)

    Xs   = eq["Xt"][t]/((1 + (ωx/(1-ωx))*(eq["Agt"][t]/eq["Ast"][t])^(εx)))
    push!(eq["Xst"],Xs)
    Xg   = X - Xs
    push!(eq["Xgt"],Xg)

    Cs = eq["est"][t]*(eq["Et"][t]/eq["Pst"][t])
    Cg = eq["egt"][t]*(eq["Et"][t]/eq["Pgt"][t])
    push!(eq["Cst"],Cs)
    push!(eq["Cgt"],Cg)

    push!(eq["Rt"],θ*eq["calAxt"][t]*eq["Kt"][t]^(θ-1))
    push!(eq["Wt"],(1-θ)*eq["calAxt"][t]*eq["Kt"][t]^θ)

    Ls = (eq["Xt"][t]/eq["Yt"][t])*(1/(((eq["Pgt"][t]*eq["Xgt"][t])/(eq["Pst"][t]*eq["Xst"][t]))+1))+(eq["Et"][t]/eq["Yt"][t])*(1/(((eq["Pgt"][t]*eq["Cgt"][t])  /(eq["Pst"][t]*eq["Cst"][t]))    +1))   
    Lg = 1-Ls

    push!(eq["Lst"],Ls)
    push!(eq["Lgt"],Lg)

    Ks = eq["Kt"][t]*Ls
    Kg = eq["Kt"][t]*Lg

    # Break loop here
    if t == 38
        break
    end

    push!(eq["Pst"], Ps[t+1])
    push!(eq["Pgt"], Pg[t+1])

    push!(eq["calAxt"],calAx[t+1])
    push!(eq["Agt"], Ag[t+1])
    push!(eq["Ast"], As[t+1])
    #push!(eq["Axt"], Ax[t+1])

    g_E  = dE_E(eq["calAxt"][t],eq["Kt"][t],θ,ρ,δ,χ,  (eq["Pst"][t+1] - eq["Pst"][t])/eq["Pst"][t])  
    E    = eq["Et"][t]*exp(g_E)
    push!(eq["Et"],E)

    g_K  = dK_dt(eq["calAxt"][t],eq["Kt"][t],eq["Et"][t], δ, θ)/eq["Kt"][t]
    K    = eq["Kt"][t]*exp(g_K)
    push!(eq["Kt"],K)
end
#------------------------------------------------------------------------------

gY = log.(eq["Yt"][2:end]) .- log.(eq["Yt"][1:end-1])
gK = log.(eq["Kt"][2:end]) .- log.(eq["Kt"][1:end-1])
gE = log.(eq["Et"][2:end]) .- log.(eq["Et"][1:end-1])
gX = log.(eq["Xt"][2:end]) .- log.(eq["Xt"][1:end-1])
eq["Rt"]


year = 1980:2017
plot(year,eq["egt"],label="Model",xlabel="Year", ylabel="Share",tilte="Share of Goods in Consumption Expenditure")
plot!(year,VAC_GOOD_SHARE,label="Data")

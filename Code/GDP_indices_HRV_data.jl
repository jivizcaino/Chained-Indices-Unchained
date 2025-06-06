cd("C:\\Users\\lezjv\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\Code")
#cd("C:\\Users\\Nacho\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\Code")
pwd()

using XLSX, DataFrames
using Plots; plotlyjs()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Import the HRV Data
file_path  = "C:\\Users\\lezjv\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\HRV2021_Data.xlsx"
#file_path  = "C:\\Users\\Nacho\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\HRV2021_Data.xlsx"
sheet_name = "Data"
data_range = "A1:CU72"

# Read the data from the Excel file
HRV_data = XLSX.readtable(file_path, sheet_name)
HRV_df   = DataFrame(HRV_data)
#HRV_df = DataFrame(HRV_data...)


#Get the data
HRV_df[,"VA_TOT"]
VA_TOT     = HRV_df[HRV_df.year .>= 1947, "VA_TOT"]      #GDP, current prices, millions of dollars
VA_TOT_P   = HRV_df[HRV_df.year .>= 1947, "VA_TOT_P"]    #Value added, price index, aggregate, 1947 = 1

VA_GOOD    = HRV_df[HRV_df.year .>= 1947, "VA_GOOD"]     #Value added, current prices, goods, millions of dollars
VA_SERV    = HRV_df[HRV_df.year .>= 1947, "VA_SERV"]     #Value added, current prices, services, millions of dollars

VAX_GOOD   = HRV_df[HRV_df.year .>= 1947, "VAX_GOOD"]    #Investment value added requirement, services, millions of dollars
VAC_GOOD   = HRV_df[HRV_df.year .>= 1947, "VAC_GOOD"]    #Consumption value added requirement, services, millions of dollars

VAX_SERV   = HRV_df[HRV_df.year .>= 1947, "VAX_SERV"]    #Investment value added requirement, services, millions of dollars
VAC_SERV   = HRV_df[HRV_df.year .>= 1947, "VAC_SERV"]    #Consumption value added requirement, services, millions of dollars

X_TOT      = HRV_df[HRV_df.year .>= 1947, "X_TOT"]       #Fixed investment expenditures, millions of current dollars
X_TOT_P    = HRV_df[HRV_df.year .>= 1947, "X_TOT_P"]     #Chain-type price indexes for fixed investment, 1947=1
X_TOT_QI   = HRV_df[HRV_df.year .>= 1947, "X_TOT_QI"]    #Chain-type quantity indexes for fixed investment, 1947=1

C_TOT      = HRV_df[HRV_df.year .>= 1947, "C_TOT"]       #Consumption defined as GDP-X, current prices, millions of dollars
C_TOT_P    = HRV_df[HRV_df.year .>= 1947, "C_TOT_P"]     #Chain-type price indexes for consumption, 1947=1
C_TOT_QI   = HRV_df[HRV_df.year .>= 1947, "C_TOT_QI"]    #Chain-type quantity indexes for consumption, 1947=1

C_GOOD_QI  = HRV_df[HRV_df.year .>= 1947, "C_GOOD_QI"]   #Chain-type quantity indexes for goods consumption, 1947=1
C_GOOD_P   = HRV_df[HRV_df.year .>= 1947, "C_GOOD_P"]    #Chain-type price indexes for goods consumption, 1947=1

C_SERV_QI  = HRV_df[HRV_df.year .>= 1947, "C_SERV_QI"]   #Chain-type quantity indexes for service consumption, 1947=1
C_SERV_P   = HRV_df[HRV_df.year .>= 1947, "C_SERV_P"]    #Chain-type price indexes for service consumption, 1947=1

#Get the Population
Nt             = HRV_df[HRV_df.year .>= 1980, "POP"]
g_n            = (log(Nt[end]/Nt[1]))/(2017-1980)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Some checks to contrast with the model
log((VA_TOT[end]/VA_TOT_P[end])/(VA_TOT[1]/VA_TOT_P[1]))

rel_P_goods = C_GOOD_P./X_TOT_P 
rel_P_serv  = C_SERV_P./X_TOT_P 

rel_P_goods[end]/rel_P_goods[1] #vs 1.0492755776360312 in the model 
rel_P_serv[end]/rel_P_serv[1]   #vs 1.5474369784534927 in the model
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Now we build a Laspayres, Paasche and Fisher index for the GDP
#Laspeyres
τ_base = 1980
τ      = τ_base - 1980 + 1  
T_base = 2017
T      = T_base - 1980 + 1  

E_Paasche_QI = ( C_GOOD_P[T].*C_GOOD_QI + C_SERV_P[T].*C_SERV_QI)./( C_SERV_P[T]*C_SERV_QI[τ] + C_GOOD_P[T]*C_GOOD_QI[τ])
P_Paasche    = (C_TOT./C_TOT[1])./E_Paasche_QI

E_Laspeyres_QI = ( C_SERV_P[τ].*C_SERV_QI .+ C_GOOD_P[τ].*C_GOOD_QI)./( C_SERV_P[τ].*C_SERV_QI[τ] .+ C_GOOD_P[τ].*C_GOOD_QI[τ])
P_Laspeyres    = (C_TOT./C_TOT[1])./E_Laspeyres_QI 

X_TOT_QI       = X_TOT_QI./X_TOT_QI[1]

y_Laspeyres = (P_Laspeyres[τ].*E_Laspeyres_QI .+ X_TOT_QI  )./(P_Laspeyres[τ].*E_Laspeyres[τ] .+ X_TOT_QI[τ] )
y_Paasche   = (P_Paasche[T].*E_Paasche_QI .+ X_TOT_QI  )./(P_Paasche[T].*E_Paasche_QI[τ] .+ X_TOT_QI[τ] )

y_FS        = sqrt.(y_Laspeyres.*y_Paasche)

log.(y_FS) .+ g_n.*( (1980 .- 1980 .+1))

plot(1980:2017, log.(y_Laspeyres) .+ g_n.*( (1980 .- 1980 .+1) : (2017 .- 1980 .+1) ), 
     label = "Laspeyres", 
     title = "GDP indices",
     xlabel = "Year",
     ylabel = "Log GDP index",
     legend = (0.15, 0.85),
     grid = true,
     lw = 2,
     color = :blue)

plot!(1980:2017, log.(y_Paasche) .+ g_n.*( (1980 .- 1980 .+1) : (2017 .- 1980 .+1) ), 
     label = "Paasche", 
     lw = 2,
     color = :red)

plot!(1980:2017, log.(y_FS) .+ g_n.*( (1980 .- 1980 .+1) : (2017 .- 1980 .+1) ), 
     label = "Fisher-Ideal", 
     lw = 2,
     color = :green)
cd("C:\\Users\\lezjv\\Dropbox\\Omar&Nacho\\Baqaee_Burstein\\Note on Structural Change\\Nacho")
pwd()

using XLSX, DataFrames
using Plots; plotlyjs()

xlsx     = XLSX.readxlsx("BEA_NIPA_Data.xlsx")

#------------------------------------------------------------------------------
#Define the Functions That Compute the Quantity Indices
function paasche_quantity_index_series(quant_df, price_df, items::Vector{String}, initial_year::String, base_year::String, comp_years::Vector{String})
    # Filter dataframes to only include the selected items
    p_base = price_df[in.(price_df.Item, Ref(items)), base_year]
    q_initial = quant_df[in.(quant_df.Item, Ref(items)), initial_year]
    result = Float64[]
    for year in comp_years
        q_comp = quant_df[in.(quant_df.Item, Ref(items)), year]
        numerator = sum(p_base .* q_comp)
        denominator = sum(p_base .* q_initial)
        push!(result, numerator / denominator)
    end
    return result
end

function laspeyres_quantity_index_series(quant_df, price_df, items::Vector{String}, base_year::String, comp_years::Vector{String})
    # Filter dataframes to only include the selected items
    p_base = price_df[in.(price_df.Item, Ref(items)), base_year]
    q_base = quant_df[in.(quant_df.Item, Ref(items)), base_year]
    result = Float64[]
    for year in comp_years
        q_comp = quant_df[in.(quant_df.Item, Ref(items)), year]
        numerator = sum(p_base .* q_comp)
        denominator = sum(p_base .* q_base)
        push!(result, numerator / denominator)
    end
    return result
end

function value_index_series(nom_df, items::Vector{String}, base_year::String, comp_years::Vector{String})
    # Filter dataframe to only include the selected items
    nom_base = nom_df[in.(nom_df.Item, Ref(items)), base_year]
    value_base = sum(nom_base)
    result = Float64[]
    for year in comp_years
        nom_comp = nom_df[in.(nom_df.Item, Ref(items)), year]
        value_comp = sum(nom_comp)
        push!(result, value_comp / value_base)
    end
    return result
end

# Example usage:
years = string.(1947:2024)
items = ["Gross domestic product", "Personal consumption expenditures"]
value_series = value_index_series(Nom_df, items, "1947", years)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Import the BEA Data: Quantity Indices
Q_sheet  = xlsx["T10103-A"]
colnames = Q_sheet["D8:CU8"][1, :]
rownames = vec(Q_sheet["B9:B23"])
Q_values = Q_sheet["D9:CU23"]

Q_df     = DataFrame(Q_values, Symbol.(colnames); makeunique=true)
insertcols!(Q_df, 1, :Item => strip.(string.(rownames)))


#------------------------------------------------------------------------------
#Import the BEA Data: Price Indices
P_sheet  = xlsx["T10104-A"]
colnames = P_sheet["D8:CU8"][1, :]
rownames = vec(P_sheet["B9:B22"])
P_values = P_sheet["D9:CU22"]

P_df     = DataFrame(P_values, Symbol.(colnames); makeunique=true)
insertcols!(P_df, 1, :Item => strip.(string.(rownames)))

#------------------------------------------------------------------------------
#Import the BEA Data: Nominal Values
Nom_sheet  = xlsx["T10105-A"]
colnames   = Nom_sheet["D8:CU8"][1, :]
rownames   = vec(Nom_sheet["B9:B22"])
Nom_values = Nom_sheet["D9:CU22"]

Nom_df     = DataFrame(Nom_values, Symbol.(colnames); makeunique=true)
insertcols!(Nom_df, 1, :Item => strip.(string.(rownames)))

#------------------------------------------------------------------------------
#Import the BEA Data: Price Indices
P_sheet  = xlsx["T10104-A"]
colnames = P_sheet["D8:CU8"][1, :]
rownames = vec(P_sheet["B9:B22"])
P_values = P_sheet["D9:CU22"]

P_df     = DataFrame(P_values, Symbol.(colnames); makeunique=true)
insertcols!(P_df, 1, :Item => strip.(string.(rownames)))

years             = string.(1947:2024)

items             = ["Nondurable goods", "Services"]
Paasche_Q_Index_1 = paasche_quantity_index_series(Q_df, P_df, items, "1947", "2024", years)
log.(Paasche_Q_Index_1)

items             = ["Durable goods","Equipment","Intellectual property products" ]
Paasche_Q_Index_2 = paasche_quantity_index_series(Q_df, P_df, items, "1947", "2024", years)
log.(Paasche_Q_Index_2)

items             = ["Structures","Residential"]
Paasche_Q_Index_3 = paasche_quantity_index_series(Q_df, P_df, items, "1947", "2024", years)
log.(Paasche_Q_Index_3)




years             = string.(1947:2024)

items             = ["Nondurable goods", "Services"]
Lasp_Q_Index_1    = laspeyres_quantity_index_series(Q_df, P_df, items, "1947", years)
log.(Lasp_Q_Index_1)

items             = ["Durable goods","Equipment","Intellectual property products" ]
Lasp_Q_Index_2    = laspeyres_quantity_index_series(Q_df, P_df, items, "1947", years)
log.(Lasp_Q_Index_2)

items             = ["Structures","Residential"]
Lasp_Q_Index_3    = laspeyres_quantity_index_series(Q_df, P_df, items, "1947", years)
log.(Lasp_Q_Index_3)



years             = string.(1947:2024)
items             = ["Nondurable goods", "Services"]
Value_Index_1     = value_index_series(Nom_df, items, "1947", years)

items             = ["Durable goods","Equipment","Intellectual property products" ]
Value_Index_2     = value_index_series(Nom_df, items, "1947", years)

items             = ["Structures","Residential"]
Value_Index_3     = value_index_series(Nom_df, items, "1947", years)
#-------------------------------------------------------------------------------

Paasche_P_Index_1 = Value_Index_1./Lasp_Q_Index_1
Paasche_P_Index_2 = Value_Index_2./Lasp_Q_Index_2
Paasche_P_Index_3 = Value_Index_3./Lasp_Q_Index_3

Lasp_P_Index_1    = Value_Index_1./Paasche_Q_Index_1
Lasp_P_Index_2    = Value_Index_2./Paasche_Q_Index_2
Lasp_P_Index_3    = Value_Index_3./Paasche_Q_Index_3

function paasche_quantity_index_agg(prices, quantities, initial_year_idx, base_year_idx)
    n_groups = length(prices)
    n_years = length(prices[1])
    idx = Float64[]
    denominator = sum(prices[g][base_year_idx] * quantities[g][initial_year_idx] for g in 1:n_groups)
    for t in 1:n_years
        numerator = sum(prices[g][base_year_idx] * quantities[g][t] for g in 1:n_groups)
        push!(idx, numerator / denominator)
    end
    return idx
end

function laspeyres_quantity_index_agg(prices, quantities, initial_year_idx)
    n_groups = length(prices)
    n_years = length(prices[1])
    idx = Float64[]
    denominator = sum(prices[g][initial_year_idx] * quantities[g][initial_year_idx] for g in 1:n_groups)
    for t in 1:n_years
        numerator = sum(prices[g][initial_year_idx] * quantities[g][t] for g in 1:n_groups)
        push!(idx, numerator / denominator)
    end
    return idx
end

quantities = [Lasp_Q_Index_1   , Lasp_Q_Index_2  , Lasp_Q_Index_3]
prices     = [Paasche_P_Index_1,Paasche_P_Index_2,Paasche_P_Index_3]

# Set the indices for the initial and base years
initial_year_idx = 1                # e.g., 1 for the first year (1947)
base_year_idx = length(years)       # e.g., last year (2024), or set as needed

paasche = paasche_quantity_index_agg(prices, quantities, 1, base_year_idx)
log.(paasche)

quantities = [Paasche_Q_Index_1,Paasche_Q_Index_2, Paasche_Q_Index_3]
prices     = [Lasp_P_Index_1   , Lasp_P_Index_2  , Lasp_P_Index_3]

laspeyres = laspeyres_quantity_index_agg(prices, quantities, initial_year_idx)
log.(laspeyres)

fisher_shell = sqrt.(paasche .* laspeyres)
log.(fisher_shell)

plot(1947:2024, log.(paasche), 
     label = "Paasche", 
     title = "GDP indices",
     xlabel = "Year",
     ylabel = "Log Index",
     legend = :topleft,
     grid = true)

plot!(1947:2024, log.(laspeyres),
     label = "Laspeyres", 
     lw = 2,
     color = :red)

plot!(1947:2024, log.(fisher_shell),
     label = "Fisher-Ideal", 
     lw = 2,
     color = :green)

#-----------------------------------------------------------------
#Now Suppose that we First Aggregate the Data Using the Fisher-Ideal Indices and Compute Paasche,Laspayres 
#and Fisher-Ideal indices afterwards


Fischer_Q_Index_1 = sqrt.(Paasche_Q_Index_1 .* Lasp_Q_Index_1)
Fischer_Q_Index_2 = sqrt.(Paasche_Q_Index_2 .* Lasp_Q_Index_2)
Fischer_Q_Index_3 = sqrt.(Paasche_Q_Index_3 .* Lasp_Q_Index_3)

Fischer_P_Index_1    = Value_Index_1./Fischer_Q_Index_1
Fischer_P_Index_2    = Value_Index_2./Fischer_Q_Index_2
Fischer_P_Index_3    = Value_Index_3./Fischer_Q_Index_3

quantities = [Fischer_Q_Index_1,Fischer_Q_Index_2,Fischer_Q_Index_3]
prices     = [Fischer_P_Index_1,Fischer_P_Index_2,Fischer_P_Index_3]

initial_year_idx = 1                # e.g., 1 for the first year (1947)
base_year_idx = length(years)       # e.g., last year (2024), or set as needed

paasche = paasche_quantity_index_agg(prices, quantities, 1, base_year_idx)
log.(paasche)

laspeyres = laspeyres_quantity_index_agg(prices, quantities, initial_year_idx)
log.(laspeyres)

fisher_shell = sqrt.(paasche .* laspeyres)
log.(fisher_shell)

plot(1947:2024, log.(paasche), 
     label = "Paasche", 
     title = "GDP indices",
     xlabel = "Year",
     ylabel = "Log Index",
     legend = :topleft,
     grid = true)

plot!(1947:2024, log.(laspeyres),
     label = "Laspeyres", 
     lw = 2,
     color = :red)

plot!(1947:2024, log.(fisher_shell),
     label = "Fisher-Ideal", 
     lw = 2,
     color = :green)

#-----------------------------------------------------------------
years             = string.(1947:2024)
items             = ["Goods","Services"]

Paasche_Q_Index_11 = paasche_quantity_index_series(Q_df, P_df, items, "1947", "2024", years)
Lasp_Q_Index_11    = laspeyres_quantity_index_series(Q_df, P_df, items, "1947", years)

Value_Index_11     = value_index_series(Nom_df, items, "1947", years)
Paasche_P_Index_11 = Value_Index_11./Lasp_Q_Index_11 
Lasp_P_Index_11    = Value_Index_11./Paasche_Q_Index_11 

items              = ["Gross private domestic investment"]
Paasche_Q_Index_21 = paasche_quantity_index_series(Q_df, P_df, items, "1947", "2024", years)
Lasp_Q_Index_21    = laspeyres_quantity_index_series(Q_df, P_df, items, "1947", years)
log.(Lasp_Q_Index_21)

Value_Index_21     = value_index_series(Nom_df, items, "1947", years)
Paasche_P_Index_21 = Value_Index_21./Lasp_Q_Index_21 
Lasp_P_Index_21    = Value_Index_21./Paasche_Q_Index_21 


quantities = [Paasche_Q_Index_11,Paasche_Q_Index_21]
prices     = [Lasp_P_Index_11   , Lasp_P_Index_21  ]

initial_year_idx = 1                # e.g., 1 for the first year (1947)
base_year_idx = length(years)       # e.g., last year (2024), or set as needed

paasche = paasche_quantity_index_agg(prices, quantities, 1, base_year_idx)
log.(paasche)


quantities = [Lasp_Q_Index_11    ,Lasp_Q_Index_21    ] 
prices     = [Paasche_P_Index_11 ,Paasche_P_Index_21 ]

laspeyres = laspeyres_quantity_index_agg(prices, quantities, initial_year_idx)
log.(laspeyres)

fisher_shell = sqrt.(paasche .* laspeyres)
log.(fisher_shell)

plot(1947:2024, log.(paasche), 
     label = "Paasche", 
     title = "GDP indices",
     xlabel = "Year",
     ylabel = "Log Index",
     legend = :topleft,
     grid = true)

plot!(1947:2024, log.(laspeyres),
     label = "Laspeyres", 
     lw = 2,
     color = :red)

plot!(1947:2024, log.(fisher_shell),
     label = "Fisher-Ideal", 
     lw = 2,
     color = :green)


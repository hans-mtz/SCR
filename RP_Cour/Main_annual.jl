# ## Loading modules
# using Distributed
# using DataFrames
# using CSV
#
# addprocs(6)
#
# ## Load modules in parallel
# @everywhere include("RP_Cournot/modules/MyModules.jl")
#
# @everywhere begin
#     using .MyModules
#     using Random
#     using SharedArrays
#     ## Setting seed
#     Random.seed!(66636)
# end
#
# ## Settin up directories
#
# dir = @__DIR__ ;
# datadir=dir*"/data/oildata.csv" ;
# resultsdir=dir*"/results/"
#
# ## Loading data
#
# pt, qt = load_data_select(datadir)

## Initializing array for results
table=SharedArray{Float64,2}(zeros(36,3))

## Testing for each year for the

@time MyModules.shared_fill_table!(table, pt, qt)

## Saving results to csv

annual=DataFrame(Year= collect(1973:2008),
                 TS=table[:,1],
                 pval=table[:,2],
                 stars=mystars.(table[:,2]))
# CSV.write(resultsdir*"annual.csv", annual)
CSV.write(resultsdir*"annual_sub_$set.csv", annual)

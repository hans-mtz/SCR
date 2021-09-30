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

mchain=SharedArray{Float64,3}(zeros(Integer(size(pt,1)/6),size(qt,2),Param.nsim))
decade=SharedArray{Float64,2}(zeros(4,2))
acc_r_t=SharedArray{Float64,1}(zeros(5))
TS, pvalue = nothing, nothing
## Testing

# Markov Chain
@time MyModules.shared_fill_chain!(mchain, pt, qt, fill_chain_sem_algo!)
# Joint test
@time TS, pvalue, acc_r_t[1] = MyModules.RP_Cournot_test_par(mchain)
# Testing by decade
@time begin
    dec = [0, 16, 36, 56, 72]
    for i=1:4
        m = mchain[dec[i]+1:dec[i+1],:,:]
        decade[i,1], decade[i,2], acc_r_t[i+1] = MyModules.RP_Cournot_test_par(m)
    end
end

## Saving results to csv

joint_sem_algo=DataFrame(Time="1980-2008",
                 TS=TS,
                 pvalue=pvalue,
                 stars=mystars.(pvalue))
CSV.write(resultsdir*"joint_ak_$set.csv", joint_sem_algo)

decade_sem_algo=DataFrame(Time=["1973-1980", "1981-1990","1991-2000","2001-2008"],
                 TS=decade[:,1],
                 pvalue=decade[:,2],
                 stars=mystars.(decade[:,2]))
CSV.write(resultsdir*"decade_ak_$set.csv", decade_sem_algo)

acc_rates_algo=DataFrame(Ar=acc_r_t)
CSV.write(resultsdir*"accp_r_ak_$set.csv",acc_rates_algo)

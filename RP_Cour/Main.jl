## Loading modules
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
println("Initializing vars ")
mchain=SharedArray{Float64,3}(zeros(Integer(size(pt,1)/6),size(qt,2),Param.nsim))
decade=SharedArray{Float64,2}(zeros(4,2))
acc_r_t=SharedArray{Float64,1}(zeros(5))

## Testing

# Markov Chain
println("Generating MCMC")
@time MyModules.shared_fill_chain!(mchain, pt, qt, fill_chain_sem!)
# Joint test
println("Testing Jointly")
@time TS, pvalue, acc_r_t[1] = MyModules.RP_Cournot_test_par(mchain)
# Testing by decade
println("Testing by decade")
@time begin
    dec = [0, 16, 36, 56, 72]
    for i=1:4
        m = mchain[dec[i]+1:dec[i+1],:,:]
        decade[i,1], decade[i,2], acc_r_t[i+1] = MyModules.RP_Cournot_test_par(m)
    end
end

## Saving results to csv
println("Saving results")
joint_sem=DataFrame(Time="1973-2008",
                 TS=TS,
                 pvalue=pvalue,
                 stars=mystars.(pvalue))
CSV.write(resultsdir*"joint_$set.csv", joint_sem)

decade_sem=DataFrame(Time=["1973-1980", "1981-1990","1991-2000","2001-2008"],
                 TS=decade[:,1],
                 pvalue=decade[:,2],
                 stars=mystars.(decade[:,2]))
CSV.write(resultsdir*"decade_$set.csv", decade_sem)

acc_rates=DataFrame(Ar=acc_r_t)
CSV.write(resultsdir*"accp_r_$set.csv",acc_rates)

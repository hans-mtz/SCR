## Loading modules
using Distributed
using Plots

addprocs(6)

# Load modules in parallel
@everywhere begin
    include("RP_Cournot/modules/CourErr.jl")
    include("RP_Cournot/modules/MarkovChain.jl")
    include("RP_Cournot/modules/Testing.jl")
    include("modules/Parallel_mod.jl")
    using Random
    ## Setting seed
    Random.seed!(66636)
end


## Settin up directories

dir = @__DIR__ ;
datadir=dir*"/data/oildata.csv" ;
resultsdir=dir*"/results/"


## CourErr

# p, q = load_data(datadir)
pt, qt = MyModules.load_data_select(datadir)

## Markov chain

# Testing
p, q = pt[1:6], qt[1:6,:]
# eδ, eq, status = warm_start(p,q)
# w, δsim, qsim = jump(p,q,eδ,eq)

MyModules.jump_algo()
mchain, status = MyModules.gchain(p,q)
malg, salg =  MyModules.gchain_algo(p,q)

histogram(1:500000,malg[1,:])
histogram(1:500000,mchain[1,:])

mchain[:,1].==mchain[:,rand(1:500000)]
malg[:,1].==malg[:,rand(1:500000)]
rand(1:500000)
## Testing

objfun(collect(1:6), m=malg[:,1:1000])
W = objfun_bbo(collect(1:6), m=malg[:,1:1000])
@time f, d = my_dist(mchain[:,1:1000], collect(1:6))
@time my_dist(mchain[:,rand(1:500000,1000)], ones(6))
@time my_dist(malg[:,rand(1:500000,1000)], ones(6))
# f = ecdf(d)
# cdf(f,W)
f(W*Param.nfast)
1-f(W*Param.nfast)
# # TS_b=T*obj_val
quantile(d,[0.05, 0.5, 0.95])

quantile(d,0.05)
quantile(d,0.95)
pval = f(W)

# Testing the functions

# p, q = pt[1:12], qt[1:12,:]
# mchain, status = gchain(p,q)

@time TS, pval = RP_Cournot_test(mchain)

## Testing

# mchain_op = mchain[:,:,1:20000]
#
# objfun(collect(1:6))
#
# @time objfun_par(collect(1:6))
# @time objfun_s(collect(1:6))
#
# @time objfun_s_bbo(collect(1:6))
# @time objfun_bbo_par(collect(1:6))
# #
# @time TS, pval = RP_Cournot_test_sem(mchain)
# @time TS, pval = RP_Cournot_test_par(mchain)

## Rejection rate
acc_rate_table=SharedArray{Float64,2}(zeros(T*2))
@distributed for i=1:72
    my_dist(mchain[i,:,:], ones(6))
end

##

MyModules.acc_rate(mchain,ones(6),f=1)

##
MyModules.sim_err(malg,ones(6)*100)

@time RP_Cournot_test(malg)

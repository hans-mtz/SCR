
using Statistics
using BlackBoxOptim
using NLopt
using Distributions

## Optimization
# Generate random index
indfast=rand(1:Param.nsim, Param.nfast)
indfast[1]=1 # Keeping the first element fixed.

# schain=SharedArray{Float64,2}(mchain[indfast,:])
mchain_op=mchain[indfast,:]

## Objective function

function objfun(gamma0...; Par::Par=Param, mchain=mchain_op)
    @unpack nsim, nfast = Par
    gamma = gamma0
    N = size(mchain[1,:])
    val=0
    logunif=log.(rand(nfast))
    geta=mchain[1,:]
    hM=zeros(N)
    for r=1:nfast
        gtry=mchain[r,:]
        for j=1:size(gtry,2)
            val=+gtry[:,j].*gamma[j]-geta[:,j].*gamma[j]
        end
        geta=logunif[r] < val[1] ? gtry : geta
        hM=+geta./nfast
    end
    mhM=mean(hM)
    var=Statistics.var(hM)

    obj=(mhM^2)/var

    return obj
end

function objfun_bbo(gamma0; Par::Par=Param, mchain=mchain_op)
    @unpack nsim, nfast = Par
    gamma = gamma0
    N = size(mchain[1,:])
    val=0
    logunif=log.(rand(nfast))
    geta=mchain[1,:]
    hM=zeros(N)
    for r=1:nfast
        gtry=mchain[r,:]
        for j=1:size(gtry,2)
            val=+gtry[:,j].*gamma[j]-geta[:,j].*gamma[j]
        end
        geta=logunif[r] < val[1] ? gtry : geta
        hM=+geta./nfast
    end
    mhM=mean(hM)
    var=Statistics.var(hM)

    obj=(mhM^2)/var

    return obj
end
## testing obj fun

objfun(collect(1:6))

objfun_bbo(repeat([200],6))
objfun_bbo(collect(1:6))
#
fop=nothing
## Solvin BBO
res = bboptimize(x->objfun_bbo(x); Sea = (-10e300,10e300), NumDimensions = 6)


## Getting results

minbb = best_fitness(res)
γ0 = best_candidate(res)
TS_bb = minbb*6
p_val_bb = ccdf(Chisq(4),TS_bb)
## Optimizing

fop=nothing
fop=Model(NLopt.Optimizer)
set_optimizer_attribute(fop, "algorithm", :LN_BOBYQA)
set_optimizer_attribute(fop, "xtol_rel", 1e-8)
@variable(fop, gamma0[1:6])

register(fop, :objfun, 6, objfun; autodiff=true )
@NLobjective(fop, Min, objfun(gamma0...))

for i=1:6
    set_start_value(gamma0[i],γ0[i])
end

## Solving in NLopt

JuMP.optimize!(fop)

## Getting Results

termination_status(fop)
primal_status(fop)

value.(gamma0)
obj_v=objective_value(fop)
# n=size(q,2)
TS=6*obj_v
p_val = ccdf(Chisq(4),TS)

## Testing function 

function RP_Cournot_test(;Par::Par=Param, mchain=mchain_op)
    N = size(mchain,2)
    # First optimization BlackBox
    res = bboptimize(x->objfun_bbo(x); Sea = (-10e300,10e300), NumDimensions = 6)
    ## Getting results
    minbb = best_fitness(res)
    γ0 = best_candidate(res)
    # TS_bb = minbb*N
    # p_val_bb = ccdf(Chisq(N),TS_bb)

    ## Second optimization NLopt
    fop=nothing
    fop=Model(NLopt.Optimizer)
    set_optimizer_attribute(fop, "algorithm", :LN_BOBYQA)
    set_optimizer_attribute(fop, "xtol_rel", 1e-8)
    @variable(fop, gamma0[1:N])

    register(fop, :objfun, N, objfun; autodiff=true )
    @NLobjective(fop, Min, objfun(gamma0...))

    for i=1:N
        set_start_value(gamma0[i],γ0[i])
    end

    ## Solving in NLopt
    JuMP.optimize!(fop)

    ## Getting Results
    termination_status(fop)
    primal_status(fop)

    γ1=value.(gamma0)
    obj_v=objective_value(fop)
    TS=N*obj_v
    p_val = ccdf(Chisq(N),TS)

    ## Second optimization for refinement
    for i=1:N
        set_start_value(gamma0[i],γ1[i])
    end
    ## Solving in NLopt
    JuMP.optimize!(fop)
    obj_v=objective_value(fop)
    TS=N*obj_v
    p_val = ccdf(Chisq(N),TS)

    fop=nothing

    return TS, p_val
end

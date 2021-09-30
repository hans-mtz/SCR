using Statistics
using BlackBoxOptim
using NLopt
using Distributions
using StatsBase

## Objective functions year

function objfun(g...; m=mchain_op)
    # @unpack nsim, nfast = Par
    gamma = g
    N, nfast = size(m)

    logunif=log.(rand(nfast))
    geta=m[:,1]
    hM=zeros(N)
    h=zeros(N,nfast)
    # ar=0
    for r=1:nfast
        val=0
        gtry=m[:,r]
        for j=1:N
            val=+gtry[j].*gamma[j]-geta[j].*gamma[j]
        end
        geta=logunif[r] < val[1] ? gtry : geta
        # ar=+ logunif[r] < val[1] ? 1 : 0
        hM=+geta
        h[:,r]=geta
    end
    # acc_r=round((ar/nfast)*100, digits=2)
    # println("Acceptance rate is : ",acc_r)
    # hM=hM./nfast
    # var=(h.-hM)*(h.-hM)'
    #
    # lambda, QM = eigen(var)
    # sel = lambda.>0
    # An = QM[:,sel]
    # dvs = An'*(h.-hM)
    # vs = An'*var*An
    # omega = inv(vs)
    # obj=dvs'*omega*dvs
    #
    # # obj=hM'*inv(var)*hM
    #
    # return obj[1]
    hM=hM./nfast
    # var=(h.-hM)*(h.-hM)'
    bn=Integer(round(sqrt(N)))
    bt=Integer(N-bn+1)
    d=zeros(N,bt)
    for b=1:bt
        ind=b:bn+b-1
        hMb=sum(h[:,ind], dims=2)./bn
        d[:,b] = (hMb.-hM)
    end
    dv= d*d'
    vhat=(bn/(bt-1)).*dv
    lambda, QM = eigen(vhat)
    sel = lambda.>0
    An = QM[:,sel]
    dvs = An'*(hM)
    vs = An'*vhat*An
    omega = inv(vs)
    t=dvs'*omega*dvs
    t[1]
end

function objfun_bbo(gamma0; m=mchain_op)
    # @unpack nfast = Par
    gamma = gamma0
    N, nfast = size(m)
    # nfast = size(m,1)

    logunif=log.(rand(nfast))
    geta=m[:,1]
    hM=zeros(N)
    h=zeros(N,nfast)
    # ar=0
    for r=1:nfast
        val=0
        gtry=m[:,r]
        for j=1:N
            val=+gtry[j].*gamma[j]-geta[j].*gamma[j]
        end
        geta=logunif[r] < val[1] ? gtry : geta
        # ar=+ logunif[r] < val[1] ? 1 : 0
        h[:,r]=geta
        hM=+geta
    end
    # acc_r=round((ar/nfast)*100, digits=2)
    # println("Acceptance rate is : ",acc_r)
    # hM=hM./nfast
    # var=(h.-hM)*(h.-hM)'
    #
    # lambda, QM = eigen(var)
    # sel = lambda.>0
    # An = QM[:,sel]
    # dvs = An'*(h.-hM)
    # vs = An'*var*An
    # omega = inv(vs)
    # obj=dvs'*omega*dvs
    #
    # # obj=hM'*inv(var)*hM
    #
    # return obj[1]

    hM=hM./nfast
    # var=(h.-hM)*(h.-hM)'
    bn=Integer(round(sqrt(N)))
    bt=Integer(N-bn+1)
    d=zeros(N,bt)
    for b=1:bt
        ind=b:bn+b-1
        hMb=sum(h[:,ind], dims=2)./bn
        d[:,b] = (hMb.-hM)
    end
    dv= d*d'
    vhat=(bn/(bt-1)).*dv
    lambda, QM = eigen(vhat)
    sel = lambda.>0
    An = QM[:,sel]
    dvs = An'*(hM)
    vs = An'*vhat*An
    omega = inv(vs)
    t=dvs'*omega*dvs
    t[1]
end
## Bootstrap testing

function my_dist(m, g; P::Par=Param)
    @unpack nboost = P
    # nboost = 20000
    N, nfast = size(m)
    gamma=g
    logunif=log.(rand(nfast))
    geta=m[:,1]
    hM=zeros(N)
    h=zeros(N,nfast)
    d=zeros(nboost)
    ar=0.0
    for r=1:nfast
        val=0.0
        gtry=m[:,r]
        for j=1:N
            val=+gtry[j].*gamma[j]-geta[j].*gamma[j]
        end
        geta=logunif[r] < val[1] ? gtry : geta

        logunif[r] < val[1] && (ar=+1.0)

        h[:,r]=geta
        hM=+geta
    end
    acc_r=100.0*(ar/nfast)
    # println("Acceptance rate is : ", acc_r)
    @show (acc_r)

    hM=hM./nfast
    var=(h.-hM)*(h.-hM)'

    for b=1:nboost
        ind=rand(1:nfast,nfast)
        hMb=sum(h[:,ind], dims=2)./nfast
        varb=(h[:,ind].-hMb)*(h[:,ind].-hMb)'
        lambda, QM = eigen(varb)
        sel = lambda.>0
        An = QM[:,sel]
        dvs = An'*(hMb.-hM)
        vs = An'*varb*An
        omega = inv(vs)
        t=dvs'*omega*dvs
        d[b]=nfast.*t[1]
    end
    # d=[nfast*(h[:,i]-hM)'*inv(var)*(h[:,i]-hM) for i=1:nfast]
    f=ecdf(d)
    f, d
end

function sim_err(m, g; P::Par=Param)
    @unpack nboost = P
    # nboost = 20000
    N, nfast = size(m)
    gamma=g
    logunif=log.(rand(nfast))
    geta=m[:,1]
    hM=zeros(N)
    h=zeros(N,nfast)
    ar=0.0
    for r=1:nfast
        val=0.0
        gtry=m[:,r]
        for j=1:N
            val=+gtry[j].*gamma[j]-geta[j].*gamma[j]
        end
        geta=logunif[r] < val[1] ? gtry : geta

        # logunif[r] < val[1] && (ar=+1.0)

        h[:,r]=geta
        hM=+geta
    end
    # acc_r=100.0*(ar/nfast)
    # # println("Acceptance rate is : ", acc_r)
    # @show (acc_r)

    hM=hM./nfast
    # var=(h.-hM)*(h.-hM)'
    bn=Integer(round(sqrt(N)))
    bt=Integer(N-bn+1)
    d=zeros(N,bt)
    for b=1:bt
        ind=b:bn+b-1
        hMb=sum(h[:,ind], dims=2)./bn
        d[:,b] = (hMb.-hM)
    end
    dv= d*d'
    vhat=(bn/(bt-1)).*dv
    lambda, QM = eigen(vhat)
    sel = lambda.>0
    An = QM[:,sel]
    dvs = An'*(hM)
    vs = An'*vhat*An
    omega = inv(vs)
    t=dvs'*omega*dvs
    t[1]
end

function acc_rate(m, g; P::Par=Param, f=0)
    # @unpack nboost = P
    # nboost = 20000
    d=ndims(m)
    if d==3
        T, N, nfast = size(m)
        geta=m[:,:,1]
    else
        T=1
        N, nfast = size(m)
        geta=m[:,1]
    end

    gamma=g
    unif=rand(nfast)

    # hM=zeros(N)
    # h=zeros(N,nfast)
    # d=zeros(nboost)
    ar=0.0
    for t = 1:T, r=1:nfast
        val=0.0
        d==3 ? gtry=m[t,:,r] : gtry=m[:,r]
        for j=1:N
            val=+gtry[j].*gamma[j]-geta[j].*gamma[j]
        end
        # geta=logunif[r] < val[1] ? gtry : geta
        if f==0
            # log.(unif[r]) < val[1] && (ar=+1.0)
            ars=log(unif[r]) < val[1] ? 1.0 : 0.0
        else
            ars=unif[r] < 1.4^(val[1]) ? 1.0 : 0.0
        end
        ar=+ars
        # ar=+logunif[r] < val[1] ? 1.0 : 0.0
        # h[:,r]=geta
        # hM=+geta
    end
    acc_r=100.0*(ar/nfast)
    # println("Acceptance rate is : ", acc_r)
    @show (acc_r)
    acc_r
end

## Testing function

function RP_Cournot_test(m ; Par::Par=Param)
    @unpack nfast = Par
    N, nsim = size(m)
    # Generate random index
    indfast=rand(1:nsim, nfast)
    indfast[1]=1 # Keeping the first element fixed.
    mchain_op=m[:,indfast]

    # First optimization BlackBox
    res = bboptimize(x->objfun_bbo(x, m=mchain_op); SearchRange = (-10e300,10e300), NumDimensions = N, TraceMode=:silent, MaxTime=400.0)
    ## Getting results
    # minbb = best_fitness(res)
    γ0 = best_candidate(res)

    # mycdf, dist = my_dist(mchain_op,γ0)
    # TS = minbb*nfast
    # p_val = mycdf(TS)
    println("First optimization from ", myid()," done.")
    # Second optimization NLopt
    fop=nothing
    fop=Model(NLopt.Optimizer)
    set_optimizer_attribute(fop, "algorithm", :LN_BOBYQA)
    set_optimizer_attribute(fop, "xtol_rel", 1e-6)
    # set_optimizer_attribute(fop, "maxtime", 300)
    @variable(fop, gamma0[1:N])

    register(fop, :objfun, N, (x...)->objfun(x..., m=mchain_op); autodiff=true )
    @NLobjective(fop, Min, objfun(gamma0...))

    for i=1:N
        set_start_value(gamma0[i],γ0[i])
    end

    ## Solving in NLopt
    JuMP.optimize!(fop)

    ## Getting Results
    # termination_status(fop)
    # primal_status(fop)

    γ1=value.(gamma0)
    # obj_v=objective_value(fop)
    # TS=obj_v*nfast
    # mycdf, dist = my_dist(mchain_op,γ1)
    # p_val = mycdf(TS)
    println("Second optimization from ", myid()," done.")
    ## Second optimization for refinement
    for i=1:N
        set_start_value(gamma0[i],γ1[i])
    end
    ## Solving in NLopt
    JuMP.optimize!(fop)

    ## Getting Results
    γ2=value.(gamma0)
    obj_v=objective_value(fop)
    TS=obj_v
    # mycdf, dist = my_dist(mchain_op,γ2)
    # p_val = 1-mycdf(TS)
    p_val = ccdf(Chisq(N),TS)
    println("Third optimization from ", myid()," done.")
    # fop=nothing
    acc_r = MyModules.acc_rate(mchain_op,γ2)
    return TS, p_val , acc_r
end


## Objective functions semester

# function objfun_s(gamma0...; mchain=mchain_op)
#     # @unpack nsim, nfast = Par
#     gamma = Vector(gamma0[1])
#     T,N,nfast = size(mchain)
#
#     logunif=log.(rand(T,nfast))
#     geta=mchain[:,:,1]
#     hM=zeros(T,N)
#     for t=1:T, r=1:nfast
#         val=0
#         gtry=mchain[t,:,r]
#         for j=1:N
#             val=+gtry[j].*gamma[j]-geta[t,j].*gamma[j]
#         end
#         geta[t,:]=logunif[t,r] < val[1] ? gtry : geta[t,:]
#         hM[t,:]=+geta[t,:]./nfast
#     end
#     mhM=sum(hM, dims=1)./T
#     omega=(hM.-mhM)'*(hM.-mhM)
#
#     obj=mhM*inv(omega)*mhM'
#
#     return obj[1]
# end
#
# function objfun_s_bbo(gamma0; mchain=mchain_op)
#     # @unpack nfast = Par
#     gamma = gamma0
#     T,N,nfast = size(mchain)
#
#     logunif=log.(rand(T,nfast))
#     geta=mchain[:,:,1]
#     hM=zeros(T,N)
#     for t=1:T, r=1:nfast
#         val=0
#         gtry=mchain[t,:,r]
#         for j=1:N
#             val=+gtry[j].*gamma[j]-geta[t,j].*gamma[j]
#         end
#         geta[t,:]=logunif[t,r] < val[1] ? gtry : geta[t,:]
#         hM[t,:]=+geta[t,:]./nfast
#     end
#     mhM=sum(hM, dims=1)./T
#     omega=(hM.-mhM)'*(hM.-mhM)
#
#     obj=mhM*inv(omega)*mhM'
#
#     return obj[1]
# end

# Testing function
#
# function RP_Cournot_test_sem(mchain ; Par::Par=Param)
#     @unpack nfast = Par
#     T, N, nsim = size(mchain)
#     # Generate random index
#     indfast=rand(1:nsim, nfast)
#     indfast[1]=1 # Keeping the first element fixed.
#     mchain_op=mchain[:,:,indfast]
#
#     # First optimization BlackBox
#     res = bboptimize(x->objfun_s_bbo(x, mchain=mchain_op); SearRange = (-10e300,10e300), NumDimensions = N, TraceMode=:silent)
#     ## Getting results
#     minbb = best_fitness(res)
#     γ0 = best_candidate(res)
#     # TS_bb = minbb*N
#     # p_val_bb = ccdf(Chisq(N),TS_bb)
#     println("First optimization from ", myid()," done.")
#     ## Second optimization NLopt
#     fop=nothing
#     fop=Model(NLopt.Optimizer)
#     set_optimizer_attribute(fop, "algorithm", :LN_BOBYQA)
#     set_optimizer_attribute(fop, "xtol_rel", 1e-8)
#     @variable(fop, gamma0[1:N])
#
#     register(fop, :objfun_s, N, (x...)->objfun_s(x..., mchain=mchain_op); autodiff=true )
#     @NLobjective(fop, Min, objfun_s(gamma0...))
#
#     for i=1:N
#         set_start_value(gamma0[i],γ0[i])
#     end
#
#     ## Solving in NLopt
#     JuMP.optimize!(fop)
#
#     ## Getting Results
#     # termination_status(fop)
#     # primal_status(fop)
#
#     γ1=value.(gamma0)
#     obj_v=objective_value(fop)
#     TS=T*obj_v
#     p_val = ccdf(Chisq(N),TS)
#     println("Second optimization from ", myid()," done.")
#     ## Second optimization for refinement
#     for i=1:N
#         set_start_value(gamma0[i],γ1[i])
#     end
#     ## Solving in NLopt
#     JuMP.optimize!(fop)
#
#     ## Getting Results
#     obj_v=objective_value(fop)
#     TS=T*obj_v
#     p_val = ccdf(Chisq(N),TS)
#     println("Third optimization from ", myid()," done.")
#     fop=nothing
#
#     return TS, p_val
# end
##

function RP_Cournot_test_par(mchain ; Par::Par=Param)
    @unpack nfast = Par
    T, N, nsim = size(mchain)
    # Generate random index
    indfast=rand(1:nsim, nfast)
    indfast[1]=1 # Keeping the first element fixed.
    mchain_op=mchain[:,:,indfast]

    # First optimization BlackBox
    res = bboptimize(x->objfun_bbo_par(x, mchain=mchain_op); SearchRange = (-10e300,10e300), NumDimensions = N, TraceMode=:silent, MaxTime=400.0)
    ## Getting results
    minbb = best_fitness(res)
    γ0 = best_candidate(res)
    # TS_bb = minbb*N
    # p_val_bb = ccdf(Chisq(N),TS_bb)
    println("First optimization from ", myid()," done.")
    ## Second optimization NLopt
    fop=nothing
    fop=Model(NLopt.Optimizer)
    set_optimizer_attribute(fop, "algorithm", :LN_BOBYQA)
    set_optimizer_attribute(fop, "xtol_rel", 1e-6)
    set_optimizer_attribute(fop, "maxtime", 300)
    @variable(fop, gamma0[1:N])

    register(fop, :objfun_par, N, (x...)->objfun_par(x..., mchain=mchain_op); autodiff=true )
    @NLobjective(fop, Min, objfun_par(gamma0...))

    for i=1:N
        set_start_value(gamma0[i],γ0[i])
    end

    ## Solving in NLopt
    JuMP.optimize!(fop)

    ## Getting Results
    # termination_status(fop)
    # primal_status(fop)

    γ1=value.(gamma0)
    obj_v=objective_value(fop)
    TS=T*obj_v
    p_val = ccdf(Chisq(N),TS)
    println("Second optimization from ", myid()," done.")
    ## Second optimization for refinement
    for i=1:N
        set_start_value(gamma0[i],γ1[i])
    end
    ## Solving in NLopt
    JuMP.optimize!(fop)

    ## Getting Results
    γ2=value.(gamma0)
    obj_v=objective_value(fop)
    TS=T*obj_v
    p_val = ccdf(Chisq(N),TS)
    println("Third optimization from ", myid()," done.")
    fop=nothing
    acc_r = MyModules.acc_rate(mchain_op, γ2)
    return TS, p_val, acc_r
end

function RP_Cournot_test_par_alt(mchain ; Par::Par=Param)
    @unpack nfast = Par
    T, N, nsim = size(mchain)
    # Generate random index
    indfast=rand(1:nsim, nfast)
    indfast[1]=1 # Keeping the first element fixed.
    mchain_op=mchain[:,:,indfast]

    # First optimization BlackBox
    res = bboptimize(x->objfun_bbo_par_alt(x, mchain=mchain_op); SearchRange = (-10e300,10e300), NumDimensions = N, TraceMode=:silent)
    ## Getting results
    minbb = best_fitness(res)
    γ0 = best_candidate(res)
    # TS_bb = minbb*N
    # p_val_bb = ccdf(Chisq(N),TS_bb)
    println("First optimization from ", myid()," done.")
    ## Second optimization NLopt
    fop=nothing
    fop=Model(NLopt.Optimizer)
    set_optimizer_attribute(fop, "algorithm", :LN_BOBYQA)
    set_optimizer_attribute(fop, "xtol_rel", 1e-6)
    @variable(fop, gamma0[1:N])

    register(fop, :objfun_par_alt, N, (x...)->objfun_par_alt(x..., mchain=mchain_op); autodiff=true )
    @NLobjective(fop, Min, objfun_par_alt(gamma0...))

    for i=1:N
        set_start_value(gamma0[i],γ0[i])
    end

    ## Solving in NLopt
    JuMP.optimize!(fop)

    ## Getting Results
    # termination_status(fop)
    # primal_status(fop)

    γ1=value.(gamma0)
    obj_v=objective_value(fop)
    TS=T*obj_v
    p_val = ccdf(Chisq(N),TS)
    println("Second optimization from ", myid()," done.")
    ## Second optimization for refinement
    for i=1:N
        set_start_value(gamma0[i],γ1[i])
    end
    ## Solving in NLopt
    JuMP.optimize!(fop)

    ## Getting Results
    γ2=value.(gamma0)
    obj_v=objective_value(fop)
    TS=T*obj_v
    p_val = ccdf(Chisq(N),TS)
    println("Third optimization from ", myid()," done.")
    fop=nothing
    acc_r = MyModules.acc_rate(mchain_op, γ2,f=1)
    return TS, p_val, acc_r
end

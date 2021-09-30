using LinearAlgebra
using JuMP
using Ipopt
using Parameters

## Setting parameters

@with_kw struct Par
    burn::Int64=0 ;
    nsim::Int64=500000 ;
    nfast::Int64=20000 ;
    nboost::Int64=20000;
end

Param=Par()

## Road map

# 1) Warm start
# 2) Chain filling
# 3) Optimization

## Warm start
# min ||g(x,e)|| st. Inequalities

function warm_start(p,q)
    T,N=size(q)
    wmst=Model(Ipopt.Optimizer)
    set_optimizer_attribute(wmst, "print_level", 0)

    # Declaring variables
    @variable(wmst, d[1:T,1:N] >= 0) # δ_it Convex cost function c'(q)
    @variable(wmst, 1.0e-16 <= qs[1:T,1:N] <= 1.0) # unobserved true q_it

    # Non-negative markups
    # @constraint(wmst,mup[t=1:T,j=1:N],d[t,j]<=p[t])
    for j=1:N
        for t=1:T
            set_upper_bound(d[t,j],p[t])
        end
    end

    # Common ratio constraings (p_t-δ_it)/q_it=(p_t-δ_jt)/q_jt
    @NLconstraint(wmst,eq[t=1:T,j=2:N],(p[t]-d[t,1])/qs[t,1]-(p[t]-d[t,j])/qs[t,j]==0)

    # Co-monotone constraint (δ_it'-δ_it)(q_it'-q_it)>=0
    for j=1:N
        @constraint(wmst,[t=1:T,s=1:T],(d[t,j]-d[s,j])*(qs[t,j]-qs[s,j])>=0)
    end

    # Setting objective function
    @NLexpression(wmst, obj, sum((p[t]*(qs[t,j] - q[t,j]))^2 for t=1:T, j=1:N))
    @NLobjective(wmst,Min,obj)

    ## Solving the Model
    JuMP.optimize!(wmst)

    # Checking convergence
    # JuMP.primal_status(wmst)==MOI.FEASIBLE_POINT || error("Unfeasible point in warm start")
    status=JuMP.primal_status(wmst)==MOI.FEASIBLE_POINT

    # Saving variable values
    dsim=value.(d)
    qsim=value.(qs)

    # Memory cleaning
    wmst=nothing
    println("Warm star done.")
    return dsim, qsim, status
end

## Build chain

# Getting random direction vector from uniform distribution on the unit sphere

# Compute interval Α to draw uniformly to get α≥0
# Interval Α has to satisfy all constraints
# 1) Sign constraints
# 2) Markup constraints
# 3) Common ratio property
# 4) Co-monotone property

## jump function

function jump(p, q, eδ, eq; n::Int64=0)
    n < 5 || return q.-eq, eδ, eq
    T, N = size(eδ)
    # Getting random direction vector from uniform distribution on the unit sphere
    v=rand(2,T,N)
    v=v.*2 .-1 # comes from unif*(1-(-1))+(-1)=unif*(ub-lb)+lb
    ξ=v./norm(v) # ξ( δ ,t, i) ξ(q*,t,i)
    unif=rand(N)

    # Initializing matrices
    αmax=ones(N).*1.0e6
    αmin=ones(N).*-1.0e6
    α=zeros(N)

    # Sign constraints for δit
    for i=1:N, t=1:T
        if ξ[1,t,i]<0
            αmax[i]=αmax[i] < -(eδ[t,i]/ξ[1,t,i]) ? αmax[i] : -(eδ[t,i]/ξ[1,t,i])
            αmin[i]=αmin[i] > (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmin[i] : (p[t]-eδ[t,i]/ξ[1,t,i])
        else
            αmax[i]=αmax[i] < (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmax[i] : (p[t]-eδ[t,i]/ξ[1,t,i])
            αmin[i]=αmin[i] > -(eδ[t,i]/ξ[1,t,i]) ? αmin[i] : -(eδ[t,i]/ξ[1,t,i])
        end
    end

    # Sign constraints for q^*it
    for i=1:N, t=1:T
        if ξ[2,t,i]<0
            αmax[i]=αmax[i] < -(eq[t,i]/ξ[2,t,i]) ? αmax[i] : -(eq[t,i]/ξ[2,t,i])
        else
            αmin[i]=αmin[i] > -(eq[t,i]/ξ[2,t,i]) ? αmin[i] : -(eq[t,i]/ξ[2,t,i])
        end
    end

    # Co-monotone constraints
    for i=1:N, t=1:T, s=1:T
        t==s && continue
        num=eq[s,i]-eq[t,i]
        den=ξ[2,s,i]-ξ[2,t,i]

        if  den > 0
            # αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den)
            if num > 0
                αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den)
            else
                αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den) # If negative
            end
        else
            # αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den)
            if num < 0
                αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den) # if negative
            else
                αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den)
            end
        end

        num=eδ[s,i]-eδ[t,i]
        den=ξ[1,s,i]-ξ[1,t,i]

        if  den > 0
            # αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den)
            if num > 0
                αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den)
            else
                αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den) # if negative
            end
        else
            # αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den)
            if num < 0
                αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den) #if negative
            else
                αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den)
            end
        end
    end

    # Draw uniformily from Α
    for i=1:N
        α[i]=αmax[i]>αmin[i] ? unif[i]*(αmax[i]-αmin[i])+αmin[i] : 0.0
    end

    # Generate random δt ∀i in t=0
    δsim=eδ.+α'.*ξ[1,:,:]

    # Common ratio constraint
    qtemp=eq.+α'.*ξ[2,:,:]
    # Pick a country i at random for every t
    pick=rand(collect(1:N),T)
    # estimate ratio for every t
    r=[(p[t]-δsim[t,i])/qtemp[t,i] for (t,i) ∈ zip(collect(1:T), pick)]
    # Use Common ratio property to estimate q^* for every i and t
    qsim=(p.-δsim)./r

    # Estimate W
    w=q.-qsim

    ## Check inequalities
    check=[(δsim[t,i]-δsim[s,i])*(qsim[t,i]-qsim[s,i]) for t=1:T, s=1:T, i=1:N]
    # minimum(check)>=0 || error("Co-monotone inequalities are not satisfied!")
    minimum(check)>=0 || jump(p,q,eδ,eq,n=n+1)

    return w, δsim, qsim
end

function jump_algo(p, q, eδ, eq; n::Int64=0)
    n < 5 || return q.-eq, eδ, eq
    T, N = size(eδ)
    # Getting random direction vector from uniform distribution on the unit sphere
    v=rand(2,T,N)
    v=v.*2 .-1 # comes from unif*(1-(-1))+(-1)=unif*(ub-lb)+lb
    ξ=v./norm(v) # ξ( δ ,t, i) ξ(q*,t,i)
    unif=rand(N)

    # Initializing matrices
    αmax=ones(N).*1.0e6
    αmin=ones(N).*-1.0e6
    βmax=ones(N).*1.0e6
    βmin=ones(N).*-1.0e6
    α=zeros(N)
    β=zeros(N)

    # Sign constraints for δit
    for i=1:N, t=1:T
        if ξ[1,t,i]<0
            αmax[i]=αmax[i] < -(eδ[t,i]/ξ[1,t,i]) ? αmax[i] : -(eδ[t,i]/ξ[1,t,i])
            αmin[i]=αmin[i] > (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmin[i] : (p[t]-eδ[t,i]/ξ[1,t,i])
        else
            αmax[i]=αmax[i] < (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmax[i] : (p[t]-eδ[t,i]/ξ[1,t,i])
            αmin[i]=αmin[i] > -(eδ[t,i]/ξ[1,t,i]) ? αmin[i] : -(eδ[t,i]/ξ[1,t,i])
        end
    end

    # Sign constraints for q^*it
    for i=1:N, t=1:T
        if ξ[2,t,i]<0
            βmax[i]=βmax[i] < -(eq[t,i]/ξ[2,t,i]) ? βmax[i] : -(eq[t,i]/ξ[2,t,i])
        else
            βmin[i]=βmin[i] > -(eq[t,i]/ξ[2,t,i]) ? βmin[i] : -(eq[t,i]/ξ[2,t,i])
        end
    end

    # Co-monotone constraints
    for i=1:N, t=1:T, s=1:T
        t==s && continue
        # For q
        num=eq[s,i]-eq[t,i]
        den=ξ[2,s,i]-ξ[2,t,i]

        if  den > 0
            # αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den)
            if num > 0
                βmin[i]=βmin[i] > -(num/den) ? βmin[i] : -(num/den)
            else
                βmax[i]=βmax[i] < -(num/den) ? βmax[i] : -(num/den) # If negative
            end
        else
            # βmax[i]=βmax[i] < -(num/den) ? βmax[i] : -(num/den)
            if num < 0
                βmin[i]=βmin[i] > -(num/den) ? βmin[i] : -(num/den) # if negative
            else
                βmax[i]=βmax[i] < -(num/den) ? βmax[i] : -(num/den)
            end
        end

        num=eδ[s,i]-eδ[t,i]
        den=ξ[1,s,i]-ξ[1,t,i]

        if  den > 0
            # αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den)
            if num > 0
                αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den)
            else
                αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den) # if negative
            end
        else
            # αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den)
            if num < 0
                αmin[i]=αmin[i] > -(num/den) ? αmin[i] : -(num/den) #if negative
            else
                αmax[i]=αmax[i] < -(num/den) ? αmax[i] : -(num/den)
            end
        end
    end

    # Draw uniformily from Α
    for i=1:N
        α[i]=αmax[i]>αmin[i] ? unif[i]*(αmax[i]-αmin[i])+αmin[i] : 0.0
    end

    # Generate random δt ∀i in t=0
    δsim=eδ.+α'.*ξ[1,:,:]

    # Common ratio constraint
    # qsim=zeros(size(eq))
    # for t=1:T, s=1:T
    #     t==s && continue
    #     pick = rand(1:N) # pick a country at random
    #
    #     if δsim[s,pick]-δsim[t,pick] >= 0
    #         # t>=2 & (minimum(qsim[s,:]-qsim[t,:]) >= 0) && continue
    #         diff = 1.0-eq[t,pick]
    #         qs= diff > 0 ? rand()*(diff)+eq[t,pick] : eq[t,pick]
    #         r=(p[s]-δsim[s,pick])/qs
    #         qsim[s,:]=(p[s].-δsim[s,:])./r
    #     else
    #         # t>=2 & (maximum(qsim[s,:]-qsim[t,:]) <= 0) && continue
    #         diff = eq[t,pick]-1e-8
    #         qs = diff > 0 ? rand()*diff+0.0 : eq[t,pick]
    #         r=(p[s]-δsim[s,pick])/qs
    #         qsim[s,:]=(p[s].-δsim[s,:])./r
    #     end
    # end

    # Draw uniformily from B
    for i=1:N
        β[i]=βmax[i]>βmin[i] ? rand()*(βmax[i]-βmin[i])+βmin[i] : 0.0
    end

    # Generate random q^*it ∀i in t=0
    qtemp=eq.+β'.*ξ[2,:,:]
    # Pick a country i at random for every t
    pick=rand(collect(1:N),T)

    # estimate ratio for every t
    r=[(p[t]-δsim[t,i])/qtemp[t,i] for (t,i) ∈ zip(collect(1:T), pick)]
    # Use Common ratio property to estimate q^* for every i and t
    qsim=(p.-δsim)./r

    # Estimate W
    w=q.-qsim

    ## Check inequalities
    check=[(δsim[t,i]-δsim[s,i])*(qsim[t,i]-qsim[s,i]) for t=1:T, s=1:T, i=1:N]
    # minimum(check)>=0 || error("Co-monotone inequalities are not satisfied!")
    minimum(check)>=0 || jump(p,q,eδ,eq,n=n+1)

    return w, δsim, qsim
end


## Fill in MCMC Markov Chain Monte Carlo integration

function gchain(p, q; Par::Par=Param)
    @unpack burn, nsim = Par
    T, N = size(q)
    eδ, eq, status = warm_start(p,q)
    ew = eq-q
    mchain=zeros(N,(nsim-burn))
    status || return mchain, status # If Unfeasible solution in warmstart don't do loop
    for r=burn+1:nsim
        W, δsim, qsim = jump(p,q,eδ,eq)
        logtrydens = norm(p.*W)-norm(p.*ew)

        if logtrydens.>log.(rand())
            ew=W
            eδ=δsim
            eq=qsim
        end

        mchain[:,r-burn]=gM(p,ew)
    end
    println("Chain is ready.")
    return mchain, status
end

function gchain_algo(p, q; Par::Par=Param)
    @unpack burn, nsim = Par
    T, N = size(q)
    eδ, eq, status = warm_start(p,q)
    ew = eq-q
    mchain=zeros(N,(nsim-burn))
    status || return mchain, status # If Unfeasible solution in warmstart don't do loop
    for r=burn+1:nsim
        W, δsim, qsim = jump_algo(p,q,eδ,eq)
        logtrydens = norm(p.*W)-norm(p.*ew)

        if logtrydens.>log.(rand())
            ew=W
            eδ=δsim
            eq=qsim
        end

        mchain[:,r-burn]=gM(p,ew)
    end
    println("Chain is ready.")
    return mchain, status
end

function gchain_alt(p, q; Par::Par=Param)
    @unpack burn, nsim = Par
    T, N = size(q)
    eδ, eq, status = warm_start(p,q)
    ew = eq-q
    mchain=zeros(N,(nsim-burn))
    status || return mchain, status # If Unfeasible solution in warmstart don't do loop
    for r=burn+1:nsim
        W, δsim, qsim = jump(p,q,eδ,eq)
        logtrydens = norm(p.*W)-norm(p.*ew)

        if rand() < 1.4^(logtrydens[1])
            ew=W
            eδ=δsim
            eq=qsim
        end

        mchain[:,r-burn]=gM(p,ew)
    end
    println("Chain is ready.")
    return mchain, status
end

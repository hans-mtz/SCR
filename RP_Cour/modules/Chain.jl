using LinearAlgebra
# using Convex, ECOS, LinearAlgebra
using JuMP
using Ipopt
using Parameters

# using NLopt
## Setting parameters

@with_kw struct Par
    burn::Int64=0 ;
    nsim::Int64=500000 ;
    nfast::Int64=20000 ;
end

Param=Par()

# 1) Warm start ✓
# 2) Chain filling
# 3) Optimization

## Gurobi set up
# ENV["GUROBI_HOME"] = "/Library/gurobi911/mac64"
# ENV["GRB_LICENSE_FILE"]="/Library/gurobi911/mac64/gurobi/gurobi.lic"




## Warm start
# min ||g(x,e)|| st. Inequalities

function warm_start(p,q)
    T,N=size(q)
    wmst=Model(Ipopt.Optimizer)
    set_optimizer_attribute(wmst, "print_level", 0)
    # set_optimizer_attribute(wmst, "algorithm", :LD_SLSQP)
    # Declaring variables
    @variable(wmst, d[1:T,1:N] >= 0) # δ_it Convex cost function c'(q)
    @variable(wmst, 1.0e-16 <= qs[1:T,1:N] <= 1) # unobserved true q_it

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
    # println("Primal status = ", JuMP.primal_status(wmst))
    JuMP.primal_status(wmst)==MOI.FEASIBLE_POINT || error("Unfeasible point in warm start")
    # Saving variable values
    dsim=value.(d)
    qsim=value.(qs)

    # Memory cleaning
    wmst=nothing
    println("Warm star done.")
    return dsim, qsim
end

## test function

# eδ, eq = warm_start(p=p[1:12],q=q[1:12,1:6])

# it works!

## Build chain

# Getting random direction vector from uniform distribution on the unit sphere

# Compute interval Α to draw uniformly to get α≥0
# Interval Α has to satisfy all constraints
# 1) Sign constraints
# 2) Markup constraints
# 3) Common ratio property
# 4) Co-monotone property


# PLAN B; algorithm
# For time t=0
# 1) Generate random δit ∀i ∈ I at time t=1
# 2) Get random q^*it for a random i and ratio  rt=(pt-δit)/q^*it
# 3) Get q^*jt ∀j≠i (pt-δjt)/rt=q^*jt
# 4) Get ranking
# For time t+1
# 1) Generate random δit sequentially, according to ranking at t=0
# 2) Generate random q^*it for a random i
# 3) get q^*jt ∀j≠i


# For time t=0
# 1) Generate random δit ∀i ∈ I at time t=1

# 1,1) compute interval Α box constraints
function jump_init(p, eδ, eq)
    T, N = size(eδ)
    # Getting random direction vector from uniform distribution on the unit sphere
    v=rand(2,T,N).*2 .-1 # comes from unif*(1-(-1))+(-1)=unif*(ub-lb)+lb
    ξ=v./norm(v) # ξ( δ ,t, i) ξ(q*,t,i)

    αmax=ones(T,N).*1.0e6
    αmin=ones(T,N).*-1.0e6
    α=zeros(N)
    unif=rand(N)

    # Sign constraints for δit
    for i=1:N, t=1:T
        if ξ[1,t,i]<0
            αmax[1,i]=αmax[1,i] < -(eδ[t,i]/ξ[1,t,i]) ? αmax[1,i] : -(eδ[t,i]/ξ[1,t,i])
            αmin[1,i]=αmin[1,i] > (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmin[1,i] : (p[t]-eδ[t,i]/ξ[1,t,i])
        else
            αmax[1,i]=αmax[1,i] < (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmax[1,i] : (p[t]-eδ[t,i]/ξ[1,t,i])
            αmin[1,i]=αmin[1,i] > -(eδ[t,i]/ξ[1,t,i]) ? αmin[1,i] : -(eδ[t,i]/ξ[1,t,i])
        end
    end

    # Sign constraints for q*it
    for i=1:N, t=1:T
        if ξ[2,t,i]<0
            αmax[1,i]=αmax[1,i] < -(eq[t,i]/ξ[2,t,i]) ? αmax[1,i] : -(eq[t,i]/ξ[2,t,i])
        else
            αmin[1,i]=αmin[1,i] > -(eq[t,i]/ξ[2,t,i]) ? αmin[1,i] : -(eq[t,i]/ξ[2,t,i])
        end
    end

    # Co-monotone constraints
    for i=1:N, t=1:T, s=1:T
        t==s && continue
        num=eq[s,i]-eq[t,i]
        den=ξ[2,s,i]-ξ[2,t,i]

        if  den > 0
            αmin[1,i]=αmin[1,i] > -(num/den) ? αmin[1,i] : -(num/den)
        else
            αmax[1,i]=αmax[1,i] < -(num/den) ? αmax[1,i] : -(num/den)
        end

        num=eδ[s,i]-eδ[t,i]
        den=ξ[1,s,i]-ξ[1,t,i]

        if  den > 0
            αmin[1,i]=αmin[1,i] > -(num/den) ? αmin[1,i] : -(num/den)
        else
            αmax[1,i]=αmax[1,i] < -(num/den) ? αmax[1,i] : -(num/den)
        end
    end

    # Draw uniformily from Α
    for i=1:N
        α[i]=αmax[1,i]>αmin[1,i] ? unif[1,i]*(αmax[1,i]-αmin[1,i])+αmin[1,i] : 0.0
    end

    # Generate random δt ∀i in t=0
    δsim=eδ[1,:].+α.*ξ[1,1,:]

    # Common ratio constraint
    qtemp=eq[1,:].+α.*ξ[2,1,:]
    # Pick a country i at random for
    pick=rand(collect(1:N))
    r=(p[1]-δsim[pick])/qtemp[pick]

    qsim=(p[1].-δsim)./r

    ## Check inequalities
    # check=[(δsim[t,i]-δsim[s,i])*(qsim[t,i]-qsim[s,i]) for t=1:T, s=1:T, i=1:N]
    # minimum(check)>=0 || error("Init: Co-monotone inequalities are not satisfied!")
    #

    return δsim, qsim, αmax, αmin, ξ
end
#
# ## Testing get_delta0 function
#
# δ0, q0=get_delta0()
#
#
# (p[1].-δ0)./q0 # Worked!
## jump function

@everywhere function jump(;p=p[1:12], q=q[1:12,1:6], eδ=eδ, eq=eq)
    T, N = size(eδ)
    # Getting random direction vector from uniform distribution on the unit sphere
    v=rand(2,T,N)
    v=v.*2 .-1 # comes from unif*(1-(-1))+(-1)=unif*(ub-lb)+lb
    ξ=v./norm(v) # ξ( δ ,t, i) ξ(q*,t,i)

    αmax=ones(T,N).*1.0e6
    αmin=ones(T,N).*-1.0e6
    α=zeros(N)
    unif=rand(N)

    # Sign constraints for δit
    for i=1:N, t=1:T
        if ξ[1,t,i]<0
            αmax[1,i]=αmax[1,i] < -(eδ[t,i]/ξ[1,t,i]) ? αmax[1,i] : -(eδ[t,i]/ξ[1,t,i])
            αmin[1,i]=αmin[1,i] > (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmin[1,i] : (p[t]-eδ[t,i]/ξ[1,t,i])
        else
            αmax[1,i]=αmax[1,i] < (p[t]-eδ[t,i]/ξ[1,t,i]) ? αmax[1,i] : (p[t]-eδ[t,i]/ξ[1,t,i])
            αmin[1,i]=αmin[1,i] > -(eδ[t,i]/ξ[1,t,i]) ? αmin[1,i] : -(eδ[t,i]/ξ[1,t,i])
        end
    end

    # Sign constraints for q*it
    for i=1:N, t=1:T
        if ξ[2,t,i]<0
            αmax[1,i]=αmax[1,i] < -(eq[t,i]/ξ[2,t,i]) ? αmax[1,i] : -(eq[t,i]/ξ[2,t,i])
        else
            αmin[1,i]=αmin[1,i] > -(eq[t,i]/ξ[2,t,i]) ? αmin[1,i] : -(eq[t,i]/ξ[2,t,i])
        end
    end

    # Co-monotone constraints
    for i=1:N, t=1:T, s=1:T
        t==s && continue
        num=eq[s,i]-eq[t,i]
        den=ξ[2,s,i]-ξ[2,t,i]

        if  den > 0
            αmin[1,i]=αmin[1,i] > -(num/den) ? αmin[1,i] : -(num/den)
        else
            αmax[1,i]=αmax[1,i] < -(num/den) ? αmax[1,i] : -(num/den)
        end

        num=eδ[s,i]-eδ[t,i]
        den=ξ[1,s,i]-ξ[1,t,i]

        if  den > 0
            αmin[1,i]=αmin[1,i] > -(num/den) ? αmin[1,i] : -(num/den)
        else
            αmax[1,i]=αmax[1,i] < -(num/den) ? αmax[1,i] : -(num/den)
        end
    end

    # Draw uniformily from Α
    for i=1:N
        α[i]=αmax[1,i]>αmin[1,i] ? unif[1,i]*(αmax[1,i]-αmin[1,i])+αmin[1,i] : 0.0
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

    ## Check inequalities
    check=[(δsim[t,i]-δsim[s,i])*(qsim[t,i]-qsim[s,i]) for t=1:T, s=1:T, i=1:N]
    minimum(check)>=0 || error("Init: Co-monotone inequalities are not satisfied!")

    # Estimate W
    w=q.-qsim

    return w, δsim, qsim
end

## testing

# w, δ, qs = jump()
# w1, δ1, qs1 = jump()
#
# w.==w1

##

function ejump(p, q, eδ, eq)
    T, N = size(eδ)
    δ0, q0, αmax, αmin, ξ = jump_init(p, eδ, eq)
    # Initializing arrays
    δsim=zeros(T,N)
    qsim=zeros(T,N)
    α=zeros(N)
    unif=rand(T,N)
    δsim[1,:]=δ0
    qsim[1,:]=q0

    # 4) Get ranking
    rank=indexin(sort(δ0, rev=true),δ0)
    # For time t+1
    # 1) Generate random δit sequentially, according to ranking at t=0

    # random draw for δ1
    b=rank[1] # pick the best firm/country from t=0
    α[1]=αmax[1,b]>αmin[1,b] ? unif[1,b]*(αmax[1,b]-αmin[1,b])+αmin[1,b] : 0.0
    δsim[:,b]=eδ[:,b]+α[1]*ξ[1,:,b]


    for h=2:T # For t=1,...,T
        for i=2:N
            j, k = rank[i-1], rank[i]

            for t=1:T
                # Sign constraints for δit
                if ξ[1,t,k]<0
                    αmax[h,k]=αmax[h,k] < (δsim[h,j]-eδ[t,k])/ξ[1,t,k] ? αmax[h,k] : (δsim[h,j]-eδ[t,k])/ξ[1,t,k]
                    αmin[h,k]=αmin[h,k] > (p[t]-eδ[t,k]/ξ[1,t,k]) ? αmin[h,k] : (p[t]-eδ[t,k]/ξ[1,t,k])
                else
                    αmax[h,k]=αmax[h,k] < (p[t]-eδ[t,k]/ξ[1,t,k]) ? αmax[h,k] : (p[t]-eδ[t,k]/ξ[1,t,k])
                    αmin[h,k]=αmin[h,k] > (δsim[h,j]-eδ[t,k])/ξ[1,t,k] ? αmin[h,k] : (δsim[h,j]-eδ[t,k])/ξ[1,t,k]
                end

                # Sign constraints for q^*it
                if ξ[2,t,k]<0
                    αmax[h,k]=αmax[1,k] < -(eq[t,k]/ξ[2,t,k]) ? αmax[1,k] : -(eq[t,k]/ξ[2,t,k])
                    # αmin[h,k]=αmin[h,k] > (qsim[t,j]-eq[t,k])/ξ[2,t,k] ? αmin[h,k] : (qsim[t,j]-eq[t,k])/ξ[2,t,k]
                else
                    αmin[h,k]=αmin[1,k] > -(eq[t,k]/ξ[2,t,k]) ? αmin[1,k] : -(eq[t,k]/ξ[2,t,k])
                    # αmax[h,k]=αmax[h,k] < (qsim[t,j]-eq[t,k])/ξ[2,t,k] ? αmax[h,k] : (qsim[t,j]-eq[t,k])/ξ[2,t,k]
                end

                # Co-monotone constraints
                for s=1:T
                    t==s && continue
                    # q
                    num=eq[s,k]-eq[t,k]
                    den=ξ[2,s,k]-ξ[2,t,k]

                    if  den > 0
                        αmin[h,k]=αmin[h,k] > -(num/den) ? αmin[h,k] : -(num/den)
                    else
                        αmax[h,k]=αmax[h,k] < -(num/den) ? αmax[h,k] : -(num/den)
                    end

                    # δ
                    num=eδ[s,k]-eδ[t,k]
                    den=ξ[1,s,k]-ξ[1,t,k]

                    if  den > 0
                        αmin[h,k]=αmin[h,k] > -(num/den) ? αmin[h,k] : -(num/den)
                    else
                        αmax[h,k]=αmax[h,k] < -(num/den) ? αmax[h,k] : -(num/den)
                    end
                end
            end

            # Draw uniformily from Α
            α[k]=αmax[h,k]>αmin[h,k] ? unif[h,k]*(αmax[h,k]-αmin[h,k])+αmin[h,k] : 0.0
            # Generate random δt ∀i in t
            δsim[h,k]=eq[h,k]+α[k]*ξ[1,h,k]
        end

        # Common ratio constraint
        # 2) Generate random q^*it for a random i
        qtemp=eq[h,:].+α.*ξ[2,h,:]
        # Pick a country i at random to generate ratio
        pick=rand(collect(1:N))
        r=(p[h]-δsim[h,pick])/qtemp[pick]
        # 3) get q^*jt ∀j≠i
        # Generate q^*jt for j≠i
        qsim[h,:]=(p[h].-δsim[h,:])./r
    end

    ## Check inequalities
    check=[(δsim[t,i]-δsim[s,i])*(qsim[t,i]-qsim[s,i]) for t=1:T, s=1:T, i=1:N]
    minimum(check)>=0 || error("Co-monotone inequalities are not satisfied!")

    W=q-qsim

    return  W, δsim, qsim
end

## Test

# W, dsim, qsim = get_deltat()

# eδ, eq = warm_start(p[1:12],q[1:12,1:6])

## Fill in MCMC Markov Chain Monte Carlo integration

function gchain(p, q; Par::Par=Param, eδ=eδ, eq=eq)
    @unpack burn, nsim = Par
    T, N = size(q)
    # eδ, eq = warm_start()
    ew=eq-q
    mchain=zeros((nsim-burn),N)
    for r=burn+1:nsim
        # W, δsim, qsim = ejump(p, q, eδ, eq)
        W, δsim, qsim = jump(p=p, q=q, eδ=eδ, eq=eq)
        logtrydens = norm(p.*W)-norm(p.*ew)

        if logtrydens.>log.(rand())
            ew=W
            eδ=δsim
            eq=qsim
        end

        mchain[r-burn,:]=gM(p,ew)

    end
    println("Chain is ready.")
    return mchain
end

## Testing chain

# mchain = gchain(p[1:12],q[1:12,1:6])

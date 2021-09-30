using Distributed
@everywhere using SharedArrays, LinearAlgebra

addprocs(4)

## Markov Chain
@everywhere function myrange_burn(q::SharedArray, burn)
    idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
       return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in range(burn, stop=size(q,1), length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end

@everywhere function fill_chain!(m::SharedArray, p, q, eδ, eq, burn)
    ew=eq-q
    range = myrange_burn(m, burn)
    for r in range
        W, δsim, qsim = jump(p=p, q=q, eδ=eδ, eq=eq)
        logtrydens = norm(p.*W)-norm(p.*ew)

        if logtrydens.>log.(rand())
            ew=W
            eδ=δsim
            eq=qsim
        end

        m[r,:]=gM(p,ew)
    end
    m
end

function shared_fill!(m, p, q, eδ, eq, burn)
    @sync begin
        for w in procs(m)
            @async remotecall_wait(fill_chain!, w, m, p, q, eδ, eq, burn)
        end
    end
    m
end;

function gchain_par(p, q; Par::Par=Param, eδ=eδ, eq=eq)
    @unpack burn, nsim = Par
    T, N = size(q)
    # eδ, eq = warm_start()

    ew=eq-q

    mchain=zeros((nsim-burn),N)
    mchain_par=SharedArray{Float64,2}(mchain)
    shared_fill!(mchain_par, p, q, eδ, eq, burn)

    println("Chain is ready.")
    return mchain_par
end

## Testin Parallel
@fetchfrom 2 myrange_burn(mchain_shared, 1000)
@fetchfrom 3 myrange_burn(mchain_shared, 1000)
mchain_par = gchain_par(p[1:12],q[1:12,1:6])

@time gchain(p[1:12],q[1:12,1:6])
@time gchain_par(p[1:12],q[1:12,1:6])

@fetchfrom 2 fill_chain!(mchain_shared,p[1:12],q[1:12,1:6], $eδ, $eq, 1000 )

## Optmization

@everywhere function myrange(q::SharedArray)
    idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
       return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in range(0, stop=size(q,1), length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end

@everywhere function hM_chunk!(hM, mchain, gamma)
   # @show (irange, jrange, trange)  # display so we can see what's happening

   # hM=zeros(size(mchain,2))
   T, N, nfast = size(mchain)
   geta=mchain[:,:,1]
   range = myrange(hM)
   for t in range, r=1:nfast
       val=0.0
       gtry=mchain[t,:,r]
       for j=1:N
           val=+gtry[j].*gamma[j]-geta[t,j].*gamma[j]
       end
       geta[t,:]=log.(rand()) < val[1] ? gtry : geta[t,:]
       hM[t,:]=+geta[t,:]./nfast
       # return hM
   end
   hM
end

function shared_hM_fill!(hM, m, g)
    @sync begin
        for w in procs(hM)
            @async remotecall_wait(hM_chunk!, w, hM, m, g)
        end
    end
    hM
end;


@everywhere function shared_hM_chunck!(hM, m, g)
    hM = hM_chunk!(hM, m, g, myrange(hM))
    return hM
end

function shared_dist_hM!(m, gamma) # display so we can see what's happening
    # hM = zeros(size(m,2))
    hM = @distributed (+) for p in workers()
        remotecall_fetch(shared_hM_chunck!, p, m, gamma)
    end
    return hM
end;


function objfun_par(gamma0::Vector...; Par::Par=Param, mchain=mchain_op)
    @unpack nsim, nfast = Par
    @everywhere gamma = $gamma0
    N = size(mchain[1,:])

    mchain_shared = SharedArray{Float64,2}(mchain)
    hM = shared_dist_hM!(mchain_shared, gamma)

    hM=hM./nfast
    mhM=mean(hM)
    var=Statistics.var(hM)

    obj=(mhM^2)/var

    return obj
end

##
@time objfun_par(collect(1:6))
@time objfun(collect(1:6))
@time objfun_bbo(collect(1:6))

@fetchfrom 2 myrange(mchain_shared)

@fetchfrom 4 myrange(mchain_shared)

myrange(mchain_shared)

## parallel obj funs

@everywhere function myrange(q::SharedArray)
    idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
       return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in range(0, stop=size(q,1), length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end

@everywhere function hM_chunk!(hM, mchain, gamma)
   # @show (irange, jrange, trange)  # display so we can see what's happening

   # hM=zeros(size(mchain,2))
   T, N, nfast = size(mchain)
   geta=mchain[:,:,1]
   range = myrange(hM)
   for t in range, r=1:nfast
       val=0.0
       gtry=mchain[t,:,r]
       for j=1:N
           val=+gtry[j].*gamma[j]-geta[t,j].*gamma[j]
       end
       geta[t,:]=log.(rand()) < val[1] ? gtry : geta[t,:]
       hM[t,:]=+geta[t,:]./nfast
       # return hM
   end
   hM
end

function shared_hM_fill!(hM, m, g)
    @sync begin
        for w in procs(hM)
            @async remotecall_wait(hM_chunk!, w, hM, m, g)
        end
    end
    hM
end;

function objfun_par(gamma0::Vector ...; mchain=mchain_op)
    # @unpack nsim, nfast = Par
    gamma = gamma0[1]
    T, N, nfast = size(mchain)

    hM = SharedArray{Float64,2}(zeros(T,N))
    shared_hM_fill!(hM, mchain, gamma)

    mhM=sum(hM, dims=1)./T
    omega=(hM.-mhM)'*(hM.-mhM)

    obj=mhM*inv(omega)*mhM'

    return obj[1]
end

function objfun_bbo_par(gamma0; mchain=mchain_op)

    gamma = gamma0
    T,N,nfast = size(mchain)

    hM = SharedArray{Float64,2}(zeros(T,N))
    shared_hM_fill!(hM, mchain, gamma)

    mhM=sum(hM, dims=1)./T
    omega=(hM.-mhM)'*(hM.-mhM)

    obj=mhM*inv(omega)*mhM'

    return obj[1]
end

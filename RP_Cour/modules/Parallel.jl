@everywhere using SharedArrays

## fill table in parallel using shared arrays

@everywhere function myrange_y(t::SharedArray,p)
    idx = indexpids(t)
    if idx == 0 # This worker is not assigned a piece
       return 1:0, 1:0
    end
    nchunks = length(procs(t))
    splits = [round(Int, s) for s in range(0, stop=size(t,1), length=nchunks+1)]
    years = [round(Int, y) for y in range(0, stop=size(p,1), length=nchunks+1)]
    splits[idx]+1:splits[idx+1], years[idx]+1:12:years[idx+1]
end

@everywhere function fill_table!(t::SharedArray, pt, qt)
    split, years = myrange_y(t,pt)
    for (i,y) in zip(split, years)
        t[i,1] == 0.0 || continue # don't do previously done evals
        @show (myid(),i ,y)
        p = pt[y:y+11]
        q = qt[y:y+11,:]
        mchain, status = gchain(p,q)
        status || continue # don't do if infeasible point
        t[i,1], t[i,2] = RP_Cournot_test(mchain)
    end
    t
end

function shared_fill_table!(t, p, q)
    @sync begin
        for w in procs(t)
            @async remotecall_wait(fill_table!, w, t, p, q)
        end
    end
    t
end;

## parallel chain

# Semesters
@everywhere function myrange_sem(t::SharedArray,p)
    idx = indexpids(t)
    if idx == 0 # This worker is not assigned a piece
       return 1:0, 1:0
    end
    nchunks = length(procs(t))
    splits = [round(Int, s) for s in range(0, stop=size(t,1), length=nchunks+1)]
    years = [round(Int, y) for y in range(0, stop=size(p,1), length=nchunks+1)]
    splits[idx]+1:splits[idx+1], years[idx]+1:6:years[idx+1]
end

@everywhere function fill_chain_sem!(t::SharedArray, pt, qt)
    split, years = myrange_sem(t,pt)
    for (i,y) in zip(split, years)
        t[i,1,1] == 0.0 || continue # don't do previously done evals
        @show (myid(),i ,y)
        p = pt[y:y+5]
        q = qt[y:y+5,:]
        t[i,:,:], status = gchain(p,q)
        # status || continue # don't do if infeasible point
        # t[i,1], t[i,2] = RP_Cournot_test(mchain)
    end
    t
end

function shared_fill_chain!(t, p, q,f::Function)
    @sync begin
        for w in procs(t)
            @async remotecall_wait(f, w, t, p, q)
        end
    end
    t
end;

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

# @everywhere function hM_dist_chunk!(hM, mchain, gamma)
# Tried to use @distributed but it actually took more time
#    T, N, nfast = size(mchain)
#    geta=mchain[:,:,1]
#    range = myrange(hM)
#    for t in range
#        hM[t,:] = @distributed (+) for r=1:nfast
#            val=0.0
#            gtry=mchain[t,:,r]
#            for j=1:N
#                val=+gtry[j].*gamma[j]-geta[t,j].*gamma[j]
#            end
#            geta[t,:]=log.(rand()) < val[1] ? gtry : geta[t,:]
#            hM[t,:]=+geta[t,:]
#        end
#    end
#    hM./nfast
# end

@everywhere function shared_hM_fill!(hM, m, g, fun::Function)
    @sync begin
        for w in procs(hM)
            @async remotecall_wait(fun, w, hM, m, g)
        end
    end
    hM
end;

##

@everywhere function objfun_par(gamma0...; mchain=mchain_op)
    # @unpack nsim, nfast = Par
    gamma = gamma0
    T, N, nfast = size(mchain)

    hM = SharedArray{Float64,2}(zeros(T,N))
    shared_hM_fill!(hM, mchain, gamma, hM_chunk!)

    mhM=sum(hM, dims=1)./T
    omega=(hM.-mhM)'*(hM.-mhM)

    obj=mhM*inv(omega)*mhM'

    return obj[1]
end

@everywhere function objfun_bbo_par(gamma0; mchain=mchain_op)

    gamma = gamma0
    T,N,nfast = size(mchain)

    hM = SharedArray{Float64,2}(zeros(T,N))
    shared_hM_fill!(hM, mchain, gamma, hM_chunk!)

    mhM=sum(hM, dims=1)./T
    omega=(hM.-mhM)'*(hM.-mhM)

    obj=mhM*inv(omega)*mhM'

    return obj[1]
end

## fill table in parallel using shared arrays
# Use mchain by semesters, test by decade

@everywhere function myrange_dec(t::SharedArray)
    idx = indexpids(t)
    if idx == 0 # This worker is not assigned a piece
       return 1:0, 1:0
    end
    nchunks = length(procs(t))
    splits = collect(1:size(t,1))
    dec = [0, 16, 36, 56, 72]
    if idx ≤ 4
        splits[idx], dec[idx]+1:dec[idx+1]
    else
        1:0, 1:0
    end
end

@everywhere function fill_table_dec!(t::SharedArray, m)
    i, split = myrange_dec(t)
    @show (myid(), i, split)
    mchain = m[split, :, :]
    t[i,1], t[i,2] = RP_Cournot_test_par(mchain)
    return t
end

function shared_fill_table_dec!(t, m)
    @sync begin
        for w in procs(t)
            @async remotecall_wait(fill_table_dec!, w, t, m)
        end
    end
    t
end;

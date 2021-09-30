#my centering fun
using DataFramesMeta

## My function

ld=stack(data, Between(:algeria, :uk))
#Getting market shares
@transform!(ld,q=:value ./ :totalworld)
#Transforming data into wide form to get market shares only
# df=unstack(ld, [:year,:variable], :month, :q, allowduplicates=true)
df=unstack(ld, [:year,:month], :variable, :q, allowduplicates=true)

sum(q[:,1,:]'.*p,dims=2)

g=ones(36,19)

pq=p.*q
j=1
for i=1:12:size(pq,1)-1
    g[j,:]=sum(pq[i:i+11,:], dims=1)
    j+=1
end

g=nothing

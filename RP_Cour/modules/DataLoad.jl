using CSV, DataFrames, DataFramesMeta


## Loading data
datadir="/Volumes/SSD Hans/Github/Research_ideas/Testing Market Power/data/Carvajal et al 2013/CDFQ_ECMA_Files/oildata.csv"

data=CSV.read(datadir, DataFrame)
data[:, Between(:algeria, :uk)]
data.id=1:size(data, 1)
ld=stack(data, Between(:algeria, :uk))
country=names(data, Between(:algeria, :uk))
@transform!(ld,q=:value ./ :totalworld)

df=unstack(ld, [:year,:variable], :month, :q, allowduplicates=true)

p=unstack(data[data.year.<2009,:], :year, :month, :oilpricereal2009)

# p array: 36 years x 12 months
p=Array(p[:,Not(:year)])

df[df.year.<2009, Not([:year, :variable])]

# q array, 12 months x 19 countries x 36 years
q=reshape(Array(ld[ld.year.<2009,:q]),12,19,36)

q[:,1,1].==@linq ld |>
    where(:year.==1973, :variable.=="algeria") |>
    select(:q)

@linq ld |>
        where(:year.==1973, :variable.=="algeria") |>
        select(:q)

q[1,1,:]

p=nothing
q=nothing

##
p=@linq data |>
    where(:year.<2009) |>
    select(:oilpricereal2009)


q = @linq data |>
        where(:year.<2009)

q = q[:,Between(:algeria, :uk)]

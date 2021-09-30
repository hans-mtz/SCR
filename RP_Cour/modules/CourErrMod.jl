# module CourErr
#
# export load_data, gM

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra
using Parameters

## Loading Data

function load_data(file)

    # Reading data
    data=CSV.read(file, DataFrame)

    # # Transforming data to wide form to get Prices
    # p=unstack(data[data.year.<2009,:], :year, :month, :oilpricereal2009)
    # # Price array p: 36 years x 12 months
    # p=Array(p[:,Not(:year)]./100)

    p = @linq data |>
         where(:year.<2009) |>
         select(:oilpricereal2009)
    p=p./100 #
    #Transforming data into long form to get Quantities
    ld=stack(data, Between(:algeria, :uk))
    #Getting market shares
    @transform!(ld,q=:value ./ :totalworld)
    #Transforming data into wide form to get market shares only
    # df=unstack(ld, [:year,:variable], :month, :q, allowduplicates=true)
    df=unstack(ld, [:year,:month], :variable, :q, allowduplicates=true)
    # Quantities array, q: 12 months x 19 countries x 36 years
    # q=reshape(Array(ld[ld.year.<2009,:q]),M,N,T)

    q = @linq df |>
            where(:year.<2009)

    q = q[:,Between(:algeria, :uk)]

    println("Data is ready.")
    return Array(p), Array(q)
end

function load_data_select(file;select=[:saudiarabia, :unitedstates, :iran, :china, :mexico, :venezuela])

    # Reading data
    data=CSV.read(file, DataFrame)

    # # Transforming data to wide form to get Prices
    # p=unstack(data[data.year.<2009,:], :year, :month, :oilpricereal2009)
    # # Price array p: 36 years x 12 months
    # p=Array(p[:,Not(:year)]./100)

    p = @linq data |>
         where(:year.<2009) |>
         select(:oilpricereal2009)
    p=p./100 #
    #Transforming data into long form to get Quantities
    ld=stack(data, Between(:algeria, :uk))
    #Getting market shares
    @transform!(ld,q=:value ./ :totalworld)
    #Transforming data into wide form to get market shares only
    # df=unstack(ld, [:year,:variable], :month, :q, allowduplicates=true)
    df=unstack(ld, [:year,:month], :variable, :q, allowduplicates=true)
    # Quantities array, q: 12 months x 19 countries x 36 years
    # q=reshape(Array(ld[ld.year.<2009,:q]),M,N,T)

    q = @linq df |>
            where(:year.<2009)

    q = q[:,select]

    println("Data is ready.")
    return Array(p), Array(q)
end

## Defining centering moment condition

@everywhere function gM(p, W)
    T, N = size(W)
    # g=ones(T, N)
    g=ones(N)
    # pq=p.*W
    # j=1

    # for i=1:12:size(pq,1)
    #     g[j,:]=sum(pq[i:i+11,:], dims=1)
    #     j+=1
    # end

    for i=1:N
        g[i]=p⋅W[:,i]
    end

    return g
end

## include warm start file
# include("warm_s.jl")


# end  # module CourErr

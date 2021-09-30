
using CSV, DataFrames, DataFramesMeta
using LinearAlgebra
using Parameters

## Loading Data

function load_data(file)

    # Reading data
    data=CSV.read(file, DataFrame)
    # Getting prices from data
    p = @linq data |>
         where(:year.<2009) |>
         select(:oilpricereal2009)
    p=p./100 # Transforming price to price per 100 barrels

    #Transforming data into long form to get Quantities
    ld=stack(data, Between(:algeria, :uk))
    #Getting market shares
    @transform!(ld,q=:value ./ :totalworld)
    # Transforming data into wide form keeping quantities and time variables
    df=unstack(ld, [:year,:month], :variable, :q, allowduplicates=true)

    # Keeping only full years
    q = @linq df |>
            where(:year.<2009)

    # Dropping time variables
    q = q[:,Between(:algeria, :uk)]

    return Array(p), Array(q)
    println("Data is ready.")
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

function gM(p, W)
    T, N = size(W)

    g=[p⋅W[:,i] for i=1:N]

    return g
end

## mystars
function mystars(x)
    if x < 0.001
        y = "***"
    elseif x < 0.05
        y = "**"
    elseif x < 0.10
        y =  "*"
    else
        y = ""
    end
end

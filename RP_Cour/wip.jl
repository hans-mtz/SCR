## Loading modules
using Distributed
using DataFrames
using CSV
addprocs(6)

## Settin up directories

@everywhere dir = @__DIR__ ;
datadir=dir*"/data/oildata.csv" ;
resultsdir=dir*"/results/"

## Load modules in parallel
@everywhere include(dir*"/modules/MyModules.jl")

@everywhere begin
    using .MyModules
    using Random
    using SharedArrays
    ## Setting seed
    Random.seed!(66636)
end

@everywhere Param=Par(nfast=1000)


## Big 6

pt, qt = load_data_select(datadir)
set="Big6"

# Run tests

# Run base case
include("Main.jl")

# A) Test alternative transition kernel

include("RP_Cour_ME_sem_algo.jl")

# B) Test alternative acceptance rule with current
# transition kernel

# include("RP_Cour_ME_sem_alt.jl")

# C) Annual test

include("Main_annual.jl")


## OPEC

OPEC=[:saudiarabia, :iran, :venezuela, :uae, :nigeria, :kuwait]
pt, qt = load_data_select(datadir, select=OPEC)
set="OPEC"

# Run tests

# Run base case
include("Main.jl")
GC.gc()
# A) Test alternative transition kernel
include("RP_Cour_ME_sem_algo.jl")
GC.gc()
# C) Annual test
include("Main_annual.jl")
GC.gc()
## Non-OPEC
NonOPEC=[:unitedstates, :china, :mexico, :uk, :canada, :norway]
pt, qt = load_data_select(datadir, select=NonOPEC)
set="NonOPEC"

# Run tests

# Run base case
include("Main.jl")
GC.gc()
# A) Test alternative transition kernel
include("RP_Cour_ME_sem_algo.jl")
GC.gc()
# C) Annual test
include("Main_annual.jl")
GC.gc()
## OPEC no Saudia Arabia
OPEC_NSA=[:iran, :venezuela, :uae, :nigeria, :kuwait, :iraq]
pt, qt = load_data_select(datadir, select=OPEC_NSA)
set="OPEC_NSA"

# Run tests

# Run base case
include("Main.jl")
GC.gc()
# A) Test alternative transition kernel
include("RP_Cour_ME_sem_algo.jl")
GC.gc()
# C) Annual test
include("Main_annual.jl")
GC.gc()

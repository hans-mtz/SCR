
## Reproducibility instructions

# 1) Open Julia Pro/Atom
# 2) Juno: Work in folder (right click on folder)
# 3) Juno: Activate environment in folder (right click on folder)

## Create environments

# 1) set in Projects folder
# 2) julia -e 'Pkg.generate("Project_Name")'
# 3) set in project folder
# 4) julia --project=. / Reproducibility instructions

## Saving package versions in Project.toml and Manifest.toml
Pkg.add("Distributed", preserve=PRESERVE_DIRECT)
Pkg.add(["DataFrames","CSV",
         "DataFramesMeta", "LinearAlgebra",
         "Parameters", "JuMP", "Ipopt", "SharedArrays",
         "Statistics","BlackBoxOptim", "NLopt", "Distributions",
         "StatsBase"], preserve=PRESERVE_DIRECT)

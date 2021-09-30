module MyModules
# CourErr.jl functions
export load_data, load_data_select, gM, mystars
# MarkovChain.jl functions
export Par, Param, warm_start, jump, gchain
# Parallel_mod.jl functions
export shared_fill_chain!, fill_chain_sem!, myrange_sem
export shared_hM_fill!, hM_chunk!, myrange, objfun_par, objfun_bbo_par
export myrange_y, fill_table!, shared_fill_table!
# Testing.jl functions
export objfun, objfun_bbo, my_dist, RP_Cournot_test, acc_rate
export RP_Cournot_test_par

# Testing alternative transition kernel and alternative acceptance rule
export objfun_par_alt, objfun_bbo_par_alt, hM_chunk_alt!, RP_Cournot_test_par_alt
export gchain_algo, jump_algo, fill_chain_sem_algo!, fill_chain_sem_alt!, gchain_alt

include("CourErr.jl")
include("MarkovChain.jl")
include("Parallel_mod.jl")
include("Testing.jl")

end  # module MyModules

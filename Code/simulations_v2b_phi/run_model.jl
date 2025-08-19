# load packages
using Distributed, Random, DelimitedFiles, DataFrames, CSV, Statistics
@everywhere using Distributed, Random, DelimitedFiles, DataFrames, CSV, Statistics

# load functions
@everywhere include("Code/simulations_v2b_phi/functions_master.jl")
@everywhere include("Code/simulations_v2b_phi/functions_setup.jl")
@everywhere include("Code/simulations_v2b_phi/functions_extinction_colonisation.jl")
@everywhere include("Code/simulations_v2b_phi/functions_coevolution_evolution.jl")

# define function for running replicas and specify model inputs
@everywhere function run_parallel(network, forestCover, pvalue)
  
  # simulation parameters
  tmax = 1000   # maximum number of timesteps

  # extinction and colonisation probabilities
  e_r = 1    # extinction - resource
  e_c = 0.50    # extinction - consumer
  c_r = 1    # colonisation - resource
  c_c = 0.60    # colonisation - consumer

  # coevolution parameters
  phi_r = 0.1   # selection gradient resource
  phi_c = 0.5   # selection gradient consumer
  alpha = 0.2   # sensitivity to trait matching
  m_r = 0.5     # strength of coevolutionary selection - resource
  m_c = 0.7     # strength of coevolutionary selection - consumer
  epsilon = 5   # critical trait mismatch

  # run simulation (see functions_master)
  df_out = dynamics(network, forestCover, pvalue, tmax, e_r, e_c, c_r, c_c, phi_c, phi_r, alpha, m_r, m_c, epsilon)

  # create output directory
  isdir("Output_v2b_phi") || mkdir("Output_v2b_phi")

  # write out simulation output
  CSV.write(string("Output_v2b_phi/",network,"_fc",forestCover,"_p",pvalue,"_er",e_r,"_ec",e_c,"_cr",c_r,"_cc",c_c,"_mr",m_r,"_mc",m_c,"_phi_c",phi_c,"phi_r",phi_r,"_alpha",alpha,"_eps",epsilon,".csv"), df_out)

end

# number of simulations to run in parallel
n_sim = 4

# run simulations in parallel
# pmap(run_parallel, [network], [forestCover], [pvalue])

# eg: to run n_sim=3 different forest covers, with network="metanetwork", and pvalue=0.5, do the following:
#pmap(run_parallel, repeat(["metanetwork"], n_sim), repeat([0.1,0.3,0.5], inner=1), repeat([0.01], outer=n_sim))
pmap(run_parallel, repeat(["metanetwork"], n_sim), repeat([0.1,0.3,0.5, 0.7], inner=1), repeat([0.5], outer=n_sim))
pmap(run_parallel, repeat(["metanetwork"], n_sim), repeat([0.1,0.3,0.5, 0.7], inner=1), repeat([0.25], outer=n_sim))
pmap(run_parallel, repeat(["metanetwork"], n_sim), repeat([0.1,0.3,0.5, 0.7], inner=1), repeat([0.01], outer=n_sim))


#################
# TO RUN THIS SCRIPT
# DO THE FOLLOWING IN THE COMMAND LINE
# CHANGE "3" TO NUMBER OF SIMULATIONS TO BE RUN IN PARALLEL
# julia -p 3 run_model.jl
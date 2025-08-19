function setup_neighborhood(i, j, n)

  # Von Neumann neighborhood
	neigh = [i+0 j-1
           i-1 j+0
           i+1 j+0
           i+0 j+1]

	# Moore's neighborhood
	#neigh = [i-1 j-1
  #         i+0 j-1
  #         i+1 j-1
  #         i-1 j+0
	#         i+1 j+0
  #  	      i-1 j+1
  #         i+0 j+1
  #         i+1 j+1]

	# remove non-existent neighbours
	neigh = neigh[(neigh[:,1].>=1) .& (neigh[:,2].>=1) .& (neigh[:,1].<=n) .&(neigh[:,2].<=n),:]

	return neigh
end


function update_patches_sequentially(x_state, x_r, x_c, z_r, z_c, x_r_old, x_c_old, z_r_old, z_c_old, n_r, n_c, pf_r, e_r, e_c, c_r, c_c, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi_c, phi_r, n, M_inc)

  # update each patch sequentially
  for i = 1:n
    for j = 1:n

      # find patch neighbourhood
      neigh = setup_neighborhood(i, j, n)
    
      # resource extinctions and colonisations (see functions_extinction_colonisation)
      x_r, z_r = resource_extinctions_and_colonisations(x_state, x_r, z_r, x_r_old, z_r_old, n_r, pf_r, e_r, c_r, neigh, i, j)
    
      # consumer extinctions and colonisations (see functions_extinction_colonisation)
      x_c, z_c = consumer_extinctions_and_colonisations(x_c, z_c, x_r_old, z_r_old, x_c_old, z_c_old, n_c, e_c, c_c, alpha, neigh, M_inc, i, j)

      #sum(isnan.(skipmissing(z_r)))
      #sum(isnan.(skipmissing(z_c)))

      # species present in current patch
      r_current_sp = findall(x_r[i,j,:] .== 1)
      c_current_sp = findall(x_c[i,j,:] .== 1)
   
      # check if at least one species present (if not, skip to next patch)
      if length(r_current_sp) == 0 && length(c_current_sp) == 0
        continue
      end
    
      # simulate coevolution or evolution (see functions_coevolution_evolution)
      z_r, z_c, x_c = simulate_coevolution_or_evolution(z_r, z_c, x_c, r_current_sp, c_current_sp, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi_c, phi_r, M_inc, i, j)

      #if sum(isnan.(skipmissing(z_r[i,j,:]))) > 0 || sum(isnan.(skipmissing(z_c[i,j,:]))) > 0
      #  println(string(g," ",i,"-",j," ",sum(isnan.(skipmissing(z_r)))," ",sum(isnan.(skipmissing(z_c)))))
      #end

    end # j loop
  end # i loop

return x_r, x_c, z_r, z_c
end

function calculate_global_abundance_and_store(x, z, n, n_sp, df_dt, g)

	# calculate global abundance
	ab = sum(x, dims=[1,2])[1,1,:] ./ (n*n)

  # calculate trait mean and sd
  z_reshaped = reshape(z, ((n*n), n_sp))
  z_mean = mean.(skipmissing.(eachcol(z_reshaped)))
  z_sd = std.(skipmissing.(eachcol(z_reshaped)))

	# store results
	df_dt[(df_dt.t.==g),["abundance"]] .= ab
  df_dt[(df_dt.t.==g),["z_mean"]] .= z_mean
  df_dt[(df_dt.t.==g),["z_sd"]] .= z_sd

	return df_dt
end

function iterate_model_through_time(x_r, x_c, z_r, z_c, x_state, n_r, n_c, pf_r, e_r, e_c, c_r, c_c, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi_c, phi_r, n, tmax, M_inc)
#function iterate_model_through_time(x_r, x_c, z_r, z_c, x_state, n_r, n_c, pf_r, e_r, e_c, c_r, c_c, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi, n, tmax, M_inc, df_dt_r, df_dt_c)

  # iterate until tmax
  for g = 2:tmax

    # make copies of occupancy and trait arrays
    x_r_old = copy(x_r)
    x_c_old = copy(x_c)
    z_r_old = copy(z_r)
    z_c_old = copy(z_c)

    # update each patch sequentially
    x_r, x_c, z_r, z_c = update_patches_sequentially(x_state, x_r, x_c, z_r, z_c, x_r_old, x_c_old, z_r_old, z_c_old, n_r, n_c, pf_r, e_r, e_c, c_r, c_c, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi_c, phi_r, n, M_inc)

    # calculate global abundance and store results
    #df_dt_r = calculate_global_abundance_and_store(x_r, z_r, n, n_r, df_dt_r, g)
    #df_dt_c = calculate_global_abundance_and_store(x_c, z_c, n, n_c, df_dt_c, g)

    #df_dt_r[df_dt_r.t.==g,:]
    
    #sum(isnan.(skipmissing(z_r)))
    #sum(isnan.(skipmissing(z_c)))

  end # g loop

  # combine global abundance dataframes
  #df_dt_r[:,"guild"] .= "resources"
  #df_dt_c[:,"guild"] .= "consumers"
  #df_dt = vcat(df_dt_r, df_dt_c)
  
  return x_r, x_c, z_r, z_c
  #return x_r, x_c, z_r, z_c, df_dt
end


function build_output_data(n_patch, n_r, n_c, z_r, z_c)

  # create dataframes for resources and conusmers
  df_out_r = DataFrame(patch_no = repeat(1:n_patch, outer=n_r),
                       guild = fill("resources", (n_patch*n_r)),
                       species = repeat(1:n_r, inner=n_patch),
                       z = reshape(z_r, (n_patch*n_r)))
                       
  df_out_c = DataFrame(patch_no = repeat(1:n_patch, outer=n_c),
                       guild = fill("consumers", (n_patch*n_c)),
                       species = repeat(1:n_c, inner=n_patch),
                       z = reshape(z_c, (n_patch*n_c)))

  # combine resources and consumers dataframes
  df_out = vcat(df_out_r, df_out_c)

  # keep present species only
  df_out = dropmissing(df_out)
  
  return df_out
end


function dynamics(network, forestCover, pvalue, tmax, e_r, e_c, c_r, c_c, phi_c, phi_r, alpha, m_r, m_c, epsilon)

    # setup generated landscapes
    # n - landscape size (in x and y directions)
    # n_patch - number of patches (n*n)
    # x_state - array of patch states for each fragmented landscape (0=destroyed, 1=forest; dimensions n x n)
    n, n_patch, x_state = setup_landscapes(forestCover, pvalue)
    
    # get interaction incidence matrix, number of resources and number of consumers, and
    # forest occurences resources (see functions_setup)
    M_inc, n_r, n_c, pf_r = get_interaction_network(network)

    # get theta parameters for resources and consumers (see functions_setup)
    theta_r, theta_c = get_theta(network)

    # setup initial patch occupancy and traits (x - species presencce/absence; z - species trait values)
    # (see functions_setup)
    x_r, x_c, z_r, z_c = setup_grids(n, n_r, n_c, theta_r, theta_c)

    #using Random
    #Random.seed!(1)

    # setup dataframe for storing transient results (see funtions_setup)
    #df_dt_r = initialise_dataframes_store_results(tmax, n_r)  # resources
    #df_dt_c = initialise_dataframes_store_results(tmax, n_c)  # consumers
    #df_dt = DataFrame(t=Int[], species=Int[], abundance=Float32[], guild=String[])

    # iterate colonisation/extinction and coevolution dynamics through timesteps
    x_r, x_c, z_r, z_c = iterate_model_through_time(x_r, x_c, z_r, z_c, x_state, n_r, n_c, pf_r, e_r, e_c, c_r, c_c, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi_c, phi_r, n, tmax, M_inc)
    #x_r, x_c, z_r, z_c, df_dt = iterate_model_through_time(x_r, x_c, z_r, z_c, x_state, n_r, n_c, pf_r, e_r, e_c, c_r, c_c, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi, n, tmax, M_inc, df_dt_r, df_dt_c)
    
    # write out simulation output
    #CSV.write(string("../../Output_v2a/dt_",network,"_fc",forestCover,"_p",pvalue,"_er",e_r,"_ec",e_c,"_cr",c_r,"_cc",c_c,"_mr",m_r,"_mc",m_c,"_alpha",alpha,"_eps",epsilon,".csv"), df_dt)

    # store final timestep results as dataframe with trait values of species in patches
    df_out = build_output_data(n_patch, n_r, n_c, z_r, z_c)

  return df_out
end
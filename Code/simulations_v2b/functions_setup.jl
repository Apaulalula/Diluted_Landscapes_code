function setup_landscapes(forestCover, pvalue)

	# get generated landscape
	landscape = readdlm(string("Data/landscapes/forestcover_",forestCover,"_pvalue_",pvalue,"_edge.csv"), Int)

	# number of patches
	n_patch = size(landscape, 1)

	# landscape size
	n = Int(sqrt(n_patch))

	# grid of patch states (0=destroyed, 1=forest)
	x_state = reshape(landscape, (n,n))

	return n, n_patch, x_state
end


function get_interaction_network(network)

	# get network incidence matrix
	M_inc = readdlm(string("Data/",network,"/Minc.csv"), ' ', Int)

	# number of resources
	n_r = size(M_inc, 1)

	# number of consumers
	n_c = size(M_inc, 2)

	# get host forest occurrences
	pf_r = readdlm(string("Data/",network,"/host_forest_occurrence.csv"))

	return M_inc, n_r, n_c, pf_r
end


function get_theta(network)

	# theta of resources
	theta_r = readdlm(string("Data/",network,"/theta_r.csv"), Float32)

	# theta of consumers
	theta_c = readdlm(string("Data/",network,"/theta_c.csv"), Float32)

	return theta_r, theta_c
end


function setup_grids(n, n_r, n_c, theta_r, theta_c)

	# grid of resources (1=present, 0=absent) for each species
	x_r = ones(Int, n, n, n_r)

	# grid of consumers (1=present, 0=absent) for each species
	x_c = ones(Int, n, n, n_c)

	# grid of resource trait values (initial trait = theta)
	z_r = Array{Union{Missing, Float32}}(repeat(transpose(theta_r), outer=(n*n)))
	z_r = reshape(z_r, (n,n,n_r))
  
	# grid of consumer trait values (initial trait = theta)
	z_c = Array{Union{Missing, Float32}}(repeat(transpose(theta_c), outer=(n*n)))
	z_c = reshape(z_c, (n,n,n_c))

	return x_r, x_c, z_r, z_c
end


function initialise_dataframes_store_results(tmax, n)

	# initialise dataframes for storing results
	df_dt = DataFrame(t = repeat(1:tmax, inner=n),
	                  species = repeat(1:n, outer=tmax),
					  abundance = Vector{Union{Missing, Float32}}(missing, tmax*n),
					  z_mean = Vector{Union{Missing, Float32}}(missing, tmax*n),
					  z_sd = Vector{Union{Missing, Float32}}(missing, tmax*n))

	return df_dt
end
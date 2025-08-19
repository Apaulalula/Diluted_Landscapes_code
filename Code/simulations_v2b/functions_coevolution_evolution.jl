function get_local_network(r_current_sp, c_current_sp, M_inc)

    # current incidence matrix
    current_inc = M_inc[r_current_sp, c_current_sp]
  
    # non-interacting species
    r_nonint_sp = Array{Int64}(undef, 0)
    c_nonint_sp = Array{Int64}(undef, 0)
  
    if sum(current_inc) == 0
      r_nonint_sp = vcat(r_nonint_sp, r_current_sp)
      c_nonint_sp = vcat(c_nonint_sp, c_current_sp)
    end
    
    return current_inc, r_nonint_sp, c_nonint_sp
end


function setup_coevolution(current_inc, r_current_sp, c_current_sp, r_nonint_sp, c_nonint_sp, z_r, z_c, theta_r, theta_c, i, j)

    # find species with interactions
    r_int = vec(mapslices(col -> any(col .!= 0), current_inc, dims = 2))
    c_int = vec(mapslices(row -> any(row .!= 0), current_inc, dims = 1))
    r_int_sp = r_current_sp[r_int]
    c_int_sp = c_current_sp[c_int]

    # find non-interacting species
    r_nonint_sp = vcat(r_nonint_sp, setdiff(r_current_sp, r_int_sp))
    c_nonint_sp = vcat(c_nonint_sp, setdiff(c_current_sp, c_int_sp))

    # current incidence matrix without non-interacting species
    current_inc_coev = current_inc[r_int.==1,c_int.==1]

    # build adjacency matrix
    f = vcat(hcat(zeros(Int, size(current_inc_coev,1),size(current_inc_coev,1)), current_inc_coev), hcat(transpose(current_inc_coev), zeros(Int,size(current_inc_coev,2),size(current_inc_coev,2))))

    # trait values
    z = vcat(z_r[i,j,r_int_sp], z_c[i,j,c_int_sp])

    # theta values
    theta = vcat(theta_r[r_int_sp], theta_c[c_int_sp])

	return r_int_sp, c_int_sp, r_nonint_sp, c_nonint_sp, f, z, theta
end


function simulate_coevolution(f, z, phi, alpha, m_r, m_c, epsilon, theta, z_r, z_c, r_int_sp, c_int_sp, i, j)

    # number of resources, consumers, and total number of species
    n_r_current = length(r_int_sp)
    n_c_current = length(c_int_sp)
    n_sp_current = n_r_current + n_c_current

    # vector of m values for all species
    m = vcat(repeat([m_r], n_r_current), repeat([m_c], n_c_current))

	# trait change due to coevolution
	#z_dif = transpose(f.*z) - f.*z
    z_dif = Float64.(transpose(f.*z) - f.*z)
    q = f.*(exp.(-alpha.*(z_dif.^2)))
    q_n = q./sum(q,dims = 2)
    q_m = q_n.* m
    z_subset_og = z_dif[1:n_r_current,(n_r_current+1):n_sp_current]
    u = 1 .* (broadcast(abs, z_subset_og) .<= epsilon)
    z_subset_mod = copy(z_subset_og)
    z_subset_mod[findall(z_subset_og .>= 0)] = z_subset_og[findall(z_subset_og .>= 0)] .- epsilon
    z_subset_mod[findall(z_subset_og .< 0)] = z_subset_og[findall(z_subset_og .< 0)] .+ epsilon
    z_dif[1:n_r_current,(n_r_current+1):n_sp_current] = u .* z_subset_mod
    sel_dif = q_m .* z_dif
    r_mut = phi .* sum(sel_dif,dims=2)
    r_env = phi .* (1 .- m) .* (theta .- z)
    z_new = z .+ r_mut .+ r_env

    #sum(isnan.(q_n))
    #sum(isnan.(r_mut))
    #sum(isnan.(r_env))
    #sum(isnan.(z_new))

    #writedlm( "z_dif.csv",  z_dif, ',')
    #writedlm( "f.csv",  f, ',')
    #writedlm( "q.csv",  q, ',')
    #writedlm( "q_n.csv",  q_n, ',')

    # update trait values arrays
    z_r[i,j,r_int_sp] .= z_new[1:n_r_current]
    z_c[i,j,c_int_sp] .= z_new[(n_r_current+1):end]

	return z_r, z_c
end


function simulate_evolution(nonint_sp, z_n, phi, theta_n, i, j)

    # theta values
    theta = theta_n[nonint_sp]

    # trait values
    z = z_n[i,j,nonint_sp]

    # trait change due to environmental selection only
    r_env = phi .* (theta .- z)
    z_new = z .+ r_env

    # update trait values arrays
    z_n[i,j,nonint_sp] .= z_new

    return z_n
end


function simulate_coevolution_or_evolution(z_r, z_c, x_c, r_current_sp, c_current_sp, theta_r, theta_c, alpha, m_r, m_c, epsilon, phi, M_inc, i, j)

    # get local network incidence matrix and noninteracting species
    current_inc, r_nonint_sp, c_nonint_sp = get_local_network(r_current_sp, c_current_sp, M_inc)

    # check if at least one interaction present -> coevolution
    if sum(current_inc) != 0
  
        # setup parameters for coevolution
        r_int_sp, c_int_sp, r_nonint_sp, c_nonint_sp, f, z, theta = setup_coevolution(current_inc, r_current_sp, c_current_sp, r_nonint_sp, c_nonint_sp, z_r, z_c, theta_r, theta_c, i, j)
  
        # simulate coevolution (see functions_coevolution_mutualistic/antagonistic)
        z_r, z_c = simulate_coevolution(f, z, phi, alpha, m_r, m_c, epsilon, theta, z_r, z_c, r_int_sp, c_int_sp, i, j)
  
    end

    # if non-interacting resources present -> environmental selection only
    #if length(r_nonint_sp) != 0 || length(c_nonint_sp) != 0
    if length(r_nonint_sp) != 0
  
        # simulate evolution (only environmental selection)
        z_r = simulate_evolution(r_nonint_sp, z_r, phi, theta_r, i, j)
        #z_c = simulate_evolution(c_nonint_sp, z_c, phi, theta_c, i, j)
  
    end

    # if non-interacting consumers present -> extinct
    if length(c_nonint_sp) != 0
  
        # non interacting parasites extinct
        z_c[i,j,c_nonint_sp] .= missing
        x_c[i,j,c_nonint_sp] .= 0
  
    end
  
    return z_r, z_c, x_c
end
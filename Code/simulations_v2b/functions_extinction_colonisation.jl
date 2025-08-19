# RESOURCE FUNCTIONS

function resource_extinctions(x_r, z_r, n_r, e_r, pf_r, x_state, i, j)

	# extinction probabilities
	pe_r = e_r .* (x_state[i,j]/5 .- pf_r.*(2/5*x_state[i,j]-1))

	# extinct species (if random number < extinction probability)
    sp_ext = findall((rand(n_r) .< pe_r)[:,1])

    # update species presence and traits
    x_r[i,j,sp_ext] .= 0
    z_r[i,j,sp_ext] .= missing

	return x_r, z_r
end


function resource_colonisations(neigh, x_r_old, z_r_old, x_r, z_r, c_r, p, x_state, i, j, h)

	# check neighbours
    for f = 1:size(neigh, 1)

		i_n = neigh[f,1]
        j_n = neigh[f,2]

    	# if resource h absent from neighbouring patch -> skip
    	if x_r_old[i_n,j_n,h] == 0
    		continue

    	# if resource h present in neighbouring patch
		else

			# colonisation probability
			pc_r = c_r * (p.*(2/5*x_state[i,j]-1) .- x_state[i,j]/5 .+ 1)

            # colonisation if random number < colonisation probability
    		if rand() < pc_r
        		x_r[i,j,h] = 1
         	   	z_r[i,j,h] = z_r_old[i_n,j_n,h]
        	    # stop colonisations
				break
       		end

    	end

    end # f loop

	return x_r, z_r
end


function resource_extinctions_and_colonisations(x_state, x_r, z_r, x_r_old, z_r_old, n_r, pf_r, e_r, c_r, neigh, i, j)

	# resource extinctions
	x_r, z_r = resource_extinctions(x_r, z_r, n_r, e_r, pf_r, x_state, i, j)

	# resource absent from patch i,j
    sp_absent = findall(x_r_old[i,j,:] .== 0)

	# resource colonisations - loop through absent resource
    for h = sp_absent

		# resource colonisations
		x_r, z_r = resource_colonisations(neigh, x_r_old, z_r_old, x_r, z_r, c_r, pf_r[h], x_state, i, j, h)

    end # h loop

	return x_r, z_r
end



# CONSUMER FUNCTIONS

function get_resource_partners(M_inc, x_r_old, i, j, h)

    # get resource partners of consumer h
	r_part = findall(M_inc[:,h].==1)

	# get resources present in current patch
	r_patch_all = findall(x_r_old[i,j,:].==1)

	# get resource partners of consumer h which are present in current patch
	r_part_sp = r_patch_all[findall(x-> x in r_part, r_patch_all)]

	return r_part_sp
end


function get_trait_matching_patch(z_r_part, z_c_old, alpha, i, j, h)

	# get trait value of consumer
	z_c_s = z_c_old[i,j,h]

	# compute trait differences
	z_dif = z_c_s .- z_r_part

	# compute trait matching
	p_match = exp.(-alpha.*((z_dif).^2))

	return p_match
end


function compute_consumer_extinction_prob(e_c, p_match)

	# initialise extinction probability
    pe_c = 1
    j_count = 0

    for m=1:length(p_match)
        j_count = j_count+1
        prod = (1 - e_c / j_count) * (1 - p_match[j_count] * e_c)
        pe_c = pe_c * prod
    end

	return pe_c
end


function consumer_extinctions(x_c, z_c, z_r_part, z_c_old, e_c, alpha, M_inc, i, j, h)

	# if no resources present -> consumer extinct
    if length(z_r_part) == 0
        x_c[i,j,h] = 0
    	z_c[i,j,h] = missing
	end

    # if at least one resource partner present in current patch
    if length(z_r_part) > 0

		# get trait matching between interacting species in current patch
		p_match = get_trait_matching_patch(z_r_part, z_c_old, alpha, i, j, h)

		# compute extinction probability
        pe_c = compute_consumer_extinction_prob(e_c, p_match)

		# extinction if random number < extinction probability
        if rand() < pe_c
            x_c[i,j,h] = 0
			z_c[i,j,h] = missing
        end
    end

	return x_c, z_c
end


function compute_consumer_colonisation_prob(z_c_old, z_r_part, c_c, alpha, i_n, j_n, h)

    # trait matching between consumer and resources
	p_match = get_trait_matching_patch(z_r_part, z_c_old, alpha, i_n, j_n, h)

    # initialise "not colonised" probability
    pc_c = 1
    j_count = 0

    for m=1:length(p_match)
		j_count = j_count+1
		prod = (1 - c_c / j_count) * (1 - p_match[j_count] * c_c)
		pc_c = pc_c * prod
    end

    return pc_c
end


function consumer_colonisations(x_c, z_c, x_c_old, z_c_old, z_r_part, c_c, alpha, neigh, i, j, h)

	# check neighbours
	for f=1:size(neigh,1)
		
		i_n = neigh[f,1]
		j_n = neigh[f,2]
		
		# if consumer h absent from neighbouring patch -> next neighbour
		if x_c_old[i_n,j_n,h] == 0
			continue
		# if consumer h present in neighbouring patch
		else
			# compute "not colonised" probability
			pc_c = compute_consumer_colonisation_prob(z_c_old, z_r_part, c_c, alpha, i_n, j_n, h)

			# colonisation if random number > "not colonised" probability
			if rand() > pc_c
				x_c[i,j,h] = 1
				z_c[i,j,h] = z_c_old[i_n,j_n,h]
				# stop colonisations
				break
			end

        end

	end # f loop

	return x_c, z_c
end


function consumer_extinctions_and_colonisations(x_c, z_c, x_r_old, z_r_old, x_c_old, z_c_old, n_c, e_c, c_c, alpha, neigh, M_inc, i, j)

	# loop through consumers
    for h=1:n_c

		# find resource partners of consumer h which are present in patch i
		r_part_sp = get_resource_partners(M_inc, x_r_old, i, j, h)
		
		# get traits of resource partners
		z_r_part = z_r_old[i,j,r_part_sp]

    	# if species present -> extinctions
        if x_c_old[i,j,h] == 1
			x_c, z_c = consumer_extinctions(x_c, z_c, z_r_part, z_c_old, e_c, alpha, M_inc, i, j, h)
            continue
        end

        # if species absent & at least one resource present -> colonisations
        if x_c_old[i,j,h] == 0 && length(r_part_sp) > 0
            x_c, z_c = consumer_colonisations(x_c, z_c, x_c_old, z_c_old, z_r_part, c_c, alpha, neigh, i, j, h)
        end

    end # h loop

	return x_c, z_c
end
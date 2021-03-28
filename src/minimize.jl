function mod_dijk(g::VariationalGraph{U, T}, v0::Array{T,1}, alpha::T) where{T<:Real, U<:Integer}
    # -------------------------------------------------------------------------
    # Intialization
    num_verts = length(v0) # save number of verts
    finished = fill(false, num_verts) # array of bools if vertex is finished
    Tv = intersect(findall(x->(x<Inf), v0), g.verts) # all terminal vertices w.r.t. v in the current subgraph
    parents = Array{U, 1}(undef, g.num_verts) # saves parents for heap structure
    v = fill(Inf64, g.num_verts) # the values we will return
    #
    loc_verts = g.glob_to_loc # useful for subgraph assignment
    # set parents and values of terminal nodes
    for i in Tv
        parents[loc_verts[i]] = -1
        v[loc_verts[i]] = v0[i]
    end
    # create a binary heap for sorting purpose
    h = DataStructures.MutableBinaryMinHeap(v)
    # -------------------------------------------------------------------------
    # Assign all values for v
    while length(h) > 0
        v_new, i = DataStructures.top_with_handle(h) # get smallest element
        DataStructures.delete!(h, i) # take this element out of the heap
        v[i] = v_new # update value, might not be set yet!
        finished[i] = true # mark as visited
        #
        # loop over all neigbors of i
        for j in g.edges_list[i]
            it = 1 # count for neighbor
            if !finished[loc_verts[j]] # node is not finished
                ind = h.node_map[loc_verts[j]] # get index of key in heap
                new_val = v_new + alpha * g.weights_list[i][it] # potentially better value
                if h.nodes[ind].value > new_val # if true update the value
                    DataStructures.update!(h, loc_verts[j], new_val)
                    parents[loc_verts[j]] = i
                end
            end
            it = it+1
        end
    end
   return v, parents
end

# Different name for mod_dijk for consistency in the naming scheme
function comp_v_low(g::VariationalGraph{U, T}, v0::Array{T,1}, alpha::T) where{T<:Real, U<:Integer}
    return mod_dijk(g, v0, alpha)
end

function comp_v_high(g::VariationalGraph{U, T}, v0::Array{T,1}, alpha::T) where{T<:Real, U<:Integer}
    Tv = findall(x->(x<Inf), v0) # all terminal vertices w.r.t. v
    v1 = fill(Inf, length(v0)) # initialize new Array
    v1[Tv] .= -v0[Tv] # set all terminal values to the neagative value
    v2, parents = mod_dijk(g, v1, alpha) # apply mod dijkstra
    return -v2, parents
end

function comp_lip_extension(g::VariationalGraph{U, T}, v0::Array{T,1}, alpha::T) where{T<:Real, U<:Integer}
    Tv0 = findall(x->(x<Inf), v0)
    v, _ = mod_dijk(g, v0, alpha)
    v[Tv0] = v0[Tv0] # reset constraint if it got lost for small alpha
    return v
end 

function comp_high_pressure(g::VariationalGraph{U, T}, v0::Array{T,1}, alpha::T) where{T<:Real, U<:Integer}
    # compute both lipschitz extension w.r.t alpha
    v_low, low_par = comp_v_low(g, v0, alpha)
    v_high, high_par = comp_v_high(g, v0, alpha)
    # we are intersted in the subgraph were vHigh > v_low
    verts = [u for u in g.verts if round(v_high[g.glob_to_loc[u]] - v_low[g.glob_to_loc[u]], digits=5) > 0]
    return subgraph(g, verts)
end

function steepest_path_star(g::VariationalGraph{U, T}, Tv0::Array{U, 1}, v0::Array{T,1}, x::U) where{T<:Real, U<:Integer}
    @assert (length(Tv0) > 0) "There is a subgraph with no terminal vertices, cannot continue!"
    # -------------------------------------------
    x_loc = g.glob_to_loc[x]
    t_1 = Tv0[rand(1:length(Tv0))]
    dist_t_1 = g.all_weights[x_loc, g.glob_to_loc[t_1]]
    if dist_t_1 == typemax(T) # no path from x_loc to t_1
        return steepest_path_star(g, setdiff(Tv0, t_1), v0, x)
    end
    # determine the steepest terminal path from t_1
    t_2 = t_1 # if we dont find any better terminal we are the winner ourselves!
    max_grad = 0.0
    for t in Tv0
        new_grad = abs(v0[t_1] - v0[t])/(g.all_weights[x_loc, g.glob_to_loc[t]] + dist_t_1)
        if max_grad < new_grad # found a bigger gradient
            max_grad = new_grad
            t_2 = t
        end
    end
    # compute Whitney and Shane extension
    v_low = minimum(v0[t] +  max_grad * g.all_weights[x_loc, g.glob_to_loc[t]] for t in Tv0)
    v_high = maximum(v0[t] -  max_grad * g.all_weights[x_loc, g.glob_to_loc[t]] for t in Tv0)
    # determin the repective terminal verts 
    Tv0_low = filter(t-> 0.0 > round(v_low + max_grad * g.all_weights[x_loc, g.glob_to_loc[t]] - v0[t], digits = 4), Tv0)
    Tv0_high = filter(t-> 0.0 < round(v_high - max_grad * g.all_weights[x_loc, g.glob_to_loc[t]] - v0[t], digits = 4), Tv0)
    # new terminal set as union
    T_new = union(Tv0_low, Tv0_high)
    if length(T_new) == 0 # we found the steepset terminal path
        if v0[t_1] >= v0[t_2]
            return t_1, t_2
        else 
            return t_2, t_1
        end
    else
        return steepest_path_star(g, T_new, v0, x)
    end
end

function vertex_steepest_path(g::VariationalGraph{U, T}, v0::Array{T,1}, s_paths::Array{U, 1}, x::U) where{T<:Real, U<:Integer}
    x_loc = g.glob_to_loc[x]
    for i in 1:g.num_verts s_paths[i]=-1 end # reset the paths
    dijk_shortest_path!(g, s_paths, x) # find shortest path fom source x
    Tv0 = findall(x->x<Inf, v0) # all terminal values
    Tv0 = filter(t->g.glob_to_loc[t]!=0, Tv0) # filter all terminal nodes that ae not in the subgraph
    if v0[x] < Inf # x is terminal node
        max_grad = 0.0
        y_max = x
        for y in Tv0
            # determin gradient of this terminal path
            new_grad = abs(v0[x] - v0[y])/g.all_weights[x_loc, g.glob_to_loc[y]] 
            if max_grad < new_grad # found a bigger gradient
                max_grad = new_grad
                y_max = y
            end
        end
        if v0[x] >= v0[y_max]
            path =  shortest_path_from_mat(g, s_paths, x, y_max)
        else
            path = shortest_path_from_mat_rev(g, s_paths, x, y_max)
        end
        len = g.all_weights[x_loc, g.glob_to_loc[y_max]]
    else # x is not terminal
        # we need to find the steepest terminal path that includes x
        t_1, t_2 = steepest_path_star(g, Tv0, v0, x)
        path_from_1_to_x = shortest_path_from_mat_rev(g, s_paths, x, t_1)
        path_from_x_to_2 = shortest_path_from_mat(g, s_paths, x, t_2)
        path = vcat(path_from_1_to_x, path_from_x_to_2[2:end])
        len = g.all_weights[x_loc, g.glob_to_loc[t_1]] + g.all_weights[x_loc, g.glob_to_loc[t_2]]
    end
    return path, len

end

function steepest_path(g::VariationalGraph{U, T}, v0::Array{T,1}) where{T<:Real, U<:Integer}
    e = rand(1:g.num_edges) # random edge
    x_1_loc, x_2 = find_edge_ind(g, e) # local indices of this edge
    x_3_loc = rand(1:g.num_verts) # random vertex
    x = [g.verts[x_1_loc], x_2, g.verts[x_3_loc]]
    s_paths = fill(-1, g.num_verts) # Array that stores shortest paths
    # steepest paths
    max_grad = typemin(T)
    P_max = U[]
    len_max = 0.0
    # Assign the three potential steepest paths and compute its gradients 
    for i = 1:3
        P, len = vertex_steepest_path(g, v0, s_paths, x[i])
        new_grad = (v0[P[1]] - v0[P[end]])/len
        if new_grad > max_grad
            max_grad = new_grad
            P_max = P
            len_max = len
        end
    end
    h = comp_high_pressure(g, v0, max_grad)
    if h.num_verts == 0
        return P_max, len_max, max_grad
    else
        return steepest_path(h, v0)
    end
end

function comp_inf_min(g_in::VariationalGraph{U, T}, v0::Array{T,1}) where{T<:Real, U<:Integer}
    # initialize fields for construction of new graph from the input
    edges_list = deepcopy(g_in.edges_list) # the graph will be modified
    weights_list =  deepcopy(g_in.weights_list)
    num_edges = g_in.num_edges
    #
    Tv0 = findall(x->x<Inf, v0) # all terminal nodes
    alpha = 0.0 # value for the maximal gradient
    # -------------------------------------------
    # first maximize over all terminal edges
    for t_1_loc = 1:g_in.num_verts
        if g_in.verts[t_1_loc] in Tv0 # terminal vertex
            s = 1
            while s <= length(edges_list[t_1_loc])
                t_2 = edges_list[t_1_loc][s]
                if t_2 in Tv0 # found terminal edge
                    new_alpha = abs(v0[t_1_loc] - v0[t_2])/weights_list[t_1_loc][s]
                    alpha = max(alpha, new_alpha) # update max gradient value
                    # now delete this terminal-terminal edge
                    deleteat!(edges_list[t_1_loc], s)
                    deleteat!(weights_list[t_1_loc], s)
                    num_edges-=1
                    ind_sym = findfirst(x->x==g_in.verts[t_1_loc], edges_list[g_in.glob_to_loc[t_2]]) 
                    if ind_sym != nothing # if edge is symmetric
                        deleteat!(edges_list[g_in.glob_to_loc[t_2]], ind_sym)
                        deleteat!(weights_list[g_in.glob_to_loc[t_2]], ind_sym)
                        num_edges-=1
                    end
                else # go to next edge
                    s+=1
                end
            end
        end
    end
    # construct new graph
    g = VariationalGraph(g_in.verts, g_in.glob_to_loc, g_in.num_verts, num_edges, edges_list, weights_list, 
                            fill(type_max(T), g_in.num_verts, g_in.num_verts), childgraphdef(g_in))
    if g.num_edges == 0
        return v0
    end
    # -------------------------------------------
    # start calculations on new graph    
    P, _, grad = steepest_path(g, v0) # find steepest path with max grad
    alpha = max(alpha, grad) # uodate alpha if grad is bigger
    v_low, par_low = comp_v_low(g, v0, alpha)
    v_high, par_high = comp_v_high(g, v0, alpha)
    # -------------------------------------------
    v = similar(v0) # assign the output
    for x in g.verts
        if x in Tv0
            v[x] = v0[x]
        else
            v[x] = 1/2 * (v_low[x] + v_high[x])
        end
    end
    return v
end

function comp_max_grad(g::VariationalGraph{U, T}, v::Array{T,1}) where{T<:Real, U<:Integer}
    # initialize values
    max_grad = 0.0
    max_t_1 = 1
    max_t_2 = 1
    # loop over all edges
    for t_1 = 1:g.num_verts
        for t_2 = 1:length(g.edges_list[t_1])
            new_grad = abs(v[g.verts[t_1]] - v[g.edges_list[t_1][t_2]])/g.weights_list[t_1][t_2]
            if new_grad > max_grad
                max_t_1 = t_1
                max_t_2 = t_2
                max_grad = new_grad
            end
        end
    end
    return max_grad, g.verts[max_t_1], g.edges_list[max_t_1][max_t_2]
end

function comp_inf_min_approx(g::VariationalGraph{U, T}, v0::Array{T,1}, rel_tol::T) where{T<:Real, U<:Integer}
    Tv0 = findall(x->x<Inf, v0) # all terminal nodes
    alpha = 0.0 # value for the maximal gradient
    # -------------------------------------------
    # find lower bound for alpha
    v_max = maximum(v0[Tv0])
    v_min = minimum(v0[Tv0])
    d_max = 0.0
    for u in 1:length(g.edges_list)
        for v in 1:length(g.edges_list[u])
            d_max = max(d_max, g.weights_list[u][v])
        end
    end
    alpha = abs(v_max - v_min)/d_max
    #--------------------------------------------
    # values for iteration
    v_new = similar(v0)
    alpha_real = alpha
    alpha_step = 0.5
    while alpha_step > rel_tol
        v_new, _ = comp_lip_extension(g, v0, alpha)
        v_new[Tv0] = v0[Tv0] # reset constraint
        alpha_real, _, _ = comp_max_grad(g, v_new)
        if alpha_real <= alpha
            while alpha <= alpha_step
                alpha_step/=2
            end 
            alpha = alpha - alpha_step
        else
            alpha = alpha + alpha_step # reset the step
            alpha_step/=2
            alpha = alpha - alpha_step # update with halfed step
        end
    end
    return v_new, alpha_real
end

function fix_steepest!(g::VariationalGraph{U, T}, v0::Array{T,1}, Tv0::Array{U, 1}, path::Array{U, 1}, len::T) where{T<:Real, U<:Integer}
    if length(path) <= 2 return end # not a free terminal path
    dist = 0 # distance to point in graph
    u_old = path[1] # starting vetex
    v_start = v0[path[1]] # we will need this value a lot
    grad_path = (v_start - v0[path[end]])/len # gradient of this steepest path
    # iterate over the path to set the new assignment
    for u in setdiff(path[2:end-1])
        v_ind = findfirst(x->x==u, g.edges_list[g.glob_to_loc[u_old]]) # find the edge u-u_old
        dist += g.weights_list[g.glob_to_loc[u_old]][v_ind] # update the distance to this point
        if !(u in Tv0)
            v0[u] = v_start - grad_path * dist # set the new value for v0
        end
        u_old = u # reset the previous node
    end
end 

function comp_lex_min(g_in::VariationalGraph{U, T}, v0_in::Array{T,1}) where{T<:Real, U<:Integer}
    v0 = copy(v0_in) # v0 will be modified in this context
    # initialize fields for construction of new graph from the input
    edges_list = deepcopy(g_in.edges_list)
    weights_list =  deepcopy(g_in.weights_list)
    num_edges = g_in.num_edges
    num_verts = g_in.num_verts
    verts = g_in.verts
    #
    Tv0 = findall(x->x<Inf, v0) # all terminal nodes
    # -------------------------------------------
    # loop while not all vertices are terminal for v0
    all_weights = fill(typemax(T), g_in.num_verts, g_in.num_verts)
    path = Tv0
    while length(Tv0) < g_in.num_verts
        # delete all terminal terminal edges
        for t_1 in path
            for t_2 in Tv0
                ind = findfirst(x->x==t_2, edges_list[g_in.glob_to_loc[t_1]])
                ind2 = findfirst(x->x==t_1, edges_list[g_in.glob_to_loc[t_2]])
                if ind != nothing
                    deleteat!(edges_list[g_in.glob_to_loc[t_1]], ind) # delete all terminal terminal
                    deleteat!(weights_list[g_in.glob_to_loc[t_1]], ind)
                    num_edges-=1
                end
                if ind2 != nothing
                    if ind==nothing
                        display(t_1)
                        display(t_2)
                    end
                    deleteat!(edges_list[g_in.glob_to_loc[t_2]], ind2) # delete all terminal terminal
                    deleteat!(weights_list[g_in.glob_to_loc[t_2]], ind2)
                    num_edges-=1
                end
            end
        end
        # e = rand(1:g.num_edges) # random edge
        # x_1_loc, x_2 = find_edge_ind(g, e) # local indices of this edge
        # x_3_loc = rand(1:g.num_verts) # random vertex
        # x = [g.verts[x_1_loc], x_2, g.verts[x_3_loc]]
        # s_paths = fill(-1, g.num_verts) # Array that stores shortest paths
        # # steepest paths
        # max_grad = typemin(T)
        # P_max = U[]
        # len_max = 0.0
        # # Assign the three potential steepest paths and compute its gradients 
        # for i = 1:3
        #     P, len = vertex_steepest_path(g, v0, s_paths, x[i])
        #     new_grad = (v0[P[1]] - v0[P[end]])/len
        #     if new_grad > max_grad
        #         max_grad = new_grad
        #         P_max = P
        #         len_max = len
        #     end
        # end
        # h = comp_high_pressure(g, v0, max_grad)


        # construct the new geph without terminal-terminal edges
        g = VariationalGraph(verts, g_in.glob_to_loc, num_verts, num_edges, edges_list, weights_list, 
                             fill!(all_weights, typemax(T)), childgraphdef(g_in))
        path, len, grad = steepest_path(g, v0) # steepest path in this graph
        fix_steepest!(g, v0, Tv0, path, len) # assignment via new steepest path
        Tv0 = findall(x->x<Inf, v0) # update terminal nodes
    end
    return v0
end


function comp_inf_lap(g_in::VariationalGraph{U, T}, v0_in::Array{T,1}) where{T<:Real, U<:Integer}
    v = copy(v0_in) # v will be the solution
    M = 1
    # initialize fields for construction of new graph from the input
    edges_list = deepcopy(g_in.edges_list)
    weights_list =  deepcopy(g_in.weights_list)
    num_edges = g_in.num_edges
    num_verts = g_in.num_verts
    verts = g_in.verts
    # loop until terminaion
    Tv0 = findall(x->x<Inf, v0_in) # all terminal nodes
    for i=1:5
        # -------------------------------------------
        # loop over all vertices
        for t_1 = 1:g_in.num_verts
            max_grad = typemax(T)
            min_grad = typemin(T)
            # loop over all neighbors
            for t_2 = 1:length(g_in.edges_list[t_1])
                new_grad = (v[g_in.edges_list[t_1][t_2]] - v[g_in.verts[t_1]])/g_in.weights_list[t_1][t_2]
                if new_grad > max_grad
                    max_grad = new_grad
                end
                if new_grad < min_grad
                    min_grad = new_grad
                end
            end
            # we found the min max value update v locally at t_1
            v[g_in.verts[t_1]] += 1/(2*M) *(max_grad + min_grad)
        end
        # reset the Dircihlet condition
        v[Tv0]=v0_in
    end
    return v
end
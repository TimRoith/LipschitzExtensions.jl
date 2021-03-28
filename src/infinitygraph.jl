###################################################################################################
"""
```julia
    GraphConfig
```
The struct used to dispatch different graph configurations.
"""
abstract type GraphConfig end
struct vargraphdef <: GraphConfig end
###################################################################################################
# The struct representing a Graph and outer constructor functions.
"""
    VariationalGraph{T, U}
A type representing an directed weighted Graph graph.
"""
struct VariationalGraph{T <: Integer, U <: Real} <: AbstractVariationalGraph{T, U}
    verts::Array{T, 1}
    glob_to_loc::Array{T, 1}
    #
    num_verts::Int
    num_edges::Int
    #
    edges_list::Array{Array{T, 1}, 1}
    weights_list::Array{Array{U, 1}, 1}
    all_weights::Array{U,2}
    #
    config::GraphConfig
end
#------------------------------------------------
# config not defined
VariationalGraph(verts::Array{T, 1}, glob_to_loc::Array{Union{Nothing, T},1}, num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1},
                 weight_list::Array{Array{U, 1}, 1}) where {T <: Integer, U <: Real} =
VariationalGraph(verts, glob_to_loc, num_verts, num_edges, edge_list, weight_list, vargraphdef())
#------------------------------------------------
VariationalGraph(num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1},
                 weight_list::Array{Array{U, 1}, 1}, config::GraphConfig) where {T <: Integer, U <: Real} =
VariationalGraph(collect(1:num_verts), collect(1:num_verts), num_verts, num_edges, edge_list, weight_list, 
                 fill(typemax(U), num_verts, num_verts), config)
VariationalGraph(num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1},
                 weight_list::Array{Array{U, 1}, 1}) where {T <: Integer, U <: Real} = 
VariationalGraph(collect(1:num_verts), collect(1:num_verts), num_verts, num_edges, edge_list, weight_list, 
                 fill(typemax(U), num_verts, num_verts), vargraphdef())              
#------------------------------------------------
# function VariationalGraph(num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1}, 
#                           weight_list::Array{Array{U, 1}, 1}, config::GraphConfig) where {T <: Integer, U <: Real}
#     edge_mat, weight_mat = list2mat(num_edges, edge_list, weight_list)
#     return VariationalGraph(num_verts, num_edges, edge_list, edge_mat, weight_list, weight_mat, config)
# end
#------------------------------------------------                
VariationalGraph(num_edges::T, edges::Array{Array{T, 1}, 1}, weights::Array{Array{U, 1}, 1}) where {T <: Integer, U <: Real} = 
VariationalGraph(length(edges), num_edges, edges, weights)
#------------------------------------------------
###################################################################################################

struct childgraphdef{T<:Integer, U<:Real} <: GraphConfig
    parent::VariationalGraph{T, U}
end

"""
```julia
    subgraph()
```
Extract a subgraph by a given subset of vertices. Only the edges that are contained in 
the new maximal edgesubset are kept.
"""
function subgraph(g::VariationalGraph{U, T}, in_verts::Array{U,1}) where {U<:Integer, T<:Real}
    verts = intersect(in_verts, g.verts) # vertices of the new graph
    # initialize new graph
    edges_list = Array{Array{U, 1}, 1}(undef, 0)
    weights_list = Array{Array{T, 1}, 1}(undef, 0)
    num_edges = 0
    num_verts = length(verts)
    # add new edges
    for u in verts
        ind = filter!(x->x!=nothing, indexin(verts, g.edges_list[g.glob_to_loc[u]]))
        push!(edges_list, g.edges_list[g.glob_to_loc[u]][ind])
        push!(weights_list, g.weights_list[g.glob_to_loc[u]][ind])
        num_edges += length(ind)
    end
    # Create the array for the local index conversion
    glob_to_loc = replace(indexin(collect(1:length(g.glob_to_loc)), verts), nothing=>0)
    return VariationalGraph(verts, glob_to_loc, num_verts, num_edges, edges_list, weights_list, 
                            fill(Inf, num_verts, num_verts), childgraphdef(g))
end


function find_edge_ind(g::VariationalGraph{U, T}, e::U) where{T<:Real, U<:Integer}
    u = 0
    it = 0
    while it < e
        u+=1
        it += length(g.edges_list[u])
    end
    if it == e
        v = length(g.edges_list[u])
    else
        v = e - (it - length(g.edges_list[u]))
    end
    return u, g.edges_list[u][v]
end

function dijk_shortest_path!(g::VariationalGraph{U, T}, s_paths::Array{U, 1}, src::U) where{T<:Real, U<:Integer}
    q = U[]
    src_loc = g.glob_to_loc[src] # local source index
    for v in g.verts # we consider all verices that are in the graph
        # if g.all_weights[src_loc, g.glob_to_loc[v]] == typemax(T) # this distance has not been set yet
        #     push!(q, v)
        # end
        push!(q, v)
    end
    g.all_weights[src_loc, src_loc] = 0
    # itarate over the list q
    while !isempty(q)
        u = argmin(g.all_weights[src_loc, g.glob_to_loc[q]]) # get pos of minimal distance value
        u = q[u] # set index of minimal value
        filter!(x->x!=u, q) # delete this vertix
        u_loc = g.glob_to_loc[u] # set u to the local index
        for v_loc = 1:length(g.edges_list[u_loc])
            if g.edges_list[u_loc][v_loc] in q
                new_len = g.weights_list[u_loc][v_loc] + g.all_weights[src_loc, u_loc]
                edge_loc = g.glob_to_loc[g.edges_list[u_loc][v_loc]]
                if g.all_weights[src_loc, edge_loc] > new_len
                    g.all_weights[src_loc, edge_loc] = new_len
                    g.all_weights[edge_loc, src_loc] = new_len
                    s_paths[edge_loc] = u
                elseif round(g.all_weights[src_loc, edge_loc] - new_len, digits = 5) == 0 # set spath anyway
                    s_paths[edge_loc] = u
                end
            end
        end
    end
end

function shortest_path_from_mat_rev(g::VariationalGraph{T, U}, s_paths::Array{T, 1}, src::T, dist::T) where{T<:Integer, U<:Real}
    it = 0 # safty, if there is no path
    path = T[] # the path we want tp return
    cur_vert = dist # start point
    while  cur_vert != src && it <= g.num_verts
        if cur_vert == -1
            return T[]
        end
        push!(path, cur_vert)
        old_vert = cur_vert # stor previous vert
        cur_vert = s_paths[g.glob_to_loc[cur_vert]]
        it+=1
    end
    push!(path, src)
    return path
end

function shortest_path_from_mat(g::VariationalGraph{T, U}, s_paths::Array{T, 1}, src::T, dist::T) where{T<:Integer, U<:Real}
    return reverse!(shortest_path_from_mat_rev(g, s_paths, src, dist))
end

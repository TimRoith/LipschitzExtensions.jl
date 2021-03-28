
###################################################################################################
# type for kernel function
abstract type kernelType end # for epsilon Graphs
#------------------------------------------------
# all edges have the same weight
struct uniform{T<:Real} <: kernelType 
    weight::T
end
uniform() = uniform(1.0) # set default weight to 1
# repective weighting function
function weight_fct(p::Array{T, 1}, conf::uniform{T}) where {T<:Real}
    return conf.weight
end
###################################################################################################
# epsBall Graph has fixed ball of radius epsilon where every vertix in this ball is considered 
# a neighbor
struct epsBall{T<:Real, U<:kernelType} <: GraphConfig
    points::Array{T,2}
    epsilon::T
    kernel::U
end
# the respectiv constructor
function VariationalGraph(gconf::epsBall{T, U}) where{T<:Real, U<:kernelType}
    num_verts = size(gconf.points, 2)
    weights = [Vector{Float64}() for _ in 1:num_verts] # Initialize Array for Edgeweights
    kdtree = BallTree(gconf.points); #Generate the ballTree from given data
    edges = NearestNeighbors.inrange(kdtree, gconf.points, gconf.epsilon, false)
    num_edges = 0
    # Modify the obtained edges to fit in our structure
    for u = 1:num_verts
        setdiff!(edges[u],u) # delete the index of u itself in the inds array
        pt = gconf.points[:, u] # set edges and weights accordingly
        for v = 1:length(edges[u])
            # Set weights via outer and inner Norm and weight function
            push!(weights[u], weight_fct(pt - gconf.points[:,edges[u][v]], gconf.kernel))
            num_edges += 1
        end
    end
    g = VariationalGraph(num_verts, num_edges, edges, weights, gconf)
end
module InfinityGraphs

using Distances
using LightGraphs
using NearestNeighbors
using Printf
using SparseArrays
using DataStructures
using PyPlot

import LightGraphs:
DiGraph, AbstractGraph, add_edge!, rem_edge!, edges,
strongly_connected_components, connected_components

import NearestNeighbors:
KDTree, BallTree, BruteTree, knn, inrange

import Distances:
euclidean

###############################################################################

abstract type AbstractVariationalGraph{T<:Integer, U<:Real} end

###############################################################################

include("infinitygraph.jl")
include("constructgraph.jl")
include("minimize.jl")
include("visualization.jl")

end # module

include("setupProblem.jl")
println(string(line, "\n Finsihed preprocessing.\n", line))
v = InfinityGraphs.comp_inf_lap(g, v0)
#
println(string(line, "\n finsihed calc, starting plot\n", line))
#
PyPlot.close("all")
InfinityGraphs.plotgraph(g, v0, v)
#InfinityGraphs.comp_inf_min(g, v)
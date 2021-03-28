include("setupProblem.jl")
println(string(line, "\n Finsihed preprocessing.\n", line))

g, v0 = setUpProblem()
v = InfinityGraphs.comp_lex_min(g, v0)

PyPlot.close("all")
InfinityGraphs.plotgraph(g, v0, v)
alpha = InfinityGraphs.comp_max_grad(g, v)
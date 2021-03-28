using InfinityGraphs
using PyPlot
#
line = repeat("=",90)
println(string(line, "\n Finsihed loading packages. \n", line))
println(string(line, "\n Which setup do you want to use? \n"))
setup = readline()
println(string(line, "\n Constructing the graph. \n", line))
# Graph setup
function setUpProblem()
    if setup == "A"
        dim_domain = 2
        num_pts = 5000
        points = hcat(rand(dim_domain, num_pts), [0 1; 0.5 0.5])
        gconf = InfinityGraphs.epsBall(points, 0.05, InfinityGraphs.uniform())
        g = InfinityGraphs.VariationalGraph(gconf)
        v0 = fill(Inf, num_pts + 2)
        v0[end - 1] = 0
        v0[end] = 1
    elseif setup == "B"
        dim_domain = 2
        num_pts = 1000
        points = hcat(rand(dim_domain, num_pts), [0 1; 0.5 0.5])
        gconf = InfinityGraphs.epsBall(points, 0.4, InfinityGraphs.uniform())
        g = InfinityGraphs.VariationalGraph(gconf)
        v0 = fill(Inf, num_pts + 2)
        v0[end - 1] = 0
        v0[end] = 1
    elseif setup == "cust"
        println(string(line, "\n dim_domain:"))
        dim_domain = parse(Int64, readline())
        println(string(line, "\n num_pts:"))
        num_pts = parse(Int64, readline())
        points = hcat(rand(dim_domain, num_pts), [0 1; 0.5 0.5])
        println(string(line, "\n epsilon:"))
        epsilon = parse(Float64, readline())
        gconf = InfinityGraphs.epsBall(points, epsilon, InfinityGraphs.uniform())
        g = InfinityGraphs.VariationalGraph(gconf)
        v0 = fill(Inf, num_pts + 2)
        v0[end - 1] = 0
        v0[end] = 1
    end
    return g, v0
end
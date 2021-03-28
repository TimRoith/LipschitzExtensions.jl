using InfinityGraphs
using PyPlot
using GraphRecipes
using Plots



dim_domain = 2
num_pts = 40
num_x = floor(Int, sqrt(num_pts))
h = 1/num_x
# const points = rand(dim_domain, num_pts)


const points = [ i*h + j*h for i=1:num_x, j=1:num_x]

gconf = InfinityGraphs.discKernel(0.2)

g = InfinityGraphs.VariationalGraph(points, gconf)

# Visulaization
colors = ["b" "r" "c" "m" "k" "g"]
for u = 1:g.num_verts
    c = "0.7"
    plot(points[1,u], points[2,u], color = c, marker = "o")
    #annotate(string(u), (points[1,u], points[2,u]))
    for i = 1:length(g.edges_list[u])
        v = g.edges_list[u][i]
        plot((points[1, u], points[1, v]), (points[2, u], points[2, v]), color = c)
    end
end

# graphplot(g.edges_list)

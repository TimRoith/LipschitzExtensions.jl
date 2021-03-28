function plotgraph(g::VariationalGraph{T,U}, v0::Array{U, 1}, v::Array{U, 1}) where{T<:Integer, U<:Real}
    #
    PyPlot.using3D()
    fig = PyPlot.figure()
    ax = fig.gca(projection="3d")
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")
    #
    Tv0 = intersect(findall(x->x<Inf, v0), g.verts)
    Iv0 = intersect(findall(x->x==Inf, v0), g.verts)
    ax.scatter(g.config.points[1,Tv0], g.config.points[2,Tv0], v0[Tv0], marker="o", color="r")
    ax.scatter(g.config.points[1,Iv0], g.config.points[2,Iv0], v[Iv0], marker=".", color="b")
    # for u in 1:g.num_verts
    #     label = string(u)
    #     ax.text(points[1,u], points[2,u], v[u], label)
    #     #annotate(string(u), (points[1,u], points[2,u]))
    #     for i = 1:length(g.edges_list[u])
    #         vv = g.edges_list[u][i]
    #         ax.plot((points[1, g.verts[u]], points[1, vv]), (points[2, g.verts[u]], points[2, vv]), (v[g.verts[u]], v[vv]), color = "k")
    #     end
    # end
end
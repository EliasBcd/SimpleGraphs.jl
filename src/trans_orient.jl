
# This code by Tara Abrishami

export transitive_orientation

"""
`transitive_orientation(G)` finds a transitive orientation of the
simple graph `G`. The result is a `SimpleDigraph`. An error is raised
if `G` does not have a transitive orientation.
"""
function transitive_orientation(G::SimpleGraph)
  err_msg = "This graph does not have a transitive orientation"
  vertices = deepcopy(vlist(G))
  edges = deepcopy(elist(G))
  V = vertex_type(G)
  diredges = Tuple{V,V}[]
  D = SimpleDigraph{V}()
  while length(diredges) != 0 || length(edges) != 0
    if length(diredges) == 0
      e = shift!(edges)
      add!(D, e[1], e[2])
    else
      e = shift!(diredges)
    end
    v1 = e[1]
    v2 = e[2]
    vs = V[]
    append!(vs,G[v1])
    append!(vs,G[v2])
    for v in vs
      if has(G, v1, v) && !has(D, v1, v) && !has(G, v, v2)
        if has(D, v, v1)
          error(err_msg)
        end
        add!(D, v1, v)
        edg = (v1, v)
        swapedg = (v, v1)
        unshift!(diredges, edg)
        edges = filter(x -> x != edg && x != swapedg,edges)
      elseif has(G, v2, v) && !has(D, v, v2) && !has(G, v, v1)
        if has(D, v2, v)
          error(err_msg)
        end
        add!(D, v, v2)
        edg = (v, v2)
        swapedg = (v2, v)
        unshift!(diredges, edg)
        edges = filter(x -> x != edg && x != swapedg,edges)
      end
    end
  end
  if length(elist(D)) == length(elist(G))
    return D
  else
    error(err_msg)
  end
end

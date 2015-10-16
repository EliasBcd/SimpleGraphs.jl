# Various functions regarding graph connectivity

export components, num_components, is_connected, spanning_forest
export find_path, dist, diam, is_cut_edge

# Find the components of the graph as a Set of subsets of the vertices
function components{T}(G::SimpleGraph{T})
    VL = vlist(G)

    # create the partition of the vertex set inside a "DisjointSets"
    # Julia structure
    parts = DisjointSets{T}(VL)
    for e in G.E
        (u,v) = e
        union!(parts,u,v)
        if num_groups(parts)==1
            break
        end
    end
    return set_of_sets(parts)
end

# count the number of connected components
function num_components{T}(G::SimpleGraph{T})
    if NV(G)<2
        return NV(G)
    end
    parts = DisjointSets{T}(vlist(G))
    for e in G.E
        (u,v) = e
        union!(parts,u,v)
        if num_groups(parts)==1
            return 1
        end
    end
    return num_groups(parts)
end

# determine if the graph is connected
function is_connected{T}(G::SimpleGraph{T})
    return num_components(G) <= 1
end

# create a spanning forest of G, i.e., a maximal acyclic spanning
# subgraph. If G is connected, this is a tree.
function spanning_forest{T}(G::SimpleGraph{T})
    H = SimpleGraph{T}()
    if NV(G) == 0
        return H
    end

    for v in G.V
        add!(H,v)
    end

    parts = DisjointSets{T}(vlist(G))

    for e in G.E
        (u,v) = e
        if in_same_set(parts,u,v)
            continue
        end
        add!(H,u,v)
        union!(parts,u,v)
        if num_groups(parts)==1
            break
        end
    end
    return H
end

# find a shortest path between vertices
function find_path(G::AbstractSimpleGraph,s,t)
    T = vertex_type(G)
    if ~has(G,s) || ~has(G,t)
        error("Source and/or target vertex is not in this graph")
    end
    if s==t
        return [s]
    end

    # set up a queue for vertex exploration
    Q = Queue(T)
    enqueue!(Q,s)

    # set up trace-back dictionary
    tracer = Dict{T,T}()
    tracer[s] = s

    while length(Q) > 0
        v = dequeue!(Q)
        Nv = G.N[v]
        for w in Nv
            if haskey(tracer,w)
                continue
            end
            tracer[w] = v
            enqueue!(Q,w)

            if w==t  # success!
                path = Array(T,1)
                path[1] = t
                while path[1] != s
                    v = tracer[path[1]]
                    prepend!(path,[v])
                end
                return path

            end
        end
    end
    return T[]   # return empty array if no path found
end

# functions for finding distances between vertices in a graph. the
# distance between two vertices is the number of edges in a shortest
# path between them. if there is no such path, the distance is
# undefined (or infinite) but since we want the return value to be an
# Int, we use -1 to signal this.

# find the distance between specified vertices
function dist(G::AbstractSimpleGraph,u,v)
    return length(find_path(G,u,v))-1
end

# find all distances from a given vertex
function dist(G::AbstractSimpleGraph, v)
    T = vertex_type(G)
    d = Dict{T,Int}()
    if !has(G,v)
        error("Given vertex is not in this graph")
    end

    d[v] = 0
    Q = Queue(T)
    enqueue!(Q,v)

    while length(Q)>0
        w = dequeue!(Q)  # get 1st vertex in the queue
        Nw = G.N[w]
        for x in Nw
            if !haskey(d,x)
                d[x] = d[w]+1
                enqueue!(Q,x)
            end
        end
    end

    # record -1 for any vertices we missed
    for x in G.V
        if !haskey(d,x)
            d[x] = -1
        end
    end

    return d
end

# find all distances between all vertices
function dist(G::AbstractSimpleGraph)
    T = vertex_type(G)
    dd = Dict{Tuple{T,T},Int}()
    vtcs = vlist(G)

    for v in vtcs
        d = dist(G,v)
        for w in vtcs
            dd[(v,w)] = d[w]
        end
    end

    return dd
end

# Calculate the diameter of a graph, but return -1 if not connected.
function diam(G::SimpleGraph)
    if is_connected(G)
        return maximum(values(dist(G)))
    end
    return -1
end

# Determine if a given edge in a graph is a cut edge. If there is no
# such edge in the graph, an error is raised.
function is_cut_edge(G::SimpleGraph, u, v)
    if !has(G,u,v)
        error("No such edge in this graph")
    end

    delete!(G,u,v)
    P = find_path(G,u,v)
    result = false
    if length(P)==0
        result = true
    end
    add!(G,u,v)
    return result
end

# When called as is_cut_edge(G,e), we assume e is a tuple or list
# whose first two entries are the end points of the edge
function is_cut_edge(G::SimpleGraph, e)
    return is_cut_edge(G,e[1],e[2])
end

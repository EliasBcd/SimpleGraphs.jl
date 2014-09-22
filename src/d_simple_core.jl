# Core definitions for directed graphs


export SimpleDigraph, IntDigraph, StringDigraph
export add!, delete!, has
export is_looped, allow_loops!, forbid_loops!, remove_loops!, loops
export show, NV, NE
export out_deg, in_deg, deg
export in_neighbors, out_neighbors
export vlist, elist

type SimpleDigraph{T}
    V::Set{T}              # vertex set of this graph
    N::Dict{T,Set{T}}      # map vertices to out-neighbors
    NN::Dict{T,Set{T}}     # map vertices to in-neighbors
    looped::Bool           # flag to indicate if loops are allowed
    function SimpleDigraph()
        V = Set{T}()
        N = Dict{T,Set{T}}()
        NN = Dict{T,Set{T}}()
        G = new(V,N,NN,true)
    end
end

SimpleDigraph() = SimpleDigraph{Any}()
IntDigraph() = SimpleDigraph{Int}()
StringDigraph() = SimpleDigraph{ASCIIString}()

vertex_type{T}(G::SimpleDigraph{T}) = T


function IntDigraph(n::Int)
    G = IntDigraph()
    for v=1:n
        add!(G,v)
    end
    return G
end

# Simple way to print a digraph
function show(io::IO, G::SimpleDigraph)
    print(io,"$(typeof(G)) ($(NV(G)) vertices)")
end

# Do we allow loops?
is_looped(G::SimpleDigraph) = G.looped

# Grant permission for loops
function allow_loops!(G::SimpleDigraph)
    G.looped = true
    nothing
end

# Remove all loops from this digraph (but don't change loop
# permission)
function remove_loops!(G::SimpleDigraph)
    if !G.looped
        return nothing
    end
    for v in G.V
        delete!(G,v,v)
    end
    nothing
end

# Forbid loops (and delete any that we might have)
function forbid_loops!(G::SimpleDigraph)
    remove_loops!(G)
    G.looped = false
    nothing
end

# List all the loops in this digraph
function loops{T}(G::SimpleDigraph{T})
    if !is_looped(G)
        return T[]
    end
    loop_set = Set{T}()
    for v in G.V
        if has(G,v,v)
            push!(loop_set,v)
        end
    end
    loop_list = collect(loop_set)
    try
        sort!(loop_list)
    end
    return loop_list
end

# Out-degree of a vertex
out_deg(G::SimpleDigraph, v) = length(G.N[v])

# In-degree of a vertex
in_deg(G::SimpleDigraph, v) = length(G.NN[v])

# The degree of a vertex is the sum of in and out degrees
deg(G::SimpleDigraph, v) = in_deg(G,v) + out_deg(G,v)

# out neighbors of a vertex
function out_neighbors(G::SimpleDigraph, v)
    result = collect(G.N[v])
    try
        sort!(result)
    end
    return result
end

# in neighbors of a vertex
function in_neighbors(G::SimpleDigraph, v)
    result = collect(G.NN[v])
    try
        sort!(result)
    end
    return result
end

# Number of vertices
NV(G::SimpleDigraph) = length(G.V)

# Number of edges
function NE(G::SimpleDigraph)
    total::Int = 0
    for v in G.V
        total += out_deg(G,v)
    end
    return total
end

# Check if this digraph has a given vertex
has(G::SimpleDigraph, v) = in(v,G.V)

# Check if this digraph has a given edge
has(G::SimpleDigraph, v, w) = has(G,v) && in(w,G.N[v])

# Add a vertex to a digraph
function add!{T}(G::SimpleDigraph{T}, v)
    if has(G,v)
        return false
    end
    push!(G.V, v)
    G.N[v] = Set{T}()
    G.NN[v] = Set{T}()
    return true
end

# Add an edge to a digraph
function add!{T}(G::SimpleDigraph{T}, v, w)
    if !G.looped && v==w
        return false
    end
    if has(G,v,w)
        return false
    end
    if !has(G,v)
        add!(G,v)
    end
    if !has(G,w)
        add!(G,w)
    end
    push!(G.N[v],w)
    push!(G.NN[w],v)
    return true
end

# Delete an edge from this digraph
function delete!(G::SimpleDigraph, v, w)
    if !has(G,v,w)
        return false
    end
    delete!(G.N[v],w)
    delete!(G.NN[w],v)
    return true
end

# Delete a vertex from this digraph
function delete!(G::SimpleDigraph, v)
    if !has(G,v)
        return false
    end
    for w in G.N[v]
        delete!(G,v,w)
    end
    for u in G.NN[v]
        delete!(G,u,v)
    end
    delete!(G.V,v)
    delete!(G.N,v)
    delete!(G.NN,v)
    return true
end

# Create a list of all vertices in the digraph
function vlist(G::SimpleDigraph)
    result = collect(G.V)
    try
        sort!(result)
    end
    return result
end

# Create a list of all edges in the digraph
function elist{T}(G::SimpleDigraph{T})
    E = Set{(T,T)}()
    for v in G.V
        for w in G.N[v]
            push!(E, (v,w))
        end
    end
    result = collect(E)
    try
        sort!(result)
    end
    return result
end




import Base.show

export SimpleGraph, IntGraph, StringGraph
export show, NV, NE, has, vertex_type, fastN!
export vlist, elist, neighbors, getindex, deg


type SimpleGraph{T}
    V::Set{T}          # Vertex set
    E::Set{(T,T)}      # Edge set
    N::Dict{T,Set{T}}  # Optional neighbor sets
    Nflag::Bool        # Tells if N is used or not (default on)
    function SimpleGraph(Nflag::Bool=true)
        V = Set{T}()
        E = Set{(T,T)}()
        N = Dict{T,Set{T}}()
        G = new(V,E,N,Nflag)
    end
end

function show(io::IO, G::SimpleGraph)
    print(io,"$(typeof(G)) ($(NV(G)) vertices, $(NE(G)) edges)")
end

# Default constructor uses Any type vertices
SimpleGraph(Nflag::Bool = true) = SimpleGraph{Any}(Nflag)

# A StringGraph has vertices of type ASCIIString.
StringGraph() = SimpleGraph{ASCIIString}()

# Create a new StringGraph by reading a file.
# Each line of the file should contain one or two tokens separated by
# space. If the line contains a single token, we add that token as a
# vertex. If the line contains two (or more) tokens, then the first
# two tokens on the line are added as an edge (assuming they are
# different). Any extra tokens on the line are ignored.
# Also: If the line begins with a # character, the line is ignored.
function StringGraph(file::String)
    G = StringGraph()
    load!(G,file)
    return G
end

# Helper function for StrinGraph(file), and can be used to add
# vertices and edges to a graph (assuming its vertex type can
# accomodate strings).
function load!(G::SimpleGraph, file::String)
    f = open(file, "r")
    while(~eof(f))
        line = chomp(readline(f))
        tokens = split(line)

        if (length(tokens) == 0)
            continue
        end

        if (tokens[1][1] == '#')
            continue
        end

        add!(G,tokens[1])
        if length(tokens) > 1
            add!(G,tokens[1],tokens[2])
        end
    end
end

# An IntGraph has integer vertices. The basic form creates an empty
# SimpleGraph{int}.
IntGraph() = SimpleGraph{Int}()

# With a postive integer argument, adds 1:n as vertex set, but no
# edges.
function IntGraph(n::Int)
    G = IntGraph()
    for v=1:n
        add!(G,v)
    end
    return G
end

# Determine the type of vertex this graph holds
vertex_type{T}(G::SimpleGraph{T}) = T


# number of vertices and edges
NV(G::SimpleGraph) = length(G.V)
NE(G::SimpleGraph) = length(G.E)

# check for membership of a vertex or edge
has(G::SimpleGraph, v) = in(v,G.V)
has(G::SimpleGraph, v, w) = in((v,w), G.E) || in((w,v), G.E)

# fastN(G,true) creates an additional data structure to speed up
# neighborhood lookups.
function fastN!{T}(G::SimpleGraph{T},flg::Bool=true)
    # if no change, do nothing
    if flg == G.Nflag
        return flg
    end

    # if setting the flag to false, erase the G.N data structure;
    # otherwise build it.
    if flg == false
        G.N = Dict{T, Set{T}}()
        sizehint(G.N, NV(G))
    else
        # build the G.N structure.
        # start with empty sets for each vertex
        for v in G.V
            G.N[v] = Set{T}()
        end
        # now iterate over the edge set and populate G.N sets
        for e in G.E
            v,w = e
            push!(G.N[v],w)
            push!(G.N[w],v)
        end
    end
    G.Nflag = flg
    return flg
end

# Create a mapping between G.V and 1:n. This is not exposed outside
# this module; it's a helper function used by other functions. This
# has been crafted to work with either SimpleGraph or SimpleDigraph
# arguments.
function vertex2idx(G)
    T = vertex_type(G)
    d = Dict{T,Int}()
    V = vlist(G)
    n = NV(G)

    for k=1:n
        d[V[k]] = k
    end

    return d
end

# get the vertices as a (sorted if possible) list
function vlist(G::SimpleGraph)
    result = collect(G.V)
    try
        sort!(result)
    end
    return result
end

# get the edges as a (sorted if possible) list
function elist(G::SimpleGraph)
    result = collect(G.E)
    try
        sort!(result)
    end
    return result
end

# Get the neighbors of a vertex
function neighbors{T}(G::SimpleGraph{T}, v)
    if ~has(G,v)
        error("Graph does not contain requested vertex")
    end

    if G.Nflag
        N = collect(G.N[v])

    else
        N = T[]
        for w in G.V
            if has(G,v,w)
                append!(N,[w])
            end
        end
    end
    try
        sort!(N)
    end
    return N
end

# Here is another way to access the neighbors of a vertex: G[v]
getindex(G::SimpleGraph,v) = neighbors(G,v)

# And here's a getindex way to check for edges: G[u,v] is a shortcut
# for has(G,u,v).
getindex(G::SimpleGraph,v,w) = has(G,v,w)

# Degree of a vertex
function deg(G::SimpleGraph, v)
    if ~has(G,v)
        error("Graph does not contain requested vertex")
    end
    if G.Nflag
        return length(G.N[v])
    end
    return length(G[v])
end

# Degree sequence
function deg{T}(G::SimpleGraph{T})
    if G.Nflag
        ds = [deg(G,v) for v in G.V]
    else
        dd = Dict{T,Int}()
        for v in G.V
            dd[v] = 0
        end
        for e in G.E
            v,w = e
            dd[v] += 1
            dd[w] += 1
        end
        ds = collect(values(dd))
    end
    sort!(ds, lt = >)
    return ds
end

# Module written by Ed Scheinerman, ers@jhu.edu
# distributed under terms of the MIT license

module SimpleGraphs
using DataStructures


import Base.show, Base.delete!, Base.isequal, Base.union
import Base.complement, Base.complement!


export SimpleGraph, IntGraph, StringGraph
export show, NV, NE, has, vertex_type, fastN!
export vlist, elist, neighbors, getindex, deg
export add!, delete!, contract!
export adjacency, laplace, incidence, copy, induce

export Complete, Path, Cycle, RandomGraph, line_graph
export RandomTree, code_to_tree
export Grid, Wheel, Cube, BuckyBall

export complement, complement!, ctranspose
export Petersen, Kneser, subsets
export cartesian, relabel, trim
export disjoint_union, union, join
export euler

export components, num_components, is_connected, spanning_forest
export find_path, dist, diam, dist_matrix, is_cut_edge
export bipartition, two_color, greedy_color, random_greedy_color

export export_simple

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

# check for graph equality
function isequal(G::SimpleGraph, H::SimpleGraph)
    if G.V != H.V || NE(G) != NE(H)
        return false
    end
    try  # check vertex sets for sortability
        sort(collect(G.V))
        sort(collect(H.V))
    catch  #if not sortable, resort to a slow method
        for e in G.E
            if ! has(H,e[1],e[2])
                return false
            end
        end
        return true
    end
    # otherwise faster to check edge sets this way
    return G.E==H.E
end

function ==(G::SimpleGraph, H::SimpleGraph)
    return isequal(G,H)
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
    sizehint(G.V,n)  # speed up this way by pre-allocating
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

# adding vertices
function add!{T}(G::SimpleGraph{T}, v::T)
    if has(G,v)
        return false
    end
    push!(G.V, v)
    if G.Nflag
        G.N[v] = Set{T}()
    end
    return true
end

# adding edges
function add!{T}(G::SimpleGraph{T}, v::T, w::T)
    if v==w
        return false
    end

    try
        if v > w
            v,w = w,v
        end
    end

    if ~has(G,v)
        add!(G,v)
    end
    if ~has(G,w)
        add!(G,w)
    end
    if has(G,v,w)
        return false
    end
    push!(G.E, (v,w))
    if G.Nflag
        push!(G.N[v], w)
        push!(G.N[w], v)
    end
    return true
end

# Create a mapping between G.V and 1:n. This is not exposed outside
# this module; it's a helper function used by other functions.
function vertex2idx{T}(G::SimpleGraph{T})
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
function neighbors{T}(G::SimpleGraph{T}, v::T)
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

# Edge deletion
function delete!(G::SimpleGraph, v, w)
    flag = false
    if has(G,v,w)
        flag = true
        delete!(G.E,(v,w))
        delete!(G.E,(w,v))
        if G.Nflag
            delete!(G.N[v],w)
            delete!(G.N[w],v)
        end
    end
    return flag
end

# Vertex deletion
function delete!(G::SimpleGraph, v)
    flag = false
    if has(G,v)
        flag = true
        Nv = G[v]
        for w in Nv
            delete!(G,v,w)
        end
        delete!(G.V, v)
        if G.Nflag
            delete!(G.N,v)
        end
    end
    return flag
end

# Contract an edge in a graph. If uv is an edge of G, we add all of
# v's neighbors to u's neighborhood and then delete v. This mutates
# the graph. Usually, vertices u and v are adjacent, but we don't
# require that. If either u or v is not a vertex of this graph, or if
# u==v, we return false. Otherwise we return true to indicate success.
function contract!(G::SimpleGraph, u, v)
    if !has(G,u) || !has(G,v) || u==v
        return false
    end

    Nv = G[v]
    for x in Nv
        add!(G,u,x)
    end
    delete!(G,v)
    return true
end

# Create an independent copy of a graph
function copy{T}(G::SimpleGraph{T})
    H = SimpleGraph{T}()
    for v in G.V
        add!(H,v)
    end
    for e in G.E
        add!(H,e[1],e[2])
    end
    return H
end

# Adjaceny Matrix
function adjacency(G::SimpleGraph)
    n = NV(G)
    A = zeros(Int,(n,n))

    # create a table from V to 1:n
    d = vertex2idx(G)

    for e in G.E
        i = d[e[1]]
        j = d[e[2]]
        A[i,j]=1
        A[j,i]=1
    end

    return A
end

# Laplace matrix
function laplace(G::SimpleGraph)
    A = adjacency(G)
    d = collect(sum(A,1))
    D = Base.diagm(d)
    L = D-A
    return L
end

# incidence matrix
function incidence(G::SimpleGraph, signed::Bool = true)
    n = NV(G)
    m = NE(G)
    M = spzeros(Int,n,m)
    d = vertex2idx(G)
    E = elist(G)
    a = 1
    b = signed ? -1 : 1

    idx = 0
    for e in E
        i = d[e[1]]
        j = d[e[2]]
        idx += 1
        M[i,idx] = a
        M[j,idx] = b
    end

    return M
end

# Given a simple graph G and a set of vertices A, form the induced
# subgraph G[A]. Note that A must be a subset of V(G).
function induce{T}(G::SimpleGraph{T}, A::Set{T})
    # Check that A is a subset of V(G)
    for v in A
        if ~has(G,v)
            error("The set A must be a subset of V(G)")
        end
    end
    H = SimpleGraph{T}()  # place to hold the answer

    # add all the vertices in A to H
    for v in A
        add!(H,v)
    end

    # The method we choose depends on the size of A. For small A, best
    # to iterate over pairs of elements of A. For large A, best to
    # iterate over G.E
    a = length(A)

    if a<2
        return H
    end

    # case of small A, iterate over pairs from A
    if a*(a-1) < 2*NE(G)
        alist = collect(A)
        for i=1:a-1
            u = alist[i]
            for j=i+1:a
                v = alist[j]
                if has(G,u,v)
                    add!(H,u,v)
                end
            end
        end

    else  # case of large A, iterate over G.E
        for e in G.E
            if in(e[1],A) && in(e[2],A)
                add!(H,e[1],e[2])
            end
        end
    end
    return H
end

# Create a complete graph
function Complete(n::Int)
    G = IntGraph(n)

    for k=1:n-1
        for j=k+1:n
            add!(G,j,k)
        end
    end
    return G
end

# Create a complete bipartite graph
function Complete(n::Int, m::Int)
    G = IntGraph(n+m)
    for u=1:n
        for v=n+1:n+m
            add!(G,u,v)
        end
    end
    return G
end

# Create the complete multipartite graph with given part sizes
function Complete(parts::Array{Int,1})
    # check all part sizes are positive
    for p in parts
        if p < 1
            error("All part sizes must be positive")
        end
    end

    n = sum(parts)
    G = IntGraph(n)

    np = length(parts)
    if np < 2
        return G
    end

    # create table of part ranges
    ranges = Array(Int,np,2)
    ranges[1,1] = 1
    ranges[1,2] = parts[1]
    for k = 2:np
        ranges[k,1] = ranges[k-1,2] + 1
        ranges[k,2] = ranges[k,1] + parts[k] - 1
    end

    # Add all edges between all parts
    for i=1:np-1
        for j=i+1:np
            for u=ranges[i,1]:ranges[i,2]
                for v=ranges[j,1]:ranges[j,2]
                    add!(G,u,v)
                end
            end
        end
    end

    return G
end

# Create a path graph on n vertices
function Path(n::Int)
    G = IntGraph(n)
    for v = 1:n-1
        add!(G,v,v+1)
    end
    return G
end

# Create a path graph from a list of vertices
function Path{T}(verts::Array{T})
    G = SimpleGraph{T}()
    n = length(verts)

    if n==1
        add!(G,verts[1])
    end
    for k = 1:n-1
        add!(G,verts[k],verts[k+1])
    end
    return G
end

# Create a cycle graph on n vertices
function Cycle(n::Int)
    if n<3
        error("Cycle requires 3 or more vertices")
    end
    G = Path(n)
    add!(G,1,n)
    return G
end

# Create the wheel graph on n vertices: a cycle on n-1 vertices plus
# an additional vertex adjacent to all the vertices on the wheel.
function Wheel(n::Int)
    if n < 4
        error("Wheel graphs must have at least 4 vertices")
    end
    G = Cycle(n-1)
    for k=1:n-1
        add!(G,k,n)
    end
    return G
end

# Create a grid graph
function Grid(n::Int, m::Int)
    G = SimpleGraph{(Int,Int)}()

    # add the vertices
    for u=1:n
        for v=1:m
            add!(G,(u,v))
        end
    end

    #horizontal edges
    for u=1:n
        for v=1:m-1
            add!(G,(u,v),(u,v+1))
        end
    end

    # vertical edges
    for v=1:m
        for u=1:n-1
            add!(G,(u,v),(u+1,v))
        end
    end
    return G
end

# Create an Erdos-Renyi random graph
function RandomGraph(n::Int, p::Real=0.5)
    G = IntGraph(n)
    for v=1:n-1
        for w=v+1:n
            if (rand() < p)
                add!(G,v,w)
            end
        end
    end
    return G
end

# Generate a random tree on vertex set 1:n. All n^(n-2) trees are
# equally likely.
function RandomTree(n::Int)

    if n<0   # but we allow n==0 to give empty graph
        error("Number of vertices cannot be negative")
    end

    if n<2
        return IntGraph(n)
    end

    code = [ mod(rand(Int),n)+1 for _ in 1:n-2 ]
    return code_to_tree(code)
end

# This is a helper function for RandomTree that converts a Prufer code
# to a tree. No checks are done on the array fed into this function.
function code_to_tree(code::Array{Int,1})
    n = length(code)+2
    G = IntGraph(n)
    degree = ones(Int,n)  # initially all 1s

    #every time a vertex appears in code[], up its degree by 1
    for c in code
        degree[c]+=1
    end

    for u in code
        for v in 1:n
            if degree[v]==1
                add!(G,u,v)
                degree[u] -= 1
                degree[v] -= 1
                break
            end
        end
    end

    last = find(degree)
    add!(G,last[1],last[2])

    return G
end

# Create the line graph of a given graph
function line_graph{T}(G::SimpleGraph{T})
    H = SimpleGraph{(T,T)}()

    m = NE(G)
    E = elist(G)
    for e in E
        add!(H,e)
    end

    for i=1:m-1
        e = E[i]
        for j=i+1:m
            f = E[j]
            if e[1]==f[1] || e[2]==f[1] || e[1]==f[2] || e[2]==f[2]
                add!(H,e,f)
            end
        end
    end
    return H
end

# Create the complement of a graph.
function complement{T}(G::SimpleGraph{T})
    H = SimpleGraph{T}()
    V = vlist(G)
    n = NV(G)

    for i=1:n
        add!(H,V[i])
    end

    for i=1:n-1
        v = V[i]
        for j=i+1:n
            w = V[j]
            if ~has(G,v,w)
                add!(H,v,w)
            end
        end
    end
    return H
end

# We can use G' to mean complement(G)
ctranspose(G::SimpleGraph) = complement(G)

# complement in place. Returns None.
function complement!(G::SimpleGraph)
    n = NV(G)
    V = vlist(G)
    for i=1:n-1
        v = V[i]
        for j=i+1:n
            w = V[j]
            if has(G,v,w)
                delete!(G,v,w)
            else
                add!(G,v,w)
            end
        end
    end
    return None
end

# Create the Cube graph with 2^n vertices
function Cube(n)
    G = StringGraph()
    for u=0:2^n-1
        for shift=0:n-1
            v = (1<<shift) $ u
            add!(G,bin(u,n), bin(v,n))
        end
    end
    return G
end

# Create the BuckyBall graph
function BuckyBall()
    G = IntGraph()
    edges = [(1,3), (1,49), (1,60), (2,4), (2,10), (2,59),
	     (3,4), (3,37), (4,18), (5,7), (5,9), (5,13),
	     (6,8), (6,10), (6,17), (7,8), (7,21), (8,22),
	     (9,10), (9,57), (11,12), (11,13), (11,21), (12,28),
	     (12,48), (13,14), (14,47), (14,55), (15,16), (15,17),
	     (15,22), (16,26), (16,42), (17,18), (18,41), (19,20),
	     (19,21), (19,27), (20,22), (20,25), (23,24), (23,32),
	     (23,35), (24,26), (24,39), (25,26), (25,31), (27,28),
	     (27,31), (28,30), (29,30), (29,32), (29,36), (30,45),
	     (31,32), (33,35), (33,40), (33,51), (34,36), (34,46),
	     (34,52), (35,36), (37,38), (37,41), (38,40), (38,53),
	     (39,40), (39,42), (41,42), (43,44), (43,47), (43,56),
	     (44,46), (44,54), (45,46), (45,48), (47,48), (49,50),
	     (49,53), (50,54), (50,58), (51,52), (51,53), (52,54),
	     (55,56), (55,57), (56,58), (57,59), (58,60), (59,60),
	     ]
    for e in edges
        add!(G,e[1],e[2])
    end
    return G
end

# Create the cartesian product of two graphs
function cartesian{S,T}(G::SimpleGraph{S}, H::SimpleGraph{T})
    K = SimpleGraph{(S,T)}()
    for v in G.V
        for w in H.V
            add!(K,(v,w))
        end
    end

    for e in G.E
        (u,v) = e
        for w in H.V
            add!(K, (u,w), (v,w))
        end
    end

    for e in H.E
        (u,v) = e
        for w in G.V
            add!(K, (w,u), (w,v))
        end
    end

    return K
end

# Use G*H for Cartesian product
function *{S,T}(G::SimpleGraph{S},H::SimpleGraph{T})
    return cartesian(G,H)
end

# Unexposed helper function for subsets() function
function array2set{T}(A::Array{T,1})
    S = Set{T}()
    for a in A
        push!(S,a)
    end
    return S
end

# Create the set of all subsets of size k of a given set
function subsets{T}(A::Set{T}, k::Int)
    # create list of lists
    L = Base.combinations(collect(A),k)
    B = Set{Set{T}}()
    for x in L
        push!(B, array2set(x))
    end
    return B
end

# The Kneser graph Kneser(n,k) has C(n,k) vertices that are the
# k-element subsets of 1:n in which two vertices are adjacent if (as
# sets) they are disjoint. The Petersen graph is Kneser(5,2).
function Kneser(n::Int,k::Int)
    A = array2set(collect(1:n))
    vtcs = collect(subsets(A,k))
    G = SimpleGraph{Set{Int}}()

    for v in vtcs
        add!(G,v)
    end

    n = length(vtcs)
    for i=1:n-1
        u = vtcs[i]
        for j=i+1:n
            v = vtcs[j]
            if length(intersect(u,v))==0
                add!(G,u,v)
            end
        end
    end

    return G
end

# Create the Petersen graph.
Petersen() = Kneser(5,2)

# Relabel the vertics of a graph based on a dictionary mapping old
# vertex names to new
function relabel{S,T}(G::SimpleGraph{S}, label::Dict{S,T})
    H = SimpleGraph{T}()
    for v in G.V
        add!(H,label[v])
    end

    for e in G.E
        u = label[e[1]]
        v = label[e[2]]
        add!(H,u,v)
    end
    return H
end

# Relabel the vertices with the integers 1:n
function relabel{S}(G::SimpleGraph{S})
    verts = vlist(G)
    n = length(verts)
    label = Dict{S,Int}()
    sizehint(label,n)

    for idx = 1:n
        label[verts[idx]] = idx
    end

    return relabel(G,label)
end

# This is an unexposed helper function that takes a graph and creates
# a new graph in which the name of each vertex has an integer
# appended. For example, if the vertex type is String in the original
# graph, the new vertices are type (String, Int).
function label_append{S}(G::SimpleGraph{S}, a::Int)
    mapper = Dict{S,(S,Int)}()
    for v in G.V
        mapper[v] = (v,a)
    end
    return relabel(G,mapper)
end

# The disjoint_union of two graphs, G and H, is a new graph consisting of
# disjoint copies of G and H.
function disjoint_union{S,T}(G::SimpleGraph{S}, H::SimpleGraph{T})
    GG = label_append(G,1)
    HH = label_append(H,2)
    if S==T
        K = SimpleGraph{(S,Int)}()
    else
        K = SimpleGraph{Any}()
    end

    for v in GG.V
        add!(K,v)
    end
    for v in HH.V
        add!(K,v)
    end

    for e in GG.E
        add!(K,e[1],e[2])
    end
    for e in HH.E
        add!(K,e[1],e[2])
    end
    return K
end

# The join of two graphs is formed by taking disjoint copies of the
# graphs and all possible edges between the two.
function join{S,T}(G::SimpleGraph{S}, H::SimpleGraph{T})
    K = disjoint_union(G,H)
    for v in G.V
        for w in H.V
            add!(K, (v,1), (w,2) )
        end
    end
    return K
end

# Create the union of two graphs. If they have vertices or edges in
# common, that's OK.
function union{S,T}(G::SimpleGraph{S}, H::SimpleGraph{T})
    if S==T
        K = SimpleGraph{S}()
    else
        K = SimpleGraph{Any}()
    end

    for v in G.V
        add!(K,v)
    end
    for v in H.V
        add!(K,v)
    end

    for e in G.E
        add!(K,e[1],e[2])
    end
    for e in H.E
        add!(K,e[1],e[2])
    end
    return K
end

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

# repeatedly remove vertices with the given degree or less
function trim(G::SimpleGraph, d::Int = 0)
    H = copy(G)
    while NV(H) > 0 && minimum(deg(H)) <= d
        for v in H.V
            if deg(H,v) <= d
                delete!(H,v)
            end
        end
    end
    return H
end

# find a shortest path between vertices
function find_path{T}(G::SimpleGraph{T},s,t)
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
        Nv = G[v]
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
function dist{T}(G::SimpleGraph{T},u::T,v::T)
    return length(find_path(G,u,v))-1
end

# find all distances from a given vertex
function dist{T}(G::SimpleGraph{T}, v::T)
    d = Dict{T,Int}()
    if !has(G,v)
        error("Given vertex is not in this graph")
    end

    d[v] = 0
    Q = Queue(T)
    enqueue!(Q,v)

    while length(Q)>0
        w = dequeue!(Q)  # get 1st vertex in the queue
        Nw = G[w]
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
function dist{T}(G::SimpleGraph{T})
    dd = Dict{(T,T),Int}()
    vtcs = vlist(G)

    for v in vtcs
        d = dist(G,v)
        for w in vtcs
            dd[(v,w)] = d[w]
        end
    end

    return dd
end

# Create the n-by-n distance matrix
function dist_matrix(G::SimpleGraph)
    vtcs = vlist(G)
    n = length(vtcs)
    dd = dist(G)

    A = zeros(Int,n,n)

    for i = 1:n
        u = vtcs[i]
        for j = 1:n
            v = vtcs[j]
            A[i,j] = dd[(u,v)]
        end
    end

    return A
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
function is_cut_edge{T}(G::SimpleGraph{T}, u::T, v::T)
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

function is_cut_edge{T}(G::SimpleGraph{T}, e::(T,T))
    return is_cut_edge(G,e[1],e[2])
end

# This function finds an Euler trail/tour given initial
# vertices. Returns the trail if it exists of [] if it does not.

function euler{T}(G::SimpleGraph{T}, u::T, v::T)
    notrail = T[]

    # perform basic checks:
    if ! (has(G,u) && has(G,v))
        error("One or both of these vertices is not in this graph")
    end

    # special case: if all vertices have degree zero
    if all( [deg(G,x)==0 for x in G.V] )
        if u==v
            return [u]
        else
            return notrail
        end
    end

    # if either vertex has degree zero, we're doomed
    if deg(G,u)==0 || deg(G,v)==0
        return notrail
    end

    # vertex degree checks
    if u==v
        for x in G.V
            if deg(G,x)%2 == 1
                return notrail
            end
        end
    else        # u != v
        for x in G.V
            if x==u || x==v
                if deg(G,x)%2==0
                    return notrail
                end
            else
                if deg(G,x)%2==1
                    return notrail
                end
            end
        end
    end

    # Remove isolates and check for connectivity
    GG = trim(G)
    if !is_connected(GG)
        return notrail
    end

    # all tests have been satisfied. Now find the trail in GG using a
    # helper function.
    return euler_work!(GG,u)
end

# special case: find an Euler tour from a specified vertex
function euler{T}(G::SimpleGraph{T},u::T)
    return euler(G,u,u)
end

# special case: find any Euler tour. If the graph is connected, any
# vertex will do but if there are isolated vertices, we don't want to
# pick one of those!
function euler{T}(G::SimpleGraph{T})
    if NV(G)==0
        return T[]
    end

    verts = vlist(G)

    # search for a vertex that isn't isolated
    for u in verts
        if deg(G,u) > 0
            return euler(G,u,u)
        end
    end

    # if we reach here, the graph only has isolated vertics. Let it
    # work from any isolated vertex (and likely fail).
    u = verts[1]
    return euler(G,u,u)

end

# private helper function for euler()
function euler_work!{T}(G::SimpleGraph{T}, u::T)
    trail = T[]

    while true
        # if last vertex
        if NV(G)==1
            append!(trail,[u])
            return trail
        end

        # get the neighbors of u
        Nu = G[u]

        ## println("Current vertex is ", u, " with neighbors ", Nu)
        ## println("Trail so far:\t", trail')
        ## println()
        ## sleep(0.5)

        # if only one neighbor delete and move on
        if length(Nu)==1
            w = Nu[1]
            delete!(G,u)
            append!(trail,[u])
            u = w
        else
            for w in Nu
                if ! is_cut_edge(G,u,w)
                    delete!(G,u,w)
                    append!(trail,[u])
                    u = w
                    break
                end
            end
        end
    end
    error("This can't happen")
end

# Create a two-coloring of a graph or die trying. Returns a map from
# the vertex set to the set {1,2}, or error if no such mapping exists.
function two_color{T}(G::SimpleGraph{T})
    f = Dict{T,Int}()
    for A in components(G)
        a = first(A)
        f[a] = 1
        Q = Deque{T}()
        push!(Q,a)
        while length(Q)>0
            v = pop!(Q)
            Nv = G[v]
            for w in Nv
                if haskey(f,w)
                    if f[w]==f[v]
                        error("Graph is not bipartite")
                    end
                else
                    f[w] = 3-f[v]
                    push!(Q,w)
                end
            end
        end
    end
    return f
end

# Create a bipartition of a graph or die trying. Returns a set {X,Y}
# that is a bipartition of the vertex set of G.
function bipartition{T}(G::SimpleGraph{T})
    f::Dict{T,Int} = two_color(G)
    X = Set{T}()
    Y = Set{T}()
    for v in G.V
        if f[v]==1
            push!(X,v)
        else
            push!(Y,v)
        end
    end

    B = Set{Set{T}}()
    push!(B,X)
    push!(B,Y)
    return B
end


# Color a graph by the greedy algorithm in the sequence specified by
# seq. The array seq must be a permutation of G.V. We don't check
# that's true!
function greedy_color{T}(G::SimpleGraph{T}, seq::Array{T,1})
    f = Dict{T,Int}()  # this is the mapping from V to colors
    maxf::Int = 0      # largest color used

    for v in seq
        colors_used = falses(maxf)  # array if color is used by N[v]
        for w in G[v]
            if haskey(f,w)  # w already colored
                colors_used[f[w]]=true  # mark that color is used
            end
        end
        # give first unused color to v
        for k in 1:maxf
            if colors_used[k] == false
                f[v] = k
                break
            end
        end
        # but if that fails, extend the number of colors available
        if !haskey(f,v)
            maxf += 1
            f[v] = maxf
        end
    end
    return f
end

# This function returns a list of the vertices of G in descending
# order by degree. The order of vertices of the same degree is
# indeterminate. NOTE: This is not exported from this module. Should
# it be?
function deg_sorted_vlist(G::SimpleGraph)
    bye = x -> -x[1]
    list = [ (deg(G,v) , v) for v in G.V ]
    sort!(list, by=bye)
    outlist = [ item[2] for item in list ]
    return outlist
end

# Apply greedy_color to the graph visiting the vertices in decreasing
# order of degree.
function greedy_color{T}(G::SimpleGraph{T})
    seq = deg_sorted_vlist(G)
    return greedy_color(G,seq)
end

# Generate multiple random orders of the vertex set and apply
# greedy_color; return one that uses the fewest colors. This do as
# well as or better than greedy_color on some decreasing order of
# degree.
function random_greedy_color{T}(G::SimpleGraph{T}, reps::Int=1)
    n = NV(G)
    bestf = greedy_color(G)  # degree order default start
    best  = maximum(values(bestf))
    seq = vlist(G)
    println("Initial coloring uses ", best, " colors")

    for k in 1:reps
        shuffle!(seq)
        f = greedy_color(G,seq)
        mx = maximum(values(f))
        if mx < best
            bestf = f
            best = mx
            println("Reduced to ", best, " colors")
        end
    end
    return bestf
end



######################################################################
#
# Convert a graph into a "simple_graph" (from Graphs.jl). We return a
# triple: the simple_graph H with vertex set 1:n, a dictionary d mapping
# the vertices of G to the vertices of H, and the inverse dictionary
# dinv from V(H) to V(G).
import Graphs.simple_graph, Graphs.add_edge!
function export_simple{T}(G::SimpleGraph{T})
    n = NV(G)

    d = vertex2idx(G)
    dinv = Dict{Int,T}()
    for k in keys(d)
        v = d[k]
        dinv[v] = k
    end

    H = simple_graph(n,is_directed=false)
    for e in G.E
        u = d[e[1]]
        v = d[e[2]]
        add_edge!(H,u,v)
    end
    return (H,d,dinv)
end


#########################################################
# These are helper functions for the DisjointSets type. #
#########################################################

# This takes a DisjointSets object and returns its ground set.
function ground_set{T}(DS::DisjointSets{T})
    G = Set{T}()
    for item in keys(DS.intmap)
        push!(G,item)
    end
    return G
end

# Create a set of sets containing the parts of this partition.
function set_of_sets{T}(DS::DisjointSets{T})
    n = num_groups(DS)
    GS = ground_set(DS)

    # Get root elements
    roots = Set{Int}()
    for item in GS
        r = find_root(DS,item)
        push!(roots,r)
    end
    roots = collect(roots)

    # Map root numbers to [1:n]
    rootmap = Dict{Int,Int}()
    for k=1:n
        rootmap[roots[k]] = k
    end

    # Create an array of sets to hold the parts
    parts = Array(Set{T},n)
    for k=1:n
        parts[k] = Set{T}()
    end

    # place each ground set item in its appropriate part
    for item in GS
        index = rootmap[find_root(DS,item)]
        push!(parts[index], item)
    end

    # now take those parts and pack them into a set
    P = Set{Set{T}}() # This is a set of sets
    for item in parts
        push!(P,item)
    end

    return P
end



end # module SimpleGraphs

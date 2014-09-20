import Base.isequal
import Base.delete!, Base.union
import Base.complement, Base.complement!
import Base.join, Base.copy 

export add!, delete!, contract!, copy, induce
export line_graph, complement, complement!, ctranspose
export cartesian, relabel, trim
export disjoint_union, union, join

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

# edge deletion
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

# vertex deletion
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

# Repeatedly remove vertices with the given degree or less until there
# are no such vertices remaining in the graph. The default trim(G)
# simply removes all isolated vertices.
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

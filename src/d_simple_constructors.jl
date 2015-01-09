export DirectedPath, DirectedCycle, DirectedComplete
export RandomDigraph, RandomTournament

# Create a directed path on n vertices
function DirectedPath(n::Int)
    if n < 1
        error("n must be positive")
    end
    G = IntDigraph(n)
    for u=1:(n-1)
        add!(G,u,u+1)
    end
    return G
end

# Create a directed cycle on n vertices
function DirectedCycle(n::Int)
    G = DirectedPath(n)
    add!(G,n,1)
    return G
end

# Create a complete digraph (all possible edges)
function DirectedComplete(n::Int, with_loops::Bool=true)
    G = IntDigraph(n)
    if !with_loops
        forbid_loops!(G)
    end
    for u=1:n
        for v=1:n
            add!(G,u,v)
        end
    end
    return G
end

# Create a random digraph (Erdos-Renyi style)
function RandomDigraph(n::Int, p::Real=0.5, with_loops=false)
    G = IntDigraph(n)
    if !with_loops
        forbid_loops!(G)
    end
    for u=1:n
        for v=1:n
            if rand() < p
                add!(G,u,v)
            end
        end
    end
    return G
end

# Create a random tournament (no loops!)
function RandomTournament(n::Int)
    G = IntDigraph()
    for u=1:n-1
        for v=u+1:n
            if rand() < 0.5
                add!(G,u,v)
            else
                add!(G,v,u)
            end
        end
    end
    return G
end

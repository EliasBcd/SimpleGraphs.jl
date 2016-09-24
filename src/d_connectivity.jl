export stronglyconnectedcomponents, isacyclic
export transitiveclosure!
export maxcycles, simplecycles, getcycles, countcycles, itercycles

"""
Supertype Visitor. (future proof)
"""
abstract Visitor

###############################################################################
### Strong connectivity

"""
```type TarjanVisitor{T} <: Visitor
    indexmap::Dict{T, Int}
    stack::Vector{T}
    index::Int
    lowlink::Dict{T, Int}
    onstack::Dict{T, Bool}
    visited::Dict{T, Bool}
    components::Vector{Vector{T}}
end```

Concrete type used for the Tarjan's algorithm for strongly connected components.

# Arguments:
* indexmap: Save the order visit
* stack: The vertices that have been visited, but are not yet part of a strongly connected
 component
* index: Number keeping the iteration we are in
* lowlink: keep the lowest index reachable from each node
* onstack: whether the each vertex is on the stack or not
* visited: whether  the vertex has been visited or not (avoid unnecessary explorations)
* components: the strongly connected components.
"""
type TarjanVisitor{T} <: Visitor
    indexmap::Dict{T, Int}
    stack::Vector{T}
    index::Int
    lowlink::Dict{T, Int}
    onstack::Dict{T, Bool}
    visited::Dict{T, Bool}
    components::Vector{Vector{T}}
end

"""
Constructor of the above type adapted to the algorithm.
"""
TarjanVisitor{T}(dg::SimpleDigraph{T}) = TarjanVisitor{T}(
Dict([v => 0 for v in dg.V]),
Vector{T}(),
1,
Dict([v => 0 for v in dg.V]),
Dict([v => false for v in dg.V]),
Dict([v => false for v in dg.V]),
Vector{Vector{T}}()
)

"""
```components{T}(source::T, vis::TarjanVisitor{T})```

Get the strongly connected components from the stack and remove necessary nodes.

# Arguments
* source: the source vertex of the strongly connected component
* visitor: the data of the visit (including the stack)

# Return
* scc: the strongly connected component
"""
function components{T}(source::T, vis::TarjanVisitor{T})
    scc = Vector{T}()
    w = pop!(vis.stack)
    while w != source
        push!(scc, w)
        vis.onstack[w] = false
        w = pop!(vis.stack)
    end
    push!(scc, w)
    vis.onstack[w] = false
    return scc
end

"""
```visit!{T}(v::T, vis::TarjanVisitor{T})```

Update the visitor with the proper values at each visit.

# Arguments
* v: the current vertex.
* vis: the information needed to get the stronglyconnected component

# Return
* vis.index +1 : update the indexing for the next node visited.
"""
function visit!{T}(v::T, vis::TarjanVisitor{T})
    vis.indexmap[v] = vis.index
    vis.visited[v] = true
    vis.lowlink[v] = vis.index
    push!(vis.stack, v)
    vis.onstack[v] = true
    return vis.index + 1
end

"""
```strongconnect!{T}(dg::SimpleDigraph{T}, v::T, visitor::TarjanVisitor{T})```

One visit in Tarjan's algorithm.

# Arguments
* dg: the digraph we wish to decompose in strongly connected components
* v: the vertex currently visited
* visitor: the information needed to get the stronglyconnected component

# Return
* visitor.index: the iteration we have finished.
"""
function strongconnect!{T}(dg::SimpleDigraph{T}, v::T, visitor::TarjanVisitor{T})
   #Set the depth index for v to the smallest unused index
    visitor.index = visit!(v, visitor)

    # Consider successors of v
    for w in dg.N[v]
        if !visitor.visited[w]
        # Successor w has not yet been visited; recurse on it
            visitor.index = strongconnect!(dg, w, visitor)
            visitor.lowlink[v]  = min(visitor.lowlink[v], visitor.lowlink[w])
        elseif visitor.onstack[w]
        # Successor w is in stack S and hence in the current SCC
            visitor.lowlink[v] = min(visitor.lowlink[v], visitor.indexmap[w])
        end
    end 
    # If v is a root node, pop the stack and generate an SCC
    if visitor.lowlink[v] == visitor.indexmap[v]
        push!(visitor.components, components(v, visitor))
    end
    return visitor.index
end

"""
```stronglyconnectedcomponents{T}(dg::SimpleDigraph{T})```

Return all strongly connected components of a graph using the recursive 
[Tarjan's algorithm](https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm)

# Arguments
* dg: the digraph we want the strongly connected components.

# Returns
* visitor.components: all the strongly connected components.
"""
function stronglyconnectedcomponents{T}(dg::SimpleDigraph{T})
    visitor = TarjanVisitor(dg)
    for  v in dg.V 
        if visitor.indexmap[v] == 0
            strongconnect!(dg, v, visitor)
        end
    end
    return visitor.components
end


"""
```isacyclic{T}(dg::SimpleDigraph{T})```

Return true if the directed graph is acyclic.
"""
function isacyclic{T}(dg::SimpleDigraph{T})
    acyclic = (length(stronglyconnectedcomponents(dg)) == NV(dg))
    return acyclic
end

###############################################################################
####  Shortest paths

"""
```shortestpath{T}(dg::SimpleDigraph{T})```

Compute the shortest path between two vertices of a given directed graph, 
using a Floyd-Warshall algorithm. 

If the distance is higher than the number of vertices in the graph, 
they are not connected.
"""
function shortestpath{T}(dg::SimpleDigraph{T})
    n = NV(dg) + 1
    unitdist = Dict{T, Int}([v=> n for v in dg.V])
    dist = Dict{T, Dict{T, Int}}([v=>copy(unitdist) for v in dg.V])
    for v in dg.V
        dist[v][v] = 0
    end
    for (u, v) in elist(dg)
        dist[u][v] = 1
    end
    for k in dg.V
        for i in dg.V
            for j in dg.V
                if dist[i][j] > (dist[i][k] + dist[k][j])
                    dist[i][j] = dist[i][k] + dist[k][j]
                end
            end
        end
    end
    return dist
end


"""
```transitiveclosure!{T}(dg::SimpleDigraph{T})```

Compute the transitive closure of a directed graph, using a Floyd-Warshall
algorithm.
"""
function transitiveclosure!{T}(dg::SimpleDigraph{T})
    for k in dg.V
        for i in dg.V
            for j in dg.V
                if (has(dg, i, k) & has(dg, k, j))
                    add!(dg, i, j)
                end
            end
        end
    end
    return dg
end

#= Not working.
"""
```
type TransitiveVisitor{T} <: Visitor
    visited::Dict{T, Bool}
    stack::Vector{T}
end
```

DFS visitor used for the transitive closure.

# Arguments
* visited: whether a vertex has been visited or not
* stack: ordered list of vertex
"""
type TransitiveVisitor{T} <: Visitor
    visited::Dict{T, Bool}
    stack::Vector{T}
end

"""
```
TransitiveVisitor{T}(dg::SimpleDigraph{T}) = TransitiveVisitor{T}(
Dict([v =>false for v in dg.V]),
Vector{T}()
)
```

Constructor of the eponym type, using the knowledge of the number of 
vertices in the directed graph.
"""
TransitiveVisitor{T}(dg::SimpleDigraph{T}) = TransitiveVisitor{T}(
Dict([v =>false for v in dg.V]),
Vector{T}()
)

"""
```closure!{T}(v::T, dg::SimpleDigraph{T}, vis::TransitiveVisitor{T})```

Add the edges to the graph when necessary. 

If the directed graph is not looped, the `add!` will not add the self-loops.
"""
function closure!{T}(v::T, dg::SimpleDigraph{T}, vis::TransitiveVisitor{T})
    for w in vis.stack
        for x in out_neighbors(dg, v)
            add!(dg, w, x)
        end
    end
end        


"""
```recursiveclosure!{T}(v::T, dg::SimpleDigraph{T}, vis::TransitiveVisitor{T})```

One iteration of the transitive closure.
"""
function recursiveclosure!{T}(v::T, dg::SimpleDigraph{T}, vis::TransitiveVisitor{T})
    vis.visited[v] = true
    push!(vis.stack, v)
    #println("iter $v")
    for w in out_neighbors(dg, v)       
        closure!(w, dg, vis)
        if !vis.visited[w]
            recursiveclosure!(w, dg, vis)
            pop!(vis.stack)
        end
    end
end

"""
```transitiveclosure!{T}(dg::SimpleDigraph{T})```

Construct the transitive closure of a directed graph on the same graph, using a DFS visit.
"""
function transitiveclosure!{T}(dg::SimpleDigraph{T})
    visitor = TransitiveVisitor(dg)
    for v in vlist(dg)
        visitor.visited[v] && continue
        recursiveclosure!(v, dg, visitor)
        visitor.stack = Vector{T}()
    end
    return dg
end            
=#
###############################################################################
#### Cycles



"""
```ncycles_n_i(n::Integer, i::Integer)```

Compute the theoretical maximum number of cycles of size `i` in a directed graph of `n`
 vertices.
"""
function ncycles_n_i(n::Integer, i::Integer) 
    return binomial(big(n), big(n-i+1)) * factorial(big(n-i))
end

"""
```maxcycles(n::Integer)```

Compute the theoretical maximum number of cycles in a directed graph of n vertices, 
assuming there are no self-loops. 
Formula coming from [Johnson's paper](http://epubs.siam.org/doi/abs/10.1137/0204007).

# Arguments
* n: the number of vertices.
"""
function maxcycles(n::Integer)
    return sum(x -> ncycles_n_i(n, x), 1:(n-1))
end

"""
``` maxcycles{T}(dg::SimpleDigraph{T}, byscc::Bool = true)```

Compute the theoretical maximum number of cycles in the directed graph, 
taking into account the self loops if necessary. 

The computation can be performed assuming the graph is complete or taking into account the 
decomposition in strongly connected components. The formula coming from 
[Johnson's paper](http://epubs.siam.org/doi/abs/10.1137/0204007).

# Arguments:
* dg: The directed graph to be considered.
* byscc: whether it should be computed knowing the strongly connected components of the
 directed or not (default:yes)

# Return
* c: the theoretical maximum number of cycles.

Note: A more efficient version is possible.
"""    
function maxcycles{T}(dg::SimpleDigraph{T}, byscc::Bool = true) 
    c = 0
    n = length(dg.V)
    if !byscc        
        c = maxcycles(n)
    else
        for scc in stronglyconnectedcomponents(dg)
            length(scc) > 1 && (c += maxcycles(length(scc)))
        end     
    end
    dg.looped && (c += n)
    return c
end

"""
```type JohnsonVisitor{T} <: Visitor
    stack::Vector{T}
    blocked::Dict{T, Bool}
    blockedmap::Dict{T, Set{T}}
end```

Composite type that regroups the information needed for Johnson's algorithm.

# Arguments
* stack: the stack of visited vertices
* blocked: boolean for each vertex, tell whether it is blocked or not
* blockedmap: tell which vertices to unblock if the key vertex is unblocked.
"""
type JohnsonVisitor{T} <: Visitor
    stack::Vector{T}
    blocked::Dict{T, Bool}
    blockedmap::Dict{T, Set{T}}
end

"""
Constructor of the visitor, using the directed graph information.
"""
JohnsonVisitor{T}(dg::SimpleDigraph{T}) = JohnsonVisitor(
Vector{T}(),
Dict([v => false for v in dg.V]),
Dict([v => Set{T}() for v in dg.V])
)

"""
```unblock!{T}(v::T, blocked::Dict{T, Bool}, B::Dict{T, Set{T}})```

Unblock the vertices recursively.

# Arguments
* v: the vertex to unblock
* blocked: tell wether a vertex is blocked or not
* B: the map that tells if the unblocking of one vertex should unblock other vertices
"""
function unblock!{T}(v::T, blocked::Dict{T, Bool}, B::Dict{T, Set{T}})
    blocked[v] = false
    for w in B[v]
        delete!(B[v], w)
        blocked[w] && unblock!(w, blocked, B)
    end
end

"""
```circuit{T}(v::T, dg::SimpleDigraph{T}, vis::JohnsonVisitor{T}, 
allcycles::Vector{Vector{T}}, startnode = v)```

One step of the recursive version of simple cycle detection, using a DFS algorithm.

The CIRCUIT function from [Johnson's algorithm](http://epubs.siam.org/doi/abs/10.1137/0204007), 
recursive version. Modify the vector of cycles, when needed.

# Arguments
* v: the vertex considered in this iteration of the DFS
* dg: the digraph from which cycles are computed
* visitor: Informations needed for the cycle computation, contains:
    * stack: the stack of parent vertices
    * blocked: tells whether a vertex has already been explored or not
    * blockedmap: mapping of the blocking/ unblocking consequences
* allcycles: output containing the cycles already detected
* startnode = v: optional argument giving the starting node. In the first iteration,
 the same as v, otherwise it should be passed.

# Returns
* done: tells whether a circuit has been found in the current exploration.
"""
function circuit{T}(v::T, dg::SimpleDigraph{T}, vis::JohnsonVisitor{T}, 
allcycles::Vector{Vector{T}}, startnode = v)
    done = false
    push!(vis.stack, v)
    vis.blocked[v] = true
    for w in dg.N[v]
        if w == startnode
            push!(allcycles, deepcopy(vis.stack))           
            done = true
            elseif !vis.blocked[w]
            circuit(w, dg, vis, allcycles, startnode) && (done = true)
        end
    end
    if done
        unblock!(v, vis.blocked, vis.blockedmap)
    else
        for w in dg.N[v]
            in(v, vis.blockedmap[w]) || push!(vis.blockedmap[w], v)
        end
    end
    pop!(vis.stack)
    return done
end


"""
```simplecycles{T}(dg::SimpleDigraph{T})```

Compute all cycles of the given directed graph, using 
[Johnson's algorithm](http://epubs.siam.org/doi/abs/10.1137/0204007). 

/!\ The number of cycles grow more than exponentially with the number of vertices, 
you might want to use the algorithm with a ceiling  on large directed graphs 
(slightly slower). If you want to have an idea of the possible number of cycles, 
look at function ```maxcycles{T}(dg::SimpleDigraph{T}, byscc::Bool = true)```.

# Arguments:
* dg: the directed graph we want to explore

# Returns:
* cycles: all the cycles of the directed graph.
"""
function simplecycles{T}(dg::SimpleDigraph{T})
    sccs = stronglyconnectedcomponents(dg)
    cycles = Vector{Vector{T}}()
    is_looped(dg) ? (min = 0) : (min = 1)
    for scc in sccs
        while length(scc) > min
            wdg = subgraph(dg, scc)
            startnode = shift!(scc)
            visitor = JohnsonVisitor(wdg)
            circuit(startnode, wdg, visitor, cycles)
        end
    end
    return cycles
end


##########################################################
#### Iterative version, using Tasks, of the previous algorithms.
"""
```circuit{T}(v::T, dg::SimpleDigraph{T}, vis::JohnsonVisitor{T}, startnode = v)```

One step of the recursive version of simple cycle detection, using a DFS algorithm.

The CIRCUIT function from [Johnson's algorithm](http://epubs.siam.org/doi/abs/10.1137/0204007),
 recursive and iterative version. Produce a cycle when needed, can be used only inside a 
 Task.

# Arguments
* v: the vertex considered in this iteration of the DFS
* dg: the digraph from which cycles are computed
* visitor: Informations needed for the cycle computation, contains:
    * stack: the stack of parent vertices
    * blocked: tells whether a vertex has already been explored or not
    * blockedmap: mapping of the blocking/ unblocking consequences
* startnode = v: optional argument giving the starting node. In the first iteration, 
the same as v, otherwise it should be passed.

# Returns
* done: tells whether a circuit has been found in the current exploration.
"""
function circuit{T}(v::T, dg::SimpleDigraph{T}, vis::JohnsonVisitor{T}, startnode = v)
    done = false
    push!(vis.stack, v)
    vis.blocked[v] = true
    for w in dg.N[v]
        if w == startnode
            produce(copy(vis.stack))        
            done = true
        elseif !vis.blocked[w]
            circuit(w, dg, vis, startnode) && (done = true)
        end
    end
    if done
        unblock!(v, vis.blocked, vis.blockedmap)
            else
        for w in dg.N[v]
            in(v, vis.blockedmap[w]) || push!(vis.blockedmap[w], v)
        end
    end
    pop!(vis.stack)
    return done
end


"""
```itercycles{T}(dg::SimpleDigraph{T})```

Compute all cycles of the given directed graph, using 
[Johnson's algorithm](http://epubs.siam.org/doi/abs/10.1137/0204007). 
It should not be visible to the end user.

/!\ The number of cycles grow more than exponentially with the number of vertices, 
you might want to use the algorithm with a ceiling  on large directed graphs 
(slightly slower). If you want to have an idea of the possible number of cycles, 
look at function ```maxcycles{T}(dg::SimpleDigraph{T}, byscc::Bool = true)```.

# Arguments:
* dg: the directed graph we want to explore
"""
function itercycles{T}(dg::SimpleDigraph{T})
    sccs = stronglyconnectedcomponents(dg)
    is_looped(dg) ? (min = 0) : (min = 1)
    for scc in sccs
        while length(scc) > min
            wdg = subgraph(dg, scc)
            startnode = shift!(scc)
            visitor = JohnsonVisitor(wdg)
            circuit(startnode, wdg, visitor)
        end
    end
end

"""
```countcycles{T}(dg::SimpleDigraph{T}, ceiling = 10^6)```

Count the number of cycles in a directed graph, using 
[Johnson's algorithm](http://epubs.siam.org/doi/abs/10.1137/0204007).

The ceiling is here to avoir memory overload if there are a lot of cycles in the graph. 
Default value is 10^6, but it can be higher or lower. You can use the function 
```maxcycles{T}(dg::SimpleDigraph{T}, byscc::Bool = true)``` to get an idea of the 
theoretical maximum number or cycles.

# Arguments
* dg: the directed graph we are interested in.
* ceiling = 10^6: the maximum number of cycles we want to search for.

# Returns
* i: the number of cycles if below the ceiling, the ceiling otherwise.
"""
function countcycles{T}(dg::SimpleDigraph{T}, ceiling = 10^6)
    t = Task(() -> itercycles(dg))
    i = -1 
    while (t.state == :runnable) & (i < ceiling)
        consume(t)
        i += 1
    end
    return i
end

"""
```getcycles{T}(dg::SimpleDigraph{T}, ceiling = 10^6)```

Search all cycles of the given directed graph, using 
[Johnson's algorithm](http://epubs.siam.org/doi/abs/10.1137/0204007), 
up to the ceiling (avoid memory overload).

If the graph is small, the ceiling will not bite and 
``itercycles{T}(dg::SimpleDigraph{T})`` is more efficient. it avoids the overhead of the 
counting and testing if the ceiling is reached.

To get an idea of the possible number of cycles, using function 
```maxcycles{T}(dg::SimpleDigraph{T}, byscc::Bool = true)``` on the directed graph. 

# Arguments:
* dg: the directed graph to explore
* ceiling = 10^6: Maximum number of cycles that will be computed. It can be raised or 
lowered.

# Returns:
* cycles: all the cycles of the directed graph.
"""
function getcycles{T}(dg::SimpleDigraph{T}, ceiling = 10^6)
    t = Task(() -> itercycles(dg))
    i = -1
    cycles = Vector{Vector{T}}()
    while (t.state == :runnable) & (i < ceiling - 1)
        cycle = consume(t)
        i += 1
        cycle == nothing || push!(cycles, cycle)
    end
    return cycles
end

"""
```getcycleslength{T}(dg::SimpleDigraph{T})```

Search all cycles of the given directed graph, using 
[Johnson's algorithm](http://epubs.siam.org/doi/abs/10.1137/0204007), 
and return their length.

To get an idea of the possible number of cycles, using function 
```maxcycles{T}(dg::SimpleDigraph{T}, byscc::Bool = true)``` on the directed graph. 

# Arguments:
* dg: the directed graph to explore

# Returns:
* cyclelengths: the lengths of all cycles, the index in the array is the length. 
"""
function getcycleslength{T}(dg::SimpleDigraph{T})
    t = Task(() -> itercycles(dg))
    cyclelength = zeros(Int, NV(dg))
    while (t.state == :runnable) 
        cycle = consume(t)
        cycle == nothing || (cyclelength[length(cycle)] +=1)
    end
    return cyclelength
end
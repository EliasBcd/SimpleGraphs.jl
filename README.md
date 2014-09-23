#### Note: This documentation is a work in progress

To get a sense of all the functionality provided, you may need to look
in the source files directly. Of course, our goal is to have this
documentation be complete. Your help is welcome.

# The `SimpleGraphs` Julia Module

This module defines two data types for working with graphs:

+ The `SimpleGraph` type represents *undirected* graphs without loops
  or multiple edges.
+ The `SimpleDigraph` type represents *directed* graphs in which there
  may be at most one directed edge `(u,v)` from a vertex `u` to a
  vertex `v`. There may also be a directed edge in the opposite
  direction, `(v,u)`.


To get started, get the package like this:
```julia
Pkg.clone("https://github.com/scheinerman/SimpleGraphs.jl.git")
```
Then, in any Julia session or program, give the command
`using SimpleGraphs`.

## Graph types

The `SimpleGraphs` module defines two data types:

+ `G = SimpleGraph()` creates a new *simple* graph.
+ `D = SimpleDigraph()` creates a new *directed* graph.

In both cases, the vertices of the graph may be of `Any` type. More
often, it is useful to restrict the type of the vertices to be of a
given type. This can be done with `SimpleGraph{T}()` or
`SimpleDigraph{T}()` where `T` is a Julia type.

Two common choices the vertex type are integers and strings. For this
purpose, we provide these special constructors:

+ `IntGraph()` and `IntDigraph()` creates graphs with `Int`
  vertices. These can be called with an optional positive integer
  argument like this: `IntGraph(n)` or `IntDigraph(n)`. This
  prepopulates the vertex set with `n` vertices named `1` through `n`.
+ `StringGraph()` and `StringDigraph()` creates graphs whose vertex
  type is `ASCIIString`.

## Fundamental operations

### Adding/deleting edges/vertices

The most basic operations for graphs are adding and deleting vertices and edges. These are done with `add!` and `delete!`. In general, if `G` is a graph or a digraph, then we have the following:

+ `add!(G,v)` adds the vertex `v` to the graph.
+ `add!(G,v,w)` adds the edge `(v,w)` to the graph. If `G` is an
  undirected graph, then `(v,w)` and `(w,v)` are the same thing. See
  also the discussion concerning loops just below.
+ `delete!(G,v)` deletes the vertex `v` from the graph as well as any
  edges that might be incident with `v`.
+ `delete!(G,v,w)` deletes the edge `(v,w)`.

The `has` function may be used to determine if a given vertex or edge is present in the graph:

+ `has(G,v)` returns `true` if `v` is a vertex of `G`.
+ `has(G,v,w)` returns `true` if `(v,w)` is an edge of `G`. In the
  case of undirected graphs, this is the same as `has(G,w,v)`.

### Vertex and edge sets

The functions `NV(G)` and `NE(G)` return the number of vertices and
edges in `G`, respectively.

The functions `vlist(G)` and `elist(G)` return a list of the vertices
and edges in `G`. These lists are sorted if possible.

### Neighborhoods and degrees

If `G` is a `SimpleGraph` and `v` is one of its vertices, then
`neighbors(G,v)` as a (possibly sorted) list. This may also be found
with `G[v]`. The degree of `v` is returned by `deg(G,v)`. A call to
`deg(G)` returns a sorted list of the degrees of all the vertices.

The situation is more involved if `G` is a `SimpleDigraph`. In this
case we have the following.
 
+ `out_neighbors(G,v)` returns a list of all vertices `w` such that
  `(v,w)` is an edge of `G`. This includes `v` itself if there's a
  loop at `v`.
+ `in_neighbors(G,v)` returns a list of all vertices `u` such that
  `(u,v)` is an edge of `G`.
+ `out_deg(G,v)` is the number of vertices in `out_neighbors(G,v)`.
+ `in_deg(G,v)` is the number of vertices in `in_neighbors(G,v)`.
+ `deg(G,v)` is the sum of `out_deg(G,v)` and `in_deg(G,v)`.

#### Technical note on neighborhoods in undirected graphs

The `SimpleGraph` data structure provides a mechanism for speedy
neighborhood lookup. This comes at a cost of holding considerable
redundant information. If the graph is very large, this may cause a
problem.

The function `fastN!` may be used to supress (or restore) this
additional data structure.

+ `fastN!(G,false)` destroys the redundant data
  structure. Neighborhood lookup is now much slower (as are any
  algorithms that might depend on it).
+ `fastN!(G,true)` restores the additional data structure. Note that
  `fastN!(G)` has the same effect.

To check if a graph possesses this data structure, the value of
`G.Nflag` may be inspected *but in no circumstance should this value
be changed except through the use of* `fastN!`.


### Loops

A `SimpleGraph` may never contain a loop (an edge whose end points are
the same). However, the ability for a `SimpleDigraph` to have loops is
at the user's discretion. The default constructor `SimpleDigraph{T}()`
creates a digraph that may have loops; use `SimpleDigraph{T}(false)`
to create a digraph that is incapable of containing loops.

In addition, we provide the following functions to inspect and
manipulate loops:

+ `is_looped(D)` tests if `D` may hold loops; it returns `true` if so
  and `false` otherwise.
+ `loops(D)` returns a list of vertices in `D` at which there is a
  loop present.
+ `allow_loops(D)` enables the digraph to have loops.
+ `remove_loops(D)` deletes all loops from `D` (if any) but does not
  alter its ability to have loops.
+ `forbid_loops(D)` deletes all loops (if any) and prevents `D` from
  having loops.
 
## Constructors

We provide a variety of functions for generating certain standard graphs. 

### Constructors for undirected graphs

In addition to the basic `SimpleGraph`, `IntGraph`, and `StringGraph`
constructors, we have the following:

+ `Complete(n)` creates a complete graph with `n` vertices numbered
  `1` through `n`.
+ `Complete(n,m) ` creates a complete bipartite graph with `n`
  vertices in one part and `m` vertices in the other.
+ `Complete([n1,n2,n3,...,nt])` creates a complete multipartite graph
  with part sizes as given in the array of integers
  `[n1,n2,n3,...,nt]`.
+ `Cube(n)` creates an `n`-dimensional cube graph. The vertices are
   `n`-long character strings of 0s and 1s. Two vertices are adjacent
   iff they differ in exactly one location.
+ `Path(n)` creates a path graph with `n` vertices. Alternatively, if
   `list` is a 1-dimensional array then this creates a graph with the
   members of `list` as vertices and with edges of the form
   `(list[k],list[k+1])`. The elements of `list` should be distinct
   *but this is not checked*.
+ `Cycle(n)` creates a cycle on `n` vertices. This requires `n` to be
   at least 3.
+ `Wheel(n)` creates a wheel graph with `n` vertices. This is formed
  from an `n-1`-cycle and additional vertex adjacent to all the
  vertices on the cycle.
+ `Grid(n,m)` creates an `n`-by-`m` grid graph.
+ `BuckyBall()` creates the 60-vertex graph that is the 1-skellaton of
  a truncated icosahedron.
+ `RandomGraph(n,p)` creates an Erdos-Renyi random graph with `n`
  vertices and edge probability `p`. The parameter `p` may be omitted,
  in which case the value 1/2 is used.
+ `RandomTree(n)` creates a random tree with `n` vertices. There are
  `n^(n-2)` trees with vertices set `{1,2,...,n}` and they are all
  equally likely.
+ `Knesser(n,k)` creates the Knesser graph. The vertices of this graph
  are the `k`-element subsets of `{1,2,...,n}`. Two vertices are
  adjacent iff they correspond to disjoint sets.
+ `Petersen()` creates the Petersen graph as `Knesser(5,2)`.

### Constructors for directed graphs

In addition to `SimpleDigraph`, `IntDigraph`, and `StringDigraph`
constructors, we have the following:

+ `DirectedComplete(n)` creates an `n`-vertex digraph with all `n*n`
  possible edges (including loops). To prevent loops, use
  `DirectedComplete(n,false)`.
+ `DirectedPath(n)` creates an `n`-vertex directed path with the `n-1`
  edges `(1,2)`, `(2,3)`, ..., `(n-1,n)`.
+ `DirectedCycle(n)` creates an `n`-vertex directed cycle with the `n`
  edges `(1,2)`, `(2,3)`, ..., `(n-1,n)`, `(n,1)`.
+ `RandomDigraph(n,p)` creates a random directed graph in which each
  of the `n*n` possible edges is present with probability `p` (whose
  default value is 1/2). Use `RandomDigraphs(n,p,false)` to prevent
  the formation of loops.
+ `RandomTournament(n)` creates a random digraph on `n` vertices. For
  each pair of distinct vertices `u` and `v` we have exactly one of
  the edges `(u,v)` or `(v,u)` with probability 1/2 each
  (independently for all pairs of vertices). This graph has no loops.

### Converting a directed graph into an undirected graph

If `D` is a `SimpleDigraph` then `simplify(D)` creates a new
`SimpleGraph` that's formed by removing directions (and loops). That
is, the new graph has the same vertices as `D` and an edge between
distinct vertices `u` and `v` if and only if `(u,v)` or `(v,u)` (or
both) is an edge of `D`.


## Graph operations

Undirected graphs only at this time. To be documented:

+ `isequal`, `G==H`
+ `copy`
+ `complement`, `complement!`, `G'`
+ `induce`
+ `contract!`
+ `cartesian`, `G*H`
+ `disjoint_union`
+ `union`
+ `join`
+ `trim`
+ `relabel`

## Paths and connectivity 

Undirected graphs only at this time. To be documented:

+ `is_connected`
+ `num_components`
+ `components`
+ `find_path`
+ `dist`
+ `dist_matrix`
+ `diam`
+ `is_cut_edge`
+ `spanning_forest`

## Graph matrices

We generate the following kinds of matrices for graphs

+ `adjacency` for both undirected and directed graphs.
+ `laplace` for undirected graphs only.
+ `incidence` for both undirected and directed graphs.
+ `dist_matrix` for undirected graphs only.


## Algorithms

Undirected graphs only at this time. To be documented:

+ `euler`
+ `bipartition`
+ `two_color`
+ `greedy_color`
+ `random_greedy_color`

## Interface to `Graphs.jl`

We provide a `convert_simple` function that takes a `SimpleGraph` as
input and returns a Julia `Graphs.simple_graph` representation of the
same graph (together with dictionaries to match up the vertex sets).

# Please Help

This is very much a work in process with a lot of more features that
can/should be added. If you're interested in contributing, please
contact me. I'm especially interested in JHU undergraduates getting
involved.

Ed Scheinerman (ers@jhu.edu)

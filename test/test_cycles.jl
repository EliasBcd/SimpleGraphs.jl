using SimpleGraphs

import Base.Test

dg1 = loadgraph("data\\largegraph.csv", separator = ',', directed = true)
dg2 = loadgraph("data\\largegraph.csv", separator = ',', directed = true, vertextype = Int)
dg3 = DirectedComplete(4, true)
dg4 = DirectedComplete(4, false)

@test typeof(getcycles(dg1)) == Vector{Vector{ASCIIString}}
@test typeof(getcycles(dg2)) == Vector{Vector{Int}}
@test countcycles(dg1) == 12532
@test countcycles(dg2) == 12532
@test countcycles(dg2, 10) == 10
@test countcycles(dg1, 10) == 10
@test maxcycles(dg4, true) == 20
@test maxcycles(dg4, false) == 20
@test maxcycles(dg3, true) == 24
@test maxcycles(dg3, false) == 24

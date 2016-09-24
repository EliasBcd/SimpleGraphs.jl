using SimpleGraphs

import Base.Test

function readedges(T::DataType)
    RPrelations = readdlm("data\\largegraph.csv", ',', Int)
    RP = Vector{Tuple{T, T}}(0)
    if T <: AbstractString
        for i in 1:size(RPrelations, 1)
            push!(RP, (string(RPrelations[i, 1]), string(RPrelations[i, 2])))
        end
        elseif T <: Integer
        for i in 1:size(RPrelations, 1)
            push!(RP, (RPrelations[i, 1], RPrelations[i, 2]))
        end
    end
    dg = SimpleDigraph{T}()
    addedges!(dg, RP)
    return dg
end

dg1 = loadgraph("data\\largegraph.csv", separator = ',', vertextype = Int, directed = true)
dg2 = deepcopy(dg1)
transitiveclosure!(dg2)
@show dg1, dg2


dg3 = loadgraph("data\\largegraph.csv", separator = ',', directed = true)
dg4 = deepcopy(dg3)
transitiveclosure!(dg4)
@show dg3, dg4

@test length(TarjanSCC(dg1)) == length(TarjanSCC(dg2))
@test length(TarjanSCC(dg3)) == length(TarjanSCC(dg4))

@test NE(dg1) == NE(dg3)
@test NE(dg2) == NE(dg4)
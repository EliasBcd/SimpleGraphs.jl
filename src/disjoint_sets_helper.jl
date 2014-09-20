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

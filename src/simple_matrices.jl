# Functions to create standard graph matrices

export adjacency, laplace, incidence, dist_matrix

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


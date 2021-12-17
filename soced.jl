using Arpack, SparseArrays, LinearAlgebra

function Hsoc(ksoc::Float64)

    # define functions used here
    include("ades.jl")
    include("acre.jl")

    # define a matrix for the Hamiltonian

    # diagonal terms
    for jj = 1:Msize

end










# defining vectors
n = 10
a = Matrix(I, n, n)

b = [1 0 0; 0 1 0; 0 0 1]
c = zeros(10,1)
d = ones(10,1)

e = sparse(a)
eigs(e)

# defining functions
function f(x)
    return x + 1
end

# g(1) will fail, but g(1.0) will work
function g(x::Float64)
    return x+1
end

# h(1) will work, but h("1") will not
function h(x::T) where {T <: Number}
    return x+1
end

function main(x)
    f(x)
    #g(x)
    h(x)

    return nothing
end

main(1)

@time main(1)

using BenchmarkTools

@benchmark main(1)

using Plots

a = rand(10,10)

heatmap(a)

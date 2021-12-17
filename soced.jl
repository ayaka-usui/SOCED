using Arpack, SparseArrays, LinearAlgebra

function Hsoc(Msize0::Int64,Np::Int64,ksoc::Float64,Omega::Float64)

    # define functions used here
    include("ades.jl")
    include("acre.jl")
    include("pascaltriangle.jl")
    include("in2b.jl")
    include("b2in.jl")
    include("nop.jl")

    Msize = Msize0*2
    matp = pascaltriangle(Msize,Np) # the size is Msize+1 times Np+1
    maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    # defines vectors and matrices
    vecmb = sparse(zeros(Int64,Msize+1));
    Hsoc = sparse(zeros(Float64,maxmatp,maxmatp));

    # define a matrix for the Hamiltonian

    # diagonal terms
    for nn = 1:maxmatp

        vecmbnn = in2b(nn,Msize,Np)

        for jj = 1:Msize

            vecmbnnj = ades(jj,vecmbnn)
            if vecmbnnj[Msize+1] == 0
               continue
            end

            for ii = 1:Msize

                vecmbnnij = acre(ii,vecmbnnj)
                energy0 = epsilon(ii,jj,Msize0,ksoc,Omega)
                if energy0 == 0
                   continue
                end

            end

        end

        for mm = 1:maxmatp

            vecmbmm = in2b(mm,Msize,Np)
            Hsoc[mm,nn] = (vecmbmm[1:Msize]' * vecmbnnij[1:Msize])*sqrt(vecmbnnij[Msize+1])*epsilon

        end

    end

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

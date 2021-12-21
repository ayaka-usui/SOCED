include("pascaltriangle.jl")
include("in2b.jl")
include("b2in.jl")

function cutMsizeEne(Msize0::Int64, Np::Int64, Ene0::Float64)

    Msize = Msize0*2

    matp = zeros(Int,Msize+1,Np+1);
    pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    for nn = 1:maxmatp

        vecmbnn = in2b(nn,Msize,Np)

        Enenn = 

    end





end

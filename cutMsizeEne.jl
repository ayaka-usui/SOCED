# include("pascaltriangle.jl")
# include("in2b.jl")
# include("b2in.jl")
include("in2bEne.jl")

function cutMsizeEne(Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Ene0minumhalf::Int64)

    # this function puts the cut-off in the number of states by a given energy
    # Ene0minumhalf equals to "the energy cut-off" - 1/2 so that Ene0minumhalf is integer. 

    Msize = Msize0*2

    # matp = zeros(Int64,Msize+1,Np+1)
    # pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    indvec = zeros(Int64,maxmatp)
    mm = 0

    for nn = 1:maxmatp

        Enenn = in2bEne(nn,Msize,Np,matp)

        if Enenn <= Ene0minumhalf
           mm = mm + 1
           indvec[mm] = nn
        end

    end

    indvec1 = copy(indvec[1:mm])

    return indvec1

end

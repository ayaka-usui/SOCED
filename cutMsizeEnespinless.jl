include("in2bEnespinless.jl")

function cutMsizeEnespinless(Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Enecutoff::Float64)

    # this function puts the cut-off in the number of states by a given energy

    maxmatp = matp[Msize0+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    indvec = zeros(Int64,maxmatp)
    mm = 0

    for nn = 1:maxmatp

        Enenn = in2bEnespinless(nn,Msize0,Np,matp)

        if Enenn < Enecutoff || isapprox(Enenn,Enecutoff) # Enenn <= Enecutoff
           mm = mm + 1
           indvec[mm] = nn
        end

    end

    indvec1 = copy(indvec[1:mm])

    return indvec1

end

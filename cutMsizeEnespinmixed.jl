include("in2bEnespinless.jl")

function cutMsizeEnespinmixed(Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Enecutoff::Float64)

    # this function puts the cut-off in the number of states by a given energy

    # maxmatp = matp[Msize0+1,Np+1] # the indices are m+1 and n+1 for N^m_ns
    maxmatp = matp[Msize0+1,Np]

    indvec = zeros(Int64,maxmatp*maxmatp)
    mm = 0

    for nn = 1:maxmatp*maxmatp

        indket2 = mod(nn,Msize0)
        if indket2 != 0
           indket1 = div(nn,Msize0)+1
        else
           indket1 = div(nn,Msize0)
           indket2 = Msize0
        end

        Enenn = in2bEnespinless(indket1,Msize0,Np-1,matp)
        Enenn += in2bEnespinless(indket2,Msize0,Np-1,matp)

        if Enenn < Enecutoff || isapprox(Enenn,Enecutoff) # Enenn <= Enecutoff
           mm = mm + 1
           indvec[mm] = nn
        end

    end

    indvec1 = copy(indvec[1:mm])

    return indvec1

end

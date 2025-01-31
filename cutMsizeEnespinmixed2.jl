include("in2bEnespinless.jl")

function cutMsizeEnespinmixed3(Msize0::Int64, Np::Int64, matp20::Matrix{Int64}, matp21::Matrix{Int64}, Enecutoff::Float64, indNp::Int64)

    # this function puts the cut-off in the number of states by a given energy

    # maxmatp = matp[Msize0+1,Np+1] # the indices are m+1 and n+1 for N^m_ns
    maxmatp20 = matp20[Msize0+1,Np-indNp+1]
    maxmatp21 = matp21[Msize0+1,indNp+1]

    indvec = zeros(Int64,maxmatp20*maxmatp21)
    mm = 0

    for nndown = 1:maxmatp20

        Enenn0 = in2bEnespinless(nndown,Msize0,Np-indNp,matp20)

        if Enenn0 < Enecutoff || isapprox(Enenn0,Enecutoff) # Enenn0 <= Enecutoff

            for nnup = 1:maxmatp21

                mm += 1
                indvec[mm] = (nndown-1)*maxmatp21 + nnup

                # Enenn = Enenn0 + in2bEnespinless(nnup,Msize0,indNp,matp21)
                #
                # if Enenn < Enecutoff || isapprox(Enenn,Enecutoff) # Enenn <= Enecutoff
                #    mm += 1
                #    indvec[mm] = (nndown-1)*maxmatp21 + nnup
                # end

            end

        end
    end

    indvec1 = copy(indvec[1:mm])

    return indvec1

end

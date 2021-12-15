function b2in(vecmb::Vector{Int},Msize::Int,Np::Int)

    # Lm = Msize+1
    Ln = Np + 1
    indfk = 1
    include("pascaltriangle.jl") # define pascaltriangle(m,n)

    for jj = 1:Msize-2 # take -2 out?
        for ii = 0:vecmb[jj]

            if ii == vecmb[jj]
               continue
            end

            sumparticles = sum(vecmb[1:jj-1]) # numnber of particles taken account so far
            indfk = indfk + pascaltriangle(Msize-jj,Np-ii-sumparticles)

        end
    end

    return indfk

end

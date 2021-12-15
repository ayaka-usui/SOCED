function b2in(vecmb::Vector{Int})

    mbsize = legnth(vecmb)-1
    Np = sum(vecmb[1:mbsize])

    # Lm = mbsize+1
    # Ln = Np + 1
    indfk = 1
    include("pascaltriangle.jl") # define pascaltriangle(m,n)

    for jj = 1:mbsize
        for ii = 0:vecmb[jj]

            if ii == vecmb[jj]
               continue
            end

            sumparticles = sum(vecmb[1:jj-1]) # numnber of particles taken account so far
            pascal = pascaltriangle(mbsize-jj+1,Np-ii-sumparticles+1) # pascalmax(m+1,n+1) = (n+m-1)!/n!(m-1)!
            indfk = indfk +

        end
    end

    return indfk

end

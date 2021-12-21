include("pascaltriangle.jl") # define pascaltriangle(m,n)
# include("pascalmax.jl") # degine pascalmax(m,n)

function b2in(vecmb::SparseVector{Int64})

    # this function return the index of many body state at the index of Fock state

    # ignore the last element since it is the coefficient
    Msize = length(vecmb)-1
    Np = sum(vecmb[1:Msize])

    # Lm = mbsize+1
    # Ln = Np + 1
    indfk = 1
    matp = zeros(Int,Msize+1,Np+1);
    pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # matp = pascaltriangle(Msize,Np) # the size is Msize+1 times Np+1
    # note the indices are m+1 and n+1 for N^m_n

    for jj = 1:Msize
        for ii = 0:vecmb[jj]

            if ii == vecmb[jj]
               continue
            end

            sumparticles = sum(vecmb[1:jj-1]) # numnber of particles taken account so far
            indfk = indfk + matp[Msize-jj+1,Np-ii-sumparticles+1] # (n+m-1)!/n!(m-1)!
            # note the indices are m+1 and n+1 for N^m_n
            # indfk = indfk + pascalmax(Msize-jj,Np-ii-sumparticles) # pascalmax(m,n) = (n+m-1)!/n!(m-1)!

        end
    end

    # unlike Ponomarev's way, I count the index of Fock state in the backwards
    # so that "1" means |Np,0,0,...>.
    indfk = matp[Msize+1,Np+1] + 1 - indfk # the indices are m+1 and n+1 for N^m_n

    return indfk

end

# include("pascaltriangle.jl") # define pascaltriangle(m,n)

function in2bEne(jj::Int64, Msize::Int64, Np::Int64, matp::Matrix{Int64})

    # this function return the index of Fock state at the energy of many body state

    indfk = jj

    # matp = zeros(Int,Msize+1,Np+1);
    # pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # note the indices are m+1 and n+1 for N^m_n

    # unlike Ponomarev's way, I count the index of Fock state in the backwards
    # so that "1" means |Np,0,0,...>.
    indfk = matp[Msize+1,Np+1] + 1 - indfk # the indices are m+1 and n+1 for N^m_n

    if indfk > matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_n
       error("indfk larger than the maximum of N^M_N.")
    end

    indfk = indfk - 1
    # vecmb = spzeros(Int64,Msize+1) #spzeros(Int64,Msize+1)

    indM = Msize-1
    indN = Np

    Ene = 0

    while true

          if indN == 0
             break
          end

          if indfk >= matp[indM+1,indN+1]

             indfk = indfk - matp[indM+1,indN+1] # the indices are m+1 and n+1 for N^m_n
             # vecmb[Msize-indM] = vecmb[Msize-indM] + 1
             indN = indN - 1

             Ene = Ene + ceil(Int64,(Msize-indM)/2)-1

          else
             indM = indM - 1
          end

    end

    Ene = Ene + Int64(Np/2)

    return Ene

end

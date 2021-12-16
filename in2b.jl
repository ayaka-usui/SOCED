function in2b(indfk::Int,Msize::Int,Np::Int)

    # this function return the index of Fock state at the index of many body state

    # indfk = indfk - 1
    vecmb = zeros(Int,Msize+1)

    include("pascaltriangle.jl") # define pascaltriangle(m,n)
    matp = pascaltriangle(Msize,Np) # the size is Msize+1 times Np+1
    # note the indices are m+1 and n+1 for N^m_n

    indN = Np
    indM = Msize-1

    while true

          if indfk == 1
             break
          end

          if indfk >= matp[indM+1,indN+1]
             indfk = indfk - matp[indM+1,indN+1]
             vecmb[Msize-indM] = vecmb[Msize-indM] + 1
             indN = indN - 1
          else
             indM = indM - 1
          end

    end

    # the last element for the coefficient
    vecmb[Msize+1] = 1

    return vecmb

end

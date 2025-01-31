function in2bEnespinless(jj::Int64, Msize0::Int64, Np::Int64, matp::Matrix{Int64})

    # this function return the index of Fock state at the energy of many body state

    # unlike Ponomarev's way, I count the index of Fock state in the backwards so that "1" means |Np,0,0,...>.
    indfk = jj
    indfk = matp[Msize0+1,Np+1] + 1 - indfk

    if indfk > matp[Msize0+1,Np+1] # the indices are m+1 and n+1 for N^m_n
       error("indfk is larger than the maximum of N^M_N.")
    end

    indfk = indfk - 1
    indM = Msize0-1
    indN = Np
    Ene = 0

    while true

          if indN == 0
             break
          end

          if indfk >= matp[indM+1,indN+1]

             indfk = indfk - matp[indM+1,indN+1] # the indices are m+1 and n+1 for N^m_n
             indN = indN - 1
             Ene = Ene + Msize0-indM - 1

          else

             indM = indM - 1

          end

    end

    Ene = Ene + Np/2
    
    return Ene

end

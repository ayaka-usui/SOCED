# include("pascaltriangle.jl") # define pascaltriangle(m,n)

function in2bindtid(jj::Int64, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, vecmbind::Vector{Int64})

    # this function return the index of Fock state at the index of many body state
    # "in2b" means "an index to a many-body state"

    # unlike Ponomarev's way, I count the index of Fock state in the backwards so that "1" means |Np,0,0,...>.
    indfk = jj
    indfk = matp[Msize0+1,Np+1] + 1 - indfk # the indices are m+1 and n+1 for N^m_n

    if indfk > matp[Msize0+1,Np+1] # the indices are m+1 and n+1 for N^m_n
       error("indfk larger than the maximum of N^M_N.")
    end

    indfk = indfk - 1
    # vecmb = spzeros(Int64,Msize0+1)
    # vecmbind = spzeros(Int64,Np)
    vecmbind .= 0
    nn = 0

    indM = Msize0-1
    indN = Np

    while true

          if indN == 0
             break
          end

          if indfk >= matp[indM+1,indN+1]
             indfk = indfk - matp[indM+1,indN+1] # the indices are m+1 and n+1 for N^m_n
             # vecmb[Msize0-indM] = vecmb[Msize0-indM] + 1
             nn = nn + 1
             vecmbind[nn] = Msize0 - indM
             indN = indN - 1
          else
             indM = indM - 1
          end

    end

    # the last element for the coefficient
    # vecmb[Msize0+1] = 1

    return vecmbind

end

using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("vijkl.jl")
# include("in2b.jl")

function Hinttest0(Msize0::Int64)

    matind = zeros(Int64,binomial(Msize0+3,4),4)
    ind0 = 0

    for n1 = 0:Msize0-1
        for n2 = 0:n1
            for n3 = 0:n2
                for n4 = 0:n3
                    ind0 += 1
                    matind[ind0,1] = n1
                    matind[ind0,2] = n2
                    matind[ind0,3] = n3
                    matind[ind0,4] = n4
                end
            end
        end
    end

    return matind

end

function indtest0(n1::Int64,n2::Int64,n3::Int64,n4::Int64)
    return binomial(n1+3,4)+binomial(n2+2,3)+binomial(n3+1,2)+binomial(n4+0,1)+1
end

function Hintfunccutoff2!(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Hintdown::ST, Hintup::ST, Hintdu::ST) where ST <: Union{SparseMatrixCSC{Float64},Array{Float64}}

    maxmatpcut = length(indvec)

    # register data of Vijkl in a vector
    ind0 = 0
    vecV = zeros(Float64,binomial(Msize0+3,4))
    for n1 = 0:Msize0-1
        for n2 = 0:n1
            for n3 = 0:n2
                for n4 = 0:n3
                    ind0 += 1
                    vecV[ind0] = vijkl(n1,n2,n3,n4)
                end
            end
        end
    end

    # defines vectors and matrices
    Hintdown .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintup .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintdu .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)

    # define a matrix for the Hamiltonian for down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        for mm = 1:nn

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            ind1 = [vecmbindnn[1]-1, vecmbindnn[2]-1, vecmbindmm[1]-1, vecmbindmm[2]-1]
            sort!(ind1,rev=true)
            ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1
            # println(ind2)
            # println(length(vecV))

            if vecmbindnn[1] == vecmbindnn[2]
               if vecmbindmm[1] == vecmbindmm[2]
                  Hintdown[mm,nn] = vecV[ind2]*factorial(Np)
               else
                  Hintdown[mm,nn] = vecV[ind2]*sqrt(factorial(Np))*2
               end
            else
               if vecmbindmm[1] == vecmbindmm[2]
                  Hintdown[mm,nn] = vecV[ind2]*sqrt(factorial(Np))*2
               else
                  Hintdown[mm,nn] = vecV[ind2]*4
               end
            end

            # if vecmbindnn[1] == vecmbindnn[2] && vecmbindmm[1] == vecmbindmm[2]
            #    Hintdown[nn,mm] = vecV[ind2]*factorial(Np)
            # elseif vecmbindnn[1] == vecmbindnn[2] || vecmbindmm[1] == vecmbindmm[2]
            #    Hintdown[nn,mm] = vecV[ind2]*sqrt(factorial(Np))*2
            # else
            #    Hintdown[nn,mm] = vecV[ind2]*4
            # end

        end

    end

    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdown .= Hintdown + Hintdown' - spdiagm(diag(Hintdown))

    # adjust the size of Hintdown to the size of total Hamiltonian
    Hintdown2 = spzeros(maxmatpcut,maxmatpcut)
    Hintdown2[1:maxmatpcut,1:maxmatpcut] = Hintdown

    # define Hintup
    Hintup2 = spzeros(maxmatpcut,maxmatpcut)
    Hintup2[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hintdown






    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdu .= Hintdu + Hintdu' - spdiagm(diag(Hintdu))

end

using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("ades.jl")
include("acre.jl")
# include("pascaltriangle.jl")
include("in2b.jl")
# include("b2in.jl")
# include("epsilon.jl")
include("Vijkl2.jl")

function Hintfunccutoff2!(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Hintdown::SparseMatrixCSC{Float64}, Hintup::SparseMatrixCSC{Float64}, Hintdu::SparseMatrixCSC{Float64})

    # construct a matrix for the interaction Hamiltonian

    Msize = Msize0*2
    maxmatpcut = length(indvec)

    # matp = zeros(Int,Msize+1,Np+1);
    # pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    # defines vectors and matrices
    Hintdown .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintup .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintdu .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);

    # calculate Vijkl in advance
    matV = zeros(Float64,Msize0,Msize0,Msize0,Msize0) #spzeros(Float64,Msize0,Msize0,Msize0,Msize0)
    for nn4 = 0:Msize0-1
        for nn3 = 0:nn4
            for nn2 = 0:nn3
                for nn1 = 0:nn2

                    if isodd(nn1+nn2+nn3+nn4)
                       continue
                    end

                    matV[nn1+1,nn2+1,nn3+1,nn4+1] = Vijkl2(nn1,nn2,nn3,nn4)

                end
            end
        end
    end

    # define a matrix for the Hamiltonian
    # Threads.@threads for nn = 1:maxmatpcut #maxmatp # parfor
    for nn = 1:maxmatpcut

        # define ket state |n>
        vecmbnn = spzeros(Int64,Msize+1)
        in2b!(indvec[nn],Msize,Np,matp,vecmbnn) #in2b(nn,Msize,Np)

        # Interactions
        for ll = 1:Msize

            # apply an anhilation operator a_ll to |n>
            vecmbnnl = spzeros(Int64,Msize+1);
            ades!(ll,vecmbnn,vecmbnnl) #vecmbnnl = ades(ll,vecmbnn)
            if vecmbnnl[Msize+1] == 0 # go back if a|0> = 0
               continue
            end

            for kk = 1:Msize

                # apply a_kk to a_ll |n>
                vecmbnnkl = spzeros(Int64,Msize+1)
                ades!(kk,vecmbnnl,vecmbnnkl) #vecmbnnkl = ades(kk,vecmbnnl)
                if vecmbnnkl[Msize+1] == 0
                   continue
                end

                for jj = 1:Msize

                    # apply a_jj^{+} to a_kk a_ll |n>
                    vecmbnnjkl = spzeros(Int64,Msize+1)
                    acre!(jj,vecmbnnkl,vecmbnnjkl) #vecmbnnjkl = acre(jj,vecmbnnkl)

                    for ii = 1:Msize

                        # apply a_ii^{+} to a_jj^{+} a_kk a_ll |n>
                        vecmbnnijkl = spzeros(Int64,Msize+1)
                        acre!(ii,vecmbnnjkl,vecmbnnijkl) #vecmbnnijkl = acre(ii,vecmbnnjkl)

                        # get the index of eigenstate of harmonic oscillator
                        n1 = ceil(Int64,ii/2)-1
                        n2 = ceil(Int64,jj/2)-1
                        n3 = ceil(Int64,kk/2)-1
                        n4 = ceil(Int64,ll/2)-1

                        # V_{ijkl}=0 if n1+n2+n3+n4 is odd
                        if isodd(n1+n2+n3+n4)
                           continue
                        end

                        # @floop for mm = 1:nn
                        for mm = 1:nn

                            vecmbmm = spzeros(Int64,Msize+1)
                            in2b!(indvec[mm],Msize,Np,matp,vecmbmm) #vecmbmm = in2b(indvec[mm],Msize,Np)

                            if vecmbnnijkl[1:Msize] == vecmbmm[1:Msize] # if <m| a_ii^{+} a_jj^{+} a_kk a_ll |n> is not zero

                               vec0 = [n1,n2,n3,n4]
                               sort!(vec0,rev=true)
                               n1 = vec0[1]
                               n2 = vec0[2]
                               n3 = vec0[3]
                               n4 = vec0[4]

                               if isodd(ii) && isodd(jj) && isodd(kk) && isodd(ll)
                                  Hintdown[mm,nn] = Hintdown[mm,nn] + matV[n1+1,n2+1,n3+1,n4+1]*sqrt(vecmbnnijkl[Msize+1])
                                  # Hintdown[mm,nn] = Hintdown[mm,nn] + Vijkl2(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               elseif iseven(ii) && iseven(jj) && iseven(kk) && iseven(ll)
                                  Hintup[mm,nn] = Hintup[mm,nn] + matV[n1+1,n2+1,n3+1,n4+1]*sqrt(vecmbnnijkl[Msize+1])
                                  # Hintup[mm,nn] = Hintup[mm,nn] + Vijkl2(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               elseif isodd(ii) && iseven(jj) && isodd(kk) && iseven(ll)
                                  Hintdu[mm,nn] = Hintdu[mm,nn] + matV[n1+1,n2+1,n3+1,n4+1]*sqrt(vecmbnnijkl[Msize+1])
                                  # Hintdu[mm,nn] = Hintdu[mm,nn] + Vijkl2(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               elseif iseven(ii) && isodd(jj) && iseven(kk) && isodd(ll)
                                  Hintdu[mm,nn] = Hintdu[mm,nn] + matV[n1+1,n2+1,n3+1,n4+1]*sqrt(vecmbnnijkl[Msize+1])
                                  # Hintdu[mm,nn] = Hintdu[mm,nn] + Vijkl2(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               end

                            end

                        end

                    end

                end

            end

        end
    end

    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdown .= Hintdown + Hintdown' - spdiagm(diag(Hintdown))
    Hintup .= Hintup + Hintup' - spdiagm(diag(Hintup))
    Hintdu .= Hintdu + Hintdu' - spdiagm(diag(Hintdu))

    # return Hintdown, Hintup, Hintdu

end

using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("ades.jl")
include("acre.jl")
include("pascaltriangle.jl")
include("in2b.jl")
include("b2in.jl")
include("epsilon.jl")
include("vijkl.jl")

function Hintfunccutoff(indvec::Vector{Int64}, Msize0::Int64, Np::Int64)

    Msize = Msize0*2
    maxmatpcut = length(indvec)

    # matp = zeros(Int,Msize+1,Np+1);
    # pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    # defines vectors and matrices
    vecmbnn = spzeros(Int64,Msize+1);
    vecmbnnj = spzeros(Int64,Msize+1);
    vecmbnnij = spzeros(Int64,Msize+1);
    # Hsoc = spzeros(ComplexF64,maxmatp,maxmatp);
    Hintdown = spzeros(Float64,maxmatp,maxmatp);
    Hintup = spzeros(Float64,maxmatp,maxmatp);
    Hintdu = spzeros(Float64,maxmatp,maxmatp);

    # define a matrix for the Hamiltonian
    for nn = 1:maxmatpcut #maxmatp

        vecmbnn = in2b(indvec[nn],Msize,Np) #in2b(nn,Msize,Np)

        # Interactions
        for ll = 1:Msize

            vecmbnnl = ades(ll,vecmbnn)
            if vecmbnnl[Msize+1] == 0
               continue
            end

            for kk = 1:Msize

                vecmbnnkl = ades(kk,vecmbnnl)
                if vecmbnnkl[Msize+1] == 0
                   continue
                end

                for jj = 1:Msize

                    vecmbnnjkl = acre(jj,vecmbnnkl)

                    for ii = 1:Msize

                        vecmbnnijkl = acre(ii,vecmbnnjkl)

                        # if isodd(ii+jj+kk+ll)
                        #    continue
                        # end

                        for mm = 1:nn

                            vecmbmm = in2b(indvec[mm],Msize,Np)

                            if vecmbnnijkl[1:Msize] == vecmbmm[1:Msize]

                               n1 = ceil(Int64,ii/2)-1
                               n2 = ceil(Int64,jj/2)-1
                               n3 = ceil(Int64,kk/2)-1
                               n4 = ceil(Int64,ll/2)-1

                               if isodd(ii) && isodd(jj) && isodd(kk) && isodd(ll)
                                  Hintdown[mm,nn] = Hintdown[mm,nn] + Vijkl(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               elseif iseven(ii) && iseven(jj) && iseven(kk) && iseven(ll)
                                  Hintup[mm,nn] = Hintup[mm,nn] + Vijkl(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               elseif isodd(ii) && iseven(jj) && isodd(kk) && iseven(ll)
                                  Hintdu[mm,nn] = Hintdu[mm,nn] + Vijkl(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               elseif iseven(ii) && isodd(jj) && iseven(kk) && isodd(ll)
                                  Hintdu[mm,nn] = Hintdu[mm,nn] + Vijkl(n1,n2,n3,n4)*sqrt(vecmbnnijkl[Msize+1])
                               end

                            end

                        end

                    end

                end

            end

        end
    end

    Hintdown = Hintdown + Hintdown' - spdiagm(diag(Hintdown))
    Hintup = Hintup + Hintup' - spdiagm(diag(Hintup))
    Hintdu = Hintdu + Hintdu' - spdiagm(diag(Hintdu))

    return Hintdown, Hintup, Hintdu

end

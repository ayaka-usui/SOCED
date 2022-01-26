using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("ades.jl")
include("acre.jl")
# include("pascaltriangle.jl")
include("in2b.jl")
# include("b2in.jl")
# include("epsilon.jl")
include("vijkl2.jl")

function Hintfunccutoff2!(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Hintdown::ST, Hintup::ST, Hintdu::ST) where ST <: Union{SparseMatrixCSC{Float64},Array{Float64}}

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

    vecmbnn = zeros(Int64,Msize+1, Threads.nthreads())
    vecmbnnl = zeros(Int64,Msize+1, Threads.nthreads());
    vecmbnnkl = zeros(Int64,Msize+1, Threads.nthreads())
    vecmbnnjkl = zeros(Int64,Msize+1, Threads.nthreads())
    vecmbnnijkl = zeros(Int64,Msize+1, Threads.nthreads())
    vecmbmm = zeros(Int64,Msize+1,Threads.nthreads())

    # pre-calculate set of indices
    #index_set = zeros(Int, Int(ceil(Msize0/2)^4),4)
    index_set = zeros(Int, div(Msize0^4,2^3),4) .- 1
    curr_index = 1
    @time for nn4 = 0:Msize0-1
        for nn3 = 0:nn4
            for nn2 = 0:nn3
                for nn1 = 0:nn2

                    if iseven(nn1+nn2+nn3+nn4)
                        index_set[curr_index,1] = nn1
                        index_set[curr_index,2] = nn2
                        index_set[curr_index,3] = nn3
                        index_set[curr_index,4] = nn4
                        curr_index += 1
                    end

                end
            end
        end
    end

    # calculate Vijkl in advance
    matV = zeros(Float64,Msize0,Msize0,Msize0,Msize0)
    matA = zeros(Rational{BigInt},2,Msize0,
                 Msize0,Msize0,Threads.nthreads())
    @time Threads.@threads for i = 1:curr_index-1
        tid = Threads.threadid()
        nn1 = index_set[i,1]
        nn2 = index_set[i,2]
        nn3 = index_set[i,3]
        nn4 = index_set[i,4]

        maxm1 = floor(Int64,nn1/2)+1
        maxm2 = floor(Int64,nn2/2)+1
        maxm3 = floor(Int64,nn3/2)+1
        maxm4 = floor(Int64,nn4/2)+1


        matV[nn1+1,nn2+1,nn3+1,nn4+1] = Vijkl2(nn1,nn2,nn3,nn4, matA[1:2,
                                                                     1:maxm2,
                                                                     1:maxm3,
                                                                     1:maxm4,
                                                                     tid])
    end

    # resetting index_set
    index_set[1:curr_index,:] .= -1
    curr_index = 1

    # define a matrix for the Hamiltonian
    println("time for for loops")
    @time Threads.@threads for nn = 1:maxmatpcut
    #@time for nn = 1:maxmatpcut
        tid = Threads.threadid()

        # define ket state |n>
        vecmbnn[:,tid] .= zeros(Int64,Msize+1)
        in2b!(indvec[nn],Msize,Np,matp,vecmbnn,tid)
        #println(vecmbnn[:,tid])

        # Interactions
        for ll = 1:Msize

            # apply an anhilation operator a_ll to |n>
            vecmbnnl[:,tid] .= zeros(Int64,Msize+1);
            ades!(ll,vecmbnn[:,tid],vecmbnnl,tid)
            if vecmbnnl[Msize+1,tid] == 0 # go back if a|0> = 0
               continue
            end

            for kk = 1:Msize

                # apply a_kk to a_ll |n>
                vecmbnnkl[:,tid] .= zeros(Int64,Msize+1)
                ades!(kk,vecmbnnl[:,tid],vecmbnnkl,tid)
                if vecmbnnkl[Msize+1,tid] == 0
                   continue
                end

                for jj = 1:Msize

                    # take interactions between (down, down), (up,up), (down,up), or (up,down) and skip other things
                    if isodd(jj+ll)
                       continue
                    end

                    # apply a_jj^{+} to a_kk a_ll |n>
                    vecmbnnjkl[:,tid] .= zeros(Int64,Msize+1)
                    acre!(jj,vecmbnnkl[:,tid],vecmbnnjkl,tid)

                    for ii = 1:Msize

                        # apply a_ii^{+} to a_jj^{+} a_kk a_ll |n>
                        vecmbnnijkl[:,tid] .= zeros(Int64,Msize+1)
                        acre!(ii,vecmbnnjkl[:,tid],vecmbnnijkl,tid)

                        # take interactions between (down, down), (up,up), (down,up), or (up,down) and skip other things
                        if isodd(ii+kk) 
                           continue
                        end

                        # get the index of eigenstate of harmonic oscillator
                        n1 = ceil(Int64,ii/2)-1
                        n2 = ceil(Int64,jj/2)-1
                        n3 = ceil(Int64,kk/2)-1
                        n4 = ceil(Int64,ll/2)-1

                        # V_{ijkl}=0 if n1+n2+n3+n4 is odd
                        if isodd(n1+n2+n3+n4)
                           continue
                        end

                        # Threads.@threads for mm = 1:nn
                        for mm = 1:nn

                            vecmbmm[:,tid] .= zeros(Int64,Msize+1)
                            #vecmbmm = zeros(Int64,Msize+1)
                            in2b!(indvec[mm],Msize,Np,matp,vecmbmm,tid)

                            vec0 = [n1,n2,n3,n4]
                            sort!(vec0)
                            n1 = vec0[1]
                            n2 = vec0[2]
                            n3 = vec0[3]
                            n4 = vec0[4]

                            if vecmbnnijkl[1:Msize,tid] == vecmbmm[1:Msize,tid]

                               if isodd(ii) && isodd(jj)
                                  Hintdown[mm,nn] = Hintdown[mm,nn] + matV[n1+1,n2+1,n3+1,n4+1]*sqrt(vecmbnnijkl[Msize+1,tid])
                               elseif iseven(ii) && iseven(jj)
                                  Hintup[mm,nn] = Hintup[mm,nn] + matV[n1+1,n2+1,n3+1,n4+1]*sqrt(vecmbnnijkl[Msize+1,tid])
                               else
                                  Hintdu[mm,nn] = Hintdu[mm,nn] + matV[n1+1,n2+1,n3+1,n4+1]*sqrt(vecmbnnijkl[Msize+1,tid])
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

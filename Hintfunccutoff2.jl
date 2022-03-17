using Arpack, SparseArrays, LinearAlgebra, SharedArrays

# define functions used here
include("vijkl.jl")
include("in2bindtid.jl")

function coefficientInttid(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    vecmbind0 = 0
    # common = 0

    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    for pp = 1:Np

        # consider to edit when parfor is inplemented
        if vecmbindnn[pp] == vecmbind0
           continue
        end

        for qq = 1:Np

            if vecmbindnn[pp] == vecmbindmm[qq]

               vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)] # ll kk
               vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)] # jj ii
               common = vecmbindnn[pp]

               element = 1.0

               # coefficients of operators for ket
               if vecmbindnn3[1] == vecmbindnn3[2] # a_{kk} a_{kk}
                  if common == vecmbindnn3[1]
                     element = sqrt(Np*(Np-1))*element
                  else
                     element = sqrt(2*1)*element
                  end
               else # vecmbindnn3[1] != vecmbindnn3[2] # a_{kk} a_{ll}
                  if common == vecmbindnn3[1] || common == vecmbindnn3[2]
                     element = sqrt(2)*element
                  # else
                     # element = 1.0*element
                  end
                  element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
               end

               # coefficients of operators for bra
               if vecmbindmm3[1] == vecmbindmm3[2] # a^+_{ii} a^+_{ii}
                  if common == vecmbindmm3[1]
                     element = sqrt(2*3)*element
                  else
                     element = sqrt(1*2)*element
                  end
               else # vecmbindmm3[1] != vecmbindmm3[2] # a^+_{ii} a^+_{jj}
                  if common == vecmbindmm3[1] || common == vecmbindmm3[2]
                     element = sqrt(2)*element
                  # else
                     # element = 1.0*element
                  end
                  element = element*2 # a^+_{ii} a^+_{jj} + a^+_{jj} a^+_{ii}
               end

               # indices for Vijkl
               ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
               sort!(ind1,rev=true)
               ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

               ind0 += 1
               vecindcoeff[ind0,1] = ind2
               vecindcoeff[ind0,2] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    vecindcoeff[1,3] = ind0

    return vecindcoeff

end

function coefficientInt2tid(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # down down up for g12

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    vecmbind0 = 0
    # common = 0

    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    for pp = 1:Np-1

        # consider to edit when parfor is implemented
        if vecmbindnn[pp] == vecmbind0
           continue
        end

        for qq = 1:Np-1

            if vecmbindnn[pp] == vecmbindmm[qq]

               vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)] # ll kk
               vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)] # jj ii
               common = vecmbindnn[pp]

               element = 1.0

               # coefficients of operators for ket
               if common == vecmbindnn3[1]
                  element = sqrt(2)*element
               # else
               #    element = 1.0*element
               end

               # coefficients of operators for bra
               if common == vecmbindmm3[1]
                  element = sqrt(2)*element
               # else
               #    element = 1.0*element
               end

               # a^+_{down}a^+_{up}a_{down}a_{up} + a^+_{up}a^+_{down}a_{up}a_{down}
               element = 2*element

               # indices for Vijkl
               ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
               sort!(ind1,rev=true)
               ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

               ind0 += 1
               vecindcoeff[ind0,1] = ind2
               vecindcoeff[ind0,2] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    vecindcoeff[1,3] = ind0

    return vecindcoeff

end

function coefficientInt3tid(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # down down up for g

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    # vecmbind0 = 0
    # common = 0

    if vecmbindnn[end] != vecmbindmm[end]
       return vecindcoeff
    end

    vecmbindnn3 .= vecmbindnn[1:2] # ll kk
    vecmbindmm3 .= vecmbindmm[1:2] # jj ii

    element = 1.0

    if vecmbindnn3[1] == vecmbindnn3[2]
       element = sqrt(2)*element
    else # vecmbindnn3[1] != vecmbindnn3[2]
       element = 2*element
    end
    if vecmbindmm3[1] == vecmbindmm3[2]
       element = sqrt(2)*element
    else # vecmbindmm3[1] != vecmbindmm3[2]
       element = 2*element
    end

    # indices for Vijkl
    ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
    sort!(ind1,rev=true)
    ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

    ind0 += 1
    vecindcoeff[ind0,1] = ind2
    vecindcoeff[ind0,2] = element
    vecindcoeff[1,3] = ind0

    return vecindcoeff

end

function in2bindmixedspins(Msize0::Int64,Np::Int64,nn::Int64,indvec2::Vector{Int64},maxmatp21::Int64,matp20::Matrix{Int64},matp21::Matrix{Int64},vecmbindnn2::Vector{Int64},vecmbindnn::Vector{Int64})
   # indvec2[nn] = (nndown-1)*maxmatp21 + nnup
   indketup = mod(indvec2[nn],maxmatp21)
   if indketup != 0
      indketdown = div(indvec2[nn],maxmatp21) + 1
   else # indketup == 0
      indketdown = div(indvec2[nn],maxmatp21)
      indketup = maxmatp21
   end
   in2bind!(indketdown,Msize0,Np-1,matp20,vecmbindnn2)
   vecmbindnn[1:Np-1] = vecmbindnn2[1:Np-1]
   in2bind!(indketup,Msize0,1,matp21,vecmbindnn2)
   vecmbindnn[Np:Np] = vecmbindnn2[1:1]
   return vecmbindnn
end

function Hintfunccutoff2!(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, Hintdown::ST, Hintup::ST, Hintdu::ST) where ST<:(SparseMatrixCSC{Float64, Ti} where Ti<:Integer)

    if Msize0 > 150
       error("The memory for the elements of sparse matrice such may be too large. I have not checked it.")
    end

    maxmatpcut = length(indvec)
    maxmatpcut2 = length(indvec2)

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
    # Hintdown .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintup .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintdu .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);

    indrow_Hint = SharedArray{Int64,1}(binomial(maxmatpcut+1,2))
    indcolumn_Hint = SharedArray{Int64,1}(binomial(maxmatpcut+1,2))
    element_Hint = SharedArray{Float64,1}(binomial(maxmatpcut+1,2))

    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)
    vecmbindnn3 = zeros(Int64,2)
    vecmbindmm3 = zeros(Int64,2)
    vecindcoeff = zeros(Float64,3,2)

    tmax = Threads.nthreads()
    vecmbindnntid = zeros(Int64,Np,tmax)
    vecmbindmmtid = zeros(Int64,Np,tmax)
    vecmbindnn3tid = zeros(Int64,2,tmax)
    vecmbindmm3tid = zeros(Int64,2,tmax)
    vecindcoefftid = zeros(Float64,3,3,tmax)

    # define a matrix for the Hamiltonian for down down down
    println("time for down down down")
    @time Threads.@threads for nn = 1:maxmatpcut # parfor
    # Threads does not work on Julia 1.6 but does on Julia 1.7
    # for nn = 1:maxmatpcut

        tid = Threads.threadid()

        # ket
        vecmbindnntid[:,tid] .= in2bindtid(indvec[nn],Msize0,Np,matp,vecmbindnntid[:,tid])

        for mm = 1:nn

            # bra
            vecmbindmmtid[:,tid] .= in2bindtid(indvec[mm],Msize0,Np,matp,vecmbindmmtid[:,tid])

            vecindcoefftid[:,:,tid] .= coefficientInttid(vecmbindnntid[:,tid],vecmbindmmtid[:,tid],vecmbindnn3tid[:,tid],vecmbindmm3tid[:,tid],vecindcoefftid[:,:,tid],Np)
            # ind0 = Int64(vecindcoefftid[1,3,tid])
            # ind2 = Int64.(vecindcoefftid[1:Int64(vecindcoefftid[1,3,tid]),1,tid])
            # Hintdown[mm,nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

            indrow_Hint[binomial(nn,2)+mm] = mm
            indcolumn_Hint[binomial(nn,2)+mm] = nn
            # element_Hint[binomial(nn,2)+mm] = sum(vecV[ind2].*vecindcoefftid[1:ind0,2,tid])
            element_Hint[binomial(nn,2)+mm] = sum(vecV[Int64.(vecindcoefftid[1:Int64(vecindcoefftid[1,3,tid]),1,tid])].*vecindcoefftid[1:Int64(vecindcoefftid[1,3,tid]),2,tid])

        end

    end

    indrow_Hint = Vector(indrow_Hint)
    indcolumn_Hint = Vector(indcolumn_Hint)
    element_Hint = Vector(element_Hint)
    deleteat!(indrow_Hint, element_Hint .== 0.0)
    deleteat!(indcolumn_Hint, element_Hint .== 0.0)
    deleteat!(element_Hint, element_Hint .== 0.0)
    Hintdown .= sparse(indrow_Hint,indcolumn_Hint,element_Hint,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)

    # define Hint for down down up
    maxmatp21 = matp21[Msize0+1,1+1]
    vecmbindnn2tid = zeros(Int64,Np,tmax)
    vecmbindmm2tid = zeros(Int64,Np,tmax)

    indrow_Hint = SharedArray{Int64,1}(binomial(maxmatpcut2+1,2))
    indcolumn_Hint = SharedArray{Int64,1}(binomial(maxmatpcut2+1,2))
    element_Hint = SharedArray{Float64,1}(binomial(maxmatpcut2+1,2))

    indrow_Hint2 = SharedArray{Int64,1}(binomial(maxmatpcut2+1,2))
    indcolumn_Hint2 = SharedArray{Int64,1}(binomial(maxmatpcut2+1,2))
    element_Hint2 = SharedArray{Float64,1}(binomial(maxmatpcut2+1,2))

    println("time for down down up")
    @time Threads.@threads for nn = 1:maxmatpcut2 # parfor # down down up for ket
    # for nn = 1:maxmatpcut2

        tid = Threads.threadid()

        # ket
        vecmbindnntid[:,tid] .= in2bindmixedspins(Msize0,Np,nn,indvec2,maxmatp21,matp20,matp21,vecmbindnn2tid[:,tid],vecmbindnntid[:,tid])

        for mm = 1:nn # down down up for bra

            # bra
            vecmbindmmtid[:,tid] .= in2bindmixedspins(Msize0,Np,mm,indvec2,maxmatp21,matp20,matp21,vecmbindmm2tid[:,tid],vecmbindmmtid[:,tid])

            vecindcoefftid[:,:,tid] .= coefficientInt2tid(vecmbindnntid[:,tid],vecmbindmmtid[:,tid],vecmbindnn3tid[:,tid],vecmbindmm3tid[:,tid],vecindcoefftid[:,:,tid],Np)
            # ind0 = Int64(vecindcoefftid[1,3,tid])
            # ind2 = Int64.(vecindcoefftid[1:ind0,1,tid])
            # Hintdu[maxmatpcut+mm,maxmatpcut+nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

            indrow_Hint2[binomial(nn,2)+mm] = mm
            indcolumn_Hint2[binomial(nn,2)+mm] = nn
            element_Hint2[binomial(nn,2)+mm] = sum(vecV[Int64.(vecindcoefftid[1:Int64(vecindcoefftid[1,3,tid]),1,tid])].*vecindcoefftid[1:Int64(vecindcoefftid[1,3,tid]),2,tid])

            vecindcoefftid[:,:,tid] .= coefficientInt3tid(vecmbindnntid[:,tid],vecmbindmmtid[:,tid],vecmbindnn3tid[:,tid],vecmbindmm3tid[:,tid],vecindcoefftid[:,:,tid],Np)
            # ind0 = Int64(vecindcoefftid[1,3,tid])
            # ind2 = Int64.(vecindcoefftid[1:ind0,1,tid])
            # Hintdown[maxmatpcut+mm,maxmatpcut+nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

            indrow_Hint[binomial(nn,2)+mm] = mm
            indcolumn_Hint[binomial(nn,2)+mm] = nn
            element_Hint[binomial(nn,2)+mm] = sum(vecV[Int64.(vecindcoefftid[1:Int64(vecindcoefftid[1,3,tid]),1,tid])].*vecindcoefftid[1:Int64(vecindcoefftid[1,3,tid]),2,tid])

        end

    end

    indrow_Hint2 = Vector(indrow_Hint2)
    indcolumn_Hint2 = Vector(indcolumn_Hint2)
    element_Hint2 = Vector(element_Hint2)
    deleteat!(indrow_Hint2, element_Hint2 .== 0.0)
    deleteat!(indcolumn_Hint2, element_Hint2 .== 0.0)
    deleteat!(element_Hint2, element_Hint2 .== 0.0)
    Hintdu[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2] .= sparse(indrow_Hint2,indcolumn_Hint2,element_Hint2,maxmatpcut2,maxmatpcut2)

    indrow_Hint = Vector(indrow_Hint)
    indcolumn_Hint = Vector(indcolumn_Hint)
    element_Hint = Vector(element_Hint)
    deleteat!(indrow_Hint, element_Hint .== 0.0)
    deleteat!(indcolumn_Hint, element_Hint .== 0.0)
    deleteat!(element_Hint, element_Hint .== 0.0)
    Hintdown[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2] .= sparse(indrow_Hint,indcolumn_Hint,element_Hint,maxmatpcut2,maxmatpcut2)

    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdown .= Hintdown + Hintdown' - spdiagm(diag(Hintdown))
    Hintdu .= Hintdu + Hintdu' - spdiagm(diag(Hintdu))

    # define Hintup
    Hintup[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hintdown[1:maxmatpcut,1:maxmatpcut]
    Hintup[end-maxmatpcut-maxmatpcut2+1:end-maxmatpcut,end-maxmatpcut-maxmatpcut2+1:end-maxmatpcut] = Hintdown[1+maxmatpcut:maxmatpcut2+maxmatpcut,1+maxmatpcut:maxmatpcut2+maxmatpcut]

    # define Hint for up up down
    Hintdu[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2] = Hintdu[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2]

end

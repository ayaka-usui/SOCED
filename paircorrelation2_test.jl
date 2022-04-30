
using Arpack, SparseArrays, LinearAlgebra
using Polynomials, SpecialPolynomials

# define functions used here
# include("vijkl.jl")

function setfunHO(xrange::LinRange{Float64},Msize0::Int64)

    # y = variable(Polynomial{Rational{Int}})
    Nx = length(xrange)
    outcome = zeros(Float64,Nx,Msize0)
    dx = abs(xrange[2]-xrange[1])

    for nn = 1:Msize0
        for jjx = 1:Nx
            outcome[jjx,nn] = basis(Hermite,nn-1)(xrange[jjx])*exp(-xrange[jjx]^2/2)
        end
        outcome[:,nn] = outcome[:,nn]/sqrt(sum(abs.(outcome[:,nn]).^2)*dx)
    end

    return outcome

end

function coefficientpair(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff0::Matrix{Int64},vecindcoeff1::Vector{Float64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff0 .= 0 #zeros(Float64,3,4)
    vecindcoeff1 .= 0.0 #zeros(Float64,3)
    ind0 = 0
    vecmbind0 = 0
    # common = 0

    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    for pp = 1:Np

        # consider to edit when parfor is implemented
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
                  # element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
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
                  # element = element*2 # a^+_{ii} a^+_{jj} + a^+_{jj} a^+_{ii}
               end

               # a^+_{ii} a^+_{jj} a_{kk} a_{ll}
               ind1 = [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[2], vecmbindnn3[1]] # ii jj kk ll
               ind0 += 1
               vecindcoeff0[ind0,1:4] = ind1
               vecindcoeff1[ind0] = element

               # a^+_{ii} a^+_{jj} a_{ll} a_{kk}
               if vecmbindnn3[1] != vecmbindnn3[2]
                  ind1 .= [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[1], vecmbindnn3[2]]
                  ind0 += 1
                  vecindcoeff0[ind0,1:4] = ind1
                  vecindcoeff1[ind0] = element
               end

               # a^+_{jj} a^+_{ii} a_{kk} a_{ll}
               if vecmbindmm3[1] != vecmbindmm3[2]
                  ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[2], vecmbindnn3[1]]
                  ind0 += 1
                  vecindcoeff0[ind0,1:4] = ind1
                  vecindcoeff1[ind0] = element
               end

               # a^+_{jj} a^+_{ii} a_{ll} a_{kk}
               if vecmbindnn3[1] != vecmbindnn3[2] && vecmbindmm3[1] != vecmbindmm3[2]
                  ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[1], vecmbindnn3[2]]
                  ind0 += 1
                  vecindcoeff0[ind0,1:4] = ind1
                  vecindcoeff1[ind0] = element
               end

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff0, vecindcoeff1, ind0

end

function coefficientpair2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff0::Matrix{Int64},vecindcoeff1::Vector{Float64},Np::Int64)

    # down down up for g12

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff0 .= 0 #zeros(Float64,3,2)
    vecindcoeff1 .= 0.0 #zeros(Float64,3)
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
               end

               # coefficients of operators for bra
               if common == vecmbindmm3[1]
                  element = sqrt(2)*element
               end

               # a^+_{down}a^+_{up}a_{down}a_{up} + a^+_{up}a^+_{down}a_{up}a_{down}
               # element = 2*element

               # a^+_{down}a^+_{up}a_{down}a_{up}
               ind1 = [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[1], vecmbindnn3[2]]
               ind0 += 1
               vecindcoeff0[ind0,1:4] = ind1
               vecindcoeff1[ind0,end] = element

               # a^+_{up}a^+_{down}a_{up}a_{down}
               ind1 .= [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[2], vecmbindnn3[1]]
               ind0 += 1
               vecindcoeff0[ind0,1:4] = ind1
               vecindcoeff1[ind0,end] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff0, vecindcoeff1, ind0

end

function coefficientpair3(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff0::Matrix{Int64},vecindcoeff1::Vector{Float64},Np::Int64)

    # down down up for g

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff0 .= 0 #zeros(Float64,3,2)
    vecindcoeff1 .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    # vecmbind0 = 0
    # common = 0

    if vecmbindnn[end] != vecmbindmm[end]
       return vecindcoeff0, vecindcoeff1, ind0
    end

    vecmbindnn3 .= vecmbindnn[1:2] # ll kk
    vecmbindmm3 .= vecmbindmm[1:2] # jj ii

    element = 1.0

    # note that the order of the indices matters if vecmbindnn3[1] != vecmbindnn3[2]
    if vecmbindnn3[1] == vecmbindnn3[2]
       element = sqrt(2)*element
    end
    if vecmbindmm3[1] == vecmbindmm3[2]
       element = sqrt(2)*element
    end

    # a^+_{ii} a^+_{jj} a_{kk} a_{ll}
    ind1 = [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[2], vecmbindnn3[1]] # ii jj kk ll
    ind0 += 1
    vecindcoeff0[ind0,1:4] = ind1 #vecmbindnn3
    vecindcoeff1[ind0,end] = element

    # a^+_{ii} a^+_{jj} a_{ll} a_{kk}
    if vecmbindnn3[1] != vecmbindnn3[2]
       ind1 .= [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[1], vecmbindnn3[2]]
       ind0 += 1
       vecindcoeff0[ind0,1:4] = ind1
       vecindcoeff1[ind0,end] = element
    end

    # a^+_{jj} a^+_{ii} a_{kk} a_{ll}
    if vecmbindmm3[1] != vecmbindmm3[2]
       ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[2], vecmbindnn3[1]]
       ind0 += 1
       vecindcoeff0[ind0,1:4] = ind1
       vecindcoeff1[ind0,end] = element
    end

    # a^+_{jj} a^+_{ii} a_{ll} a_{kk}
    if vecmbindnn3[1] != vecmbindnn3[2] && vecmbindmm3[1] != vecmbindmm3[2]
       ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[1], vecmbindnn3[2]]
       ind0 += 1
       vecindcoeff0[ind0,1:4] = ind1
       vecindcoeff1[ind0,end] = element
    end

    return vecindcoeff0, vecindcoeff1, ind0

end

function paircorrelation_fun_test(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, psi::Vector{ComplexF64}, xrange::LinRange{Float64}, yrange::LinRange{Float64})

    maxmatpcut = length(indvec)
    maxmatpcut2 = length(indvec2)

    # defines vectors and matrices
    # Mpairdown3int .= 0
    # Mpairdown2up1int .= 0
    # Mpairdown1up2int .= 0
    # Mpairup3int .= 0
    # Mpairdown3float .= 0.0
    # Mpairdown2up1float .= 0.0
    # Mpairdown1up2float .= 0.0
    # Mpairup3float .= 0.0
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)
    vecmbindnn3 = zeros(Int64,2)
    vecmbindmm3 = zeros(Int64,2)
    vecindcoeff0 = zeros(Int64,3*4,4) #zeros(Int64,3,4)
    vecindcoeff1 = zeros(Float64,3*4) #zeros(Float64,3)

    # pair correlation
    Nx = length(xrange)
    Ny = length(yrange)
    # fun_nudown = zeros(Float64,Nx,Ny)
    # fun_nudu = zeros(Float64,Nx,Ny)
    # fun_nuup = zeros(Float64,Nx,Ny)
    phiHO = setfunHO(xrange,Msize0)
    # fun_nudown1 = zeros(Float64,Nx,Ny)
    # fun_nuup1 = zeros(Float64,Nx,Ny)
    fun_nu = zeros(Float64,Nx,Ny)
    dx = abs(xrange[2]-xrange[1])

    vecmbindnn = [1, 1, 1]
    vecmbindmm = [1, 1, 1]
    vecindcoeff0, vecindcoeff1, ind0 = coefficientpair(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff0,vecindcoeff1,Np)

    # println(vecindcoeff0[1:ind0,1:4])
    # println(vecindcoeff1[1:ind0])

    for jjx = 1:Nx
        for jjy = 1:Ny
            # for jj = 1:ind0
                fun_nu[jjx,jjy] += + 1/2*3*2*phiHO[jjx,1]*phiHO[jjy,1]*phiHO[jjx,1]*phiHO[jjy,1]
                # fun_nu[jjx,jjy] = fun_nu[jjx,jjy] + abs(psi[1])^2*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*phiHO[jjy,vecindcoeff0[jj,2]]*phiHO[jjx,vecindcoeff0[jj,3]]*phiHO[jjy,vecindcoeff0[jj,4]]
            # end
        end
    end

    println(sum(fun_nu)/Np/(Np-1)*dx^2)

    vecmbindnn = [1, 1, 2]
    vecmbindmm = [1, 1, 1]
    vecindcoeff0, vecindcoeff1, ind0 = coefficientpair(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff0,vecindcoeff1,Np)

    # println(vecindcoeff0[1:ind0,1:4])
    # println(vecindcoeff1[1:ind0])

    for jjx = 1:Nx
        for jjy = 1:Ny
            # for jj = 1:ind0
                fun_nu[jjx,jjy] += + 2/2*sqrt(3*2*2)*phiHO[jjx,1]*phiHO[jjy,1]*phiHO[jjx,1]*phiHO[jjy,2]
                # fun_nu[jjx,jjy] = fun_nu[jjx,jjy] + 2*real(conj(psi[1])*psi[2])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*phiHO[jjy,vecindcoeff0[jj,2]]*phiHO[jjx,vecindcoeff0[jj,3]]*phiHO[jjy,vecindcoeff0[jj,4]]
            # end
        end
    end

    println(sum(fun_nu)/Np/(Np-1)*dx^2)

    fun_nu .= fun_nu/Np/(Np-1)

    return fun_nu

end

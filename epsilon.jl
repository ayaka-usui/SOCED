include("delta.jl")

function epsilon(ii::Int64, jj::Int64, ksoc::Float64, Omega::Float64)

    energy0 = 0 # Int64
    njj = ceil(Int64,jj/2) - 1
    nii = ceil(Int64,ii/2) - 1

    if ii == jj

       energy0 = njj + 1/2 # Float64

    elseif iseven(ii+jj)

       energy0 = 1im*ksoc/sqrt(2)*(-1)^(ii)*(sqrt(njj+1)*delta(ii-jj-2) - sqrt(njj)*delta(ii-jj+2)) # Float64
       # with ii odd, (-1)^ii=-1 and ceil(jj/2)=(jj+1)/2 since jj is odd
       # with ii even, (-1)^ii=1 and ceil(jj/2)=jj/2 since jj is even

    elseif abs(ii-jj) == 1 && nii == njj # && isodd(ii+jj)

       energy0 = Omega/2.0 # Float64

    end

    return energy0

end

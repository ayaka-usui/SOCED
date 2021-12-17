function epsilon(ii::Int,jj::Int,Msize0::Int,ksoc::Float64,Omega::Float64)

    energy0 = 0.0

    if ii == jj

       energy0 = jj + 1/2

    elseif isodd(ii) && isodd(jj)

       energy0 = 1im*ksoc/sqrt(2)*(-1)*(sqrt(ceil(jj/2)+1)*)

    elseif iseven(ii) && iseven(jj)

    elseif abs(ii-jj) == 1 # && isodd(ii*jj)

       energy0 = Omega/2

    end

    return

end

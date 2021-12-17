function epsilon(ii::Int,jj::Int,Msize0::Int,ksoc::Float64)

    energy0 = 0.0

    if ii == jj
       energy0 = energy0 + jj + 1/2
    end

    if isodd(ii) && isodd(jj)

       energy0 = energy0 + 1im*ksoc/sqrt(2)*(-1)*(sqrt(ceil(jj/2))+1)

    end

    return

end

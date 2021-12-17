function epsilon(ii::Int,jj::Int,Msize0::Int,ksoc::Float64,Omega::Float64)

    function delta(x::Int)
      if x == 0
         return 1
      else
         return 0
    end

    energy0 = 0 # Int

    if ii == jj

       energy0 = jj + 1/2 # Float64

    elseif iseven(ii*jj)

       energy0 = 1im*ksoc/sqrt(2)*(-1)^(ii)*(sqrt(ceil(jj/2)+1)*delta(ii-jj-2) - sqrt(ceil(jj/2))*delta(ii-jj+2)) # Float64
       # with ii odd, (-1)^ii=-1 and ceil(jj/2)=(jj+1)/2 since jj is odd
       # with ii even, (-1)^ii=1 and ceil(jj/2)=jj/2 since jj is even

    elseif abs(ii-jj) == 1 # && isodd(ii*jj)

       energy0 = Omega/2.0 # Float64

    end

    return energy0

end

function recurringA0(M2::Int64,m1::Int64,m2::Int64,m3::Int64,m4::Int64)
    # return (M-m1-m2-m3-m4+4)/(2M-2m1-2m2-2m3-2m4+8)/(2M-2m1-2m2-2m3-2m4+7)
    return BigFloat(1)/BigFloat(M2-2m1-2m2-2m3-2m4+7) #BigFloat(1/(M2-2m1-2m2-2m3-2m4+7))
end

function recurringA1(n1::Int64,n2::Int64,n3::Int64,n4::Int64,m1::Int64,m2::Int64,m3::Int64,m4::Int64)

    # M = Int64((n1 + n2 + n3 + n4)/2) # integer

    if m1 == m2 == m3 == 0
       coeff = (-1)*BigFloat(n4-2m4+2)*BigFloat(n4-2m4+1)/BigFloat(m4)*recurringA0(n1+n2+n3+n4,1,1,1,m4)
    elseif m1 == m2 == 0
       coeff = (-1)*BigFloat(n3-2m3+2)*BigFloat(n3-2m3+1)/BigFloat(m3)*recurringA0(n1+n2+n3+n4,1,1,m3,m4)
    elseif m1 == 0
       coeff = (-1)*BigFloat(n2-2m2+2)*BigFloat(n2-2m2+1)/BigFloat(m2)*recurringA0(n1+n2+n3+n4,1,m2,m3,m4)
    else
       coeff = (-1)*BigFloat(n1-2m1+2)*BigFloat(n1-2m1+1)/BigFloat(m1)*recurringA0(n1+n2+n3+n4,m1,m2,m3,m4)
    end

    return coeff

end

function Vijkl2(n1::Int64,n2::Int64,n3::Int64,n4::Int64)

    # function to calculate V_{ijkl}

    if isodd(n1 + n2 + n3 + n4)
       return 0.0
    end

    M = Int64((n1 + n2 + n3 + n4)/2)
    maxm1 = floor(Int64,n1/2)+1
    maxm2 = floor(Int64,n2/2)+1
    maxm3 = floor(Int64,n3/2)+1
    maxm4 = floor(Int64,n4/2)+1
    matA = zeros(BigFloat,2,maxm2,maxm3,maxm4) #zeros(Float64,2,maxm2,maxm3,maxm4)

    # m1=m2=m3=m4=1
    matA[1,1,1,1] = BigFloat(1.0)
    # matA[1,1,1,1] = exp(-2M*log(2)+sum(log.(M+1:2*M))-1/2*(sum(log.(2:n1))+sum(log.(2:n2))+sum(log.(2:n3))+sum(log.(2:n4))))

    # m1=m2=m3=1
    for m4 = 1:maxm4-1
        # matA[1,1,1,m4+1] = matA[1,1,1,m4]*(-1)*2*(n4-2m4+2)*(n4-2m4+1)/m4*(M-m4+1)/(2M-2m4+2)/(2M-2m4+1)
        matA[1,1,1,m4+1] = matA[1,1,1,m4]*recurringA1(n1,n2,n3,n4,0,0,0,m4)
    end

    # m1=m2=1
    for m3 = 1:maxm3-1
        # Threads.@threads for m4 = 1:maxm4 # parfor
        for m4 = 1:maxm4
            # matA[1,1,m3+1,m4] = matA[1,1,m3,m4]*(-1)*2*(n3-2m3+2)*(n3-2m3+1)/m3*(M-m3-m4+2)/(2M-2m3-2m4+4)/(2M-2m3-2m4+3)
            matA[1,1,m3+1,m4] = matA[1,1,m3,m4]*recurringA1(n1,n2,n3,n4,0,0,m3,m4)
        end
    end

    # m1=1
    for m2 = 1:maxm2-1
        # Threads.@threads for m34 = 1:maxm3*maxm4 # parfor
        for m34 = 1:maxm3*maxm4
            m3 = div(m34,maxm4)+1
            m4 = mod(m34,maxm4) # m34-div(m34,maxm4)*maxm4
            if m4 == 0
               m3 = m3-1 # div(m34,maxm4)
               m4 = maxm4 # m34-(div(m34,maxm4)-1)*maxm4 = m34-(m3-1)*maxm4 = maxm4
            end
            # matA[1,m2+1,m3,m4] = matA[1,m2,m3,m4]*(-1)*2*(n2-2m2+2)*(n2-2m2+1)/m2*(M-m2-m3-m4+3)/(2M-2m2-2m3-2m4+6)/(2M-2m2-2m3-2m4+5)
            matA[1,m2+1,m3,m4] = matA[1,m2,m3,m4]*recurringA1(n1,n2,n3,n4,0,m2,m3,m4)
        end
    end

    sumA = sum(matA[1,:,:,:])
    for m1 = 1:maxm1-1

        # Threads.@threads for m234 = 1:maxm2*maxm3*maxm4 # parfor
        for m234 = 1:maxm2*maxm3*maxm4
            m2 = div(m234,maxm3*maxm4)+1
            m3 = div(m234 - (m2-1)*maxm3*maxm4,maxm4)+1
            m4 = m234 - (m2-1)*maxm3*maxm4 - (m3-1)*maxm4
            if mod(m234,maxm3*maxm4) == 0
               m2 = m2-1
               m3 = maxm3
               m4 = maxm4
            elseif mod(m234 - (m2-1)*maxm3*maxm4,maxm4) == 0
               m3 = m3-1
               m4 = maxm4
            end
            # matA[2,m2,m3,m4] = matA[1,m2,m3,m4]*(-1)*2*(n1-2m1+2)*(n1-2m1+1)/m1*(M-m1-m2-m3-m4+4)/(2M-2m1-2m2-2m3-2m4+8)/(2M-2m1-2m2-2m3-2m4+7)
            matA[2,m2,m3,m4] = matA[1,m2,m3,m4]*recurringA1(n1,n2,n3,n4,m1,m2,m3,m4)
        end

        sumA = sumA + sum(matA[2,:,:,:])
        matA[1,:,:,:] = matA[2,:,:,:]
        matA[2,:,:,:] .= 0

    end

    sumA = sumA*exp(-2M*log(2)+sum(log.(M+1:2*M))-1/2*(sum(log.(2:n1))+sum(log.(2:n2))+sum(log.(2:n3))+sum(log.(2:n4))))

    # if sumA < 0
    #    sumA = -sumA
    #    sumA = log(sumA) -2M*log(2)+sum(log.(M+1:2*M))-1/2*(sum(log.(2:n1))+sum(log.(2:n2))+sum(log.(2:n3))+sum(log.(2:n4)))
    #    sumA = exp(sumA)
    #    sumA = -sumA
    # else
       # sumA = log(sumA) -2M*log(2)+sum(log.(M+1:2*M))-1/2*(sum(log.(2:n1))+sum(log.(2:n2))+sum(log.(2:n3))+sum(log.(2:n4)))
       # sumA = exp(sumA)
    # end

    sumA = sumA/sqrt(2*pi)/2

    return sumA
    # return matA

end

function A0000(n1::Int64,n2::Int64,n3::Int64,n4::Int64)
   M = Int64((n1 + n2 + n3 + n4)/2)
   A0 = -2M*log(2)+sum(log.(M+1:2*M))-1/2*(sum(log.(2:n1))+sum(log.(2:n2))+sum(log.(2:n3))+sum(log.(2:n4)))
   return exp(A0)/sqrt(2*pi)/2
end

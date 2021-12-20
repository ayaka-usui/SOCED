function Vijkl(ii::Int64,jj::Int64,kk::Int64,ll::Int64,Msize::Inf64)

    n1 = ii
    n2 = jj
    n3 = kk
    n4 = ll
    if isodd(n1 + n2 + n3 + n4)
       return 0
    end

    M = (n1 + n2 + n3 + n4)/2
    matA = zeros(2,Msize,Msize,Msize)
    maxm1 = ceil(Int64,n1/2)+1
    maxm2 = ceil(Int64,n2/2)+1
    maxm3 = ceil(Int64,n3/2)+1
    maxm4 = ceil(Int64,n4/2)+1

    # m1=m2=m3=m4=1
    matA[1,1,1,1] = exp(-2M*log(2)+sum(log.(M+1:2*M))-1/2*(log.(2:n1)+log.(2:n2)+log.(2:n3)+log.(2:n4)))

    # m1=m2=m3=1
    for m4 = 1:maxm4-1
        matA[1,1,1,m4+1] = matA[1,1,1,m4]*2*(n4-2m4+2)*(n4-2m4+1)/m4*(M-m4+1)/(2M-2m4+2)/(2M-2m4+1)
    end

    # m1=m2=1
    for m3 = 1:maxm3-1
        for m4 = 1:maxm4 #parfor
            matA[1,1,m3+1,m4] = matA[1,1,m3,m4]*2*(n3-2m3+2)*(n3-2m3+1)/m3*(M-m3-m4+2)/(2M-2m3-2m4+4)/(2M-2m3-2m4+3)
        end
    end

    # m1=1
    for m2 = 1:maxm2-1
        for m34 = 1:maxm3*maxm4 # parfor
            m3 = div(m34,maxm4)+1
            m4 = mod(m34,maxm4) # m34 - div(m34,maxm4)*maxm4
            if m4 == 0
               m3 = m3-1 #div(m34,maxm4)
               m4 = maxm4
            end
            matA[1,m2+1,m3,m4] = matA[1,m2,m3,m4]*2*(n2-2m2+2)*(n2-2m2+1)/m2*(M-m2-m3-m4+3)/(2M-2m2-2m3-2m4+6)/(2M-2m2-2m3-2m4+5)
        end
    end

    for m1 = 1:maxm1
        for m234 = 1:maxm2*maxm3*maxm4 # parfor
            m2 = div(m234,maxm3*maxm4)+1
            m3 = div(m234 - (m2-1)*maxm3*maxm4,maxm4)+1
            m4 = m234 - (m2-1)*maxm3*maxm4 - (m3-1)*maxm4
            if m4 == 0
               m3 = m3-1
               m4 = maxm4
            end
        end
    end

end

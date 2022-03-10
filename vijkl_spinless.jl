using SpecialFunctions

function vijkl_spinless(a::Int64,b::Int64,c::Int64,d::Int64) # Abel's approach

    if isodd(a+b+c+d)
       return 0.0
    end

    ind0 = [a, b, c, d]
    sort!(ind0,rev=true)
    a = ind0[1]
    b = ind0[2]
    c = ind0[3]
    d = ind0[4]

    integral = 0.0

    for r = 0:d

        logf_r = -1/2*log(2)-2*log(pi)
        logf_r += 1/2*(logfactorial(c)+logfactorial(d)-logfactorial(a)-logfactorial(b))
        logf_r += -logfactorial(d-r)-logfactorial(c-d+r)-logfactorial(r)

        loggam1 = logabsgamma((a+b-c+d+1)/2-r)
        loggam2 = logabsgamma((a-b+c-d+1)/2+r)
        loggam3 = logabsgamma((-a+b+c-d+1)/2+r)
        logf_r += loggam1[1] + loggam2[1] + loggam3[1]
        phase_r = loggam1[2] * loggam2[2] * loggam3[2]

        integral += exp(logf_r)*phase_r

   end

   return integral/2

end

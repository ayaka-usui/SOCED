function correctionint(Msize0::Int64,g0::Float64,lambda0::Float64)

    nu = (lambda0 - 1)/2 # lambda0 is the lowest energy E_{COM}+E_{rel}, where E_{COM} = 1/2 and E_{rel} = 1/2+2nu
    matcalN = floor((Msize0-1)/2)

    invgc = 1/sqrt(nu) * log((sqrt(matcalN+1)+sqrt(nu))/(sqrt(matcalN+1)-sqrt(nu))) +
            1/2/sqrt(matcalN+1)/(matcalN+1-nu) +
            1/24/sqrt((matcalN+1)^3)/(matcalN+1-nu) +
            1/12/sqrt(matcalN+1)/(matcalN+1-nu)^2
    invgc = invgc/2/sqrt(2)/pi

    return g0/(1-g0*invgc)

end

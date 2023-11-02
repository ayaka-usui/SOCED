
function test_ver0(arrayspect)

    N = length(arrayspect)

    for jj = 1:N

        if arrayspect[jj] > 0.01
           println(jj)
           break
        end

    end

end
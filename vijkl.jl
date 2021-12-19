# function gammahalf(nn::Int64)
#
#    loggamma = log(pi)/2 + sum(log.(1:2*nn)) - 2*nn*log(2) - sum(log.(1:nn))
#    return exp(loggamma)
#
# end

function Ixx(ii::Int64, jj::Int64 ,kk::Int64, ll::Int64)

   for m1 = 1:floor(Int64,ii/2)
       for m2 = 1:floor(Int64,jj/2)
           for m3 = 1:floor(Int64,kk/2)
               for m4 = 1:floor(Int64,ll/2)

                   if isodd(ii+jj+kk+ll)
                      return
                   else # iseven(ii+jj+kk+ll)
                      integral = 1/sqrt(2)/2^((ii+jj+kk+ll)/2-m1-m2-m3-m4)*gammahalf((ii+jj+kk+ll)/2-m1-m2-m3-m4)
                   end

               end
           end
       end
   end

end

function Vijkl(ii::Int64,jj::Int64,kk::Int64,ll::Int64)


end

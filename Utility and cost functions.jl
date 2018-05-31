module utilityfunctions
export u, disut, disutprime, invdisutprime
#utility from consumption

function u(c)
    log(c)
end

#disutility from search effort
function disut(e::Real)
 return(e^2)
end

#derivative of this utility function
function disutprime(e::Real)
    return(2*e)
end

#inverse of the utility function
function invdisutprime(x::Real)
    return(1/2 * x)
end

end

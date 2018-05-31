module utilityfunctions
export u, disut, disutprime, invdisutprime, weights
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

function weights(e::Real, p_e_e::Real, p_u_u::Real; λ::Real = 1)
    param = [e, p_e_e, p_u_u, λ]
    for p in param
        if p .> 1
            error("$p if higher than one")
        end
        if p .< 0
            error("$p is lower than zero")
        end
    end
    ω_1 = λ * e * p_e_e
    ω_2 = λ * e * (1 - p_e_e)
    ω_3 = (1 - λ * e) * (1 - p_u_u)
    ω_4 = (1 - λ * e) * p_u_u
    return([ω_1, ω_2, ω_3,ω_4])
end

end

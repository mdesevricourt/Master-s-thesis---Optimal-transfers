module Searcheffort
using utilityfunctions
export Optimeffort

Default = Dict("n_e" => 0.5, "p_e|e" => 0.9, "p_u|u" => 0.9, "wage" => 600.0, "home production" => 50, "Budget" => 300.0, "Market tightness" => 0.8)
Default["tax rate"] = 0.2
Default["net wage"] = Default["wage"] * (1 - Default["tax rate"])

function expectu(e::Real, T::Real = 0, τ::Real = 0; wb::Real = Default["wage"], wn::Real = Default["net wage"], δ::Real = Default["home production"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], λ::Real = Default["Market tightness"])
    ω = weights(e, p_e_e, p_u_u; λ = λ)
    return(ω[1] * u(wn + T) + ω[2] * u(wb + T + τ) + ω[3] * u(δ + T) + ω[4] * u(δ + T + τ) - disut(e))
end

function Optimeffort(T::Real, τ::Real = 0; δ::Real = Default["home production"], wb::Real = Default["wage"], wn::Real = Default["net wage"], p_e_e = Default["p_e|e"], p_u_u = Default["p_u|u"], λ::Real  = Default["Market tightness"])
    #= I figured that an analytical solution was possible. So I implemented it. There are some weird cases where a numerical solution does stlightly better, probably
    due to computational approximiations, but the difference is so small it's s barely noticeable. I decided to get rid of this.
    =#
    e1 = invdisutprime(λ * (p_e_e * u(wn + T) + (1 - p_e_e) * u(wb + T + τ) - (1 - p_u_u) * u(δ + T) - p_u_u * u(δ + T + τ)))
    #enforcing the constraint on e1
    if e1 > 1
        e1 = 1
    elseif e1 < 0
        e1 = 0
    end

    #Checking for corner solutions

    if expectu(0, T, τ, wb =  wb, wn = wn, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ) > expectu(1, T, τ,  wb =  wb, wn = wn, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ)
        if expectu(0, T, τ, wb =  wb, wn = wn, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ) > expectu(e1, T, τ, wb =  wb, wn = wn, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ)
            return(0)
        else
            return(e1)
            end
    elseif expectu(1, T, τ, wb =  wb, wn = wn, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ) > expectu(e1, T, τ, wb =  wb, wn = wn, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ)
        return(1)
    else
        return(e1)
    end
end

end

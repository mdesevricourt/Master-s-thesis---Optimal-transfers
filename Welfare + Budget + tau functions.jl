module otherfunctions
using utilityfunctions, Optimeffort
export Welfare, budgetcons, budgetplot, τofT, netwage

Default = Dict("n_e" => 0.5, "p_e|e" => 0.9, "p_u|u" => 0.9, "wage" => 600.0, "home production" => 50, "Budget" => 300.0, "Market tightness" => 0.8)
Default["tax rate"] = 0.2

function netwage(gw::Real = Default["wage"] , ρ::Real = Default["tax rate"])
    return(gw * (1- ρ))
end

function Welfare(T::Real, τ::Real; e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], wb::Real = Default["wage"], ρ::Real = Default["tax rate"] , δ::Real = Default["home production"]; λ::Real =Default["Market tightness"])
    wn = netwage(wb, ρ)
    ω = weights(e, p_e_e, p_u_u, λ = λ)
    return(ω[1] * u(wn + T) + ω[2] * u(wb + T + τ) + ω[3] * u(δ + T) + ω[4] * u(δ + T + τ) - disut(e))
end

function budget(;ρ::Real = Default["tax rate"], gw::Real = Default["wage"], e::Real = 1, p_e_e::Real = Default["p_e|e"])
    return(ρ * gw * e * p_e_e)
end

function budgetcons(T::Real, τ::Real; gw = Default["wage"], ρ::Real = Default["tax rate"], n_e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], λ = Default["Market tightness"])
    nw = netwage(gw, ρ)
    ω = weights(n_e, p_e_e, p_u_u, λ = λ)
    B = budget(ρ = ρ, gw = gw, e = n_e, p_e_e = p_e_e)
    if (T + (ω[2] + ω[4]) * τ) > B
        return(false)
    else
        return(true)
    end
end

function τofT(T::Real; gw::Real = Default["wage"], ρ::Real = Default["tax rate"], δ::Real = Default["home production"], p_e_e = Default["p_e|e"], p_u_u = Default["p_u|u"], λ::Real  = Default["Market tightness"])
    nw = netwage(gw, ρ)
    ω = weights(1, p_e_e, p_u_u, λ = λ)
    τ = 0
    counter = 0
    lb = 0
    ub = (budget(ρ = ρ, gw = gw, e = 1, p_e_e = p_e_e) - T)/ (ω[2] + ω[4] )
    e = 0
    B = 0
    while ub - lb > 1
        for t in [lb, (lb + ub)/2, ub]
            counter += 1
            e = Optimeffort(T, t; wb = gw, wn = nw , δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ )
            B = budget(ρ = ρ, gw = gw, e = e, p_e_e = p_e_e)
            ω = weights(e, p_e_e, p_u_u; λ = λ)
            if t <= (B - T)/ (ω[2] + ω[4] )
                lb = t
                continue
            else
                ub = t
                break
            end
        end
    end
    println(counter)
    for t in Array(linspace(lb,ub,11))
        counter += 1
        e = Optimeffort(T, t;  wb = gw, wn = nw , δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
        B = budget(ρ = ρ, gw = gw, e = e, p_e_e = p_e_e)
        ω = weights(e, p_e_e, p_u_u; λ = λ)
        if t <= (B - T)/ (ω[2] + ω[4] )
            τ = t
            continue
        else
            break
        end
    end
    return([τ, Optimeffort(T, τ; wb = gw, wn = nw , δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)])
end


end

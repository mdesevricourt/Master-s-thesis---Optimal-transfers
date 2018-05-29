using Optim
using Plots
using Distributions
using JLD2
using LaTeXStrings
#push!(LOAD_PATH, "C:/Users/Maxime/.julia/v0.6/")
gr()
Default = Dict("n_e" => 0.5, "p_e|e" => 0.9, "p_u|u" => 0.9, "wage" => 600.0, "home production" => 50, "Budget" => 300.0, "Market tightness" => 0.8)

function u(c)
    log(c)
end

function Welfare(T::Real, τ::Real , e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], w::Real = Default["wage"], δ::Real = Default["home production"]; λ::Real =Default["Market tightness"])
    ω = weights(e, p_e_e, p_u_u, λ = λ)
    return(ω[1] * u(w + T) + ω[2] * u(w + T + τ) + ω[3] * u(δ + T) + ω[4] * u(δ + T + τ) - disut(e))
end


function weights(e::Real, p_e_e::Real, p_u_u::Real; λ::Real = Default["Market tightness"])
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

#disutility from search effort
function disut(e::Real)
 return(e^2)
end

function disutprime(e::Real)
    return(2*e)
end
function invdisutprime(x::Real)
    return(1/2 * x)
end

function expectu(e::Real, T::Real = 0, τ::Real = 0; w::Real = Default["wage"], δ::Real = Default["home production"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], λ::Real = Default["Market tightness"])
    ω = weights(e, p_e_e, p_u_u; λ = λ)
    return(ω[1] * u(w + T) + ω[2] * u(w + T + τ) + ω[3] * u(δ + T) + ω[4] * u(δ + T + τ) - disut(e))
end



function Optimeffort(T::Real, τ::Real = 0; w::Real = Default["wage"], δ::Real = Default["home production"], p_e_e = Default["p_e|e"], p_u_u = Default["p_u|u"], λ::Real  = Default["Market tightness"])
    #= I figured that an analytical solution was possible. So I implemented it. There are some weird cases where a numerical solution does stlightly better, probably
    due to computational approximiations, but the difference is so small it's s barely noticeable. I decided to get rid of this.
    =#
    e1 = invdisutprime(λ * (p_e_e * u(w + T) + (1 - p_e_e) * u(w + T + τ) - (1 - p_u_u) * u(δ + T) - p_u_u * u(δ + T + τ)))
    #enforcing the constraint on e1
    if e1 > 1
        e1 = 1
    elseif e1 < 0
        e1 = 0
    end

    #Checking for corner solutions

    if expectu(0, T, τ, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ) > expectu(1, T, τ,  w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ)
        if expectu(0, T, τ,  w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ) > expectu(e1, T, τ, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ)
            return(0)
        else
            return(e1)
            end
    elseif expectu(1, T, τ, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ) > expectu(e1, T, τ, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = λ)
        return(1)
    else
        return(e1)
    end
end

function budgetcons(T::Real, τ::Real; n_e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], B::Real = Default["Budget"], λ = Default["Market tightness"])
    ω = weights(n_e, p_e_e, p_u_u, λ = λ)
    if (T + (ω[2] + ω[4]) * τ) > B
        return(false)
    else
        return(true)
    end
end

function τofT(T::Real; w::Real = Default["wage"], δ::Real = Default["home production"], B = Default["Budget"], p_e_e = Default["p_e|e"], p_u_u = Default["p_u|u"], λ::Real  = Default["Market tightness"])
    if T == B
        return([0, Optimeffort(T, 0; w = w , δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)])
    end
    mini = B - T
    e = 1
    ω = weights(e, p_e_e, p_u_u, λ = λ)
    UB = (B - T)/ (ω[2] + ω[4] )
    τ = mini
    counter = 0
    lb = 0
    ub = (B - T)/ (ω[2] + ω[4] )
    for i in 1:6
        for t in Array(linspace(lb,ub, 5))
            counter += 1
            e = Optimeffort(T, t; w = w , δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ )
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
    for t in Array(linspace(lb,ub, max(floor(ub - lb), 10)))
        counter += 1
        e = Optimeffort(T, t; w = w , δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ )
        ω = weights(e, p_e_e, p_u_u; λ = λ)
        if t <= (B - T)/ (ω[2] + ω[4] )
            τ = t
            continue
        else
            break
        end
    end
    return([τ, Optimeffort(T, τ; w = w , δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ )])
end



function Ω(T::Real; p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], B::Real = Default["Budget"], w::Real = Default["wage"], δ::Real = Default["home production"], λ = Default["Market tightness"])
    τ, e = τofT(T, w = w , p_e_e = p_e_e, p_u_u = p_u_u, B = B, λ = λ)
    if (T < 0) || (τ < 0) || !budgetcons(T, τ; n_e = e, p_e_e = p_e_e, p_u_u = p_u_u, B = B, λ = λ)
        return(- Inf)
    else
        return(Welfare(T, τ, e, p_e_e, p_u_u, w, δ, λ = λ))
    end
end



function Solve(; p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], B::Real = Default["Budget"], w::Real = Default["wage"], δ::Real = Default["home production"], λ = Default["Market tightness"])
    Param = Dict{String, Real}(
    "p_e|e" => p_e_e,
    "p_u|u" => p_u_u,
    "wage" => w,
    "home production" => δ,
    "Budget" => float(B),
    "Market tightness" => λ)

    function Obj(T::Real)
        return(- Ω(T,  p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], B = Param["Budget"], w = Param["wage"], δ = Param["home production"], λ = λ))
    end

    Results = optimize(Obj, 0, B, Brent())


    W1 = Ω(0, p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], B = Param["Budget"], w = Param["wage"], δ = Param["home production"], λ = λ)
    W2 = Ω(B, p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], B = Param["Budget"], w = Param["wage"], δ = Param["home production"], λ = λ)

    if W1 >= - Results.minimum
        T = 0
        W = W1
    elseif W2 >= - Results.minimum
        T = B
        W = W2
    else
        T = Results.minimizer
        W = - Results.minimum
    end

    τ, e = τofT(T, w = Param["wage"], δ = Param["home production"], B = Param["Budget"], p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], λ = Param["Market tightness"])
    Param["Optimal T"] = T
    Param["Optimal τ"] = τ
    Param["e"] = e
    Param["Total Welfare"] = W
    return(Param)
end


function Tofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["Optimal T"])
end

function tofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["Optimal τ"])
end

function eofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["e"])
end



function p_e_eplots(;lb::Real=0, ub::Real=1)
    p1 = plot(Tofp_e_e, lb, ub, label = "Universal Transfer")
    plot!(p1, tofp_e_e, lb, ub, label = "Unemployed benefit")
    p2 =  plot(eofp_e_e, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end

#p_u_u

function Tofp_u_u(p_u_u::Real)
    return(Solve(p_u_u = p_u_u)["Optimal T"])
end

function tofp_u_u(p_u_u::Real)
    return(Solve(p_u_u = p_u_u)["Optimal τ"])
end

function eofp_u_u(p_u_u::Real)
    return(Solve(p_u_u = p_u_u)["e"])
end


function p_u_uplots(lb::Real=0, ub::Real=1)

    p1 = plot(Tofp_u_u, lb, ub, label = "Universal Transfer")
    plot!(p1, tofp_u_u, lb, ub, label = "Unemployed benefit")
    p2 =  plot(eofp_u_u, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end

#Wage
function Tofwage(w::Real)
    return(Solve(w = w)["Optimal T"])
end

function tofwage(w::Real)
    return(Solve(w = w)["Optimal τ"])
end

function eofwage(w::Real)
    return(Solve(w = w)["e"])
end


function wageplot(lb::Real=0, ub::Real=5000)

    p1 = plot(Tofwage, lb, ub, label = "Universal Transfer")
    plot!(p1, tofwage, lb, ub, label = "Unemployment benefit")
    p2 =  plot(eofwage, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end


#Budget

function TofB(T::Real)
    return(Solve(B = T)["Optimal T"])
end

function tofB(T::Real)
    return(Solve(B = T)["Optimal τ"])
end

function eofB(T::Real)
    return(Solve(B = T)["e"])
end



function budgetplot(lb::Real=0, ub::Real=50)

    p1 = plot(TofB, lb, ub, label = "Universal Transfer")
    info("Universal transfer done")
    plot!(p1, tofB, lb, ub, label = "Unemployment benefit")
    info("Categorical benefit done")
    p2 =  plot(eofB, lb, ub, label = "Search Effort")
    info("Search effort done")
    p = plot(p1, p2)
    return(p)
end


function budgetplot2(lb::Real = 0,  ub::Real = 100)

    T = []
    τ = []
    e = []
    Budget = linspace(lb, ub, 1000)
    for i in Budget
        Dico = Solve(B = i)
        push!(T, Dico["Optimal T"])
        push!(τ, Dico["Optimal τ"])
        push!(e, Dico["e"])
    end
    p1 = plot(Budget, T, label ="Universal Transfer")
    plot!(p1, Budget, τ, label = "Unemployment Benefit")
    p2 = plot(Budget, e, label = "Search effort")
    plot(p1,p2)
end
#gr()
#plotlyjs()

# plot(x,y,z, xlabel = "Transfer to employed", ylabel = "Prob of working", zcolor=reverse(z), m=(10,0.8,:blues,stroke(0)),leg=false,cbar=true,w=5)
function wageplot(;lb::Real=0, ub::Real=5000)
    p1 = plot(Tofwage, lb, ub, label = "Universal Transfer")
    plot!(p1, tofwage, lb, ub)
    p2 =  plot(eofwage, lb, ub, label = "Unemployment Benefit")
    p = plot(p1, p2)
    return(p)
end

function Tofdelta(w::Real)
    return(Solve(δ = w)["Optimal T"])
end

function tofdelta(w::Real)
    return(Solve(δ= w)["Optimal τ"])
end

function eofdelta(w::Real)
    return(Solve(δ = w)["e"])
end


function homeprodplot(lb::Real=0, ub::Real=Default["wage"])

    p1 = plot(Tofdelta, lb, ub, label = "Universal Transfer")
    plot!(p1, tofdelta, lb, ub, label = "Unemployment benefit")
    p2 =  plot(eofdelta, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end

#module Baseline
#Version of the model with two transfers: one categorical and one universal.
#Importing modules
using Optim
using Plots
#using ApproxFun

export u, Welfare, weights, τofT,  Ω, Solve, Obj

Default = Dict("n_e" => 0.5, "p_e|e" => 1.0, "p_u|u" => 1.0, "wage" => 600.0, "home production" => 50.0, "Budget" => 300.0, "Market tightness" => 1)

#Utility
function u(c)
    log(c)
end

#This is changed!!!
function Welfare(T::Real, τ::Real , n_e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], w::Real = Default["wage"], δ::Real = Default["home production"]; λ::Real =1)
    ω = weights(n_e, p_e_e, p_u_u, λ = λ)
    return(ω[1] * u(w + T) + ω[2] * u(w + T + τ) + ω[3] * u(δ + T) + ω[4] * u(δ + T + τ))
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

#Assuming the budget constraint is saturated, gives tu as a function of T, taking n_e exogenous
function τofT(T::Real = 0; n_e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"],  B::Real = Default["Budget"])
    ω = weights(n_e, p_e_e, p_u_u)
    return((B - T)/ (ω[2] + ω[4]))
end

function budgetcons(T::Real, τ::Real; n_e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], B::Real = Default["Budget"])
    ω = weights(n_e, p_u_u, p_u_u)
    if (T + (ω[2] + ω[4]) * τ) > B
        return(false)
    else
        return(true)
    end
end

function Ω(T; n_e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], B::Real = Default["Budget"], w::Real = Default["wage"], δ::Real = Default["home production"])
    τ = τofT(T, n_e = n_e, p_e_e = p_e_e, p_u_u = p_u_u, B= B)
    if (T < 0) || (τ < 0) || !budgetcons(T, τ; n_e = n_e, p_e_e = p_e_e, p_u_u = p_u_u, B = B)
        return(- Inf)
    else
        return(Welfare(T, τ, n_e, p_e_e, p_u_u, w, δ))
    end
end

function Solve(;n_e::Real = Default["n_e"], p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], B::Real = Default["Budget"], w::Real = Default["wage"], δ::Real = Default["home production"], λ = Default["Market tightness"])
    Param = Dict{String, Real}(
    "e" => n_e,
    "p_e|e" => p_e_e,
    "p_u|u" => p_u_u,
    "wage" => w,
    "home production" => δ,
    "Budget" => float(B),
    "Market tightness" => λ)

    function Obj(T::Real)
        return(- Ω(T, n_e = Param["e"], p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], B = Param["Budget"], w = Param["wage"], δ = Param["home production"]))
    end

    Results = optimize(Obj, 0, B, Brent())

    Tcor = 0
    W1 = Ω(0, n_e = Param["e"], p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], B = Param["Budget"], w = Param["wage"], δ = Param["home production"])
    W2 = Ω(B, n_e = Param["e"], p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], B = Param["Budget"], w = Param["wage"], δ = Param["home production"])

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


    Param["Optimal T"] = T
    Param["Optimal τ"] = τofT(T, n_e = Param["e"], p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], B = Param["Budget"])
    Param["Total Welfare"] = - Results.minimum
    return(Param)
end

#end



function Tofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["Optimal T"])
end

function tofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["Optimal τ"])
end

function eofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["e"])
end



function p_e_eplots(lb::Real=0.5, ub::Real=1)
    p1 = plot(Tofp_e_e, lb, ub, label = "Universal Transfer")
    plot!(p1, tofp_e_e, lb, ub, label = "Unemployed benefit")
    #p2 =  plot(eofp_e_e, lb, ub, label = "search effort")
    p = plot(p1)
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
    p = plot(p1)
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
    return(Solve( w= w)["e"])
end


function wageplot(lb::Real=0, ub::Real=5000)

    p1 = plot(Tofwage, lb, ub, label = "Universal Transfer")
    plot!(p1, tofwage, lb, ub, label = "Unemployment benefit")
    #p2 =  plot(eofwage, lb, ub, label = "search effort")
    p = plot(p1)
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
    #p2 =  plot(eofB, lb, ub, label = "Search Effort")
    #info("Search effort done")
    p = plot(p1)
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

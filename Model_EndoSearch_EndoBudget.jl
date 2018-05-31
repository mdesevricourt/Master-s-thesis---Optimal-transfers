if !("C:/Users/Maxime/.julia/v0.6/" in LOAD_PATH)
    push!(LOAD_PATH, "C:/Users/Maxime/.julia/v0.6/")
end

using Optim, Distributions, Plots; gr()
#using JLD2
#using LaTeXStrings
include("C:/Users/Maxime/Documents/MasterEco/Mémoire/Optimal-Transfers/Utility and cost functions.jl")
include("C:/Users/Maxime/Documents/MasterEco/Mémoire/Optimal-Transfers/SearchEffortmodule.jl")
include("C:/Users/Maxime/Documents/MasterEco/Mémoire/Optimal-Transfers/Welfare + Budget + tau functions.jl")
using utilityfunctions
using Searcheffort
using otherfunctions #welfare, budget, and τofT functions

Default = Dict("n_e" => 0.5, "p_e|e" => 0.9, "p_u|u" => 0.9, "wage" => 600.0, "home production" => 50, "Budget" => 300.0, "Market tightness" => 0.8)
Default["tax rate"] = 0.2

function Ω(T::Real; p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], gw::Real = Default["wage"],  ρ::Real = Default["tax rate"], δ::Real = Default["home production"], λ = Default["Market tightness"])
    τ, e = τofT(T, gw = gw, ρ = ρ,  p_e_e = p_e_e, p_u_u = p_u_u, δ = δ, λ = λ)
    if (T < 0) || (τ < 0) || !budgetcons(T, τ; gw = gw, ρ = ρ, n_e = e, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
        return(- Inf)
    else
        return(Welfare(T, τ; e = e, p_e_e = p_e_e, p_u_u = p_u_u, wb = gw,ρ = ρ, δ = δ, λ = λ))
    end
end



function Solve(; p_e_e::Real = Default["p_e|e"], p_u_u::Real = Default["p_u|u"], gw::Real = Default["wage"], δ::Real = Default["home production"], ρ::Real = Default["tax rate"], λ = Default["Market tightness"])
    Param = Dict{String, Real}(
    "p_e|e" => p_e_e,
    "p_u|u" => p_u_u,
    "wage" => gw,
    "home production" => δ,
    "Market tightness" => λ,
    "tax rate" => ρ)
    Param["net wage"] = netwage(gw, ρ)

    function Obj(T::Real)
        return(- Ω(T,  gw= gw, ρ = ρ, p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], δ = Param["home production"], λ = λ))
    end

    UBT = budget(ρ = ρ, gw = gw, e = 1, p_e_e = p_e_e) #upper bound on T, when everybody works (probablement infaisable car income effect)

    Results = optimize(Obj, 0, UBT, Brent())


    W1 = Ω(0, gw= gw, ρ = ρ, p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], δ = Param["home production"], λ = λ)
    W2 = Ω(UBT, gw= gw, ρ = ρ, p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], δ = Param["home production"], λ = λ)

    if W1 >= - Results.minimum
        T = 0
        W = W1
    elseif W2 >= - Results.minimum
        T = UBT
        W = W2
    else
        T = Results.minimizer
        W = - Results.minimum
    end

    τ, e, B = τofT(T, gw = Param["wage"], ρ = ρ, δ = Param["home production"], p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], λ = Param["Market tightness"])
    Param["Optimal T"] = T
    Param["Optimal τ"] = τ
    Param["e"] = e
    Param["Total Welfare"] = W
    Param["Budget"] = B
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
    info("Universal transfer done")
    plot!(p1, tofp_e_e, lb, ub, label = "Unemployed benefit")
    info("Unemployement benefit done")
    p2 =  plot(eofp_e_e, lb, ub, label = "search effort")
    info("search effort done")
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
    info("Universal transfer done")
    plot!(p1, tofp_u_u, lb, ub, label = "Unemployment benefit")
    info("Unemployment benefit done")
    p2 =  plot(eofp_u_u, lb, ub, label = "search effort")
    info("Search effort done")
    p = plot(p1, p2)
    return(p)
end

#Wage
function Tofwage(w::Real)
    return(Solve(gw = w)["Optimal T"])
end

function tofwage(w::Real)
    return(Solve(gw = w)["Optimal τ"])
end

function eofwage(w::Real)
    return(Solve(gw = w)["e"])
end


function wageplot(lb::Real=0, ub::Real=5000)

    p1 = plot(Tofwage, lb, ub, label = "Universal Transfer")
    plot!(p1, tofwage, lb, ub, label = "Unemployment benefit")
    p2 =  plot(eofwage, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end


#Budget

function Tofρ(T::Real)
    return(Solve(ρ = T)["Optimal T"])
end

function tofρ(T::Real)
    return(Solve(ρ = T)["Optimal τ"])
end

function eofρ(T::Real)
    return(Solve(ρ = T)["e"])
end



function taxrateplot(lb::Real=0, ub::Real=1)

    p1 = plot(Tofρ, lb, ub, label = "Universal Transfer", legend=:inside)
    info("Universal transfer done")
    plot!(p1, tofρ, lb, ub, label = "Unemployment benefit")
    info("Categorical benefit done")
    p2 =  plot(eofρ, lb, ub, label = "Search Effort")
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
    info("Universal transfer done")
    plot!(p1, tofdelta, lb, ub, label = "Unemployment benefit")
    info("Unemployment benefit done")
    p2 =  plot(eofdelta, lb, ub, label = "search effort")
    info("Searc effort done")
    p = plot(p1, p2)
    return(p)
end

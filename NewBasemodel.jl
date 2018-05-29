#New Baseline
using Optim
using Plots
#using PyPlot
using ApproxFun


#Parameters of the model
P = Dict("n_e" => 0.5, "n_u" => 0.5, "p_e|e" => 1.0, "p_u|u" => 1.0, "wage" => 600.0, "home production" => 50.0, "Budget" => 300.0, "Market tightness" => 1)

#utility function, which I choose to be log for the moment
function u(c)
    log(c)
    end

#Welfare function gives the total welfare in society, taking as arguments the transfers and the parameters of the model
function Welfare(τe::Real, τu::Real, n_e::Real, p_e_e::Real, p_u_u::Real, w::Real, δ::Real ; λ::Real =1)
    ω = weights(n_e, p_e_e, p_u_u)
    return(ω[1] * u(w + τe) + ω[2] * u(w + τu) + ω[3] * u(δ + τe) + ω[4] * u(δ + τu))
end

#Function that computes the weights of the different types of population given search effort, probability of fraud given employed, probability of non take-up given unemployed, and
#the market tightness. All parameters must be between 0 and 1.

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

#Function that returns true if positivity constraint on consumption are not satisfied
function cons(τe::Real, τu::Real, wage::Real =P["wage"], delta::Real =P["home production"])
    if ((τu <= - delta) | (τe <= - delta) |  (τe <= - wage) | (τu <= - wage))
        return(false)
    else
        return(true)
        end
    end

#Function that gives τ_e as a fonction of τu given the budget constraint is saturated
function teoftu(tu::Real, n_e::Real, p_e_e::Real, p_u_u::Real, T::Real)
    ω = weights(n_e, p_u_u, p_u_u)
    return((T - (ω[2] + ω[4] )* tu)/ (ω[1] + ω[3]))
end

#Objective Function with constraints (budget constraint is saturated, if positivty constraints are not satisfied, returns -Inf)
function Ω(τu; n_e::Real = P["n_e"], p_e_e::Real = P["p_e|e"], p_u_u::Real = P["p_u|u"], T::Real = P["Budget"], w::Real = P["wage"], δ::Real = P["home production"])
    τe = teoftu(τu, n_e, p_e_e, p_u_u, T)
    if cons(τe, τu, w, δ) == false #positivity constraints for consumption
        return(- Inf)
    else
        return(Welfare(τe, τu, n_e, p_e_e, p_u_u, w, δ))
        end
end

function Solve(Param::Dict = P)
    function Obj(τu::Real)e
        return(- Ω(τu, n_e = Param["n_e"], p_e_e = Param["p_e|e"], p_u_u = Param["p_u|u"], T = Param["Budget"], w = Param["wage"], δ = Param["home production"]))
    end
    #Results = optimize(univwrapper, -Param["home production"], (Param["Budget"] + Param["home production"])/Param["n_u"], Brent())
    ω = weights(Param["n_e"], Param["p_e|e"], Param["p_u|u"], λ = Param["Market tightness"])
    Results = optimize(Obj, -Param["home production"], (Param["Budget"] +  (ω[1] + ω[3]) * Param["home production"])/(ω[2] + ω[4] ), Brent())
    Dico = Dict("Optimal τu" => Results.minimizer, "Optimal τe" =>  teoftu(Results.minimizer, Param["n_e"], Param["p_e|e"], Param["p_u|u"], Param["Budget"]), "Total Welfare" => - Results.minimum)
end

function Basemodel(;n_e = P["n_e"]::Real, p_e_e = P["p_e|e"]::Real, p_u_u = P["p_u|u"]::Real, wage = P["wage"]::Real, homeprod = P["home production"]::Real, Budget = P["Budget"]::Real, λ::Real = P["Market tightness"])
    P = Dict("n_e" => n_e, "n_u" => 1 - n_e, "p_e|e" => p_e_e, "p_u|u" => p_u_u, "wage" => wage, "home production" => homeprod, "Budget" => Budget, "Market tightness" =>  λ)
    Solve(P)
    end

Basemodel()

function tuofn_e(n::Real)
    return(Basemodel(n_e = n)["Optimal τu"])
    end


function tuofT(T::Real)
    return(Basemodel(Budget = T)["Optimal τu"])
    end


function tuofp_e_e(p::Real)
    return(Basemodel(p_e_e = p)["Optimal τu"])
    end

function tuofp_u_u(p::Real)
    return(Basemodel(p_u_u = p)["Optimal τu"])
    end

function tuofwage(w::Real)
    return(Basemodel(wage =w)["Optimal τu"])
    end


function tuofhomeprod(δ::Real)
    return(Basemodel(homeprod = δ)["Optimal τu"])
    end


#Optimal τe
function teofn_e(n::Real)
    return(Basemodel(n_e = n)["Optimal τe"])
    end

function teofT(T::Real)
    return(Basemodel(Budget = T)["Optimal τe"])
    end


function teofp_e_e(p::Real)
    return(Basemodel(p_e_e = p)["Optimal τe"])
    end

function teofp_u_u(p::Real)
    return(Basemodel(p_u_u = p)["Optimal τe"])
    end

function teofwage(w::Real)
    return(Basemodel(wage =w)["Optimal τe"])
    end

function teofhomeprod(δ::Real)
    return(Basemodel(homeprod = δ)["Optimal τe"])
    end

function plottingBasemodel()
    p1 = plot(tuofn_e, 0, 1, xlabel = "%age of employed people", label = "tu")
    plot!(p1, teofn_e, 0 , 1, label = "te")
    p2 = plot(tuofp_e_e,0.5,1, xlabel = "Pourcentage of non-fraud among the employed", label = "tu")
    plot!(p2, teofp_e_e, 0.5 , 1, label = "te")
    p3 = plot(tuofp_u_u, 0.5,1, xlabel = "Pourcentage of take-up among unemployed", label = "tu")
    plot!(p3, teofp_u_u, 0.5 , 1, label = "tu")
    p4 = plot(tuofT,0 ,1000, xlabel = "Budget of the State", label = "tu")
    plot!(p4, teofT, 0 , 1000, label = "tu")
    p5 = plot(tuofwage,0 ,1000, xlabel = "Wage of the employed" , label = "tu")
    plot!(p5, teofwage, 0 , 1000, label = "tu")
    p6 = plot(tuofhomeprod, 0 ,1000, xlabel ="Home production", label = "tu")
    plot!(p6, teofhomeprod, 0 , 1000, label = "te")

    plot(p1, p5, p3, p4, p2 , p6, layout = (3,2))
    end

plottingBasemodel()
#So far everything seems right

#Now I'm going to try to find parameters so that τe = τu. One obvious solution seems to set wage = home production but it is not
#very realistic

function Difference(n_e::Real = P["n_e"], p_u_u::Real = P["p_u|u"], p_e_e::Real = P["p_e|e"] , Par::Dict = P)
    Param = copy(Par)
    Param["n_e"] = n_e
    Param["n_u"] = 1 - n_e
    Param["p_e|e"] = p_e_e
    Param["p_u|u"] = p_u_u
    Dico = Solve(Param)
    τe = Dico["Optimal τe"]
    τu = Dico["Optimal τu"]
    return(τu - τe)
    end

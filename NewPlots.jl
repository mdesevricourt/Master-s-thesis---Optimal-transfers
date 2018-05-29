module Plottransfers


using Plots
#=
export p_e_eplots, p_u_uplots, wageplot, budgetplot, budgetplot2

using Baseline
=#
push!(LOAD_PATH, "C:/Users/Maxime/.julia/v0.6/")
function Tofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["Optimal T"])
end

function tofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["Optimal τ"])
end

function eofp_e_e(p_e_e::Real)
    return(Solve(p_e_e = p_e_e)["e"])
end



function p_e_eplots(;lb::Real=0.5, ub::Real=1)
    p1 = plot(Tofp_e_e, lb, ub, label = "T")
    plot!(p1, tofp_e_e, lb, ub, label = "t")
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
    return(Solve(wage = w)["Optimal T"])
end

function tofwage(w::Real)
    return(Solve(wage = w)["Optimal τ"])
end

function eofwage(w::Real)
    return(Solve(wage = w)["e"])
end


function wageplot(lb::Real=0, ub::Real=5000)

    p1 = plot(Tofwage, lb, ub, label = "Universal Transfer")
    plot!(p1, tofwage, lb, ub, label = "Unemployment benefit")
    p2 =  plot(eofwage, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end


#Budget

function TofT(T::Real)
    return(Solve(B = T)["Optimal T"])
end

function tofT(T::Real)
    return(Solve(B = T)["Optimal τ"])
end

function eofT(T::Real)
    return(Solve(B = T)["e"])
end



function budgetplot(lb::Real=0, ub::Real=50)

    p1 = plot(TofT, lb, ub, label = "Universal Transfer")
    plot!(p1, tofT, lb, ub, label = "Unemployment benefit")
    p2 =  plot(eofT, lb, ub, label = "Search Effort")
    p = plot(p1, p2)
    return(p)
end


function budgetplot2(lb::Real = 0,  ub::Real = 100)

    τe = []
    τu = []
    e = []
    Budget = linspace(lb, ub, 1000)
    for i in Budget
        Dico = Solve(B = i)
        push!(τe, Dico["Optimal T"])
        push!(τu, Dico["Optimal τ"])
        push!(e, Dico["e"])
    end
    p1 = plot(Budget, τe, label ="Universal Transfer")
    plot!(p1, Budget, τu, label = "Unemployment Benefit")
    p2 = plot(Budget, e, label = "Search effort")
    plot(p1,p2)
end
#gr()
#plotlyjs()
function plot3D(f; x_lb = 0, x_ub = 1000, y_lb = x_lb, y_ub = x_ub, n = 1000)
x = repeat(linspace(0,1000,n), outer = n)
y = repeat(linspace(0,1000,n), inner = n)
z = [f(x[i], y[i]) for i in 1:length(x)]
plot(x,y,z)
end
# plot(x,y,z, xlabel = "Transfer to employed", ylabel = "Prob of working", zcolor=reverse(z), m=(10,0.8,:blues,stroke(0)),leg=false,cbar=true,w=5)
function wageplot(;lb::Real=0, ub::Real=5000)
    p1 = plot(Tofwage, lb, ub, label = "Universal Transfer")
    plot!(p1, tofwage, lb, ub)
    p2 =  plot(eofwage, lb, ub, label = "Unemployment Benefit")
    p = plot(p1, p2)
    return(p)
end

end

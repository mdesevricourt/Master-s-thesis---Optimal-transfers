#Plot functions

include("NewBasemodel.jl")
include("Modellingsearcheffort.jl")

#Endogenous search effort
PlotDefault = Dict("varint" => "Optimal τu", "param_key" => :wage)

#P_e_e
function tuofp_e_e(p_e_e::Real)
    return(EndoSearch(p_e_e = p_e_e)["Optimal τu"])
end

function teofp_e_e(p_e_e::Real)
    return(EndoSearch(p_e_e = p_e_e)["Optimal τe"])
end

function eofp_e_e(p_e_e::Real)
    return(EndoSearch(p_e_e = p_e_e)["e"])
end


function p_e_eplots(;lb::Real=0.5, ub::Real=1)

    p1 = plot(tuofp_e_e, lb, ub, label = "tu")
    plot!(p1, teofp_e_e, lb, ub, label = "te")
    p2 =  plot(eofp_e_e, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end

#p_u_u

function tuofp_u_u(p_u_u::Real)
    return(EndoSearch(p_u_u = p_u_u)["Optimal τu"])
end

function teofp_u_u(p_u_u::Real)
    return(EndoSearch(p_u_u = p_u_u)["Optimal τe"])
end

function eofp_u_u(p_u_u::Real)
    return(EndoSearch(p_u_u = p_u_u)["e"])
end


function p_u_uplots(lb::Real=0, ub::Real=1)

    p1 = plot(tuofp_u_u, lb, ub, label = "tu")
    plot!(p1, teofp_u_u, lb, ub, label = "te")
    p2 =  plot(eofp_u_u, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end

#Wage
function tuofwage(w::Real)
    return(EndoSearch(wage = w)["Optimal τu"])
end

function teofwage(w::Real)
    return(EndoSearch(wage = w)["Optimal τe"])
end

function eofwage(w::Real)
    return(EndoSearch(wage = w)["e"])
end


function wageplot(lb::Real=0, ub::Real=5000)

    p1 = plot(tuofwage, lb, ub, label = "tu")
    plot!(p1, teofwage, lb, ub, label = "te")
    p2 =  plot(eofwage, lb, ub, label = "search effort")
    p = plot(p1, p2)
    return(p)
end


#Budget

function tuofT(T::Real)
    return(EndoSearch(Budget = T)["Optimal τu"])
end

function teofT(T::Real)
    return(EndoSearch(Budget = T)["Optimal τe"])
end

function eofT(T::Real)
    return(EndoSearch(Budget = T)["e"])
end

function diffofT(T::Real)
    return(EndoSearch(Budget = T)["Corner solution gain"])
    end

function IntWelofT(T::Real)
    return(EndoSearch(Budget = T)["Interior Welfare"])
    end

function CorWelofT(T::Real)
    return(EndoSearch(Budget = T)["Corner Welfare"])
    end

function budgetplot(lb::Real=0, ub::Real=50)

    p1 = plot(tuofT, lb, ub, label = "tu")
    plot!(p1, teofT, lb, ub, label = "te")
    p2 =  plot(eofT, lb, ub, label = "search effort")
    p = plot(p1, p2, p3)
    return(p)
end


function budgetplot2(lb::Real = 0,  ub::Real = 100)

    τe = []
    τu = []
    e = []
    Budget = linspace(lb, ub, 1000)
    for i in Budget
        Dico = EndoSearch(Budget = i)
        push!(τe, Dico["Optimal τe"])
        push!(τu, Dico["Optimal τu"])
        push!(e, Dico["e"])
    end
    p1 = plot(Budget, τe, label ="te")
    plot!(p1, Budget, τu, label = "tu")
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
    p1 = plot(tuofwage, lb, ub, label = "tu")
    plot!(p1, teofwage, lb, ub)
    p2 =  plot(eofwage, lb, ub, label = "te")
    p = plot(p1, p2)
    return(p)
end

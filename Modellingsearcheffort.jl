
#Modelling search effort
P = Dict("n_e" => 0.5, "n_u" => 0.5, "p_e|e" => 0.9, "p_u|u" => 0.6, "wage" => 200.0, "home production" => 50.0, "Budget" => 1000.0, "Market tightness" => 1)
#New Welfare
function EndoWelfare(τe::Real, τu::Real, e::Real, w::Real=P["wage"], δ::Real= P["home production"],  p_e_e::Real = P["p_e|e"], p_u_u::Real = P["p_u|u"]; λ::Real = 1)
    ω = weights(e, p_e_e, p_u_u; λ =  λ)
    return(ω[1] * u(w + τe) + ω[2] * u(w + τu) + ω[3] * u(δ + τe) + ω[4] * u(δ + τu))
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

#utility from a given level of effort

function expectu(e::Real, τe::Real = 0, τu::Real = 0; w::Real = P["wage"], δ::Real = P["home production"], p_e_e::Real = P["p_e|e"], p_u_u::Real = P["p_u|u"], λ::Real = P["Market tightness"])
    ω = weights(e, p_e_e, p_u_u; λ = λ)
    return(ω[1] * u(w + τe) + ω[2] * u(w + τu) + ω[3] * u(δ + τe) + ω[4] * u(δ + τu) - disut(e))
end


function Optimeffort(τe::Real, τu::Real = 0; w::Real = P["wage"], δ::Real = P["home production"], p_e_e = P["p_e|e"], p_u_u = P["p_u|u"], λ::Real  = P["Market tightness"])
    #= I figured that an analytical solution was possible. So I implemented it. There are some weird cases where a numerical solution does stlightly better, probably
    due to computational approximiations, but the difference is so small it's s barely noticeable. I decided to get rid of this.
    =#
    e1 = invdisutprime(λ * (p_e_e * u(w + τe) + (1 - p_e_e) * u(w + τu) - (1 - p_u_u) * u(δ + τe) - p_u_u * u(δ + τu)))
    #enforcing the constraint on e1
    if e1 > 1
        e1 = 1
    elseif e1 < 0
        e1 = 0
    end

    #Checking for corner solutions

    if expectu(0, τe, τu, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = 1) > expectu(1, τe, τu,  w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = 1)
        if expectu(0, τe, τu,  w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = 1) > expectu(e1, τe, τu, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = 1)
            return(0)
        else
            return(e1)
            end
    elseif expectu(1, τe, τu, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = 1) > expectu(e1, τe, τu, w =  w, δ = δ, p_e_e =  p_e_e, p_u_u = p_u_u, λ = 1)
        return(1)
    else
        return(e1)
    end
end


function budgcons(τe::Real, τu::Real, e::Real, T::Real, p_e_e::Real, p_u_u::Real; λ::Real = 1)
    ω = weights(e, p_e_e, p_u_u; λ = λ)
    if  (ω[2] + ω[4] )* τu + (ω[1] + ω[3]) * τe > T
        return(false)
    else
        return(true)
end
end

function Ω(τe::Real= 0, τu::Real= 0; T::Real = P["Budget"], w::Real = P["wage"], δ::Real=P["home production"], p_e::Real = P["p_e|e"], p_u::Real = P["p_u|u"], λ::Real = P["Market tightness"])
    if cons(τe, τu, w, δ) == false #positivity constraints for consumption
        return(NaN)
    else
        e = Optimeffort(τe, τu, w = w, δ = δ, p_e_e = p_e, p_u_u = p_u, λ = λ)
        if budgcons(τe, τu, e, T, p_e, p_u, λ = λ) == false
            return(NaN)
        else
            return(EndoWelfare(τe, τu, e, w, δ, p_e, p_u, λ = λ))
        end
    end
end



function EndoSearch(;  p_e_e::Real = P["p_e|e"], p_u_u::Real = P["p_u|u"], wage::Real = P["wage"], homeprod::Real = P["home production"], Budget::Real = P["Budget"], λ::Real = P["Market tightness"])

    Dico = Dict{String, Real}("p_e|e" => p_e_e,
    "p_u|u" => p_u_u,
    "wage" => wage,
    "home production" => homeprod,
    "Budget" => float(Budget),
    "Market tightness" => λ)

    function Obj(τ::Vector)
        return(- Ω(τ[1], τ[2], T = Dico["Budget"], w = Dico["wage"], δ = Dico["home production"], p_e = Dico["p_e|e"], p_u = Dico["p_u|u"], λ = Dico["Market tightness"]))
    end

    function g!(storage::Vector, τ::Vector)
        τe = τ[1]
        τu = τ[2]
         if cons(τe, τu, Dico["wage"], Dico["home production"]) == false #positivity constraints for consumption
             storage[1] = NaN
            storage[2] = NaN
         else
             e = Optimeffort(τe, τu, w = Dico["wage"], δ = Dico["home production"], p_e_e = Dico["p_e|e"], p_u_u = Dico["p_u|u"], λ = Dico["Market tightness"])
             if budgcons(τe, τu, e, Dico["Budget"], Dico["p_e|e"], Dico["p_u|u"], λ = Dico["Market tightness"]) == false
                 storage[1] = NaN
                 storage[2] = NaN
             else
                storage[1] = - partOmte(τe, τu; w = Dico["wage"], δ = Dico["home production"], p_e_e = Dico["p_e|e"], p_u_u = Dico["p_u|u"], λ = Dico["Market tightness"])
                storage[2] = - partOmtu(τe, τu; w = Dico["wage"], δ = Dico["home production"], p_e_e = Dico["p_e|e"], p_u_u = Dico["p_u|u"], λ = Dico["Market tightness"])
             end
         end
    end


    Init_val = []
    a = rand(Uniform(0,1), 6)
    push!(Init_val, [Dico["Budget"],0.0], [0.0,Dico["Budget"]], [Dico["Budget"]-0.1,Dico["Budget"]]-0.1, [a[1]* Dico["Budget"],  a[2]* Dico["Budget"]], [a[3]* Dico["Budget"],
    a[4]* Dico["Budget"]], [a[5]* Dico["Budget"], a[6]* Dico["Budget"]] )


    Results = optimize(Obj, [0.0,0.0])
    Best_Welfare = - Results.minimum

    for val in Init_val
            Current_Results = optimize(Obj, val)
            if - Current_Results.minimum > Best_Welfare
                Results = Current_Results
                Best_Welfare = - Results.minimum
            end
    end

    for val in Init_val
             #info(val)
            try
                Current_Results = optimize(Obj, g!, val)
                if - Current_Results.minimum > Best_Welfare
                    Results = Current_Results
                    Best_Welfare = - Results.minimum
                    info("Gradient does better")

                end
            catch
                storage = [0.0,0.0]
                g!(storage, val)
                println(storage)
                println(val)
            end

    end

    Welfare = - Results.minimum
    InteriorWelfare = - Results.minimum
    #info(Results.initial_x)
    #Checking for convergence
    # if Results.g_converged == false & Results.f_converged == false & Results.x_converged == false
    #     warn("Convergence has not taken place")
    # end
    τ =  Results.minimizer

    #Checking for corner solutions in which all the budget is spent on one group only and the other transfer is minimal

    corner = false
    # τcorner = Cornersolutions(w = wage, δ = homeprod, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ, T = Budget)
    # CornerWelfare = Ω(τcorner[1], τcorner[2], T = Dico["Budget"], w = Dico["wage"], δ = Dico["home production"], p_e = Dico["p_e|e"], p_u = Dico["p_u|u"], λ = λ)
    # Differencecorner = CornerWelfare - Welfare
    #
    # if Differencecorner >= 0
    #     #warn("Corner solution")
    #     τ = τcorner
    #     Welfare = CornerWelfare
    #     corner = true
    # end


    Dico["Optimal τu"] = τ[2]
    Dico["Optimal τe"] =  τ[1]
    Dico["Total Welfare"] =  Welfare
    Dico["e"] = Optimeffort(τ[1], τ[2], w = Dico["wage"], δ = Dico["home production"], p_e_e = Dico["p_e|e"], p_u_u = Dico["p_u|u"], λ = λ)
    Dico["Convergence"] = Results.g_converged
    Dico["Corner solution"] = corner
    #Dico["Corner solution gain"] = Differencecorner
    Dico["Interior Welfare"] = InteriorWelfare
    #Dico["Corner Welfare"] = CornerWelfare
    return(Dico)
end

#Function that looks for corner solutions, that is, solution in which the whole budget is spent on one type of agents
function Cornersolutions(; w::Real = 200, δ::Real = 50, p_e_e::Real = 1, p_u_u::Real = 1, λ::Real = 1, T = 1000, maxiter = 1000)
    #whole budget spent on employed, assuming optimal effort goes to 1. Here, what could go wrong is that with λ < 1, some people are still going to be unemployed
    #I may want to tax unemployed to give money to employed

    ω = weights(1, p_e_e, p_u_u, λ = λ)
    τumin = - min(w, δ) + 0.01
    τemax = (T - (ω[2] + ω[4]) * τumin / (ω[1] + ω[3]))
    emax =  Optimeffort(τemax, τumin, w = w, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    convergence = false
    iteration = 0
    if emax !=1
        # info("Loop emax")
        while(! convergence & (iteration <= maxiter))
         eprec = copy(emax)
         iteration += 1
         ω = weights(emax, p_e_e, p_u_u, λ = λ)
         τemax = (T - (ω[2] + ω[4]) * τumin / (ω[1] + ω[3]))
         emax =  Optimeffort(τemax, τumin, w = w, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
            if isapprox(emax, eprec, atol =  1e-2)
                convergence = true
            end
            if iteration == 10
        #        warn("10 iterations")
            end
            if ((iteration == maxiter) & (convergence == false))
        #        warn("No convergence for emax")
            end
        end

     end

    #whole budget spent on unemployed, assuming optimal effort goes to
    ω = weights(0, p_e_e, p_u_u, λ = λ)
    τemin = - min(w, δ) + 0.01
    τumax = (T - (ω[1] + ω[3]) * τemin )/ (ω[2] + ω[4])
    emin =  Optimeffort(τemin, τumax, w = w, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    convergence = false
    iteration = 0
    n = 100 #number of possible arbitrary values
    eposs = sort(linspace(0,1,n), rev =true) #possible values of e we are going to look through if we are really desperate about finding a solution

    if emin != 0
        #info("Loop emin")
        i = 1 #iteration to look for fixed point

        while(! convergence & (iteration <= maxiter)) #as long as there is no convergence and the number of iteration is not reached
            #increase number of iteration
            #println(τumax)
            #println("emin = $emin")
            eprec = copy(emin) #emin of the precedent iteration is kept in memory for comparison
            ω = weights(emin, p_e_e, p_u_u, λ = λ) #calculate population weights at current emin
            τumax = (T - (ω[1] + ω[3]) * τemin / (ω[2] + ω[4])) #calulcate max affordable τu at current emin

            if (τumax == Inf) & (i < n) #if τu affordable is so high the computer cannot deal with it, then aux grands maux, les grands remèdes
                foundnexte = false

                while foundnexte == false #while we keep finding impossible τu

                    # warn("τumax is Infinity")
                    emin = eposs[i] #loop through values of e in descending order
                    i += 1 #keep track of where we are in terms of possible values for e
                    ω = weights(emin, p_e_e, p_u_u, λ = λ) #new population weights for that e
                    #info(emin)
                    if (T - (ω[1] + ω[3]) * τemin / (ω[2] + ω[4])) < Inf #new τumax
                        foundnexte = true
                        end
                    if i > n
                        #warn("we have looped through all possible values")
                        break
                        end
                    end

                continue
                end

            emin =  Optimeffort(τemin, τumax, w = w,  p_e_e = p_e_e, p_u_u = p_u_u, λ = λ) #new optimal effort at current maximal transfer


            if isapprox(emin, eprec, atol=  1e-3)

                convergence = true
                #println("τumax = $τumax", "emin = $emin", "epred = $eprec")
                end

            if iteration == 10
                #warn("10 iterations")
                end

            if ((iteration == maxiter) & (convergence == false))
                #warn("No convergence for emin")
                ω = weights(0, p_e_e, p_u_u, λ = λ)
                τumax = (T - (ω[1] + ω[3]) * τemin )/ (ω[2] + ω[4])
                #println(ω)
                #println(τumax)
                end

            iteration += 1
            end
        end

#    println("τemax = $τemax")
#    println("τumax = $τumax")

    if Ω(τemax, τumin, T = T, w = w, δ = δ, p_e = p_e_e, p_u = p_u_u, λ = λ) > Ω(τemin, τumax, T = T, w = w, δ = δ, p_e = p_e_e, p_u = p_u_u, λ = λ)
        τ = [τemax, τumin]
    elseif Ω(τemax, τumin, T = T, w = w, δ = δ, p_e = p_e_e, p_u = p_u_u, λ = λ) < Ω(τemin, τumax, T = T, w = w, δ = δ, p_e = p_e_e, p_u = p_u_u, λ = λ)
        τ = [τemin, τumax]
    else
        warn("No better corner solution found")
        τ = [0,0]
    end
    return(τ) #returns the best corner solution
    end

#=
Problem = []
for T in linspace(0,2000, 1000)
    Dico = EndoSearch(Budget = T)
    if Dico["Convergence"] == false
        push!(Problem, T)
    end
end
=#


#Gradient Method

function uprime(x::Real)
    return(1/x)
end

function eprimete(τe::Real = 0, τu::Real = 0 ; w::Real = P["wage"], δ::Real = P["home production"], p_e_e = P["p_e|e"], p_u_u = P["p_u|u"], λ::Real  = P["Market tightness"])
    e1 = Optimeffort(τe, τu, w = w, δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    if (e1 == 1) | (e1 == 0)
        return(0)
    else
        return(0.5* λ*(p_e_e * uprime(w + τe) - (1 - p_u_u) * uprime(δ + τe) ))
        end
end
function eprimetu(τe::Real = 0, τu::Real = 0 ; w::Real = P["wage"], δ::Real = P["home production"], p_e_e = P["p_e|e"], p_u_u = P["p_u|u"], λ::Real  = P["Market tightness"])
    e1 = Optimeffort(τe, τu, w = w, δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    if (e1 == 1) | (e1 == 0)
        return(0)
    else
        return(0.5* λ*((1-p_e_e) * uprime(w + τu) -  p_u_u * uprime(δ + τu) ))
    end
end

function partOmte(τe::Real = 0, τu::Real = 0 ; w::Real = P["wage"], δ::Real = P["home production"], p_e_e = P["p_e|e"], p_u_u = P["p_u|u"], λ::Real  = P["Market tightness"])
    part = eprimete(τe, τu, w = w, δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    e = Optimeffort(τe, τu, w = w, δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    return( λ* part * (p_e_e * u(w + τe) + (1 - p_e_e) * u(w + τu)) + λ * e * p_e_e * uprime(w+ τe) - λ * part * (p_u_u * u(δ + τu) + (1 - p_u_u) *  u(δ + τe)  ) + (1 - λ*e) * (1-p_u_u) * uprime(δ + τe) )
end

function partOmtu(τe::Real = 0, τu::Real = 0 ; w::Real = P["wage"], δ::Real = P["home production"], p_e_e = P["p_e|e"], p_u_u = P["p_u|u"], λ::Real  = P["Market tightness"])
    part = eprimetu(τe, τu ; w = w, δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    e = Optimeffort(τe, τu, w = w, δ = δ, p_e_e = p_e_e, p_u_u = p_u_u, λ = λ)
    return(λ * part * (p_e_e * u(w + τe) + (1 - p_e_e) * u(w + τu)) + λ * e * (1 - p_e_e) * uprime(w + τu) - λ * part * (p_u_u * u(δ + τu) + (1 - p_u_u) *  u(δ + τe)) + (1 - λ * e) * p_u_u * uprime(δ + τu) )
end


#= used to find corner solutions in a too simple way
if Ω(0, Dico["Budget"], T = Dico["Budget"], w = Dico["wage"], δ = Dico["home production"], p_e = Dico["p_e|e"], p_u = Dico["p_u|u"], λ = λ) > Welfare
# warn("Corner solution")
τ = [0, Dico["Budget"]]
# println(Optimeffort(τ[1], τ[2], w = Dico["wage"], δ = Dico["home production"], p_e_e = Dico["p_e|e"], p_u_u = Dico["p_u|u"]))
Welfare = Ω(τ[1], τ[2], T = Dico["Budget"], w = Dico["wage"], δ = Dico["home production"], p_e = Dico["p_e|e"], p_u = Dico["p_u|u"], λ = λ)
corner = true
end

if Ω(Dico["Budget"], 0, T = Dico["Budget"], w = Dico["wage"], δ = Dico["home production"], p_e = Dico["p_e|e"], p_u = Dico["p_u|u"], λ = λ) > Welfare
# warn("Corner solution")
τ = [Dico["Budget"], 0]
# println(Optimeffort(τ[1], τ[2], w = Dico["wage"], δ = Dico["home production"],  p_e_e = Dico["p_e|e"], p_u_u = Dico["p_u|u"]))
Welfare = Ω(τ[1], τ[2], T = Dico["Budget"], w = Dico["wage"], δ = Dico["home production"], p_e = Dico["p_e|e"], p_u = Dico["p_u|u"], λ = λ)
corner = true
end
=#



function EndoTest(;  p_e_e::Real = P["p_e|e"], p_u_u::Real = P["p_u|u"], wage::Real = P["wage"], homeprod::Real = P["home production"], Budget::Real = P["Budget"], λ::Real = P["Market tightness"])

    Dico = Dict{String, Real}("p_e|e" => p_e_e,
    "p_u|u" => p_u_u,
    "wage" => wage,
    "home production" => homeprod,
    "Budget" => Budget,
    "Market tightness" => λ)

    function Obj(τ::Vector)
        return(- Ω(τ[1], τ[2], T = Dico["Budget"], w = Dico["wage"], δ = Dico["home production"], p_e = Dico["p_e|e"], p_u = Dico["p_u|u"], λ = Dico["Market tightness"]))
    end



    Init_val = []
    push!(Init_val, [0.0,0.0], [Dico["Budget"],0.0], [0.0,Dico["Budget"]])

    push!(Results,
        optimize(Obj, [0.0,0.0]).minimum,
        optimize(Obj, [Dico["Budget"],0.0]).minimum,
        optimize(Obj, [0.0,Dico["Budget"]]).minimum,
        optimize(Obj, [Dico["Budget"]/2,Dico["Budget"]/2]).minimum,
        optimize(Obj, g!, [0.0,0.0], GradientDescent()).minimum,
        optimize(Obj, g!, [Dico["Budget"],0.0], GradientDescent()).minimum,
        optimize(Obj, g!, [0.0,Dico["Budget"]], GradientDescent()).minimum,
        optimize(Obj, g!, [Dico["Budget"]/2,Dico["Budget"]/2], GradientDescent()).minimum )
    return(Results)
end

function check_nonnegativity(f::Expression,OneSol::Bool=false)
    vars=variables(f)
    dimension=length(vars)
    HomotopyContinuation.@var t 
    vars_t=Vector{Variable}()
    push!(vars_t,t)
    for v in vars
        push!(vars_t,v)
    end
    f_t=0
    exponents, coeffs=exponents_coefficients(f,vars)
    num_pos_coeffs=0
    num_neg_coeffs=0
    neg_coeffs=[]
    pos_coeffs=[]
    pos_exponents=[]
    neg_exponents=[]
    for (j,c) in enumerate(coeffs)
        if c<0
            num_neg_coeffs=num_neg_coeffs+1
            push!(neg_exponents,exponents[:,j])
            push!(neg_coeffs,c)
        else
            num_pos_coeffs=num_pos_coeffs+1
            push!(pos_exponents,exponents[:,j])
            push!(pos_coeffs,c)
        end
    end
    if OneSol==false
        #Code for the general case
        for j in 1:num_pos_coeffs
            f_t=f_t+pos_coeffs[j]*prod(vars[i]^pos_exponents[j][i] for i in 1:dimension)
        end
        for j in 1:num_neg_coeffs
            f_t=f_t+t*neg_coeffs[j]*prod(vars[i]^neg_exponents[j][i] for i in 1:dimension)
        end
        
        eqs=[f_t]
        for v in vars
            push!(eqs,v*differentiate(f_t,v))
        end
        F=System(eqs;variables=vars_t)
        solution=HomotopyContinuation.solve(F)
        cert=certificates(certify(F,solution))
        t_coords=[]
        for c in cert
            if HomotopyContinuation.is_positive(c)
                push!(t_coords,solution_candidate(c)[1])
            end
        end
        t_min=minimum(t_coords)
        if t_min<1
            return false, t_min
        else 
            return true, t_min
        end
        
    else
        #Code for when the support admits a unique positive solution

        #Construcion of starting system
        
        #Matrix to compute bary.coord. of first negative exponent vector 
        A=hcat(vcat(transpose(ones(Int,length(pos_exponents))),hcat(pos_exponents...)),-vcat(1,neg_exponents[1]))
        
        #P is the polyhedron: ker(A)\cap\R^n_{\geq 0}
        P=Oscar.polyhedron((-Matrix{Int}(I, num_pos_coeffs+1, num_pos_coeffs+1),zeros(Int,num_pos_coeffs+1)),(A,zeros(Int,dimension+1)))
        
        interior=Oscar.relative_interior_point(P)
        vec_numer=Vector{BigInt}()
        vec_denom=Vector{BigInt}()
        for v in interior
            push!(vec_numer,BigInt(numerator(v)))
            push!(vec_denom,BigInt(denominator(v)))
        end

        interior_float=vec_numer./vec_denom #barycentric coordinates of the first negative exponent vector
        
        initial_coeffs=float.(zeros(length(coeffs)))

        for j in 1:num_pos_coeffs
            initial_coeffs[j]=interior_float[j]
        end
        initial_coeffs[num_pos_coeffs+1]=-interior_float[end]

        #Setting the homotopy
        HomotopyContinuation.@var p[1:(length(coeffs))]
        
        for j in 1:num_pos_coeffs
            f_t=f_t+p[j]*prod(vars[i]^pos_exponents[j][i] for i in 1:dimension)
        end
        for j in 1:num_neg_coeffs
            f_t=f_t+p[j+num_pos_coeffs]*t*prod(vars[i]^neg_exponents[j][i] for i in 1:dimension)
        end
        
        eqs=[f_t]
        for v in vars
            push!(eqs,v*differentiate(f_t,v))
        end
        F=System(eqs;variables=vars_t,parameters=p)
        
        start_solutions=float.(ones(dimension+1))
        reordered_coeffs=vcat(pos_coeffs,neg_coeffs)
        
        solution=HomotopyContinuation.solve(F,[start_solutions];start_parameters=Float64.(initial_coeffs),target_parameters=reordered_coeffs,seed=0x68a5c2c6)

        #To certify, first we need to substitute the parameters in the system
        eqs_sub=subs(eqs,Dict(p.=>reordered_coeffs))
        F_sub=HomotopyContinuation.System(eqs_sub)
        
        cert=certificates(certify(F_sub,solutions(solution)))
        candidate=[]
        for c in cert
            if HomotopyContinuation.is_positive(c)==true
                candidate=solution_candidate(c)
            end
        end
        t_coord=real(candidate[1])
        if t_coord<1
            return false, t_coord
        else 
            return true, t_coord
        end
    end
end
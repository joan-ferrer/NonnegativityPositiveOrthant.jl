function check_nonnegativity(f::Expression,OneSol::Bool=false)
    vars=variables(f)
    dimension=length(vars)
    HomotopyContinuation.@var t 
    vars_t=variables(f)
    push!(vars_t,t)
    f_t=0
    exponents, coeffs=exponents_coefficients(f,vars)
    if OneSol==false
        #Code for the general case
        for (j,c) in enumerate(coeffs)
            if c<0
                f_t=f_t+c*t*prod(vars[i]^exponents[i,j] for i in 1:dimension) #should we allow different $h$?
            else
                f_t=f_t+c*prod(vars[i]^exponents[i,j] for i in 1:dimension)
            end
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
                push!(t_coords,real(solution_candidate(c)[end]))
            end
        end
        t_min=min(t_coords)
        if t_min<1
            return false, t_min
        else 
            return true, t_min
        end
    else
        #Code for when the support admits a unique positive solution

        #Construcion of starting system
        k=0
        num_pos_coeffs=0
        num_neg_coeffs=0
        pos_exponents=[]
        for (j,c) in enumerate(coeffs)
            if c<0
                if k==0
                    k=j
                end
                num_neg_coeffs=num_neg_coeffs+1
            else
                num_pos_coeffs=num_pos_coeffs+1
                push!(pos_exponents,exponents[:,j])
            end
        end
        
        #Matrix to compute bary.coord. of first negative exponent vector 
        A=hcat(vcat(transpose(ones(Int,length(pos_exponents))),hcat(pos_exponents...)),-vcat(1,exponents[:,k]))
        
        #P is the polyhedron: ker(A)\cap\R^n_{\geq 0}
        P=Oscar.polyhedron((-Matrix{Int}(I, num_pos_coeffs+1, num_pos_coeffs+1),zeros(Int,num_pos_coeffs+1)),(A,zeros(Int,dimension+1)))
        
        interior=Oscar.relative_interior_point(P)
        vec_numer=Vector{BigInt}()
        vec_denom=Vector{BigInt}()
        for v in interior
            push!(vec_numer,BigInt(numerator(v)))
            push!(vec_denom,BigInt(denominator(v)))
        end
        
        interior_float=vec_numer./vec_denom
        initial_coeffs=float.(zeros(length(coeffs)))
        for (j,c) in enumerate(coeffs)
            if exponents[:,j] in pos_exponents
                ind=findfirst(x->x==exponents[:,j],pos_exponents)
                initial_coeffs[j]=interior_float[ind]
            end
        end
        initial_coeffs[k]=-interior_float[end]

        #Setting the homotopy
        HomotopyContinuation.@var p[1:(length(coeffs))]

        for (j,c) in enumerate(coeffs)
            if c<0
                f_t=f_t+p[j]*t*prod(vars[i]^exponents[i,j] for i in 1:dimension) #should we allow different $h$?
            else
                f_t=f_t+p[j]*prod(vars[i]^exponents[i,j] for i in 1:dimension)
            end
        end
        
        eqs=[f_t]
        for v in vars
            push!(eqs,v*differentiate(f_t,v))
        end
        F=System(eqs;variables=vars_t,parameters=p)
        start_solutions=float.(ones(dimension+1))
        solution=HomotopyContinuation.solve(F,[start_solutions];start_parameters=initial_coeffs,target_parameters=coeffs)

        #To certify, first we need to substitute the parameters in the system
        eqs_sub=subs(eqs,Dict(p.=>coeffs))
        F_sub=HomotopyContinuation.System(eqs_sub)

        cert=certificates(certify(F_sub,solution))
        t_coords=[]
        for c in cert
            if HomotopyContinuation.is_positive(c)
                push!(t_coords,real(solution_candidate(c)[end]))
            end
        end
        t_min=minimum(t_coords)
        if t_min<1
            return false, t_min
        else 
            return true, t_min
        end
    end
end


using NonnegativityPositiveOrthant
using Test
using HomotopyContinuation
using Oscar

#=
@testset "NonnegativityPositiveOrthant.jl" begin
    HomotopyContinuation.@var x[1:4]
    f=x[1]^40+x[2]^40+x[3]^40+x[4]^40-float((10/9)^(9/10)*40^(1/10))*x[1]*x[2]*x[3]*x[4]+1
    @test check_nonnegativity(f,true)  == (true, 1.0)
    
end

#Ctrl+] a REPL (Ctr+Shift+P) per testejar coses del package


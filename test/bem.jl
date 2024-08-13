using Test
using Propeller 
using YAML
using FlightConditions

config= YAML.load_file("../config/rotor/test.yaml")
rotor = Rotor(config)
fc = FlightCondition(0.0)

@testset "Hover" begin 
    res = rpm2thrust(rotor, 3000.0, [0, 0, 0]; flight_condition=fc, lb_T=0, ub_T=100)
    @test isapprox(res[1], 23.8572669801628; rtol=1e-6)
    @test isapprox(res[2], 0.7047250606549583; rtol=1e-6)
end 

@testset "Convergence Vortex Ring" begin
    res = rpm2thrust(rotor, 3000.0, [0, 0, -20]; flight_condition=fc, lb_T=0, ub_T=100)
    @test res[3] == false
end

@testset "Forward flight" begin 
    res = rpm2thrust(rotor, 3000.0, [10, 5, 3]; flight_condition=fc, lb_T=0, ub_T=100)
    @test isapprox(res[1], 22.751350257076357; rtol=1e-6)
    @test isapprox(res[2], 0.5283360139316857; rtol=1e-6)
end
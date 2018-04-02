@testset "Set parameters" begin
    lpm = MathProgBase.LinearQuadraticModel(GLPKSolverLP(it_lim=9513, tol_bnd=4.149))
    @test lpm.param.it_lim == 9513
    @test lpm.param.tol_bnd == 4.149
    mipm = MathProgBase.LinearQuadraticModel(GLPKSolverMIP(it_lim=5910, tol_obj=1.52e-3))
    @test mipm.smplxparam.it_lim == 5910
    @test mipm.param.tol_obj == 1.52e-3
end

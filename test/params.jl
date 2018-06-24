@testset "Set parameters" begin
    lpm = MathProgBase.LinearQuadraticModel(GLPKSolverLP(it_lim=9513, tol_bnd=4.149))
    @test lpm.param.it_lim == 9513
    @test lpm.param.tol_bnd == 4.149
    mipm = MathProgBase.LinearQuadraticModel(GLPKSolverMIP(it_lim=5910, tol_obj=1.52e-3))
    @test mipm.smplxparam.it_lim == 5910
    @test mipm.param.tol_obj == 1.52e-3
end

@testset "Set MPB parameters" begin
    # setparameters! on model
    lpm = MathProgBase.LinearQuadraticModel(GLPKSolverLP())
    MathProgBase.setparameters!(lpm, TimeLimit=23.0)
    @test lpm.param.tm_lim == 23000.0

    # setparameters! on solver
    lps = GLPKSolverLP()
    MathProgBase.setparameters!(lps, TimeLimit=23.0)
    lpm2 = MathProgBase.LinearQuadraticModel(lps)
    @test lpm2.param.tm_lim == 23000.0
end

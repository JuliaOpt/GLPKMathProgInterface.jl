using GLPKMathProgInterface, MathProgBase, Base.Test

@testset "Invalid bounds" begin
    solver = GLPKSolverLP()
    m = MathProgBase.LinearQuadraticModel(solver)
    # invalid variable bounds
    @test_throws GLPK.GLPKError MathProgBase.loadproblem!(m, [1 0], [0, 0], [1, -1], [0, 0], [0], [1], :Min)
    # invalid constraint bounds
    @test_throws GLPK.GLPKError MathProgBase.loadproblem!(m, [1 0], [0, 0], [1, 1], [0, 0], [0], [-1], :Min)
    MathProgBase.loadproblem!(m, [1 0], [0, 0], [1, 1], [0, 0], [0], [1], :Min)
    # Those do not throws errors since the invalid bounds may be temporary
    MathProgBase.setvarLB!(m, [2, 0])
    MathProgBase.setvarUB!(m, [1, -1])
    MathProgBase.setconstrLB!(m, [2])
    MathProgBase.setconstrUB!(m, [-1])
end

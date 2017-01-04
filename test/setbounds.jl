using GLPKMathProgInterface, MathProgBase, Base.Test

@testset "Invalid bounds" begin
    solver = GLPKSolverLP()
    m = MathProgBase.LinearQuadraticModel(solver)
    function testmodel(sense)
        MathProgBase.optimize!(m)
        @test MathProgBase.status(m) == :Infeasible
        @test MathProgBase.getobjval(m) == (sense == :Min ? Inf : -Inf)
        @test MathProgBase.getinfeasibilityray(m) == zeros(1)
        @test_throws ErrorException MathProgBase.getsolution(m)
        @test_throws ErrorException MathProgBase.getconstrsolution(m)
        @test_throws ErrorException MathProgBase.getreducedcosts(m)
        @test_throws ErrorException MathProgBase.getconstrduals(m)
        @test_throws ErrorException MathProgBase.getunboundedray(m)
    end
    # invalid variable bounds
    MathProgBase.loadproblem!(m, [1 0], [0, 0], [1, -1], [0, 0], [0], [1], :Min)
    testmodel(:Min)
    # invalid constraint bounds
    MathProgBase.loadproblem!(m, [1 0], [0, 0], [1, 1], [0, 0], [0], [-1], :Max)
    testmodel(:Max)
    MathProgBase.loadproblem!(m, [1 0], [0, 0], [1, 1], [0, 0], [0], [1], :Min)
    # Those do not throws errors since the invalid bounds may be temporary
    MathProgBase.setvarLB!(m, [2, 0])
    MathProgBase.setvarUB!(m, [1, -1])
    MathProgBase.setconstrLB!(m, [2])
    MathProgBase.setconstrUB!(m, [-1])
    testmodel(:Min)
end

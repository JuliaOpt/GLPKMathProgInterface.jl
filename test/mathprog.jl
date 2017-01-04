@testset "MathProgBase tests" begin
    const mathprogbase_test = joinpath(Pkg.dir("MathProgBase"), "test")

    include(joinpath(mathprogbase_test, "linprog.jl"))
    linprogtest(GLPKSolverLP())

    include(joinpath(mathprogbase_test, "linproginterface.jl"))
    linprogsolvertest(GLPKSolverLP())

    include(joinpath(mathprogbase_test, "mixintprog.jl"))
    mixintprogtest(GLPKSolverMIP())

    include(joinpath(mathprogbase_test,"conicinterface.jl"))
    coniclineartest(GLPKSolverLP())
    coniclineartest(GLPKSolverMIP())
end

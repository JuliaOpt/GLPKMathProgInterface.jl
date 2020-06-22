const mathprogbase_test = joinpath(dirname(pathof(MathProgBase)), "..", "test")

@testset "MathProgBase tests" begin

    include(joinpath(mathprogbase_test, "linprog.jl"))
    linprogtest(GLPKSolverLP())

    include(joinpath(mathprogbase_test, "linproginterface.jl"))
    linprogsolvertest(GLPKSolverLP())

    include(joinpath(mathprogbase_test, "mixintprog.jl"))
    mixintprogtest(GLPKSolverMIP())

    include(joinpath(mathprogbase_test,"conicinterface.jl"))
    coniclineartest(GLPKSolverLP())

    solver = GLPKSolverMIP(msg_lev=GLPK.MSG_ALL)

    # Silent should override GLPK.MSG_ALL
    MathProgBase.setparameters!(solver, Silent=true, TimeLimit=100.0)
    coniclineartest(solver)

end

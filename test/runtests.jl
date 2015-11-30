using GLPKMathProgInterface

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

lpm = GLPKInterfaceLP.LinearQuadraticModel(GLPKSolverLP(it_lim=9513, tol_bnd=4.149))
@test lpm.param.it_lim == 9513
@test lpm.param.tol_bnd == 4.149
mipm = GLPKInterfaceMIP.LinearQuadraticModel(GLPKSolverMIP(it_lim=5910, tol_obj=1.52e-3))
@test mipm.smplxparam.it_lim == 5910
@test mipm.param.tol_obj == 1.52e-3

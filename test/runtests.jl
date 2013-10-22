using GLPKMathProgInterface

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(GLPKSolverLP())

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(GLPKSolverLP())

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
mixintprogtest(GLPKSolverMIP())

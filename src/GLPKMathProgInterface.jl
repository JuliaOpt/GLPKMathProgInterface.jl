__precompile__()
module GLPKMathProgInterface

export
    GLPKInterfaceLP,
    GLPKInterfaceMIP,
    GLPKSolverLP,
    GLPKSolverMIP

include("GLPKInterfaceBase.jl")
include("GLPKInterfaceLP.jl")
include("GLPKInterfaceMIP.jl")

const GLPKSolverLP = GLPKInterfaceLP.GLPKSolverLP
const GLPKSolverMIP = GLPKInterfaceMIP.GLPKSolverMIP

# Enables GLPK to act as a conic solver
import MathProgBase
MathProgBase.ConicModel(s::Union{GLPKSolverLP,GLPKSolverMIP}) = MathProgBase.LPQPtoConicBridge(MathProgBase.LinearQuadraticModel(s))
MathProgBase.supportedcones(::Union{GLPKSolverLP,GLPKSolverMIP}) = [:Free,:Zero,:NonNeg,:NonPos]

end

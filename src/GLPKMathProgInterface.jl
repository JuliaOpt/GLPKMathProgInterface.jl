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

end

module GLPKMathProgInterface

export
    GLPKInterfaceLP,
    GLPKInterfaceMIP,
    GLPKSolverLP,
    GLPKSolverMIP

include(joinpath(dirname(Base.source_path()), "GLPKInterfaceBase.jl"))
include(joinpath(dirname(Base.source_path()), "GLPKInterfaceLP.jl"))
include(joinpath(dirname(Base.source_path()), "GLPKInterfaceMIP.jl"))

const GLPKSolverLP = GLPKInterfaceLP.GLPKSolverLP
const GLPKSolverMIP = GLPKInterfaceMIP.GLPKSolverMIP

end

module GLPKMathProgInterface

export
    GLPKInterfaceLP,
    GLPKInterfaceMIP

include(joinpath(dirname(Base.source_path()), "GLPKInterfaceBase.jl"))
include(joinpath(dirname(Base.source_path()), "GLPKInterfaceLP.jl"))
include(joinpath(dirname(Base.source_path()), "GLPKInterfaceMIP.jl"))

end

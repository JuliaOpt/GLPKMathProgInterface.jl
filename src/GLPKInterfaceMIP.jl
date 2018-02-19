module GLPKInterfaceMIP

import GLPK
importall MathProgBase.SolverInterface
importall ..GLPKInterfaceBase

export
    GLPKSolverMIP,
    GLPKCallbackData,
    optimize!,
    loadproblem!,
    writeproblem,
    getvarLB,
    setvarLB!,
    getvarUB,
    setvarUB!,
    getconstrLB,
    setconstrLB!,
    getconstrUB,
    setconstrUB!,
    getobj,
    setobj!,
    addvar!,
    addconstr!,
    setsense!,
    getsense,
    numvar,
    numconstr,
    setvartype!,
    status,
    getobjval,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver,
    setlazycallback!,
    setcutcallback!,
    setheuristiccallback!,
    cbaddlazy!,
    cbgetlpsolution,
    cbgetmipsolution,
    cbgetstate,
    cbgetobj,
    cbgetbestbound,
    cbgetexplorednodes,
    cbaddcut!,
    cbaddsolution

mutable struct GLPKMathProgModelMIP <: GLPKMathProgModel
    inner::GLPK.Prob
    param::GLPK.IntoptParam
    smplxparam::GLPK.SimplexParam
    lazycb::Union{Function,Void}
    cutcb::Union{Function,Void}
    heuristiccb::Union{Function,Void}
    infocb::Union{Function,Void}
    objbound::Float64
    cbdata::MathProgCallbackData
    binaries::Vector{Int}
    userlimit::Bool
    function GLPKMathProgModelMIP()
        lpm = new(GLPK.Prob(), GLPK.IntoptParam(), GLPK.SimplexParam(),
                  nothing, nothing, nothing, nothing, -Inf)
        lpm.cbdata = GLPKCallbackData(lpm)
        lpm.binaries = Int[]
        lpm.userlimit = false
        return lpm
    end
end

function Base.copy(m::GLPKMathProgModelMIP)

    m2 = GLPKMathProgModelMIP()

    GLPK.copy_prob(m2.inner, m.inner, GLPK.ON)

    m2.param = deepcopy(m.param)
    m2.smplxparam = deepcopy(m.smplxparam)

    m.lazycb == nothing || Base.warn_once("Callbacks can't be copied, lazy callback ignored")
    m.cutcb == nothing || Base.warn_once("Callbacks can't be copied, cut callback ignored")
    m.heuristiccb == nothing || Base.warn_once("Callbacks can't be copied, heuristic callback ignored")
    m.infocb == nothing || Base.warn_once("Callbacks can't be copied, info callback ignored")

    m2.objbound = m.objbound

    m.cbdata == nothing || Base.warn_once("Callbacks can't be copied, callbackdata ignored")

    m2.binaries = deepcopy(m.binaries)
    m2.userlimit = m.userlimit

    return m2
end

mutable struct GLPKCallbackData <: MathProgCallbackData
    model::GLPKMathProgModelMIP
    tree::Ptr{Void}
    state::Symbol
    reason::Cint
    sol::Vector{Float64}
    vartype::Vector{Symbol}
    GLPKCallbackData(model::GLPKMathProgModelMIP) = new(model, C_NULL, :Other, -1, Float64[], Char[])
end

mutable struct GLPKSolverMIP <: AbstractMathProgSolver
    presolve::Bool
    opts
    GLPKSolverMIP(;presolve::Bool=false, opts...) = new(presolve, opts)
end

callback_abort(stat, tree) = (stat == :Exit && GLPK.ios_terminate(tree))

function _internal_callback(tree::Ptr{Void}, info::Ptr{Void})
    cb_data = unsafe_pointer_to_objref(info)::GLPKCallbackData
    lpm = cb_data.model
    cb_data.tree = tree

    reason = GLPK.ios_reason(tree)
    cb_data.reason = reason

    if reason == GLPK.ISELECT
        #println("reason=SELECT")
        cb_data.state = :Intermediate
    elseif reason == GLPK.IPREPRO
        #println("reason=PREPRO")
        cb_data.state = :MIPNode
    elseif reason == GLPK.IROWGEN
        #println("reason=ROWGEN")
        cb_data.state = :MIPNode
        # if the current solution is actually integer feasible, then
        # return MIPSol status.
        _initsolution!(cb_data)
        cbgetlpsolution(cb_data, cb_data.sol)
        # tol_int = 1e-5 by default
        # TODO: query from GLPK
        all_integer = true
        for i in 1:length(cb_data.sol)
            v::Symbol = cb_data.vartype[i]
            (v == :Int || v == :Bin) || continue
            if abs(cb_data.sol[i]-round(cb_data.sol[i])) > 1e-5
                all_integer = false
                break
            end
        end
        if all_integer
            cb_data.state = :MIPSol
        end
        fill!(cb_data.sol, NaN)

        if lpm.lazycb != nothing
            stat = lpm.lazycb(cb_data)
            callback_abort(stat,tree)
        end
    elseif reason == GLPK.IHEUR
        #println("reason=HEUR")
        cb_data.state = :MIPNode
        if lpm.heuristiccb != nothing
            stat = lpm.heuristiccb(cb_data)
            callback_abort(stat,tree)
        end
    elseif reason == GLPK.ICUTGEN
        #println("reason=CUTGEN")
        cb_data.state = :MIPNode
        if lpm.cutcb != nothing
            stat = lpm.cutcb(cb_data)
            callback_abort(stat,tree)
        end
    elseif reason == GLPK.IBRANCH
        #println("reason=BRANCH")
        cb_data.state = :MIPNode
    elseif reason == GLPK.IBINGO
        #println("reason=BINGO")
        cb_data.state = :MIPSol
    else
        error("internal library error")
    end

    # Doesn't seem like there's a natural "reason" to put this with,
    # so let's just call it everywhere for now
    if lpm.infocb != nothing
        stat = lpm.infocb(cb_data)
        callback_abort(stat,tree)
    end

    bn = GLPK.ios_best_node(tree)
    bn != 0 && (lpm.objbound = GLPK.ios_node_bound(tree, bn))

    return
end

function LinearQuadraticModel(s::GLPKSolverMIP)
    lpm = GLPKMathProgModelMIP()
    lpm.param.msg_lev = GLPK.MSG_ERR
    lpm.smplxparam.msg_lev = GLPK.MSG_ERR
    if s.presolve
        lpm.param.presolve = GLPK.ON
    end

    lpm.param.cb_func = cfunction(_internal_callback, Void, (Ptr{Void}, Ptr{Void}))
    lpm.param.cb_info = pointer_from_objref(lpm.cbdata)

    for (k,v) in s.opts
        if k in [:cb_func, :cb_info]
            warn("ignored option: $(string(k)); use the MathProgBase callback interface instead")
            continue
        end
        i = findfirst(x->x==k, fieldnames(typeof(lpm.param)))
        s = findfirst(x->x==k, fieldnames(typeof(lpm.smplxparam)))
        if (VERSION < v"0.7-" && i > 0) || (VERSION >= v"0.7-" && i !== nothing)
            t = typeof(lpm.param).types[i]
            setfield!(lpm.param, i, convert(t, v))
        elseif (VERSION < v"0.7-" && s > 0) || (VERSION >= v"0.7-" && s !== nothing)
            t = typeof(lpm.smplxparam).types[s]
            setfield!(lpm.smplxparam, s, convert(t, v))
        else
            warn("Ignored option: $(string(k))")
            continue
        end
    end

    return lpm
end

function setparameters!(s::GLPKSolverMIP; mpboptions...)
    opts = collect(Any, s.opts)
    for (optname, optval) in mpboptions
        if optname == :TimeLimit
            push!(opts, (:tm_lim,round(Int,1000*optval))) # milliseconds
        elseif optname == :Silent
            if optval == true
                push!(opts, (:msg_lev,GLPK.MSG_OFF))
            end
        else
            error("Unrecognized parameter $optname")
        end
    end
    s.opts = opts
    nothing
end

function setparameters!(m::GLPKMathProgModelMIP; mpboptions...)
    for (optname, optval) in mpboptions
        if optname == :TimeLimit
            m.param.tm_lim = round(Int,1000*optval)
        elseif optname == :Silent
            if optval == true
                m.param.msg_lev = GLPK.MSG_OFF
                m.smplxparam.msg_lev = GLPK.MSG_OFF
            end
        else
            error("Unrecognized parameter $optname")
        end
    end
end

setlazycallback!(m::GLPKMathProgModel, f::Union{Function,Void}) = (m.lazycb = f)
setcutcallback!(m::GLPKMathProgModel, f::Union{Function,Void}) = (m.cutcb = f)
setheuristiccallback!(m::GLPKMathProgModel, f::Union{Function,Void}) = (m.heuristiccb = f)
setinfocallback!(m::GLPKMathProgModel, f::Union{Function,Void}) = (m.infocb = f)

_check_tree(d::GLPKCallbackData, funcname::AbstractString) =
    (d.tree != C_NULL && d.reason != -1) || error("$funcname can only be called from within a callback")

cbgetstate(d::GLPKCallbackData) = d.state

function cbgetlpsolution(d::GLPKCallbackData, output::Vector)
    _check_tree(d, "cbgetlpsolution")
    lp = GLPK.ios_get_prob(d.tree)
    n = GLPK.get_num_cols(lp)
    length(output) >= n || error("output vector is too short")

    for c = 1:n
        output[c] = GLPK.get_col_prim(lp, c)
    end
    return output
end

function cbgetlpsolution(d::GLPKCallbackData)
    _check_tree(d, "cbgetlpsolution")
    lp = GLPK.ios_get_prob(d.tree)
    n = GLPK.get_num_cols(lp)
    output = Vector{Float64}(n)

    for c = 1:n
        output[c] = GLPK.get_col_prim(lp, c)
    end
    return output
end


function cbgetmipsolution(d::GLPKCallbackData, output::Vector)
    # assuming we're in the lazy callback where
    # the LP solution is actually integral.
    # If we add an informational callback for GLPK.IBINGO,
    # then this will need to be modified.
    return cbgetlpsolution(d, output)
end
cbgetmipsolution(d::GLPKCallbackData) = cbgetlpsolution(d)

function cbgetbestbound(d::GLPKCallbackData)
    _check_tree(d, "cbbestbound")
    lpm = d.model
    return lpm.objbound
end

function cbgetobj(d::GLPKCallbackData)
    _check_tree(d, "cbgetobj")
    lp = GLPK.ios_get_prob(d.tree)
    return GLPK.mip_obj_val(lp)
end

function cbgetexplorednodes(d::GLPKCallbackData)
    _check_tree(d, "cbgetexplorednodes")
    a, _, t = GLPK.ios_tree_size(d.tree)
    return t - a
end

function cbaddlazy!(d::GLPKCallbackData, colidx::Vector, colcoef::Vector, sense::Char, rhs::Real)
    #println("Adding lazy")
    (d.tree != C_NULL && d.reason == GLPK.IROWGEN) ||
        error("cbaddlazy! can only be called from within a lazycallback")
    length(colidx) == length(colcoef) || error("colidx and colcoef have different legths")
    if sense == '='
        bt = GLPK.FX
        rowlb = rhs
        rowub = rhs
    elseif sense == '<'
        bt = GLPK.UP
        rowlb = -Inf
        rowub = rhs
    elseif sense == '>'
        bt = GLPK.LO
        rowlb = rhs
        rowub = Inf
    else
        error("sense must be '=', '<' or '>'")
    end
    # allocating a new vector is not efficient
    solution = cbgetmipsolution(d)
    # if the cut does not exclude the current solution, ignore it
    val = dot(colcoef,solution[colidx])
    if (rowlb - 1e-8 <= val <= rowub + 1e-8)
        # would be better to use GLPK's internal tolerances
        #Base.warn_once("Ignoring lazy constraint which is already satisfied")
        return
    end

    lp = GLPK.ios_get_prob(d.tree)
    GLPK.add_rows(lp, 1)
    m = GLPK.get_num_rows(lp)
    GLPK.set_mat_row(lp, m, colidx, colcoef)
    GLPK.set_row_bnds(lp, m, bt, rowlb, rowub)
    return
end

function cbaddcut!(d::GLPKCallbackData, colidx::Vector, colcoef::Vector, sense::Char, rhs::Real)
    #println("Adding cut")
    (d.tree != C_NULL && d.reason == GLPK.ICUTGEN) ||
        error("cbaddcut! can only be called from within a cutcallback")
    if sense == '<'
        bt = GLPK.UP
    elseif sense == '>'
        bt = GLPK.LO
    elseif sense == '='
        error("unsupported sense in cut plane '='")
    else
        error("sense must be '<' or '>'")
    end
    GLPK.ios_add_row(d.tree, "", 101, colidx, colcoef, bt, rhs)
    return
end

function _initsolution!(d::GLPKCallbackData)
    lp = GLPK.ios_get_prob(d.tree)
    n = GLPK.get_num_cols(lp)
    length(d.sol) == n && return
    resize!(d.sol, n)
    fill!(d.sol, NaN)
    return
end

function _fillsolution!(d::GLPKCallbackData)
    lp = GLPK.ios_get_prob(d.tree)
    n = GLPK.get_num_cols(lp)
    sol = d.sol
    for c = 1:n
        isnan(sol[c]) || continue
        sol[c] = GLPK.mip_col_val(lp, c)
    end
end

function cbaddsolution!(d::GLPKCallbackData)
    #println("Adding sol")
    (d.tree != C_NULL && d.reason == GLPK.IHEUR) ||
        error("cbaddsolution! can only be called from within a heuristiccallback")
    _initsolution!(d)
    _fillsolution!(d)
    # test feasibility of solution, would be better if GLPK supported this
    l = getvarLB(d.model)
    u = getvarUB(d.model)
    for i in 1:length(l)
        if d.sol[i] < l[i] - 1e-6 || d.sol[i] > u[i] + 1e-6
            warn("Ignoring infeasible solution from heuristic callback")
            return
        end
    end
    A = getconstrmatrix(d.model)
    lb = getconstrLB(d.model)
    ub = getconstrUB(d.model)
    y = A*d.sol
    for i in 1:length(lb)
        if y[i] < lb[i] - 1e-6 || y[i] > ub[i] + 1e-6
            warn("Ignoring infeasible solution from heuristic callback")
            return
        end
    end
    GLPK.ios_heur_sol(d.tree, d.sol)
    fill!(d.sol, NaN)
end

function cbsetsolutionvalue!(d::GLPKCallbackData, idx::Integer, val::Real)
    _check_tree(d, "cbsetsolutionvalue!")
    _initsolution!(d)
    d.sol[idx] = val
end

function setsense!(lpm::GLPKMathProgModelMIP, sense)
    lp = lpm.inner
    if sense == :Min
        GLPK.set_obj_dir(lp, GLPK.MIN)
        lpm.objbound = -Inf
    elseif sense == :Max
        GLPK.set_obj_dir(lp, GLPK.MAX)
        lpm.objbound = Inf
    else
        error("unrecognized objective sense: $sense")
    end
end

function setvartype!(lpm::GLPKMathProgModelMIP, vartype::Vector{Symbol})
    lp = lpm.inner
    lpm.binaries = Int[]
    ncol = numvar(lpm)
    @assert length(vartype) == ncol
    for i in 1:ncol
        if vartype[i] == :Int
            coltype = GLPK.IV
        elseif vartype[i] == :Cont
            coltype = GLPK.CV
        elseif vartype[i] == :Bin
            push!(lpm.binaries, i)
            coltype = GLPK.IV
        else
            error("invalid variable type: $(vartype[i])")
        end
        GLPK.set_col_kind(lp, i, coltype)
    end
end

const vartype_map = Dict(
    GLPK.CV => :Cont,
    GLPK.IV => :Int,
    GLPK.BV => :Bin
)

function getvartype(lpm::GLPKMathProgModelMIP)
    lp = lpm.inner
    ncol = numvar(lpm)
    coltype = Array{Symbol}(ncol)
    for i in 1:ncol
        ct = GLPK.get_col_kind(lp, i)
        coltype[i] = vartype_map[ct]
        if i in lpm.binaries
            coltype[i] = :Bin
        elseif coltype[i] == :Bin # GLPK said it was binary, but we didn't tell it
            coltype[i] = :Int
        end
    end
    return coltype
end

function optimize!(lpm::GLPKMathProgModelMIP)
    vartype = getvartype(lpm)
    lb = getvarLB(lpm)
    ub = getvarUB(lpm)
    old_lb = copy(lb)
    old_ub = copy(ub)
    for c in 1:length(vartype)
        vartype[c] in [:Int,:Bin] && (lb[c] = ceil(lb[c]); ub[c] = floor(ub[c]))
        vartype[c] == :Bin && (lb[c] = max(lb[c],0.0); ub[c] = min(ub[c],1.0))
    end
    lpm.cbdata.vartype = vartype
    try
        setvarLB!(lpm, lb)
        setvarUB!(lpm, ub)
        if lpm.param.presolve == GLPK.OFF
            ret_ps = GLPK.simplex(lpm.inner, lpm.smplxparam)
            ret_ps != 0 && return ret_ps
        end
        ret = GLPK.intopt(lpm.inner, lpm.param)
        if ret == GLPK.EMIPGAP || ret == GLPK.ETMLIM || ret == GLPK.ESTOP
            lpm.userlimit = true
        end
    finally
        setvarLB!(lpm, old_lb)
        setvarUB!(lpm, old_ub)
    end
end

function status(lpm::GLPKMathProgModelMIP)
    if lpm.userlimit
        return :UserLimit
    end
    s = GLPK.mip_status(lpm.inner)
    if s == GLPK.UNDEF
        if lpm.param.presolve == GLPK.OFF && GLPK.get_status(lpm.inner) == GLPK.NOFEAS
            return :Infeasible
        else
            return :Error
        end
    end
    if s == GLPK.OPT
        return :Optimal
    elseif s == GLPK.INFEAS
        return :Infeasible
    elseif s == GLPK.UNBND
        return :Unbounded
    elseif s == GLPK.FEAS
        return :Feasible
    elseif s == GLPK.NOFEAS
        return :Infeasible
    elseif s == GLPK.UNDEF
        return :Undefined
    else
        error("internal library error")
    end
end

function getobjval(lpm::GLPKMathProgModelMIP)
    status = GLPK.mip_status(lpm.inner)
    if status == GLPK.UNDEF || status == GLPK.NOFEAS
        # no feasible solution so objective is NaN
        return NaN
    end

    return GLPK.mip_obj_val(lpm.inner)
end

function getobjbound(lpm::GLPKMathProgModelMIP)
    # This is a hack. We observed some cases where mip_status == OPT
    # and objval and objbound didn't match.
    # We can fix this case, but objbound may still be incorrect in
    # cases where the solver terminates early.
    if GLPK.mip_status(lpm.inner) == GLPK.OPT
        return GLPK.mip_obj_val(lpm.inner)
    else
        return lpm.objbound
    end
end

function getsolution(lpm::GLPKMathProgModelMIP)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    status = GLPK.mip_status(lpm.inner)
    if status == GLPK.UNDEF || status == GLPK.NOFEAS
        # no feasible solution to return
        return fill(NaN, n)
    end

    x = Array{Float64}(n)
    for c = 1:n
        x[c] = GLPK.mip_col_val(lp, c)
    end
    return x
end

function getconstrsolution(lpm::GLPKMathProgModelMIP)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    x = Array{Float64}(m)
    for r = 1:m
        x[r] = GLPK.mip_row_val(lp, r)
    end
    return x
end

end

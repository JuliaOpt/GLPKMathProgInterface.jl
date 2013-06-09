module GLPKInterfaceMIP

import GLPK
importall LinprogSolverInterface
importall ..GLPKInterfaceBase

export
    GLPKSolverMIP,
    model,
    optimize,
    loadproblem,
    writeproblem,
    getvarLB,
    setvarLB,
    getvarLB,
    setvarLB,
    getconstrLB,
    setconstrLB,
    getconstrUB,
    setconstrUB,
    getobj,
    setobj,
    addvar,
    addconstr,
    updatemodel,
    setsense,
    getsense,
    numvar,
    numconstr,
    setvartype,
    status,
    getobjval,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver

type GLPKSolverMIP <: GLPKSolver
    inner::GLPK.Prob
    param::GLPK.IntoptParam
    smplxparam::GLPK.SimplexParam
    objbound::Vector{Float64}
end

function model(;GLPKpresolve=false, kwargs...)
    if length(kwargs) != 0
        warn("Unknown option(s) to GLPK MIP solver: ", join([string(x[1]) for x in kwargs], ", "))
    end
    lpm = GLPKSolverMIP(GLPK.Prob(), GLPK.IntoptParam(), GLPK.SimplexParam(), [-Inf])
    lpm.param.msg_lev = GLPK.MSG_ERR
    lpm.smplxparam.msg_lev = GLPK.MSG_ERR
    if GLPKpresolve
        lpm.param.presolve = GLPK.ON
    end

    function cb_callback(tree::Ptr{Void}, info::Ptr{Void})
        bn = GLPK.ios_best_node(tree)
        if bn == 0
            return
        end
        ret = pointer_to_array(convert(Ptr{Float64}, info), 1, false)
        ret[1] = GLPK.ios_node_bound(tree, bn)
        return
    end
    lpm.param.cb_func = cfunction(cb_callback, Void, (Ptr{Void}, Ptr{Void}))
    lpm.param.cb_info = convert(Ptr{Void}, lpm.objbound)

    return lpm
end

function setsense(lpm::GLPKSolverMIP, sense)
    lp = lpm.inner
    if sense == :Min
        GLPK.set_obj_dir(lp, GLPK.MIN)
        lpm.objbound[1] = -Inf
    elseif sense == :Max
        GLPK.set_obj_dir(lp, GLPK.MAX)
        lpm.objbound[1] = Inf
    else
        error("Unrecognized objective sense $sense")
    end
end

function setvartype(lpm::GLPKSolverMIP, vartype)
    lp = lpm.inner
    ncol = numvar(lpm)
    @assert length(vartype) == ncol
    coltype = map(vartype) do c
        if c == 'I'
            return GLPK.IV
        elseif c == 'C'
            return GLPK.CV
        else
            error("invalid var type $c")
        end
    end
    for i in 1:ncol
        GLPK.set_col_kind(lp, i, coltype[i])
    end
end

function getvartype(lpm::GLPKSolverMIP)
    lp = lpm.inner
    ncol = numvar(m)
    coltype = Array(Char, ncol)
    for i in 1:ncol
        ct = GLPK.get_col_kind(lp, i)
        coltype[i] = (ct == GLPK.CV ? 'C' : 'I')
    end
    return coltype
end

function optimize(lpm::GLPKSolverMIP)
    if lpm.param.presolve == GLPK.OFF
        ret_ps = GLPK.simplex(lpm.inner, lpm.smplxparam)
        ret_ps != 0 && return ret_ps
    end
    GLPK.intopt(lpm.inner, lpm.param)
end

function status(lpm::GLPKSolverMIP)
   s = GLPK.mip_status(lpm.inner)
   if s == GLPK.UNDEF && lpm.param.presolve == GLPK.OFF
       s = GLPK.get_status(lpm.inner)
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
       error("Internal library error")
   end
end

getobjval(lpm::GLPKSolverMIP) = GLPK.mip_obj_val(lpm.inner)

getobjbound(lp::GLPKSolverMIP) = lp.objbound[1]

function getsolution(lpm::GLPKSolverMIP)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    x = Array(Float64, n)
    for c = 1:n
        x[c] = GLPK.mip_col_val(lp, c)
    end
    return x
end

function getconstrsolution(lpm::GLPKSolverMIP)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    x = Array(Float64, m)
    for r = 1:m
        x[r] = GLPK.mip_row_val(lp, r)
    end
    return x
end

end

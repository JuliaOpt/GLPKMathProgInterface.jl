module GLPKInterfaceLP

import GLPK
importall LinprogSolverInterface
importall ..GLPKInterfaceBase

export
    GLPKSolverLP,
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
    status,
    getobjval,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver

type GLPKSolverLP <: GLPKSolver
    inner::GLPK.Prob
    method::Symbol
    param::Union(GLPK.SimplexParam, GLPK.InteriorParam)
end

function model(;GLPKpresolve=false, GLPKmethod=:Simplex, kwargs...)
    if length(kwargs) != 0
        warn("Unknown option(s) to GLPK LP solver: ", join([string(x[1]) for x in kwargs], ", "))
    end
    if GLPKmethod == :Simplex || GLPKmethod == :Exact
        param = GLPK.SimplexParam()
        if GLPKpresolve
            param.GLPKpresolve = GLPK.ON
        end
    elseif GLPKmethod == :InteriorPoint
        param = GLPK.InteriorParam()
        if GLPKpresolve
            warn("Ignored option: GLPKpresolve")
        end
    else
        error("""
              Unknown method for GLPK LP solver: $GLPKmethod
                     Allowed methods:
                       :Simplex
                       :Exact
                       :InteriorPoint""")
    end
    param.msg_lev = GLPK.MSG_ERR
    lpm = GLPKSolverLP(GLPK.Prob(), GLPKmethod, param)
    return lpm
end

function optimize(lpm::GLPKSolverLP)
    if lpm.method == :Simplex
        solve = GLPK.simplex
    elseif lpm.method == :Exact
        solve = GLPK.exact
    elseif lpm.method == :InteriorPoint
        solve = GLPK.interior
    else
        error("bug")
    end
    return solve(lpm.inner, lpm.param)
end

function status(lpm::GLPKSolverLP)
   if lpm.method == :Simplex || lpm.method == :Exact
       get_status = GLPK.get_status
   elseif lpm.method == :InteriorPoint
       get_status = GLPK.ipt_status
   else
       error("bug")
   end
   s = get_status(lpm.inner)
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

function getobjval(lpm::GLPKSolverLP)
    if lpm.method == :Simplex || lpm.method == :Exact
        get_obj_val = GLPK.get_obj_val
    elseif lpm.method == :InteriorPoint
        get_obj_val = GLPK.ipt_obj_val
    else
        error("bug")
    end
    return get_obj_val(lpm.inner)
end

function getsolution(lpm::GLPKSolverLP)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_col_prim = GLPK.get_col_prim
    elseif lpm.method == :InteriorPoint
        get_col_prim = GLPK.ipt_col_prim
    else
        error("bug")
    end

    x = Array(Float64, n)
    for c = 1:n
        x[c] = get_col_prim(lp, c)
    end
    return x
end

function getconstrsolution(lpm::GLPKSolverLP)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_row_prim = GLPK.get_row_prim
    elseif lpm.method == :InteriorPoint
        get_row_prim = GLPK.ipt_row_prim
    else
        error("bug")
    end

    x = Array(Float64, m)
    for r = 1:m
        x[r] = get_row_prim(lp, r)
    end
    return x
end

function getreducedcosts(lpm::GLPKSolverLP)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_col_dual = GLPK.get_col_dual
    elseif lpm.method == :InteriorPoint
        get_col_dual = GLPK.ipt_col_dual
    else
        error("bug")
    end

    x = Array(Float64, n)
    for c = 1:n
        x[c] = get_col_dual(lp, c)
    end
    return x
end

function getconstrduals(lpm::GLPKSolverLP)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_row_dual = GLPK.get_row_dual
    elseif lpm.method == :InteriorPoint
        get_row_dual = GLPK.ipt_row_dual
    else
        error("bug")
    end

    x = Array(Float64, m)
    for r = 1:m
        x[r] = get_row_dual(lp, r)
    end
    return x
end

end

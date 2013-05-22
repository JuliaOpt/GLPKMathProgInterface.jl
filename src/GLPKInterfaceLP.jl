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
    param::GLPK.SimplexParam
end

function model(;kwargs...)
    if length(kwargs) != 0
        warn("GLPK LP solver does not yet support options")
    end
    lpm = GLPKSolverLP(GLPK.Prob(), GLPK.SimplexParam())
    lpm.param.msg_lev = GLPK.MSG_ERR
    #lpm.param.presolve = GLPK.ON
    return lpm
end

optimize(lpm::GLPKSolverLP) = GLPK.simplex(lpm.inner, lpm.param)

function status(lpm::GLPKSolverLP)
   s = GLPK.get_status(lpm.inner)
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

getobjval(lpm::GLPKSolverLP) = GLPK.get_obj_val(lpm.inner)

function getsolution(lpm::GLPKSolverLP)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    x = Array(Float64, n)
    for c = 1:n
        x[c] = GLPK.get_col_prim(lp, c)
    end
    return x
end

function getconstrsolution(lpm::GLPKSolverLP)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    x = Array(Float64, m)
    for r = 1:m
        x[r] = GLPK.get_row_prim(lp, r)
    end
    return x
end

function getreducedcosts(lpm::GLPKSolverLP)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    x = Array(Float64, n)
    for c = 1:n
        x[c] = GLPK.get_col_dual(lp, c)
    end
    return x
end

function getconstrduals(lpm::GLPKSolverLP)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    x = Array(Float64, m)
    for r = 1:m
        x[r] = GLPK.get_row_dual(lp, r)
    end
    return x
end


end

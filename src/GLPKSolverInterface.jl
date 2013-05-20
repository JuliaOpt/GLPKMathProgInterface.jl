module GLPKSolverInterface

import GLPK

require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
importall LinprogSolverInterface

export GLPKSolver,
    model,
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
    optimize,
    status,
    getobjval,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver


type GLPKSolver <: LinprogSolver
    inner::GLPK.Prob()
    params::GLPK.SimplexParam()
end

function model()
    lpm = GLPKSolver(GLPK.SimplexParam())
    lpm.param.msg_lev = GLPK.MSG_ERR
    lpm.param.presolve = GLPK.ON
    return lpm
end

function loadproblem(lpm::GLPKSolver, filename::String)
    if endswith(filename, ".mps") || endswith(filename, ".mps.gz")
       read_mps(lpm.inner, GLPK.MPS_FILE, filename)
   elseif endswith(filename, ".lp") || endswith(filename, ".lp.gz")
       read_lp(lpm.inner, filename)
   elseif endswith(filename, ".prob") || endswith(filename, ".prob.gz")
       read_prob(lpm.inner, filename)
   else
       error("unrecognized input format extension in $filename")
    end
end   

nonnull(x) = (x != Nothing && !isempty(x))

function loadproblem(lpm::GLPKSolver, A::AbstractMatrix, collb, colub, obj, rowlb, rowub)
    lp = lpm.inner

    m, n = size(A)

    if m == 0 || n == 0
        error("empty problem matrix A")
    end

    function checksize(x, l, str) = 
        if nonnull(x) && length(x) != l
            error("size of $str is incompatible with size of A")
        end
    end

    checksize(collb, m, "collb")
    checksize(colub, m, "colub")
    checksize(rowlb, n, "rowlb")
    checksize(rowub, n, "rowub")
    checksize(obj, n, "obj")

    (ia, ja, ar) = findnz(A)

    glp_erase_prob(lp)

    function getbounds(lb, ub, i)
        if nonnull(lb) && nonnull(ub) && lb[i] != -Inf && ub[i] != Inf
            if lb[i] != ub[i]
                return lb[i], ub[i], GLPK.DB
            else
                return lb[i], ub[i], GLPK.FX
            end
        elseif nonnull(lb) && lb[i] != -Inf
            return lb[i], Inf, GLPK.LB
        elseif nonnull(ub) && ub[i] != Inf
            return -Inf, ub[i], GLPK.UB
        else
            return -Inf, Inf, GLPK.FR
        end
    end

    GLPK.add_rows(lp, m)
    for r = 1:m
        #println("  r=$r b=$(b[r])")
        l, u, t = getbounds(rowlb, rowub, r)
        GLPK.set_row_bnds(lp, r, t, l, u)
    end

    GLPK.add_cols(lp, n)
    if nonnull(obj)
        for c = 1:n
            GLPK.set_obj_coef(lp, c, obj[c])
        end
    else
        for c = 1:n
            GLPK.set_obj_coef(lp, c, 0)
        end
    end

    GLPK.add_cols(lp, n)
    for c = 1:n
        #println("  r=$r b=$(b[r])")
        l, u, t = getbounds(collb, colub, c)
        GLPK.set_col_bnds(lp, c, t, l, u)
    end

    GLPK.load_matrix(lp, ia, ja, ar)
    return lpm
end

#writeproblem(m, filename::String)

function getvarLB(lpm::GLPKSolver)
    lp = lpm.inner
    n = get_num_cols(lp)
    lb = Array(Float64, n)
    for c = 1:n
        l = get_col_lb(lp, c)
        if l <= -realmax(Float64)
            l = -Inf
        end
        lb[c] = l
    end
    return lb
end

function setvarLB(lpm::GLPKSolver, collb)
    lp = lpm.inner
    n = get_num_cols(lp)
    if nonnull(collb) && length(collb) != n
        error("invalid size of collb")
    end
    for c = 1:n
        u = get_col_ub(lp, c)
        if u >= realmax(Float64)
            u = Inf
        end
        if nonnull(collb) && collb[c] != -Inf
            l = collb[c]
            if u < Inf
                if l != u
                    GLPK.set_col_bnds(lp, c, GLPK.DB, l, u)
                else
                    GLPK.set_col_bnds(lp, c, GLPK.FX, l, u)
                end
            else
                GLPK.set_col_bnds(lp, c, GLPK.LB, l, 0.0)
            end
        else
            if u < Inf
                GLPK.set_col_bnds(lp, c, GLPK.UB, 0.0, u)
            else
                GLPK.set_col_bnds(lp, c, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getvarUB(lpm::GLPKSolver)
    lp = lpm.inner
    n = get_num_cols(lp)
    ub = Array(Float64, n)
    for c = 1:n
        u = get_col_ub(lp, c)
        if u >= realmax(Float64)
            u = Inf
        end
        ub[c] = u
    end
    return ub
end

function setvarUB(lpm::GLPKSolver, colub)
    lp = lpm.inner
    n = get_num_cols(lp)
    if nonnull(colub) && length(colub) != n
        error("invalid size of colub")
    end
    for c = 1:n
        l = get_col_lb(lp, c)
        if l <= -realmax(Float64)
            l = -Inf
        end
        if nonnull(colub) && colub[c] != Inf
            u = colub[c]
            if l > -Inf
                if l != u
                    GLPK.set_col_bnds(lp, c, GLPK.DB, l, u)
                else
                    GLPK.set_col_bnds(lp, c, GLPK.FX, l, u)
                end
            else
                GLPK.set_col_bnds(lp, c, GLPK.UB, 0.0, u)
            end
        else
            if l > -Inf
                GLPK.set_col_bnds(lp, c, GLPK.LB, l, 0.0)
            else
                GLPK.set_col_bnds(lp, c, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getconstrLB(lpm::GLPKSolver)
    lp = lpm.inner
    m = get_num_rows(lp)
    lb = Array(Float64, m)
    for r = 1:m
        l = get_row_lb(lp, r)
        if l <= -realmax(Float64)
            l = -Inf
        end
        lb[r] = l
    end
    return lb
end

function setconstrLB(lpm::GLPKSolver, rowlb)
    lp = lpm.inner
    m = get_num_rows(lp)
    if nonnull(collb) && length(rowlb) != m
        error("invalid size of rowlb")
    end
    for r = 1:m
        u = get_row_ub(lp, r)
        if u >= realmax(Float64)
            u = Inf
        end
        if nonnull(rowlb) && rowlb[r] != -Inf
            l = rowlb[c]
            if u < Inf
                if l != u
                    GLPK.set_row_bnds(lp, r, GLPK.DB, l, u)
                else
                    GLPK.set_row_bnds(lp, r, GLPK.FX, l, u)
                end
            else
                GLPK.set_row_bnds(lp, r, GLPK.LB, l, 0.0)
            end
        else
            if u < Inf
                GLPK.set_row_bnds(lp, r, GLPK.UB, 0.0, u)
            else
                GLPK.set_row_bnds(lp, r, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getconstrUB(lpm::GLPKSolver)
    lp = lpm.inner
    m = get_num_rows(lp)
    ub = Array(Float64, m)
    for r = 1:m
        u = get_col_ub(lp, r)
        if u >= realmax(Float64)
            u = Inf
        end
        ub[r] = u
    end
    return ub
end

function setconstrUB(lpm::GLPKSolver, rowub)
    lp = lpm.inner
    m = get_num_rows(lp)
    if nonnull(rowub) && length(rowub) != m
        error("invalid size of rowub")
    end
    for r = 1:m
        l = get_row_lb(lp, r)
        if l <= -realmax(Float64)
            l = -Inf
        end
        if nonnull(rowub) && rowub[r] != Inf
            u = rowub[r]
            if l > -Inf
                if l != u
                    GLPK.set_row_bnds(lp, r, GLPK.DB, l, u)
                else
                    GLPK.set_row_bnds(lp, r, GLPK.FX, l, u)
                end
            else
                GLPK.set_row_bnds(lp, r, GLPK.UB, 0.0, u)
            end
        else
            if l > -Inf
                GLPK.set_row_bnds(lp, r, GLPK.LB, l, 0.0)
            else
                GLPK.set_row_bnds(lp, r, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getobj(lpm::GLPKSolver)
    lp = lpm.inner
    n = get_num_cols(lp)
    obj = Array(Float64, n)
    for c = 1:n
        l = get_obj_coef(lp, c)
        obj[c] = l
    end
    return obj
end

function setobj(lpm::GLPKSolver, obj)
    lp = lpm.inner
    n = get_num_cols(lp)
    if nonnull(obj) && length(obj) != n
        error("invalid size of obj")
    end
    for c = 1:n
        if nonnull(obj)
            GLPK.set_obj_coef(lp, c, obj[c])
        else
            GLPK.set_obj_coef(lp, c, 0.0)
        end
    end
end

function addvar(lpm::GLPKSolver, rowidx::Vector, rowcoef::Vector, collb::Real, colub::Real, objcoef:Real)
    if length(rowidx) != length(rowcoef)
        error("rowidx and rowcoef have different legths")
    end
    lp = lpm.inner
    GLPK.add_cols(lp, 1)
    n = GLPK.get_num_cols(lp)
    GLPK.set_mat_col(lp, n, rowidx, rowcoef)
    if collb > -Inf && colub < Inf
        if collb != colub
            bt = GLPK.DB
        else
            bt = GLPK.FX
        end
    elseif collb > -Inf
        bt = GLPK.LB
    elseif colub < Inf
        bt = GLPK.UB
    else
        bt = GLPK.FR
    end
    GLPK.set_col_bnds(lp, n, bt, collb, colub)
    GLPK.set_obj_coef(lp, n, objcoef)
    return
end

function addconstr(lpm::GLPKSolver, colidx::Vector, colcoef::Vector, rowlb::Real, rowub::Real)
    if length(colidx) != length(colcoef)
        error("colidx and colcoef have different legths")
    end
    lp = lpm.inner
    GLPK.add_rows(lp, 1)
    m = GLPK.get_num_rows(lp)
    GLPK.set_mat_row(lp, m, colidx, colcoef)
    if rowlb > -Inf && rowub < Inf
        if rowlb != rowub
            bt = GLPK.DB
        else
            bt = GLPK.FX
        end
    elseif rowlb > -Inf
        bt = GLPK.LB
    elseif rowub < Inf
        bt = GLPK.UB
    else
        bt = GLPK.FR
    end
    GLPK.set_row_bnds(lp, m, bt, rowlb, rowub)
    return
end

updatemodel(m::GLPKSolver) = nothing

function setsense(lpm::GLPKSolver, sense)
    lp = lpm.inner
    if sense == :Min
        GLPK.set_obj_dir(lp, GLPK.MIN)
    elseif sense == :Max
        GLPK.set_obj_dir(lp, GLPK.MAX)
    else
        error("Unrecognized objective sense $sense")
    end
end

function getsense(lpm::GLPKSolver)
    lp = lpm.inner
    s = GLPK.get_obj_dir(lp)
    if s == GLPK.MIN
        return :Min
    elseif s == GLPK.MAX
        return :Max
    else
        error("Internal library error")
    end
end

numvar(lpm::GLPKSolver) = GLPK.get_num_cols(lpm.inner) 
numconstr(lpm::GLPKSolver) = GLPK.get_num_rows(lpm.inner)

optimize(lpm::GLPKSolver) = GLPK.simplex(lpm.inner, lpm.params)

function status(lpm::GLPKSolver)
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
       return :NoFeasible
   elseif s == GLPK.UNDEF
       return :Undefined
   else
       error("Internal library error")
   end
end

getobjval(lpm::GLPKSolver) = GLPK.get_obj_val(lpm.inner)

# TODO

function getsolution(lpm::GLPKSolver)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    x = Array(Float64, n)
    for c = 1:n
        x[c] = GLPK.get_col_prim(lp, c)
    end
    return x
end

function getconstrsolution(lpm::GLPKSolver)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    x = Array(Float64, m)
    for r = 1:m
        x[r] = GLPK.get_row_prim(lp, r)
    end
    return x
end

function getreducedcosts(lpm::GLPKSolver)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    x = Array(Float64, n)
    for c = 1:n
        x[c] = GLPK.get_col_dual(lp, c)
    end
    return x
end

function getconstrduals(lpm::GLPKSolver)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    x = Array(Float64, m)
    for r = 1:m
        x[r] = GLPK.get_row_dual(lp, r)
    end
    return x
end

getrawsolver(lpm::GLPKSolver) = lpm.inner

end

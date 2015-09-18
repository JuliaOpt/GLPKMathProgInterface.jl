module GLPKInterfaceBase

using Compat

import GLPK

importall MathProgBase.SolverInterface

export
    GLPKMathProgModel,
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
    getconstrmatrix,
    getobj,
    setobj!,
    addvar!,
    addconstr!,
    updatemodel!,
    setsense!,
    getsense,
    numvar,
    numconstr,
    getrawsolver


abstract GLPKMathProgModel <: AbstractMathProgModel

function loadproblem!(lpm::GLPKMathProgModel, filename::AbstractString)
    if endswith(filename, ".mps") || endswith(filename, ".mps.gz")
       GLPK.read_mps(lpm.inner, GLPK.MPS_FILE, filename)
    elseif endswith(filename, ".lp") || endswith(filename, ".lp.gz")
       GLPK.read_lp(lpm.inner, filename)
    elseif endswith(filename, ".prob") || endswith(filename, ".prob.gz")
       GLPK.read_prob(lpm.inner, filename)
    else
       error("unrecognized input format extension in $filename")
    end
end

nonnull(x) = (x != nothing && !isempty(x))

function loadproblem!(lpm::GLPKMathProgModel, A::AbstractMatrix, collb, colub, obj, rowlb, rowub, sense)
    lp = lpm.inner

    m, n = size(A)

    function checksize(x, l, str)
        if nonnull(x) && length(x) != l
            error("size of $str is incompatible with size of A")
        end
    end

    checksize(collb, n, "collb")
    checksize(colub, n, "colub")
    checksize(rowlb, m, "rowlb")
    checksize(rowub, m, "rowub")
    checksize(obj, n, "obj")

    (ia, ja, ar) = findnz(A)

    GLPK.erase_prob(lp)

    function getbounds(lb, ub, i)
        if nonnull(lb) && nonnull(ub) && lb[i] != -Inf && ub[i] != Inf
            if lb[i] != ub[i]
                return lb[i], ub[i], GLPK.DB
            else
                return lb[i], ub[i], GLPK.FX
            end
        elseif nonnull(lb) && lb[i] != -Inf
            return lb[i], Inf, GLPK.LO
        elseif nonnull(ub) && ub[i] != Inf
            return -Inf, ub[i], GLPK.UP
        else
            return -Inf, Inf, GLPK.FR
        end
    end

    m > 0 && GLPK.add_rows(lp, m)
    for r = 1:m
        #println("  r=$r b=$(b[r])")
        l, u, t = getbounds(rowlb, rowub, r)
        GLPK.set_row_bnds(lp, r, t, l, u)
    end

    n > 0 && GLPK.add_cols(lp, n)
    for c = 1:n
        #println("  r=$r b=$(b[r])")
        l, u, t = getbounds(collb, colub, c)
        GLPK.set_col_bnds(lp, c, t, l, u)
    end

    if nonnull(obj)
        for c = 1:n
            GLPK.set_obj_coef(lp, c, obj[c])
        end
    else
        for c = 1:n
            GLPK.set_obj_coef(lp, c, 0)
        end
    end

    m > 0 && n > 0 && GLPK.load_matrix(lp, ia, ja, ar)
    setsense!(lpm, sense)

    return lpm
end

#writeproblem(m, filename::AbstractString)

function getvarLB(lpm::GLPKMathProgModel)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    lb = Array(Float64, n)
    for c = 1:n
        l = GLPK.get_col_lb(lp, c)
        if l <= -realmax(Float64)
            l = -Inf
        end
        lb[c] = l
    end
    return lb
end

function setvarLB!(lpm::GLPKMathProgModel, collb)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    if nonnull(collb) && length(collb) != n
        error("invalid size of collb")
    end
    for c = 1:n
        u = GLPK.get_col_ub(lp, c)
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
                GLPK.set_col_bnds(lp, c, GLPK.LO, l, 0.0)
            end
        else
            if u < Inf
                GLPK.set_col_bnds(lp, c, GLPK.UP, 0.0, u)
            else
                GLPK.set_col_bnds(lp, c, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getvarUB(lpm::GLPKMathProgModel)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    ub = Array(Float64, n)
    for c = 1:n
        u = GLPK.get_col_ub(lp, c)
        if u >= realmax(Float64)
            u = Inf
        end
        ub[c] = u
    end
    return ub
end

function setvarUB!(lpm::GLPKMathProgModel, colub)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    if nonnull(colub) && length(colub) != n
        error("invalid size of colub")
    end
    for c = 1:n
        l = GLPK.get_col_lb(lp, c)
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
                GLPK.set_col_bnds(lp, c, GLPK.UP, 0.0, u)
            end
        else
            if l > -Inf
                GLPK.set_col_bnds(lp, c, GLPK.LO, l, 0.0)
            else
                GLPK.set_col_bnds(lp, c, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getconstrLB(lpm::GLPKMathProgModel)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)
    lb = Array(Float64, m)
    for r = 1:m
        l = GLPK.get_row_lb(lp, r)
        if l <= -realmax(Float64)
            l = -Inf
        end
        lb[r] = l
    end
    return lb
end

function setconstrLB!(lpm::GLPKMathProgModel, rowlb)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)
    if nonnull(rowlb) && length(rowlb) != m
        error("invalid size of rowlb")
    end
    for r = 1:m
        u = GLPK.get_row_ub(lp, r)
        if u >= realmax(Float64)
            u = Inf
        end
        if nonnull(rowlb) && rowlb[r] != -Inf
            l = rowlb[r]
            if u < Inf
                if l != u
                    GLPK.set_row_bnds(lp, r, GLPK.DB, l, u)
                else
                    GLPK.set_row_bnds(lp, r, GLPK.FX, l, u)
                end
            else
                GLPK.set_row_bnds(lp, r, GLPK.LO, l, 0.0)
            end
        else
            if u < Inf
                GLPK.set_row_bnds(lp, r, GLPK.UP, 0.0, u)
            else
                GLPK.set_row_bnds(lp, r, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getconstrUB(lpm::GLPKMathProgModel)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)
    ub = Array(Float64, m)
    for r = 1:m
        u = GLPK.get_row_ub(lp, r)
        if u >= realmax(Float64)
            u = Inf
        end
        ub[r] = u
    end
    return ub
end

function setconstrUB!(lpm::GLPKMathProgModel, rowub)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)
    if nonnull(rowub) && length(rowub) != m
        error("invalid size of rowub")
    end
    for r = 1:m
        l = GLPK.get_row_lb(lp, r)
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
                GLPK.set_row_bnds(lp, r, GLPK.UP, 0.0, u)
            end
        else
            if l > -Inf
                GLPK.set_row_bnds(lp, r, GLPK.LO, l, 0.0)
            else
                GLPK.set_row_bnds(lp, r, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getconstrmatrix(lpm::GLPKMathProgModel)
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)
    n = GLPK.get_num_cols(lp)
    colwise = [GLPK.get_mat_col(lp,i) for i in 1:n]
    nnz = 0
    for i in 1:n
        nnz += length(colwise[i][1])
    end
    colptr = Array(Int,n+1)
    rowval = Array(Int,nnz)
    nzval  = Array(Float64,nnz)
    cur_nnz = 1
    for i in 1:n
        colptr[i] = cur_nnz
        ind,vals = colwise[i]
        p = sortperm(ind) # indices must be sorted
        rowval[cur_nnz:(cur_nnz+length(ind)-1)] = ind[p]
        nzval[cur_nnz:(cur_nnz+length(ind)-1)] = vals[p]
        cur_nnz += length(ind)
    end
    colptr[n+1] = cur_nnz
    return SparseMatrixCSC(m,n,colptr,rowval,nzval)
end


function getobj(lpm::GLPKMathProgModel)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    obj = Array(Float64, n)
    for c = 1:n
        l = GLPK.get_obj_coef(lp, c)
        obj[c] = l
    end
    return obj
end

function setobj!(lpm::GLPKMathProgModel, obj)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
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

function addvar!(lpm::GLPKMathProgModel, rowidx::Vector, rowcoef::Vector, collb::Real, colub::Real, objcoef::Real)
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
        bt = GLPK.LO
    elseif colub < Inf
        bt = GLPK.UP
    else
        bt = GLPK.FR
    end
    GLPK.set_col_bnds(lp, n, bt, collb, colub)
    GLPK.set_obj_coef(lp, n, objcoef)
    return
end

function addconstr!(lpm::GLPKMathProgModel, colidx::Vector, colcoef::Vector, rowlb::Real, rowub::Real)
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
        bt = GLPK.LO
    elseif rowub < Inf
        bt = GLPK.UP
    else
        bt = GLPK.FR
    end
    GLPK.set_row_bnds(lp, m, bt, rowlb, rowub)
    return
end

updatemodel!(m::GLPKMathProgModel) = nothing

function setsense!(lpm::GLPKMathProgModel, sense)
    lp = lpm.inner
    if sense == :Min
        GLPK.set_obj_dir(lp, GLPK.MIN)
    elseif sense == :Max
        GLPK.set_obj_dir(lp, GLPK.MAX)
    else
        error("Unrecognized objective sense $sense")
    end
end

function getsense(lpm::GLPKMathProgModel)
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

numvar(lpm::GLPKMathProgModel) = GLPK.get_num_cols(lpm.inner) 
numconstr(lpm::GLPKMathProgModel) = GLPK.get_num_rows(lpm.inner)

getrawsolver(lpm::GLPKMathProgModel) = lpm.inner

end


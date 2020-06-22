module GLPKInterfaceLP

using SparseArrays
using LinearAlgebra

import GLPK
import MathProgBase
const MPB = MathProgBase
using ..GLPKInterfaceBase

export GLPKSolverLP

mutable struct GLPKMathProgModelLP <: GLPKMathProgModel
    inner::GLPK.Prob
    method::Symbol
    param::Union{GLPK.SimplexParam, GLPK.InteriorParam}
    infeasible_bounds::Bool
end

mutable struct GLPKSolverLP <: MPB.AbstractMathProgSolver
    presolve::Bool
    method::Symbol
    opts
    function GLPKSolverLP(;presolve::Bool=false, method::Symbol=:Simplex, opts...)
        method in [:Simplex, :Exact, :InteriorPoint] ||
            error("""
                  Unknown method for GLPK LP solver: $method
                         Allowed methods:
                           :Simplex
                           :Exact
                           :InteriorPoint""")
        new(presolve, method, opts)
    end
end

function Base.copy(m::GLPKMathProgModelLP)

    m2inner = GLPK.Prob()

    GLPK.copy_prob(m2inner, m.inner, GLPK.ON)

    return GLPKMathProgModelLP(m2inner, m.method, deepcopy(m.param), m.infeasible_bounds)
end

function MPB.LinearQuadraticModel(s::GLPKSolverLP)
    if s.method == :Simplex || s.method == :Exact
        param = GLPK.SimplexParam()
        if s.presolve
            param.presolve = GLPK.ON
        end
    elseif s.method == :InteriorPoint
        param = GLPK.InteriorParam()
        if s.presolve
            @warn "Ignored option: presolve"
        end
    else
        error("This is a bug")
    end
    param.msg_lev = GLPK.MSG_ERR
    for (k,v) in s.opts
        i = findfirst(x->x==k, fieldnames(typeof(param)))
        if (VERSION < v"0.7-" && i > 0) || (VERSION >= v"0.7-" && i !== nothing)
            t = typeof(param).types[i]
            setfield!(param, i, convert(t, v))
        else
            @warn "Ignored option: $(string(k))"
        end
    end
    lpm = GLPKMathProgModelLP(GLPK.Prob(), s.method, param, false)
    return lpm
end

function MPB.setparameters!(s::GLPKSolverLP; mpboptions...)
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

function MPB.setparameters!(m::GLPKMathProgModelLP; mpboptions...)
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

function MPB.optimize!(lpm::GLPKMathProgModelLP)
    lpm.infeasible_bounds = false
    lp = lpm.inner
    for c in 1:MPB.numvar(lpm)
        if GLPK.get_col_lb(lp, c) > GLPK.get_col_ub(lp, c)
            lpm.infeasible_bounds = true
            break
        end
    end
    if !lpm.infeasible_bounds
        for r in 1:MPB.numconstr(lpm)
            if GLPK.get_row_lb(lp, r) > GLPK.get_row_ub(lp, r)
                lpm.infeasible_bounds = true
                break
            end
        end
    end
    if !lpm.infeasible_bounds
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
end

function MPB.status(lpm::GLPKMathProgModelLP)
    if lpm.infeasible_bounds
        return :Infeasible
    end
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

function MPB.getobjval(lpm::GLPKMathProgModelLP)
    if lpm.infeasible_bounds
        if GLPK.get_obj_dir(lpm.inner) == GLPK.MAX
            return -Inf
        else
            return Inf
        end
    end
    if lpm.method == :Simplex || lpm.method == :Exact
        get_obj_val = GLPK.get_obj_val
    elseif lpm.method == :InteriorPoint
        get_obj_val = GLPK.ipt_obj_val
    else
        error("bug")
    end
    return get_obj_val(lpm.inner)
end

function check_feasible_bounds(lpm::GLPKMathProgModelLP, name::String)
    if lpm.infeasible_bounds
        error("$name is not available when some constraint bounds are infeasible (lower bound > upper bound)")
    end
end


function MPB.getsolution(lpm::GLPKMathProgModelLP)
    check_feasible_bounds(lpm, "getsolution")
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_col_prim = GLPK.get_col_prim
    elseif lpm.method == :InteriorPoint
        get_col_prim = GLPK.ipt_col_prim
    else
        error("bug")
    end

    return [get_col_prim(lp, i) for i in 1:n]
end

function MPB.getconstrsolution(lpm::GLPKMathProgModelLP)
    check_feasible_bounds(lpm, "getconstrsolution")
    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_row_prim = GLPK.get_row_prim
    elseif lpm.method == :InteriorPoint
        get_row_prim = GLPK.ipt_row_prim
    else
        error("bug")
    end

    return [get_row_prim(lp, i) for i in 1:m]
end

function MPB.getreducedcosts(lpm::GLPKMathProgModelLP)
    check_feasible_bounds(lpm, "getreducedcosts")
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_col_dual = GLPK.get_col_dual
    elseif lpm.method == :InteriorPoint
        get_col_dual = GLPK.ipt_col_dual
    else
        error("bug")
    end

    return [get_col_dual(lp, i) for i in 1:n]
end

function MPB.getconstrduals(lpm::GLPKMathProgModelLP)
    check_feasible_bounds(lpm, "getconstrduals")

    lp = lpm.inner
    m = GLPK.get_num_rows(lp)

    if lpm.method == :Simplex || lpm.method == :Exact
        get_row_dual = GLPK.get_row_dual
    elseif lpm.method == :InteriorPoint
        get_row_dual = GLPK.ipt_row_dual
    else
        error("bug")
    end

    return [get_row_dual(lp, i) for i in 1:m]
end

# The functions getinfeasibilityray and getunboundedray are adapted from code
# taken from the LEMON C++ optimization library. This is the copyright notice:
#
### Copyright (C) 2003-2010
### Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
### (Egervary Research Group on Combinatorial Optimization, EGRES).
###
### Permission to use, modify and distribute this software is granted
### provided that this copyright notice appears in all copies. For
### precise terms see the accompanying LICENSE file.
###
### This software is provided "AS IS" with no warranty of any kind,
### express or implied, and with no claim as to its suitability for any
### purpose.

function MPB.getinfeasibilityray(lpm::GLPKMathProgModelLP)
    if lpm.infeasible_bounds
        # See https://github.com/JuliaOpt/GLPKMathProgInterface.jl/pull/34
        return zeros(MPB.numconstr(lpm))
    end

    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
    elseif lpm.method == :InteriorPoint
        error("getinfeasibilityray is not available when using the InteriorPoint method")
    else
        error("bug")
    end

    m = GLPK.get_num_rows(lp)

    ray = zeros(m)

    ur = GLPK.get_unbnd_ray(lp)
    if ur != 0
        if ur <= m
            k = ur
            get_stat = GLPK.get_row_stat
            get_bind = GLPK.get_row_bind
            get_prim = GLPK.get_row_prim
            get_ub = GLPK.get_row_ub
        else
            k = ur - m
            get_stat = GLPK.get_col_stat
            get_bind = GLPK.get_col_bind
            get_prim = GLPK.get_col_prim
            get_ub = GLPK.get_col_ub
        end

        get_stat(lp, k) == GLPK.BS || error("unbounded ray is primal (use getunboundedray)")

        ray[get_bind(lp, k)] = (get_prim(lp, k) > get_ub(lp, k)) ? -1 : 1

        GLPK.btran(lp, ray)
    else
        eps = 1e-7
        for i = 1:m
            idx = GLPK.get_bhead(lp, i)
            if idx <= m
                k = idx
                get_prim = GLPK.get_row_prim
                get_ub = GLPK.get_row_ub
                get_lb = GLPK.get_row_lb
            else
                k = idx - m
                get_prim = GLPK.get_col_prim
                get_ub = GLPK.get_col_ub
                get_lb = GLPK.get_col_lb
            end

            res = get_prim(lp, k)
            if res > get_ub(lp, k) + eps
                ray[i] = -1
            elseif res < get_lb(lp, k) - eps
                ray[i] = 1
            else
                continue # ray[i] == 0
            end

            if idx <= m
                ray[i] *= GLPK.get_rii(lp, k)
            else
                ray[i] /= GLPK.get_sjj(lp, k)
            end
        end

        GLPK.btran(lp, ray)

        for i = 1:m
            ray[i] /= GLPK.get_rii(lp, i)
        end
    end

    return ray
end

function MPB.getunboundedray(lpm::GLPKMathProgModelLP)
    check_feasible_bounds(lpm, "getreducedcosts")

    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
    elseif lpm.method == :InteriorPoint
        error("getunboundedray is not available when using the InteriorPoint method")
    else
        error("bug")
    end

    m = GLPK.get_num_rows(lp)
    n = GLPK.get_num_cols(lp)

    ray = zeros(n)

    ur = GLPK.get_unbnd_ray(lp)
    if ur != 0
        if ur <= m
            k = ur
            get_stat = GLPK.get_row_stat
            get_dual = GLPK.get_row_dual
        else
            k = ur - m
            get_stat = GLPK.get_col_stat
            get_dual = GLPK.get_col_dual
            ray[k] = 1
        end

        get_stat(lp, k) != GLPK.BS || error("unbounded ray is dual (use getinfeasibilityray)")

        for (ri, rv) in zip(GLPK.eval_tab_col(lp, ur)...)
            ri > m && (ray[ri - m] = rv)
        end

        if (GLPK.get_obj_dir(lp) == GLPK.MAX) ⊻ (get_dual(lp, k) > 0)
            ray .*= -1.0
        end
    else
        for i = 1:n
            ray[i] = GLPK.get_col_prim(lp, i)
        end
    end

    return ray
end

end

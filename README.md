GLPKMathProgInterface.jl
========================

Interface between the [GLPK.jl] wrapper and [MathProgBase.jl].
With this package, you can use the GLPK solver in MathProgBase.jl

Testing status: [![Build Status](https://api.travis-ci.org/JuliaOpt/GLPKMathProgInterface.jl.png)](https://travis-ci.org/JuliaOpt/GLPKMathProgInterface.jl)

The linear programming solver accepts the following optional
keywords:

 * `method`: accepts `:Simplex` (default), `:Exact`, `:InteriorPoint`
 * `presolve`: accepts `false` (default) or `true`. If set to `true`,
    status reporting will always report `:Undefined` in case of
    presolver failure (i.e. also for infeasible problems).

Example:

```Julia
using MathProgBase
using GLPKMathProgInterface
sol = linprog([-1,0], [2 1], '<', 1.5, GLPKSolverLP(method=:Exact, presolve=true)
```

The mixed integer programming solver accepts the following optional
keywords:

 * `presolve`: accepts `false` (default) or `true`. If set to `true`,
    status reporting will always report `:Undefined` in case of
    presolver failure (i.e. also for infeasible problems). If set
    to `false`, presolving will be performed via the simplex method.

Example:

```Julia
using MathProgBase
using GLPKMathProgInterface
sol = mixintprog(-[5,3,3], [1 8 2], '<', 9, 'I', 0, 1, GLPKSolverMIP(presolve=true))
```

[GLPK.jl]: https://github.com/JuliaOpt/GLPK.jl
[MathProgBase.jl]: https://github.com/JuliaOpt/MathProgBase.jl

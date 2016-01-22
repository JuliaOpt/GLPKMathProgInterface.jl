GLPKMathProgInterface.jl
========================

[![Build Status](https://travis-ci.org/JuliaOpt/GLPKMathProgInterface.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/GLPKMathProgInterface.jl)
[![GLPKMathProgInterface](http://pkg.julialang.org/badges/GLPKMathProgInterface_0.3.svg)](http://pkg.julialang.org/?pkg=GLPKMathProgInterface&ver=0.3)
[![GLPKMathProgInterface](http://pkg.julialang.org/badges/GLPKMathProgInterface_0.4.svg)](http://pkg.julialang.org/?pkg=GLPKMathProgInterface&ver=0.4)

Interface between the [GLPK.jl] wrapper and [MathProgBase.jl].
With this package, you can use the GLPK solver in MathProgBase.jl

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
sol = linprog([-1,0], [2 1], '<', 1.5, GLPKSolverLP(method=:Exact, presolve=true))
```

Additionally, all low-level parameters which can be passed to the GLPK problem objects, such as
`msg_lev`, `tol_bnd` etc., can also be passed as keyword arguments, for finer control. Refer to the
[GLPK.jl documentation] for more details on the available parameters. Invalid option
names will be ignored with a warning (note that the set of valid options is different for the
`:Simplex` and `:InteriorPoint` methods). This is a low-level feature, therefore options with invalid
arguments may cause the solver to crash.

For example, to get a more verbose output in the above case, write:

```Julia
sol = linprog([-1,0], [2 1], '<', 1.5, GLPKSolverLP(method=:Exact, presolve=true, msg_lev=GLPK.MSG_ON)
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
sol = mixintprog(-[5,3,3], [1 8 2], '<', 9, :Int, 0, 1, GLPKSolverMIP(presolve=true))
```

Additionally, as in the case of linear programming solver (see above), all low-level
parameters which can be passed to the GLPK problem objects can also be passed as keyword
arguments to the `GLPKSolverMIP` constructor for finer control. The constructor accepts
arguments both for the MIP solver and the LP simplex solver (which may be used by the MIP
solver internally), so any field of an `IntoptParam` or of a `SimplexParam` can be passed;
in case of common options (e.g. `msg_lev`) the parameters are set for both solvers.
The only exceptions are the callback pointers `cb_func` and `cb_info`, which can not
be passed in this way, since the callback interface is only accessible via MathProgBase.
As in the LP case, invalid options are ignored with a warning.
This is a low level feature, and passing invalid parameter arguments may crash the solver.

[GLPK.jl]: https://github.com/JuliaOpt/GLPK.jl
[MathProgBase.jl]: https://github.com/JuliaOpt/MathProgBase.jl
[GLPK.jl documentation]: https://gplkjl.readthedocs.org/en/latest/glpk.html

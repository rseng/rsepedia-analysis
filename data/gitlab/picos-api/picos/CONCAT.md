Changelog
=========

This file documents major changes to PICOS. The format is based on
`Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_.

.. _2.4: https://gitlab.com/picos-api/picos/compare/v2.3...v2.4
.. _2.3: https://gitlab.com/picos-api/picos/compare/v2.2...v2.3
.. _2.2: https://gitlab.com/picos-api/picos/compare/v2.1...v2.2
.. _2.1: https://gitlab.com/picos-api/picos/compare/v2.0...v2.1
.. _2.0: https://gitlab.com/picos-api/picos/compare/v1.2.0...v2.0
.. _1.2.0: https://gitlab.com/picos-api/picos/compare/v1.1.3...v1.2.0
.. _1.1.3: https://gitlab.com/picos-api/picos/compare/v1.1.2...v1.1.3
.. _1.1.2: https://gitlab.com/picos-api/picos/compare/v1.1.1...v1.1.2
.. _1.1.1: https://gitlab.com/picos-api/picos/compare/v1.1.0...v1.1.1
.. _1.1.0: https://gitlab.com/picos-api/picos/compare/v1.0.2...v1.1.0
.. _1.0.2: https://gitlab.com/picos-api/picos/compare/v1.0.1...v1.0.2
.. _1.0.1: https://gitlab.com/picos-api/picos/compare/v1.0.0...v1.0.1
.. _1.0.0: https://gitlab.com/picos-api/picos/compare/b65a05be...v1.0.0
.. _0.1.3: about:blank
.. _0.1.2: about:blank
.. _0.1.1: about:blank
.. _0.1.0: about:blank


`2.4`_ - 2022-02-12
--------------------------------------------------------------------------------

*The performance update.*

.. rubric:: Added

- Support for noncovex quadratic constraints with Gurobi 9 (or later).
- Setting :data:`UNRELIABLE_STRATEGIES <picos.settings.UNRELIABLE_STRATEGIES>`
  to enable passing of problems to solvers that nominally support them but have
  proven unreliable.
- Setting :data:`PREFER_GUROBI_MATRIX_INTERFACE
  <picos.settings.PREFER_GUROBI_MATRIX_INTERFACE>` and option
  :ref:`gurobi_matint <option_gurobi_matint>` to toggle between Gurobi's legacy
  and matrix interface.
- Option :ref:`mosek_basic_sol <option_mosek_basic_sol>` to let MOSEK
  (Optimizer) compute a basic solution for LPs.

.. rubric:: Changed

- The performance for solving problems with large data has been improved

  - drastically for CVXOPT and MOSEK (Optimizer; LPs in particular),
  - significantly for Cplex and SCIP, and
  - subtly for GLPK, Gurobi and ECOS.

  This is most notable for LPs with a dense constraint matrix where the overhead
  for data passing can be significant in relation to the search time.

- The performance of :func:`picos.sum` when summing a large number of
  (bi-)affine expressions has been improved drastically.
- When possible, Gurobi is now interfaced through its matrix interface, which is
  faster for large data. This requires Gurobi 9 (or later) and SciPy.
- By default, solving with MOSEK (Optimizer) does not return a basic LP solution
  any more. Use :ref:`mosek_basic_sol <option_mosek_basic_sol>` to control this.
- The default value of :ref:`cvxopt_kktsolver <option_cvxopt_kktsolver>` is now
  :obj:`None` and means "try the fast ``"chol"`` first and fall back to the
  reliable ``"ldl"`` on error".
- Dualization now makes use of variable bounds to reduce the number of auxiliary
  constraints.
- The Python interface used to communicate with a solver is now mentioned in
  various log messages and exceptions.

.. rubric:: Fixed

- On-the-fly loading of a data vector in a multiplication with a matrix
  expression.
- Maximization of a squared norm not being detected as a nonconvex quadratic
  objective and being passed to solvers that do not support it.


`2.3`_ - 2021-10-07
--------------------------------------------------------------------------------

*The syntactic sugar update.*

.. rubric:: Important

- When forming linear matrix inequalities with the ``<<`` or ``>>`` operator,
  if one operand is an :math:`n \times n` matrix and the other is an
  :math:`n`-dimensional vector (or a scalar), the latter is now understood as
  (respectively broadcasted along) the main diagonal of an :math:`n \times n`
  diagonal matrix. In particular ``X >> 1`` is now understood as :math:`X
  \succeq I` as opposed to :math:`X \succeq J`. If you want to express a
  constraint :math:`X \succeq \alpha J` where :math:`J` is a matrix of all ones,
  use the new :func:`picos.J`.

.. rubric:: Added

- Support for the OSQP solver.
- On-the-fly loading of :mod:`scipy.sparse` matrices. (See new note
  :ref:`numscipy`.)
- Ability to negate or scale any expression and to sum any two expressions with
  the same or with a different type. This is established through a new
  :class:`~picos.expressions.exp_wsum.WeightedSum` fallback class. Convex or
  concave weighted sums can be used as an objective or in a constraint like any
  other expression.
- Properties :attr:`~picos.valuable.Valuable.sp`,
  :attr:`~picos.valuable.Valuable.np` and :attr:`~picos.valuable.Valuable.np2d`
  to query the value of an expression as a SciPy or NumPy type. (See new class
  :class:`~picos.valuable.Valuable` for all value query options.)
- Ability to use :func:`numpy.array` directly on valued PICOS objects, returning
  a zero, one or two-dimensional array depending on the shape of the value.
- New method :meth:`~picos.modeling.problem.Problem.require` and an equivalent
  overload for ``+=`` to add constraints to a
  :meth:`~picos.modeling.problem.Problem`.
- Cached functions :func:`~picos.I`, :func:`~picos.J`, and :func:`~picos.O` that
  create, respectively, an identity matrix, a matrix of all ones, and a zero
  matrix.
- Cached properties :attr:`BiaffineExpression.rowsum
  <picos.expressions.exp_biaffine.BiaffineExpression.rowsum>` and
  :attr:`~picos.expressions.exp_biaffine.BiaffineExpression.colsum` to
  complement the existing property
  :attr:`~picos.expressions.exp_biaffine.BiaffineExpression.sum` and an argument
  ``axis`` to :func:`picos.sum` for the same purpose.
- Option to give a name to :class:`problems <picos.modeling.problem.Problem>`
  via the first initialization argument or the
  :attr:`~picos.modeling.problem.Problem.name` property.
- Ability to perform some algebraic operations on :class:`objectives
  <picos.modeling.objective.Objective>`.
- Support for solving nonconvex continuous
  quadratic programs (QPs) with CPLEX and Gurobi. Gurobi further allows convex
  quadratic constraints to be present.
- Ability to
  :meth:`reshape <picos.expressions.exp_biaffine.BiaffineExpression.reshaped>`
  affine expressions in C-order, like NumPy.
- Ability to pass constant values to :func:`picos.sum`, :func:`~picos.min` and
  :func:`~picos.max`.
- Global option :data:`settings.RETURN_SOLUTION
  <picos.settings.RETURN_SOLUTION>` that controls whether
  :meth:`~picos.modeling.problem.Problem.solve` returns a
  :class:`~picos.modeling.solution.Solution`.
- Methods :class:`Samples.shuffled <picos.expressions.samples.Samples.shuffled>`
  and :class:`~picos.expressions.samples.Samples.kfold`.
- Support for MOSEK remote optimization with the :ref:`mosek_server
  <option_mosek_server>` option.
- Option :ref:`cplex_vmconfig <option_cplex_vmconfig>` to load a virtual machine
  configuration file with CPLEX.
- Function :func:`picos.patch_scipy_array_priority` to work around `SciPy#4819
  <https://github.com/scipy/scipy/issues/4819>`__.

.. rubric:: Changed

- The performance of solving semidefinite programs with trivial linear matrix
  inequalities of the form ``X >> 0`` using MOSEK (Optimizer) has been improved
  dramatically. Depending on your problem, you might experience this speedup
  when using the :ref:`dualize <option_dualize>` option.
- :attr:`Problem.minimize <picos.modeling.problem.Problem.minimize>` and
  :attr:`Problem.maximize <picos.modeling.problem.Problem.maximize>` are now
  properties that you can assign a minimization or maximization objective to,
  respectively.
- All expression types as well as the classes
  :class:`~picos.modeling.problem.Problem` and
  :class:`~picos.modeling.objective.Objective` now share the same interface to
  query their (objective) value. In particular, the new
  :attr:`~picos.valuable.Valuable.np` property can be used on all.
- Solving with ``duals=True`` will now raise an exception when duals were
  returned by the solver but not all could be converted. Use the default of
  ``duals=None`` to accept also incomplete duals.
- The new argument ``name`` is the only optional argument to
  :class:`~picos.modeling.problem.Problem` that may be passed as a positional
  argument; the arguments ``copyOptions`` and ``useOptions`` must now be passed
  as keyword arguments.

.. rubric:: Fixed

- Running ``setup.py`` under Python 3.6 and earlier.
- Bad shebang lines; all are now properly reading ``#!/usr/bin/env python3``.
- Incorrect duals returned by MOSEK (Fusion).
- An assertion failure when multiplying some quadratic expressions with a
  negative scalar.
- A false expression being created when multiplying a
  :class:`~picos.expressions.exp_detrootn.DetRootN` with a negative scalar.
- An exception when multiplying a scalar power with a constant.
- A modify-during-iteration issue that could result in a suboptimal solver being
  chosen.
- Building piecewise affine functions from a mix of certain and random
  expressions.
- A failure when computing the convex hull of a
  :class:`ScenarioPerturbationSet <picos.uncertain.ScenarioPerturbationSet>`
  with few points.
- Detection of string groups where the variable part is at the start or end of
  the strings.
- CVXOPT reacting inconsistently to some infeasible problems.
- A potential variable clash when reformulating a
  :class:`~picos.constraints.con_matnorm.NuclearNormConstraint`.
- Grammatical issues when printing variable groups of a problem.

.. rubric:: Removed

- The deprecated functions :attr:`Problem.minimize
  <picos.modeling.problem.Problem.minimize>` and
  :attr:`Problem.maximize <picos.modeling.problem.Problem.maximize>`. See
  **Changed** for the new meaning of these names.
- The deprecated arguments ``it`` and ``indices`` to :func:`picos.sum`.


`2.2`_ - 2021-02-09
--------------------------------------------------------------------------------

*The Python 3 update.*

.. rubric:: Important

- PICOS now requires Python 3.4 or later; Python 2 support was dropped.

.. rubric:: Added

- A synopsis to the :exc:`NoStrategyFound <.strategy.NoStrategyFound>`
  exception, explaining why strategy search failed.

.. rubric:: Fixed

- Optimizing matrix :math:`(p,q)`-norms when columns of the matrix are constant.
- Refining norms over a sparse constant term to a constant affine expression.
- Gurobi printing empty lines to console when dual retrieval fails.

.. rubric:: Changed

- A bunch of Python 2 compatibility code was finally removed.
- Exception readability has been improved using Python 3's ``raise from`` syntax
  where applicable.
- The ``__version_info__`` field now contains integers instead of strings.
- :attr:`QuadraticExpression.scalar_factors
  <.exp_quadratic.QuadraticExpression.scalar_factors>` is now :obj:`None`
  instead of an empty tuple when no decomposition into scalar factors is known.

.. rubric:: Deprecated

- :attr:`QuadraticExpression.quadratic_forms
  <.exp_quadratic.QuadraticExpression.quadratic_forms>`, as write access would
  leave the expression in an inconsistent state. (At your own risk, use the
  equivalent ``_sparse_quads`` instead.)


`2.1`_ - 2020-12-29
--------------------------------------------------------------------------------

*The robust optimization update.*

.. rubric:: Important

- The sign of dual values for affine equality constraints has been fixed by
  inversion.

.. rubric:: Added

- Support for a selection of robust optimization (RO) and distributionally
  robust stochastic programming (DRO) models through a new
  :mod:`picos.uncertain` namespace. You may now solve

  - scenario-robust conic programs via :class:`ScenarioPerturbationSet
    <picos.uncertain.ScenarioPerturbationSet>`,
  - conically robust linear programs and robust conic quadratic programs under
    ellipsoidal uncertainty via :class:`ConicPerturbationSet
    <picos.uncertain.ConicPerturbationSet>` and :class:`UnitBallPerturbationSet
    <picos.uncertain.UnitBallPerturbationSet>`, and
  - least squares and piecewise linear stochastic programs where the data
    generating distribution is defined ambiguously through a Wasserstein ball or
    through bounds on its first two moments via :class:`WassersteinAmbiguitySet
    <picos.uncertain.WassersteinAmbiguitySet>` and :class:`MomentAmbiguitySet
    <picos.uncertain.MomentAmbiguitySet>`, respectively.

- New function :func:`picos.block` to create block matrices efficiently.
- New convenience class :class:`picos.Samples` for data-driven applications.
- New set class :class:`picos.Ellipsoid` (has overlap with but a different
  scope than :class:`picos.Ball`).
- Support for :meth:`matrix reshuffling
  <picos.expressions.exp_biaffine.BiaffineExpression.reshuffled>` (aka *matrix
  realignment*) used in quantum information theory.
- Ability to define cones of fixed dimensionality and :class:`product cones
  <picos.ProductCone>` thereof.
- Ability to query the :attr:`solver-reported objective value
  <.solution.Solution.reported_value>` (useful with RO and DRO objectives).
- Methods :meth:`Problem.conic_form <.problem.Problem.conic_form>` and
  :meth:`reformulated <.problem.Problem.reformulated>` for internal use and
  educational purposes.
- New module :mod:`picos.settings` defining global options that can be set
  through environment variables prefixed with ``PICOS_``. Among other things,
  you can now blacklist all proprietary solvers for an application by passing
  ``PICOS_NONFREE_SOLVERS=False`` to the Python interpreter.
- A new base class :class:`BiaffineExpression
  <.exp_biaffine.BiaffineExpression>` for all (uncertain) affine expression
  types. This gives developers extending PICOS a framework to support models
  with parameterized data.
- Support for :meth:`factoring out
  <.exp_biaffine.BiaffineExpression.factor_out>` variables and parameters
  from (bi)affine vector expression.
- Support for :meth:`replacing <.expression.Expression.replace_mutables>`
  variables and parameters with affine expressions of same shape to perform a
  change of variables in a mathematical sense.
- Support for SCIP Optimization Suite 7.
- CVXOPT-specific solution search options
  :ref:`cvxopt_kktsolver <option_cvxopt_kktsolver>` and :ref:`cvxopt_kktreg
  <option_cvxopt_kktreg>`.

.. rubric:: Fixed

- Quadratic expressions created from a squared norm failing to decompose due to
  a numerically singular quadratic form.
- Solution objects unintendedly sharing memory.
- Solution search options that take a dictionary as their argument.
- Solution search with :ref:`assume_conic <option_assume_conic>` set to
  :obj:`False`.
- The :class:`EpigraphReformulation <picos.reforms.EpigraphReformulation>`
  falsely claiming that it can reformulate any nonconvex objective.
- A division by zero that could occur when computing the solution search
  overhead.
- An exception with functions that look for short string descriptions, in
  particular with :meth:`picos.sum`.

.. rubric:: Changed

- The functions :func:`picos.max` and :func:`picos.min` can now be used to
  express the maximum over a list of convex and the minimum over a list of
  concave expressions, respectively.
- Squared norms are now implemented as a subclass of quadratic expressions
  (:class:`SquaredNorm <picos.SquaredNorm>`), skipping an unnecessary
  decomposition on constraint creation.
- Commutation matrices used internally for various algebraic tasks are now
  retrieved from a centralized cached function, improving performance.
- The string description of :class:`Problem <.problem.Problem>` instances is not
  enclosed by dashed lines any more.


`2.0`_ - 2020-03-03
--------------------------------------------------------------------------------

*The backend update.*

.. rubric:: Important

This is a major release featuring vast backend rewrites as well as interface
changes. Programs written for older versions of PICOS are expected to raise
deprecation warnings but should otherwise work as before. The following lists
notable exceptions:

- The solution returned by :meth:`~.problem.Problem.solve` is now an instance of
  the new :class:`~picos.Solution` class instead of a dictionary.
- If solution search fails to find an optimal primal solution, PICOS will now
  raise a :class:`~picos.SolutionFailure` by default. Old behavior of not
  raising an exception is achieved by setting ``primals=None`` (see
  :ref:`primals <option_primals>` and :ref:`duals <option_duals>` options).
- The definition of the :math:`L_{p,q}`-norm has changed: It no longer refers
  to the :math:`p`-norm of the :math:`q`-norms of the matrix rows but to the
  :math:`q`-norm of the :math:`p`-norms of the matrix columns. This matches
  the definition you would find `on
  Wikipedia <https://en.wikipedia.org/wiki/Matrix_norm#L2,1_and_Lp,q_norms>`_
  and should reduce confusion for new users. See :class:`~picos.Norm`.
- The signs in the Lagrange dual problem of a conic problem are now more
  consistent for all cones, see :ref:`duals`. In particular the signs of dual
  values for (rotated) second order conic constraints have changed and the
  problem obtained by :attr:`Problem.dual <.problem.Problem.dual>` (new for
  :meth:`~.problem.Problem.as_dual`) has a different (but equivalent) form.

.. rubric:: Added

- A modular problem reformulation framework. Before selecting a solver, PICOS
  now builds a map of problem types that your problem can be reformulated to
  and makes a choice based on the expected complexity of the reposed problem.
- An object oriented interface to solution search options. See
  :class:`~picos.Options`.
- Support for arbitrary objective functions via an epigraph reformulation.
- Support for MOSEK 9.
- Support for ECOS 2.0.7.
- Support for multiple subsystems with :func:`~picos.partial_trace`.
- Quick-solve functions :func:`picos.minimize` and :func:`picos.maximize`.
- Lower and upper diagonal matrix variable types.
- :class:`~picos.SecondOrderCone` and :class:`~picos.RotatedSecondOrderCone`
  sets to explicitly create the associated constraints. *(You now need to use
  these if you want to obtain a conic instead of a quadratic dual.)*
- Possibility to use :func:`picos.sum` to sum over the elements of a single
  multidimensional expression.
- Possibility to create a :class:`~picos.Ball` or :class:`~picos.Simplex` with a
  non-constant radius.
- Many new properties (postfix operations) to work with affine expressions; for
  instance ``A.vec`` is a faster and cached way to express the vectorization
  ``A[:]``.
- Options :ref:`assume_conic <option_assume_conic>` and
  :ref:`verify_prediction <option_verify_prediction>`.
- An option for every solver to manipulate the chances of it being selected
  (e.g. :ref:`penalty_cvxopt <option_penalty_cvxopt>`).
- Ability to run doctests via ``test.py``.

.. rubric:: Fixed

The following are issues that were fixed in an effort of their own. If a bug is
not listed here, it might still be fixed as a side effect of some of the large
scale code rewrites that this release ships.

- Upgrading the PyPI package via pip.
- A regression that rendered the Kronecker product unusable.
- Noisy exception handling in a sparse matrix helper function.
- Shape detection for matrices given by string.
- The :ref:`hotstart <option_hotstart>` option when solving with CPLEX.
- Low precision QCP duals from Gurobi.

.. rubric:: Changed

- All algebraic expression code has been rewritten and organized in a new
  :mod:`~picos.expressions` package. In particular, real and complex expressions
  are distinguished more clearly.
- All algebraic expressions are now immutable.
- The result of any unary operation on algebraic expressions (e.g. negation,
  transposition) is cached (only computed once per expression).
- Slicing of affine expressions is more powerful, see :ref:`slicing`.
- Loading of constant numeric data has been unified, see
  :func:`~picos.expressions.data.load_data`.
- Variables are now created independently of problems by instanciating one of
  the new :mod:`variable types <picos.expressions.variables>`.
  *(*:meth:`Problem.add_variable <.problem.Problem.add_variable>`
  *is deprecated.)*
- Constraints are added to problems as they are; any transformation is done
  transparently during solution search.
- In particular, :math:`x^2 \leq yz` is now initially a (nonconvex) quadratic
  constraint and transformation to a conic constraint is controlled by the new
  :ref:`assume_conic <option_assume_conic>` option.
- Expressions constrained to be positive semidefinite are now required to be
  symmetric/hermitian by their own definition. *(Use*
  :class:`~picos.SymmetricVariable` *or* :class:`~picos.HermitianVariable`
  *whenever applicable!)*
- Options passed to :meth:`~.problem.Problem.solve` are only used for that
  particular search.
- The default value for the :ref:`verbosity <option_verbosity>` option (formerly
  ``verbose``) is now :math:`0`.
- Available solvers are only imported when they are actually being used, which
  speeds up import of PICOS on platforms with many solvers installed.
- The code obeys PEP 8 and PEP 257 more strongly. Exceptions: D105, D203, D213,
  D401, E122, E128, E221, E271, E272, E501, E702, E741.
- Production testing code was moved out of the :mod:`picos` package.

.. rubric:: Removed

- The ``NoAppropriateSolverError`` exception that was previously raised by
  :meth:`~.problem.Problem.solve`. This is replaced by the new
  :class:`~picos.SolutionFailure` exception with error code :math:`1`.
- Some public functions in the :mod:`~picos.tools` module that were originally
  meant for internal use.

.. rubric:: Deprecated

This section lists deprecated modules, functions and options with their
respective replacement or deprecation reason on the right hand side.
Deprecated entities produce a warning and will be removed in a future release.

- The :mod:`~picos.tools` module as a whole. It previously contained both
  algebraic functions for the user as well as functions meant for internal use.
  The former group of functions can now be imported directly from the
  :mod:`picos` namespace (though some are also individually deprecated). The
  other functions were either relocated (but can still be imported from
  :mod:`~picos.tools` while it lasts) or removed.
- In the :class:`~.problem.Problem` class:

  - :meth:`~.problem.Problem.add_variable`,
    :meth:`~.problem.Problem.remove_variable`,
    :meth:`~.problem.Problem.set_var_value`
    → variables are instanciated directly and added to problems automatically
  - :meth:`~.problem.Problem.minimize` → :func:`picos.minimize`
  - :meth:`~.problem.Problem.maximize` → :func:`picos.maximize`
  - :meth:`~.problem.Problem.set_option`
    → assign to attributes or items of :attr:`Problem.options <picos.Options>`
  - :meth:`~.problem.Problem.update_options`
    → :meth:`options.update <.options.Options.update>`
  - :meth:`~.problem.Problem.set_all_options_to_default`
    → :meth:`options.reset <.options.Options.reset>`
  - :meth:`~.problem.Problem.obj_value` → :attr:`~.valuable.Valuable.value`
  - :meth:`~.problem.Problem.is_continuous`
    → :attr:`~.problem.Problem.continuous`
  - :meth:`~.problem.Problem.is_pure_integer`
    → :attr:`~.problem.Problem.pure_integer`
  - :meth:`~.problem.Problem.verbosity`
    → :ref:`options.verbosity <option_verbosity>`
  - :meth:`~.problem.Problem.as_dual` → :attr:`~.problem.Problem.dual`
  - :meth:`~.problem.Problem.countVar`,
    :meth:`~.problem.Problem.countCons`,
    :meth:`~.problem.Problem.numberOfVars`,
    :meth:`~.problem.Problem.numberLSEConstraints`,
    :meth:`~.problem.Problem.numberSDPConstraints`,
    :meth:`~.problem.Problem.numberQuadConstraints`,
    :meth:`~.problem.Problem.numberConeConstraints`
    → were meant for internal use
  - arguments ``it``, ``indices`` and ``key`` to
    :meth:`~.problem.Problem.add_list_of_constraints` → are ignored

- All expression types:

  - constraint creation via ``<`` → ``<=``
  - constraint creation via ``>`` → ``>=``
  - :meth:`~.expression.Expression.is_valued`
    → :attr:`~.valuable.Valuable.valued`
  - :meth:`~.expression.Expression.set_value`
    → assign to :attr:`~.valuable.Valuable.value`

- Affine expressions:

  - :meth:`~.exp_biaffine.BiaffineExpression.fromScalar`
    → :meth:`~.exp_biaffine.BiaffineExpression.from_constant`
    or :func:`picos.Constant`
  - :meth:`~.exp_biaffine.BiaffineExpression.fromMatrix`
    → :meth:`~.exp_biaffine.BiaffineExpression.from_constant`
    or :func:`picos.Constant`
  - :meth:`~.exp_biaffine.BiaffineExpression.hadamard` → ``^``
  - :meth:`~.exp_biaffine.BiaffineExpression.isconstant`
    → :meth:`~.expression.Expression.constant`
  - :meth:`~.exp_biaffine.BiaffineExpression.same_as`
    → :meth:`~.exp_biaffine.BiaffineExpression.equals`
  - :meth:`~.exp_biaffine.BiaffineExpression.transpose`
    → :attr:`~.exp_biaffine.BiaffineExpression.T`
  - :attr:`~.exp_biaffine.BiaffineExpression.Tx`
    → :meth:`~.exp_biaffine.BiaffineExpression.partial_transpose`
  - :meth:`~.exp_biaffine.BiaffineExpression.conjugate`
    → :attr:`~.exp_biaffine.BiaffineExpression.conj`
  - :meth:`~.exp_biaffine.BiaffineExpression.Htranspose`
    → :attr:`~.exp_biaffine.BiaffineExpression.H`
  - :meth:`~.exp_biaffine.BiaffineExpression.copy`
    → expressions are immutable
  - :meth:`~.exp_biaffine.BiaffineExpression.soft_copy`
    → expressions are immutable

- Algebraic functions and shorthands in the ``picos`` namespace:

  - :func:`~picos.tracepow` → :class:`~picos.PowerTrace`
  - :func:`~picos.new_param` → :func:`~picos.Constant`
  - :func:`~picos.flow_Constraint` → :class:`~picos.FlowConstraint`
  - :func:`~picos.diag_vect` → :func:`~picos.maindiag`
  - :func:`~picos.simplex` → :class:`~picos.Simplex`
  - :func:`~picos.truncated_simplex` → :class:`~picos.Simplex`
  - arguments ``it`` and ``indices`` to :func:`~picos.sum` → are ignored

- Solution search options:

  - ``allow_license_warnings``
    → :ref:`license_warnings <option_license_warnings>`
  - ``verbose`` → :ref:`verbosity <option_verbosity>` (takes an integer)
  - ``noprimals`` → :ref:`primals <option_primals>` (the meaning is inverted)
  - ``noduals`` → :ref:`duals <option_duals>` (the meaning is inverted)
  - ``tol`` →  ``*_fsb_tol`` and ``*_ipm_opt_tol``
  - ``gaplim`` → :ref:`rel_bnb_opt_tol <option_rel_bnb_opt_tol>`
  - ``maxit`` → :ref:`max_iterations <option_max_iterations>`
  - ``nbsol`` → :ref:`max_fsb_nodes <option_max_fsb_nodes>`
  - ``pool_relgap`` → :ref:`pool_rel_gap <option_pool_rel_gap>`
  - ``pool_absgap`` → :ref:`pool_abs_gap <option_pool_abs_gap>`
  - ``lboundlimit`` → :ref:`cplex_lwr_bnd_limit <option_cplex_lwr_bnd_limit>`
  - ``uboundlimit`` → :ref:`cplex_upr_bnd_limit <option_cplex_upr_bnd_limit>`
  - ``boundMonitor`` → :ref:`cplex_bnd_monitor <option_cplex_bnd_monitor>`
  - ``solve_via_dual`` → :ref:`dualize <option_dualize>` (may not be :obj:`None`
    any more)


`1.2.0`_ - 2019-01-11
--------------------------------------------------------------------------------

.. rubric:: Important

- :attr:`A scalar expression's value <.valuable.Valuable.value>` and
  :attr:`a scalar constraint's dual <.constraint.Constraint.dual>` are returned
  as scalar types as opposed to 1×1 matrices.
- The dual value returned for rotated second order cone constraints is now a
  proper member of the dual cone (which equals the primal cone up to a factor of
  :math:`4`). Previously, the dual of an equivalent second order cone constraint
  was returned.
- The Python 2/3 compatibility library ``six`` is no longer a dependency.

.. rubric:: Added

- Support for the ECOS solver.
- Experimental support for MOSEK's new Fusion API.
- Full support for exponential cone programming.
- A production testing framework featuring around 40 novel optimization test
  cases that allows quick selection of tests, solvers, and solver options.
- A "glyph" system that allows the user to adjust the string representations of
  future expressions and constraints. For instance, :func:`picos.latin1()
  <picos.glyphs.latin1>` disables use of unicode symbols.
- Support for symmetric variables with all solvers, even if they do not support
  semidefinite programming.

.. rubric:: Changed

- Solver implementations each have a source file of their own, and a common
  interface that makes implementing new solvers easier.
- Likewise, constraint implementations each have a source file of their own.
- The implementations of CPLEX, Gurobi, MOSEK and SCIP have been rewritten.
- Solver selection takes into account how well a problem is supported,
  distinguishing between native, secondary, experimental and limited support.
- Unsupported operations on expressions now produce meaningful exceptions.
- :meth:`add_constraint <.problem.Problem.add_constraint>` and
  :meth:`add_list_of_constraints <.problem.Problem.add_list_of_constraints>`
  always return the constraints
  passed to them.
- :meth:`add_list_of_constraints <.problem.Problem.add_list_of_constraints>`
  and :func:`picos.sum` find a short string representation automatically.

.. rubric:: Removed

- The old production testing script.
- Support for the SDPA solver.
- Support for sequential quadratic programming.
- The options ``convert_quad_to_socp_if_needed``, ``pass_simple_cons_as_bound``,
  ``return_constraints``, ``handleBarVars``, ``handleConeVars`` and
  ``smcp_feas``.
- Support for GLPK and MOSEK through CVXOPT.

.. rubric:: Fixed

- Performance issues when exporting variable bounds to CVXOPT.
- Hadamard product involving complex matrices.
- Adding constant terms to quadratic expression.
- Incorrect or redundant expression string representations.
- GLPK handling of the default ``maxit`` option.
- Miscellaneous solver-specific bugs in the solvers that were re-implemented.


`1.1.3`_ - 2018-10-05
--------------------------------------------------------------------------------

.. rubric:: Added

- Support for the solvers GLPK and SCIP.
- PICOS packages `on Anaconda Cloud <https://anaconda.org/picos/picos>`_.
- PICOS packages `in the Arch Linux User Repository
  <https://aur.archlinux.org/packages/?SeB=b&K=python-picos>`_.

.. rubric:: Changed

- The main repository has moved to
  `GitLab <https://gitlab.com/picos-api/picos>`_.
- Releases of packages and documentation changes are
  `automated <https://about.gitlab.com/features/gitlab-ci-cd/>`_ and thus more
  frequent. In particular, post release versions are available.
- Test bench execution is automated for greater code stability.
- Improved test bench output.
- Improved support for the SDPA solver.
- :func:`~picos.partial_trace` can handle rectangular subsystems.
- The documentation was restructured; examples were converted to Python 3.

.. rubric:: Fixed

- Upper bounding the norm of a complex scalar.
- Multiplication with a complex scalar.
- A couple of Python 3 specific errors, in particular when deleting constraints.
- All documentation examples are reproducible with the current state of PICOS.


`1.1.2`_ - 2016-07-04
--------------------------------------------------------------------------------

.. rubric:: Added

- Ability to dynamically add and remove constraints.
- Option ``pass_simple_cons_as_bound``, see below.

.. rubric:: Changed

- Improved efficiency when processing large expressions.
- Improved support for the SDPA solver.
- :meth:`add_constraint <.problem.Problem.add_constraint>` returns a handle to
  the constraint when the option `return_constraints` is set.
- New signature for the function :func:`~picos.partial_transpose`, which can now
  transpose arbitrary subsystems from a kronecker product.
- PICOS no longer turns constraints into variable bounds, unless the new option
  ``pass_simple_cons_as_bound`` is enabled.

.. rubric:: Fixed

- Minor bugs with complex expressions.


`1.1.1`_ - 2015-08-29
--------------------------------------------------------------------------------

.. rubric:: Added

- Support for the SDPA solver.
- Partial trace of an affine expression, see :func:`~picos.partial_trace`.

.. rubric:: Changed

- Improved PEP 8 compliance.

.. rubric:: Fixed

- Compatibility with Python 3.


`1.1.0`_ - 2015-04-15
--------------------------------------------------------------------------------

.. rubric:: Added

- Compatibility with Python 3.

.. rubric:: Changed

- The main repository has moved to `GitHub <https://github.com/gsagnol/picos>`_.


`1.0.2`_ - 2015-01-30
--------------------------------------------------------------------------------

.. rubric:: Added

- Ability to read and write problems in
  `conic benchmark format <http://cblib.zib.de/>`_.
- Support for inequalities involving the sum of the :math:`k` largest or
  smallest elements of an affine expression, see :func:`~picos.sum_k_largest`
  and :func:`~picos.sum_k_smallest`.
- Support for inequalities involving the sum of the :math:`k` largest or
  smallest eigenvalues of a symmetric matrix, see
  :func:`~picos.sum_k_largest_lambda`, :func:`~picos.sum_k_smallest_lambda`,
  :func:`~picos.lambda_max` and :func:`~picos.lambda_min`.
- Support for inequalities involving the :math:`L_{p,q}`-norm of an affine
  expression, see :func:`~picos.norm`.
- Support for equalities involving complex coefficients.
- Support for antisymmetric matrix variables.
- Set expressions that affine expressions can be constrained to be an element
  of, see :func:`~picos.ball`, :func:`~picos.simplex` and
  :func:`~picos.truncated_simplex`.
- Shorthand functions :meth:`maximize <.problem.Problem.maximize>` and
  :meth:`minimize <.problem.Problem.minimize>` to specify the objective function
  of a problem and solve it.
- Hadamard (elementwise) product of affine expression, as an overload of the
  ``^`` operator, read :ref:`the tutorial on overloads <overloads>`.
- Partial transposition of an aAffine Expression, see
  :func:`~picos.partial_transpose`.

.. rubric:: Changed

- Improved efficiency of the sparse SDPA file format writer.
- Improved efficiency of the complex to real transformation.

.. rubric:: Fixed

- Scalar product of hermitian matrices.
- Conjugate of a complex expression.


`1.0.1`_ - 2014-08-27
--------------------------------------------------------------------------------

.. rubric:: Added

- Support for semidefinite programming over the complex domain, see
  :ref:`the documentation on complex problems <complex>`.
- Helper function to input (multicommodity) graph flow problems, see
  :ref:`the tutorial on flow constraints <flowcons>`.
- Additional argument to :func:`~picos.tracepow`, to represent constraints
  of the form :math:`\operatorname{trace}(M X^p) \geq t`.

.. rubric:: Changed

- Significantly improved slicing performance for affine expressions.
- Improved performance when loading data.
- Improved performance when retrieving primal solution from CPLEX.
- The documentation received an overhaul.


`1.0.0`_ - 2013-07-19
--------------------------------------------------------------------------------

.. rubric:: Added

- Ability to express rational powers of affine expressions with the ``**``
  operator, traces of matrix powers with :func:`~picos.tracepow`,
  (generalized) p-norms with :func:`~picos.norm` and :math:`n`-th roots of a
  determinant with :func:`~picos.detrootn`.
- Ability to specify variable bounds directly rather than by adding constraints,
  see :meth:`add_variable <.problem.Problem.add_variable>`.
- Problem dualization.
- Option ``solve_via_dual`` which controls passing the dual problem to the
  solver instead of the primal problem. This can result in a significant
  speedup for certain problems.
- Semidefinite programming interface for MOSEK 7.0.
- Options ``handleBarVars`` and ``handleConeVars`` to customize how SOCPs and
  SDPs are passed to MOSEK. When these are set to ``True``, PICOS tries to
  minimize the number of variables of the MOSEK instance.

.. rubric:: Changed

- If the chosen solver supports this, updated problems will be partially
  re-solved instead of solved from scratch.

.. rubric:: Removed

- Option ``onlyChangeObjective``.


`0.1.3`_ - 2013-04-17
--------------------------------------------------------------------------------

.. rubric:: Added

- A :func:`~picos.geomean` function to construct geometric mean inequalities
  that will be cast as rotated second order cone constraints.
- Options ``uboundlimit`` and ``lboundlimit`` to tell CPLEX to stop the search
  as soon as the given threshold is reached for the upper and lower bound,
  respectively.
- Option ``boundMonitor`` to inspect the evolution of CPLEX lower and upper
  bounds.
- Ability to use the weak inequality operators as an alias for the strong ones.

.. rubric:: Changed

- The solver search time is returned in the dictionary returned by
  :meth:`solve <.problem.Problem.solve>`.

.. rubric:: Fixed

- Access to dual values of fixed variables with CPLEX.
- Evaluation of constant affine expressions with a zero coefficient.
- Number of constraints not being updated in
  :meth:`remove_constraint <.problem.Problem.remove_constraint>`.


`0.1.2`_ - 2013-01-10
--------------------------------------------------------------------------------

.. rubric:: Fixed

- Writing SDPA files. The lower triangular part of the constraint matrix was
  written instead of the upper triangular part.
- A wrongly raised :class:`IndexError` from
  :meth:`remove_constraint <.problem.Problem.remove_constraint>`.


`0.1.1`_ - 2012-12-08
--------------------------------------------------------------------------------

.. rubric:: Added

- Interface to Gurobi.
- Ability to give an initial solution to warm-start mixed integer optimizers.
- Ability to get a reference to a constraint that was added.

.. rubric:: Fixed

- Minor bugs with quadratic expressions.


`0.1.0`_ - 2012-06-22
--------------------------------------------------------------------------------

.. rubric:: Added

- Initial release of PICOS.
Contribution Guide
==================

Filing a bug report or feature request
--------------------------------------

.. rubric:: Via GitLab

If you have a GitLab account, just head to PICOS' official
`issue tracker <https://gitlab.com/picos-api/picos/issues>`_.

.. rubric:: Via mail

If you don't have a GitLab account you can still create an issue by writing to
`incoming+picos-api/picos@incoming.gitlab.com
<incoming+picos-api/picos@incoming.gitlab.com>`_. Unlike issues created directly
on GitLab, issues created by mail are *not* publicly visible.

Submitting a code change
------------------------

The canonical way to submit a code change is to

1. fork the `PICOS repository on GitLab <https://gitlab.com/picos-api/picos>`_,
2. clone your fork and make your application use it instead of your system's
   PICOS installation,
3. optionally create a local topic branch to work with,
4. modify the source and commit your changes, and lastly
5. make a pull request on GitLab so that we can test and merge your changes.

Code style
----------

Set your linter to enforce :pep:`8` and :pep:`257` except for the following
codes:

.. code::

    D105,D203,D213,D401,E122,E128,E221,E271,E272,E501,E702,E741

Our line width limit is ``80`` characters.

Release procedure
-----------------

.. rubric:: Version scheme

When installed from the git repository or from a source distribution (sdist),
PICOS versions have the format ``MAJOR.MINOR.PATCH``, where ``PATCH`` is the
commit distance to the last minor release. When installed from a source tarball
that carries no metadata from either git or setuptools, the version format is
just ``MAJOR.MINOR`` as the commit distance defining the ``PATCH`` bit is not
known. Note that the ``PATCH`` bit is not strictly incremental as not every
commit to the ``master`` branch is released individually.

.. rubric:: Bumping the version

To bump the major or minor version, run ``version.py -b MAJOR.MINOR``. This
commits that base version to the ``picos/.version`` file and creates an
annotated tag (``vMAJOR.MINOR``) for that commit. The release of version
``MAJOR.MINOR.0`` is then published by pushing the tagged commit to the top of
``origin/master``, which triggers a GitLab CI/CD pipeline for that commit. By
the same logic, the ``PATCH`` bit is bumped and the resulting version is
published automatically whenever a number of commits is pushed to the ``master``
branch between two minor versions.

Note that source distributions are aware of the ``PATCH`` bit as setuptools
writes it to the ``picos/.version`` file in the source tree.

.. rubric:: Justification

The result of this unorthodox release procedure is that bugfix releases can be
made quickly simply by pushing a commit to ``master``. On the other hand,
changes that should go into the next minor or major release must remain on topic
branches and pushed onto ``master`` together with the commit from
``version.py -b``.

Implementing a test case
------------------------

Production test sets are implemented in the files in the ``tests`` folder that
start with ``ptest_``. If you want to add to our test pool, feel free to either
extend these files or create a new set, whatever is appropriate. Make sure that
the tests you add are not too computationally expensive since they are also run
as part of our continuous integration pipeline whenever a commit is pushed to
GitLab.

Implementing a solver
---------------------

If you want to implement support for a new solver, all you have to do is update
``solvers/__init__.py`` where applicable, and add a file named
``solver_<name>.py`` in the same directory with your implementation. We
recommend that you read two or three of the existing solver implementations to
get an idea of how things are done. If you want to know exactly how PICOS
communicates with your implementation, refer to the solver base class in
``solver.py``.
Introduction
============

PICOS is a user friendly Python API to several conic and integer programming
solvers, designed to be used by both application developers and researchers as
well as instructors teaching courses on mathematical optimization. It allows you
to enter an optimization problem as a **high level model**, with painless
support for **(complex) vector and matrix variables** and **multidimensional
algebra**. Your model will be transformed to the standard form understood by an
appropriate solver that is available at runtime. This makes your application
**portable** as users have the choice between several commercial and open source
solvers.

Features
--------

PICOS supports the following solvers and problem types. To use a solver, you
need to separately install it along with the Python interface listed here.

.. _Apache-2.0: https://www.apache.org/licenses/LICENSE-2.0
.. _GPL-3: https://www.gnu.org/licenses/gpl-3.0.html
.. _MIT: https://opensource.org/licenses/MIT
.. _ZIB: https://scip.zib.de/academic.txt

.. list-table::
    :header-rows: 1

    * - | Solver
        |
      - | Python
        | interface
      - | `LP <https://en.wikipedia.org/wiki/Linear_programming>`_
        |
      - | `SOCP <https://en.wikipedia.org/wiki/Second-order_cone_programming>`_,
        | `QCQP <https://en.wikipedia.org/wiki/Quadratically_constrained_quadratic_program>`_
      - | `SDP <https://en.wikipedia.org/wiki/Semidefinite_programming>`_
        |
      - | `EXP <https://docs.mosek.com/modeling-cookbook/expo.html>`_
        |
      - | `MIP <https://en.wikipedia.org/wiki/Integer_programming>`_
        |
      - | License
        |
    * - `CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`_
      - included
      - Yes
      - Yes
      -
      -
      - Yes
      - non-free
    * - `CVXOPT <https://cvxopt.org/>`_
      - native
      - Yes
      - Yes
      - Yes
      - `GP <https://en.wikipedia.org/wiki/Geometric_programming>`_
      -
      - `GPL-3`_
    * - `ECOS <https://github.com/embotech/ecos>`_
      - `ecos-python <https://github.com/embotech/ecos-python>`_
      - Yes
      - Yes
      -
      - Yes
      - Yes
      - `GPL-3`_
    * - `GLPK <https://www.gnu.org/software/glpk/>`_
      - `swiglpk <https://github.com/biosustain/swiglpk>`_
      - Yes
      -
      -
      -
      - Yes
      - `GPL-3`_
    * - `Gurobi <http://www.gurobi.com/products/gurobi-optimizer>`_
      - `gurobipy <https://www.gurobi.com>`_
      - Yes
      - Yes
      -
      -
      - Yes
      - non-free
    * - `MOSEK <https://www.mosek.com/>`_
      - included
      - Yes
      - Yes
      - Yes
      -
      - Yes
      - non-free
    * - `OSQP <https://osqp.org>`_
      - native
      - Yes
      - `QP <https://en.wikipedia.org/wiki/Quadratic_programming>`_
      -
      -
      -
      - `Apache-2.0`_
    * - `SCIP <http://scip.zib.de/>`_
      - `PySCIPOpt <https://github.com/SCIP-Interfaces/PySCIPOpt/>`_
      - Yes
      - Yes
      -
      -
      - Yes
      - `ZIB`_/`MIT`_
    * - `SMCP <http://smcp.readthedocs.io/en/latest/>`_
      - native
      -
      -
      - Yes
      -
      -
      - `GPL-3`_

.. rubric:: Example

This is what it looks like to solve a multidimensional mixed integer program
with PICOS:

>>> import picos as pc
>>> P = pc.Problem()
>>> x = pc.IntegerVariable("x", 2)
>>> P += 2*x <= 11
>>> P.maximize = pc.sum(x)
>>> P.solve(solver="glpk")  # Optional: Use GLPK as backend.
<feasible primal solution (claimed optimal) from glpk>
>>> P.value
10.0
>>> print(x)
[ 5.00e+00]
[ 5.00e+00]

You can head to our
`quick examples <https://picos-api.gitlab.io/picos/quick.html>`_ or the
`tutorial <https://picos-api.gitlab.io/picos/tutorial.html>`_ for more.

Installation
------------

As of release 2.2, PICOS requires **Python 3.4** or later.

.. rubric:: Via pip

If you are using `pip <https://pypi.org/project/pip/>`_ you can run
``pip install picos`` to get the latest version.

.. rubric:: Via Anaconda

If you are using `Anaconda <https://anaconda.org/>`_ you can run
``conda install -c picos picos`` to get the latest version.

.. rubric:: Via your system's package manager

.. list-table::
    :header-rows: 1
    :stub-columns: 1

    * - Distribution
      - Latest major version
      - Latest version
    * - Arch Linux
      - `python-picos <https://aur.archlinux.org/packages/python-picos/>`__
      - `python-picos-git <https://aur.archlinux.org/packages/python-picos-git/>`__

If you are packaging PICOS for additional platforms, please let us know.

.. rubric:: From source

The PICOS source code can be found on `GitLab
<https://gitlab.com/picos-api/picos>`_. There are only two dependencies:

- `NumPy <https://numpy.org/>`_
- `CVXOPT`_

Documentation
-------------

The full documentation can be browsed `online
<https://picos-api.gitlab.io/picos/>`__ or downloaded `in PDF form
<https://picos-api.gitlab.io/picos/picos.pdf>`__.

Credits
-------

.. rubric:: Developers

- `Guillaume Sagnol <http://page.math.tu-berlin.de/~sagnol/>`_ has started work
  on PICOS in 2012.
- `Maximilian Stahlberg <https://orcid.org/0000-0002-0190-2693>`_ is extending
  and co-maintaining PICOS since 2017.

.. rubric:: Contributors

For an up-to-date list of all code contributors, please refer to the
`contributors page <https://gitlab.com/picos-api/picos/-/graphs/master>`_.
Should a reference from before 2019 be unclear, see also the `old contributors
page <https://github.com/gsagnol/picos/graphs/contributors>`_ on GitHub.

Citing
------

The preferred way to cite PICOS in your research is our `JOSS paper
<https://joss.theoj.org/papers/10.21105/joss.03915>`_:

.. code-block:: bibtex

  @article{PICOS,
    author  = {Guillaume Sagnol and Maximilian Stahlberg},
    journal = {Journal of Open Source Software},
    title   = {{PICOS}: A {Python} interface to conic optimization solvers},
    year    = {2022},
    issn    = {2475-9066},
    month   = feb,
    number  = {70},
    pages   = {3915},
    volume  = {7},
    doi     = {10.21105/joss.03915},
  }

If citing a specific version of PICOS is necessary, then we offer also `source
deposits on Zenodo <https://doi.org/10.5281/zenodo.6052843>`_.

License
-------

PICOS is free and open source software and available to you under the terms of
the `GNU GPL v3 <https://gitlab.com/picos-api/picos/raw/master/LICENSE.txt>`_.
.. _changelog:

.. include:: ../CHANGELOG.rst
.. _optdes:

Optimal Experimental Design
===========================

Optimal experimental design is a theory
at the interface of statistics and optimization,
which studies how to allocate some statistical trials
within a set of available design points.
The goal is to allow for the best possible
estimation of an unknown parameter :math:`\theta`.
In what follows, we assume the standard linear model with
multiresponse experiments: a trial in the :math:`i^{\textrm{th}}`
design point gives a multidimensional observation that
can be written as :math:`y_i = A_i^T \theta+\epsilon_i`,
where :math:`y_i` is of dimension :math:`l_i`,
:math:`A_i` is a :math:`m \times l_i-` matrix,
and the error vectors :math:`\epsilon_i` are i.i.d. with a unit variance.

Several optimization criteria exist, leading to different SDP, SOCP and LP
formulations.
As such, optimal experimental design problens are natural examples for problems
in conic optimization. For a review of the different formulations
and more references, see :ref:`[1] <optdes_refs>`.

The code below initializes the data used in all the examples of this page.
It should be run prior to any of the codes presented in this page.

>>> import cvxopt as cvx
>>> import picos
>>> #---------------------------------#
>>> # First generate some data :      #
>>> #       _ a list of 8 matrices A  #
>>> #       _ a vector c              #
>>> #---------------------------------#
>>> A = [cvx.matrix([[1,0,0,0,0],
...                  [0,3,0,0,0],
...                  [0,0,1,0,0]]),
...      cvx.matrix([[0,0,2,0,0],
...                  [0,1,0,0,0],
...                  [0,0,0,1,0]]),
...      cvx.matrix([[0,0,0,2,0],
...                  [4,0,0,0,0],
...                  [0,0,1,0,0]]),
...      cvx.matrix([[1,0,0,0,0],
...                  [0,0,2,0,0],
...                  [0,0,0,0,4]]),
...      cvx.matrix([[1,0,2,0,0],
...                  [0,3,0,1,2],
...                  [0,0,1,2,0]]),
...      cvx.matrix([[0,1,1,1,0],
...                  [0,3,0,1,0],
...                  [0,0,2,2,0]]),
...      cvx.matrix([[1,2,0,0,0],
...                  [0,3,3,0,5],
...                  [1,0,0,2,0]]),
...      cvx.matrix([[1,0,3,0,1],
...                  [0,3,2,0,0],
...                  [1,0,0,2,0]])
... ]
>>> c = cvx.matrix([1,2,3,4,5])

Multi-response c-optimal design (SOCP)
--------------------------------------

We compute the c-optimal design (``c=[1,2,3,4,5]``)
for the observation matrices ``A[i].T`` from the variable ``A`` defined above.
The results below suggest that we should allocate 12.8% of the
experimental effort on design point #5, and 87.2% on the design point #7.

.. rubric:: Primal problem

The SOCP for multiresponse c-optimal design is:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\substack{\mu \in \mathbb{R}^s\\
                        \forall i \in [s],\ z_i \in \mathbb{R}^{l_i}}}{\mbox{minimize}}
                      & \sum_{i=1}^s \mu_i\\
   &\mbox{subject to} & \sum_{i=1}^s A_i z_i = c\\
   &                  & \forall i \in [s],\ \Vert z_i \Vert_2 \leq \mu_i,
   \end{eqnarray*}
   \end{center}

>>> # create the problem, variables and params
>>> c_primal_SOCP = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> cc = picos.Constant('c', c)
>>> z  = [picos.RealVariable('z[{0}]'.format(i), AA[i].size[1])
...        for i in range(s)]
>>> mu = picos.RealVariable('mu', s)

>>> # define the constraints and objective function
>>> cones = c_primal_SOCP.add_list_of_constraints([abs(z[i]) <= mu[i] for i in range(s)])
>>> lin   = c_primal_SOCP.add_constraint(picos.sum([AA[i] * z[i] for i in range(s)]) == cc)
>>> c_primal_SOCP.set_objective('min', (1|mu) )
>>> print(c_primal_SOCP)
Second Order Cone Program
  minimize ∑(mu)
  over
    3×1 real variable z[i] ∀ i ∈ [0…7]
    8×1 real variable mu
  subject to
    ‖z[i]‖ ≤ mu[i] ∀ i ∈ [0…7]
    ∑(A[i]·z[i] : i ∈ [0…7]) = c

>>> #solve the problem and retrieve the optimal weights of the optimal design.
>>> solution = c_primal_SOCP.solve(solver='cvxopt')
>>> mu = mu.value
>>> w = mu / sum(mu) #normalize mu to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[...]
[...]
[ 1.28e-01]
[...]
[ 8.72e-01]
[...]

The ``[...]`` above indicate a numerical zero entry
(*i.e., which can be something like 2.84e-10*).
We use the ellipsis ``...`` instead for clarity and compatibility with **doctest**.

.. rubric:: Dual problem

This is only to check that we obtain the same solution with the dual problem,
and to provide one additional example in this tutorial:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{u \in \mathbb{R}^m}{\mbox{maximize}}
                      & c^T u\\
   &\mbox{subject to} & \forall i \in [s],\ \Vert A_i^T u \Vert_2 \leq 1
   \end{eqnarray*}
   \end{center}


>>> # create the problem, variables and params
>>> c_dual_SOCP = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> cc = picos.Constant('c',c)
>>> u  = picos.RealVariable('u',c.size)
>>> # define the constraints and objective function
>>> cones = c_dual_SOCP.add_list_of_constraints(
...         [abs(AA[i].T*u)<=1 for i in range(s)])
>>> c_dual_SOCP.set_objective('max', (cc|u) )
>>> print(c_dual_SOCP)#
Second Order Cone Program
  maximize ⟨c, u⟩
  over
    5×1 real variable u
  subject to
    ‖A[i]ᵀ·u‖ ≤ 1 ∀ i ∈ [0…7]
>>> #solve the problem and retrieve the weights of the optimal design
>>> solution = c_dual_SOCP.solve(solver='cvxopt')
>>> mu = [cons.dual[0] for cons in cones] #Lagrangian duals of the SOC constraints
>>> mu = cvx.matrix(mu)
>>> w=mu/sum(mu) #normalize mu to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[...]
[...]
[ 1.28e-01]
[...]
[ 8.72e-01]
[...]


Single-response c-optimal design (LP)
-------------------------------------

When the observation matrices are row vectors (single-response framework),
the SOCP above reduces to a simple LP, because the variables
:math:`z_i` are scalar.
We solve below the LP for the case where there are 11
available design points, corresponding to the columns of the matrices
``A[4]``, ``A[5]``, ``A[6]``, and ``A[7][:,:-1]`` defined in the preambule.

The optimal design allocates 3.37% to point #5 (2nd column of ``A[5]``),
27.9% to point #7 (1st column of ``A[6]``),
11.8% to point #8 (2nd column of ``A[6]``),
27.6% to point #9 (3rd column of ``A[6]``),
and 29.3% to point #11 (2nd column of ``A[7]``).

>>> # create the problem, variables and params
>>> c_primal_LP = picos.Problem()
>>> A1 = [cvx.sparse(a[:,i],tc='d') for i in range(3) for a in A[4:]] #12 column vectors
>>> A1 = A1[:-1] # remove the last design point (it is the same as the last-but-one)
>>> s = len(A1)
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A1)] # each AA[i].T is a 1 x 5 observation matrix
>>> cc = picos.Constant('c', c)
>>> z = [picos.RealVariable('z[{0}]'.format(i), 1) for i in range(s)]
>>> mu = picos.RealVariable('mu', s)

>>> #define the constraints and objective function
>>> abs_con = c_primal_LP.add_list_of_constraints([abs(z[i]) <= mu[i] for i in range(s)])
>>> lin_con = c_primal_LP.add_constraint(picos.sum([AA[i]*z[i] for i in range(s)]) == cc)
>>> c_primal_LP.set_objective('min', (1|mu))

Note that there are no cone constraints, because
the constraints of the form :math:`|z_i| \leq \mu_i` are handled as two
inequalities when :math:`z_i` is scalar, so the problem is a LP indeed:

>>> print(c_primal_LP)
Linear Program
  minimize ∑(mu)
  over
    1×1 real variable z[i] ∀ i ∈ [0…10]
    11×1 real variable mu
  subject to
    |z[i]| ≤ mu[i] ∀ i ∈ [0…10]
    ∑(A[i]·z[i] : i ∈ [0…10]) = c

>>> #solve the problem and retrieve the weights of the optimal design
>>> solution = c_primal_LP.solve(solver='cvxopt')
>>> mu = mu.value
>>> w = mu / sum(mu) #normalize mu to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[...]
[...]
[ 3.37e-02]
[...]
[ 2.79e-01]
[ 1.18e-01]
[ 2.76e-01]
[...]
[ 2.93e-01]

SDP formulation of c-optimal design
-----------------------------------

We give below the SDP for c-optimality, in primal and dual
form. You can observe that we obtain the same results as
with the SOCP presented earlier:
12.8% on design point #5, and 87.2% on design point #7.

.. rubric:: Primal problem

The SDP formulation of the c-optimal design problem is:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\mu \in \mathbb{R}^s}{\mbox{minimize}}
                      & \sum_{i=1}^s \mu_i\\
   &\mbox{subject to} & \sum_{i=1}^s \mu_i A_i A_i^T \succeq c c^T,\\
   &                  & \mu \geq 0.
   \end{eqnarray*}
   \end{center}

>>> # create the problem, variables and params
>>> c_primal_SDP = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> cc = picos.Constant('c', c)
>>> mu = picos.RealVariable('mu',s)
>>> # define the constraints and objective function
>>> lmi = c_primal_SDP.add_constraint(
...         picos.sum([mu[i] * AA[i] * AA[i].T for i in range(s)]) >> cc*cc.T)
>>> lin_cons = c_primal_SDP.add_constraint(mu >= 0)
>>> c_primal_SDP.set_objective('min', (1|mu) )
>>> print(c_primal_SDP)
Semidefinite Program
  minimize ∑(mu)
  over
    8×1 real variable mu
  subject to
    ∑(mu[i]·A[i]·A[i]ᵀ : i ∈ [0…7]) ≽ c·cᵀ
    mu ≥ 0

>>> #solve the problem and retrieve the weights of the optimal design
>>> solution = c_primal_SDP.solve(solver='cvxopt')
>>> w = mu.value
>>> w = w / sum(w) #normalize mu to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[...]
[...]
[ 1.28e-01]
[...]
[ 8.72e-01]
[...]

.. rubric:: Dual problem

This is only to check that we obtain the same solution with the dual problem,
and to provide one additional example in this tutorial:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{X \in \mathbb{R}^{m \times m}}{\mbox{maximize}}
                      &  c^T X c\\
   &\mbox{subject to} & \forall i \in [s],\ \langle A_i A_i^T,\ X \rangle \leq 1,\\
   &                  &  X \succeq 0.
   \end{eqnarray*}
   \end{center}

>>> #create the problem, variables and params
>>> c_dual_SDP = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> cc = picos.Constant('c', c)
>>> m  = c.size[0]
>>> X  = picos.SymmetricVariable('X',(m,m))

>>> #define the constraints and objective function
>>> lin_cons = c_dual_SDP.add_list_of_constraints(
...                  [(AA[i]*AA[i].T | X ) <= 1 for i in range(s)])
>>> psd = c_dual_SDP.add_constraint(X>>0)
>>> c_dual_SDP.set_objective('max', cc.T*X*cc)

>>> print(c_dual_SDP)
Semidefinite Program
  maximize cᵀ·X·c
  over
    5×5 symmetric variable X
  subject to
    ⟨A[i]·A[i]ᵀ, X⟩ ≤ 1 ∀ i ∈ [0…7]
    X ≽ 0

>>> # solve the problem and retrieve the weights of the optimal design
>>> solution = c_dual_SDP.solve(solver='cvxopt')
>>> mu = [cons.dual for cons in lin_cons] #Lagrangian duals of the linear constraints
>>> mu = cvx.matrix(mu)
>>> w = mu / sum(mu) #normalize mu to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[...]
[...]
[ 1.28e-01]
[...]
[ 8.72e-01]
[...]

And the optimal positive semidefinite matrix X is:

>>> print(X)
[ 5.92e-03  8.98e-03  2.82e-03 -3.48e-02 -1.43e-02]
[ 8.98e-03  1.36e-02  4.27e-03 -5.28e-02 -2.17e-02]
[ 2.82e-03  4.27e-03  1.34e-03 -1.66e-02 -6.79e-03]
[-3.48e-02 -5.28e-02 -1.66e-02  2.05e-01  8.39e-02]
[-1.43e-02 -2.17e-02 -6.79e-03  8.39e-02  3.44e-02]

A-optimality (SOCP)
-------------------

We compute the A-optimal design
for the observation matrices ``A[i].T`` defined in the preambule.
The optimal design allocates
24.9% on design point #3,
14.2% on point #4,
8.51% on point #5,
12.1% on point #6,
13.2% on point #7,
and 27.0% on point #8.

.. rubric:: Primal problem

The SOCP for the A-optimal design problem is:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\substack{\mu \in \mathbb{R}^s\\
                        \forall i \in [s],\ Z_i \in \mathbb{R}^{l_i \times m}}}{\mbox{minimize}}
                      & \sum_{i=1}^s \mu_i\\
   &\mbox{subject to} & \sum_{i=1}^s A_i Z_i = I\\
   &                  & \forall i \in [s],\ \Vert Z_i \Vert_F \leq \mu_i,
   \end{eqnarray*}
   \end{center}

>>> # create the problem, variables and params
>>> A_primal_SOCP = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> Z = [picos.RealVariable('Z[{0}]'.format(i), AA[i].T.size) for i in range(s)]
>>> mu = picos.RealVariable('mu', s)

>>> #define the constraints and objective function
>>> cone_cons = A_primal_SOCP.add_list_of_constraints(
...                     [abs(Z[i]) <= mu[i] for i in range(s)])
>>> lin_cons = A_primal_SOCP.add_constraint(
...                      picos.sum([AA[i] * Z[i] for i in range(s)]) == 'I')
>>> A_primal_SOCP.set_objective('min', (1|mu) )
>>> print(A_primal_SOCP)
Second Order Cone Program
  minimize ∑(mu)
  over
    3×5 real variable Z[i] ∀ i ∈ [0…7]
    8×1 real variable mu
  subject to
    ‖Z[i]‖ ≤ mu[i] ∀ i ∈ [0…7]
    ∑(A[i]·Z[i] : i ∈ [0…7]) = I

>>> # solve the problem and retrieve the weights of the optimal design
>>> solution = A_primal_SOCP.solve(solver='cvxopt')
>>> w = mu.value
>>> w = w / sum(w) #normalize mu to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[ 2.49e-01]
[ 1.42e-01]
[ 8.51e-02]
[ 1.21e-01]
[ 1.32e-01]
[ 2.70e-01]

.. rubric:: Dual problem

This is only to check that we obtain the same solution with the dual problem,
and to provide one additional example in this tutorial:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{U \in \mathbb{R}^{m \times m}}{\mbox{maximize}}
                      &  \mbox{trace}\ U\\
   &\mbox{subject to} & \forall i \in [s],\ \Vert A_i^T U \Vert_2 \leq 1
   \end{eqnarray*}
   \end{center}

>>> #create the problem, variables and params
>>> D_SOCPual_A=picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> m  = AA[0].size[0]
>>> U  = picos.RealVariable('U',(m,m))
>>> #define the constraints and objective function
>>> cone_cons = D_SOCPual_A.add_list_of_constraints(
...       [abs(AA[i].T*U) <= 1 for i in range(s)])
>>> D_SOCPual_A.set_objective('max', 'I'|U)
>>> print(D_SOCPual_A)
Second Order Cone Program
  maximize tr(U)
  over
    5×5 real variable U
  subject to
    ‖A[i]ᵀ·U‖ ≤ 1 ∀ i ∈ [0…7]

>>> # solve the problem and retrieve the weights of the optimal design
>>> solution = D_SOCPual_A.solve(solver='cvxopt')
>>> mu = [cons.dual[0] for cons in cone_cons] # Lagrangian duals of the SOC constraints
>>> mu = cvx.matrix(mu)
>>> w = mu / sum(mu) # normalize mu to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[ 2.49e-01]
[ 1.42e-01]
[ 8.51e-02]
[ 1.21e-01]
[ 1.32e-01]
[ 2.70e-01]

A-optimality with multiple constraints (SOCP)
---------------------------------------------

A-optimal designs can also be computed by SOCP
when the vector of weights :math:`\mathbf{w}` is subject
to several linear constraints.
To give an example, we compute the A-optimal design for
the observation matrices given in the preambule, when the weights
must satisfy: :math:`\sum_{i=0}^3 w_i \leq 0.5` and :math:`\sum_{i=4}^7 w_i \leq 0.5`.
This problem has the following SOCP formulation:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\substack{\mathbf{w} \in \mathbb{R}^s\\
                        \mu \in \mathbb{R}^s\\
                        \forall i \in [s],\ Z_i \in \mathbb{R}^{l_i \times m}}}{\mbox{minimize}}
                      & \sum_{i=1}^s \mu_i\\
   &\mbox{subject to} & \sum_{i=1}^s A_i Z_i = I\\
   &                  & \sum_{i=0}^3 w_i \leq 0.5\\
   &                  & \sum_{i=4}^7 w_i \leq 0.5\\
   &                  & \forall i \in [s],\ \Vert Z_i \Vert_F^2 \leq \mu_i w_i,
   \end{eqnarray*}
   \end{center}

The optimal solution allocates 29.7% and 20.3% to the design points #3 and #4,
and  respectively 6.54%, 11.9%, 9.02% and 22.5% to the design points #5 to #8:

>>> # create the problem, variables and params
>>> A_multiconstraints = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> mu = picos.RealVariable('mu',s)
>>> w  = picos.RealVariable('w',s)
>>> Z  = [picos.RealVariable('Z[{0}]'.format(i), AA[i].T.size)
...                          for i in range(s)]
>>> # define the constraints and objective function
>>> lin_cons0 = A_multiconstraints.add_constraint(
...         picos.sum([AA[i] * Z[i] for i in range(s)]) == 'I')
>>> lin_cons1 = A_multiconstraints.add_constraint( (1|w[:4]) <= 0.5)
>>> lin_cons2 = A_multiconstraints.add_constraint( (1|w[4:]) <= 0.5)
>>> cone_cons = A_multiconstraints.add_list_of_constraints(
...       [ abs(Z[i]) **2 <= mu[i] * w[i] for i in range(s)])
>>> A_multiconstraints.set_objective('min', (1|mu) )
>>> print(A_multiconstraints)
Quadratically Constrained Program
  minimize ∑(mu)
  over
    3×5 real variable Z[i] ∀ i ∈ [0…7]
    8×1 real variables mu, w
  subject to
    ∑(A[i]·Z[i] : i ∈ [0…7]) = I
    ∑(w[:4]) ≤ 0.5
    ∑(w[4:]) ≤ 0.5
    ‖Z[i]‖² ≤ mu[i]·w[i] ∀ i ∈ [0…7]

>>> # solve the problem and retrieve the weights of the optimal design
>>> solution = A_multiconstraints.solve(solver='cvxopt')
>>> w = w.value
>>> w = w / sum(w) # normalize w to get the optimal weights

The optimal design is:

>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[ 2.97e-01]
[ 2.03e-01]
[ 6.54e-02]
[ 1.19e-01]
[ 9.02e-02]
[ 2.25e-01]


Exact A-optimal design (MISOCP)
-------------------------------

In the exact version of A-optimality, a number :math:`N \in \mathbb{N}`
of trials is given, and the goal is to find the optimal number of times
:math:`n_i \in \mathbb{N}` that a trial on design point #i should be performed,
with :math:`\sum_i n_i =N`.

The SOCP formulation of A-optimality for constrained designs
also accept integer constraints, which results in a MISOCP for exact A-optimality:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\substack{\mathbf{t} \in \mathbb{R}^s\\
                        \mathbf{n} \in \mathbb{N}^s\\
                        \forall i \in [s],\ Z_i \in \mathbb{R}^{l_i \times m}}}{\mbox{minimize}}
                      & \sum_{i=1}^s t_i\\
   &\mbox{subject to} & \sum_{i=1}^s A_i Z_i = I\\
   &                  & \forall i \in [s],\ \Vert Z_i \Vert_F^2 \leq n_i t_i,\\
   &                  & \sum_{i=1}^s n_i = N.
   \end{eqnarray*}
   \end{center}

The exact optimal design is :math:`\mathbf{n}=[0,0,5,3,2,2,3,5]`:

>>> # create the problem, variables and params
>>> A_exact = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> m  = AA[0].size[0]
>>> N  = picos.Constant('N', 20) # number of trials allowed
>>> I = picos.Constant('I', cvx.spmatrix([1]*m,range(m),range(m),(m,m))) #identity matrix
>>> Z = [picos.RealVariable('Z[{0}]'.format(i), AA[i].T.size) for i in range(s)]
>>> n = picos.IntegerVariable('n', s)
>>> t = picos.RealVariable('t', s)

>>> # define the constraints and objective function
>>> cone_cons = A_exact.add_list_of_constraints(
...         [ abs(Z[i])**2 <= n[i] * t[i] for i in range(s)])
>>> lin_cons = A_exact.add_constraint(
...          picos.sum([AA[i]*Z[i] for i in range(s)]) == I)
>>> wgt_cons = A_exact.add_constraint( (1|n) <= N )
>>> A_exact.set_objective('min',1|t)
>>> print(A_exact)
Mixed-Integer Quadratically Constrained Program
  minimize ∑(t)
  over
    8×1 integer variable n
    3×5 real variable Z[i] ∀ i ∈ [0…7]
    8×1 real variable t
  subject to
    ‖Z[i]‖² ≤ n[i]·t[i] ∀ i ∈ [0…7]
    ∑(A[i]·Z[i] : i ∈ [0…7]) = I
    ∑(n) ≤ N

>>> #solve the problem and display the optimal design
>>> solution = A_exact.solve()# doctest:+SKIP
>>> print(n)# doctest:+SKIP
[...]
[...]
[ 5.00e+00]
[ 3.00e+00]
[ 2.00e+00]
[ 2.00e+00]
[ 3.00e+00]
[ 5.00e+00]

.. note::

    The above output is not validated as we lack an appropriate solver on
    the build server.

Approximate and exact D-optimal design ((MI)SOCP)
-------------------------------------------------

The D-optimal design problem has a SOCP formulation involving a
geometric mean in the objective function:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\substack{\mathbf{L} \in \mathbb{R}^{m \times m}\\
                        \mathbf{w} \in \mathbb{R}^s\\
                        \forall i \in [s],\ V_i \in \mathbb{R}^{l_i \times m}}}{\mbox{maximize}}
                      & \left(\prod_{i=1}^m L_{i,i}\right)^{1/m}\\
   &\mbox{subject to} & \sum_{i=1}^s A_i V_i = L,\\
   &                  & L\ \mbox{lower triangular},\\
   &                  & \Vert V_i \Vert_F \leq \sqrt{m}\ w_i,\\
   &                  & \sum_{i=1}^s w_i \leq 1.
   \end{eqnarray*}
   \end{center}

By introducing a new variable :math:`t` such that
:math:`t \leq \left(\prod_{i=1}^m L_{i,i}\right)^{1/m}`, we can pass
this problem to PICOS with the function :func:`~picos.geomean`,
which reformulates the geometric mean inequality as a set of equivalent second order cone
constraints.
The example below allocates respectively 22.7%, 3.38%, 1.65%, 5.44%, 31.8% and 35.1%
to the design points #3 to #8.

>>> #create the problem, variables and params
>>> D_SOCP = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> m  = AA[0].size[0]
>>> mm = picos.Constant('m', m)
>>> L = picos.RealVariable('L', (m,m))
>>> V = [picos.RealVariable('V['+str(i)+']', AA[i].T.size) for i in range(s)]
>>> w = picos.RealVariable('w',s)
>>> # additional variable to handle the geometric mean in the objective function
>>> t = picos.RealVariable('t',1)

>>> # define the constraints and objective function
>>> lin_cons = D_SOCP.add_constraint(picos.sum([AA[i]*V[i] for i in range(s)]) == L)
>>> # L is lower triangular
>>> lowtri_cons = D_SOCP.add_list_of_constraints( [L[i,j] == 0
...                for i in range(m)
...                for j in range(i+1,m)])
>>> cone_cons = D_SOCP.add_list_of_constraints([abs(V[i]) <= (mm**0.5)*w[i]
...                                                 for i in range(s)])
>>> wgt_cons = D_SOCP.add_constraint(1|w <= 1)
>>> geomean_cons = D_SOCP.add_constraint(t <= picos.geomean(picos.maindiag(L)))
>>> D_SOCP.set_objective('max',t)

>>> #solve the problem and display the optimal design
>>> print(D_SOCP)
Optimization Problem
  maximize t
  over
    1×1 real variable t
    3×5 real variable V[i] ∀ i ∈ [0…7]
    5×5 real variable L
    8×1 real variable w
  subject to
    L = ∑(A[i]·V[i] : i ∈ [0…7])
    L[i,j] = 0 ∀ (i,j) ∈ zip([0,0,…,2,3],[1,2,…,4,4])
    ‖V[i]‖ ≤ m^(1/2)·w[i] ∀ i ∈ [0…7]
    ∑(w) ≤ 1
    geomean(maindiag(L)) ≥ t

>>> solution = D_SOCP.solve(solver='cvxopt')
>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[...]
[...]
[ 2.27e-01]
[ 3.38e-02]
[ 1.65e-02]
[ 5.44e-02]
[ 3.18e-01]
[ 3.51e-01]


As for the A-optimal problem, there is an alternative SOCP formulation
of D-optimality :ref:`[2] <optdes_refs>`, in which integer constraints may be added.
This allows us to formulate the exact D-optimal problem as a MISOCP.
For :math:`N=20`,
we obtain the following N-exact D-optimal design:
:math:`\mathbf{n}=[0,0,5,1,0,1,6,7]`:

>>> # create the problem, variables and params
>>> D_exact = picos.Problem()
>>> L = picos.RealVariable('L',(m,m))
>>> V = [picos.RealVariable('V['+str(i)+']',AA[i].T.size) for i in range(s)]
>>> T = picos.RealVariable('T', (s,m))
>>> n = picos.IntegerVariable('n', s)
>>> N = picos.Constant('N', 20)
>>> # additional variable to handle the geomean inequality
>>> t = picos.RealVariable('t',1)

>>> # define the constraints and objective function
>>> lin_cons = D_exact.add_constraint(
...         picos.sum([AA[i]*V[i] for i in range(s)]) == L)
>>> # L is lower triangular
>>> lowtri_cons = D_exact.add_list_of_constraints( [L[i,j] == 0
...                                  for i in range(m)
...                                  for j in range(i+1,m)])
>>> cone_cons = D_exact.add_list_of_constraints([ abs(V[i][:,k])**2 <= n[i]/N*T[i,k]
...                 for i in range(s) for k in range(m)])
>>> lin_cons2 = D_exact.add_list_of_constraints([(1|T[:,k]) <= 1
...                       for k in range(m)])
>>> wgt_cons = D_exact.add_constraint(1|n <= N)
>>> geomean_cons = D_exact.add_constraint(t <= picos.geomean( picos.maindiag(L)))
>>> D_exact.set_objective('max',t)
>>> print(D_exact)
Mixed-Integer Optimization Problem
  maximize t
  over
    8×1 integer variable n
    1×1 real variable t
    3×5 real variable V[i] ∀ i ∈ [0…7]
    5×5 real variable L
    8×5 real variable T
  subject to
    L = ∑(A[i]·V[i] : i ∈ [0…7])
    L[i,j] = 0 ∀ (i,j) ∈ zip([0,0,…,2,3],[1,2,…,4,4])
    ‖V[i][:,j]‖² ≤ n[i]/N·T[i,j] ∀ (i,j) ∈
      zip([0,0,…,7,7],[0,1,…,3,4])
    ∑(T[:,i]) ≤ 1 ∀ i ∈ [0…4]
    ∑(n) ≤ N
    geomean(maindiag(L)) ≥ t

>>> #solve the problem and display the optimal design
>>> solution = D_exact.solve()# doctest:+SKIP
>>> print(n)# doctest:+SKIP
[...]
[...]
[ 5.00e+00]
[ 1.00e+00]
[...]
[ 1.00e+00]
[ 6.00e+00]
[ 7.00e+00]

.. note::

    The above output is not validated as we lack an appropriate solver on
    the build server.

Former MAXDET formulation of the D-optimal design (SDP)
-------------------------------------------------------

A so-called MAXDET Programming formulation of the D-optimal design
has been known since the late 90's :ref:`[3] <optdes_refs>`, and
can be reformulated as a SDP thanks to the :func:`~picos.detrootn` function.
The following code finds the same design as the SOCP approach presented above.

>>> # problem, variables and parameters
>>> D_MAXDET = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> m  = AA[0].size[0]
>>> w = picos.RealVariable('w', s, lower=0)
>>> t = picos.RealVariable('t', 1)
>>> # constraint and objective
>>> wgt_cons = D_MAXDET.add_constraint(1|w <= 1)
>>> Mw = picos.sum([w[i] * AA[i] * AA[i].T for i in range(s)])
>>> detrootn_cons = D_MAXDET.add_constraint(t <= picos.DetRootN(Mw))
>>> D_MAXDET.set_objective('max', t)

>>> print(D_MAXDET)
Optimization Problem
  maximize t
  over
    1×1 real variable t
    8×1 real variable w (bounded below)
  subject to
    ∑(w) ≤ 1
    det(∑(w[i]·A[i]·A[i]ᵀ : i ∈ [0…7]))^(1/5) ≥ t

>>> #solve and display
>>> solution = D_MAXDET.solve(solver='cvxopt')
>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[ ...]
[ ...]
[ 2.27e-01]
[ 3.38e-02]
[ 1.65e-02]
[ 5.44e-02]
[ 3.18e-01]
[ 3.51e-01]



General Phi_p optimal design (SDP)
----------------------------------

The A- and D-optimal design problems presented above can be obtained as special cases of the general
Kiefer :math:`\Phi_p-` optimal design problem, where :math:`p` is a real in :math:`(-\infty,1]` :

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{w \in \mathbb{R}^s}{\mbox{maximize}}
                      &\quad \left( \frac{1}{m} \operatorname{trace}\ \big(\sum_{i=1}^s w_i A_i A_i^T \big)^p \right)^{1/p} \\
   &\textrm{subject to} &\quad w\geq0,\ \sum_{i=1}^s w_i \leq 1.
   \end{eqnarray*}

   \end{center}

These problems are easy to enter in PICOS, thanks to the :func:`~picos.tracepow` function,
that automatically replaces inequalities involving trace of matrix powers as a set of equivalent linear matrix
inequalities (SDP) (cf. :ref:`[4] <optdes_refs>` ). Below are two examples with :math:`p=0.2` and :math:`p=-3`,
allocating respectively (20.6%, 0.0%, 0.0%, 0.92%, 40.8%, 37.7%), and
(24.8%, 16.6%, 10.8%, 14.1%, 7.84%, 26.0%) of the trials to the design points 3 to 8.

>>> #problems, variables and parameters
>>> P0dot2_SDP  = picos.Problem()
>>> Pminus3_SDP = picos.Problem()
>>> AA = [picos.Constant('A[{0}]'.format(i), Ai)
...       for i, Ai in enumerate(A)] # each AA[i].T is a 3 x 5 observation matrix
>>> s  = len(AA)
>>> m  = AA[0].size[0]
>>> w = picos.RealVariable('w', s, lower=0)
>>> t = picos.RealVariable('t', 1)

>>> # constraint and objective
>>> wgt02_cons = P0dot2_SDP.add_constraint(1|w <= 1)
>>> wgtm3_cons = Pminus3_SDP.add_constraint(1|w <= 1)

>>> Mw = picos.sum([w[i]*AA[i]*AA[i].T for i in range(s)])

>>> tracep02_cons = P0dot2_SDP.add_constraint(t <= picos.PowerTrace(Mw, 0.2))
>>> P0dot2_SDP.set_objective('max', t)

>>> tracepm3_cons = Pminus3_SDP.add_constraint(t >= picos.PowerTrace(Mw, -3))
>>> Pminus3_SDP.set_objective('min', t)

>>> # p=0.2
>>> print(P0dot2_SDP)
Optimization Problem
  maximize t
  over
    1×1 real variable t
    8×1 real variable w (bounded below)
  subject to
    ∑(w) ≤ 1
    tr(∑(w[i]·A[i]·A[i]ᵀ : i ∈ [0…7])^(1/5)) ≥ t

>>> #solve and display
>>> solution = P0dot2_SDP.solve(solver='cvxopt')
>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[ ...]
[ ...]
[ 2.06e-01]
[ ...]
[ ...]
[ 9.20e-03]
[ 4.08e-01]
[ 3.77e-01]

>>> # p=-3
>>> print(Pminus3_SDP)
Optimization Problem
  minimize t
  over
    1×1 real variable t
    8×1 real variable w (bounded below)
  subject to
    ∑(w) ≤ 1
    tr(∑(w[i]·A[i]·A[i]ᵀ : i ∈ [0…7])^(-3)) ≤ t
>>> solution = Pminus3_SDP.solve(solver='cvxopt')
>>> print(w)# doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
[ ...]
[ ...]
[ 2.48e-01]
[ 1.66e-01]
[ 1.08e-01]
[ 1.41e-01]
[ 7.83e-02]
[ 2.60e-01]

.. _optdes_refs:

References
----------

        1. "`Computing Optimal Designs of multiresponse Experiments reduces to
           Second-Order Cone Programming <http://arxiv.org/abs/0912.5467>`_", G. Sagnol,
           *Journal of Statistical Planning and Inference*,
           141(5), p. *1684-1708*, 2011.

        2. "`Computing exact D-optimal designs by mixed integer second order cone
           programming <http://arxiv.org/abs/1307.4953>`_",
           G. Sagnol and R. Harman, Submitted: arXiv:1307.4953.

        3. "`Determinant maximization with linear matrix inequality
           constraints <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.38.7483&rep=rep1&type=pdf>`_",
           L. Vandenberghe, S. Boyd and S.P. Wu, *SIAM journal on matrix analysis and applications*,
           19(2), 499-533, 1998.

        4. "`On the semidefinite representations of real functions applied to symmetric
           matrices <http://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/1751>`_", G. Sagnol,
           *Linear Algebra and its Applications*,
           439(10), p. *2829-2843*, 2013.
.. _duals:

Dual Values
===========

Picos typically reformulates optimization problems as
conic programs of the form

.. math::
   :nowrap:

   \begin{center}
   $\begin{array}{cclc}
   \underset{\mathbf{x} \in \mathbb{R}^n}{\mbox{minimize}}
                      & \mathbf{c}^T \mathbf{x} + \gamma & &\\
   \mbox{subject to} & A_i(\mathbf{x}) & \succeq_{K_i} \mathbf{b}_i,\ \forall i \in I,
   \end{array}$
   \end{center}

where each :math:`A_i` is a linear map from :math:`\mathbb{R}^n` to a linear
space containing the cone :math:`K_i`, and the generalized conic inequality
:math:`\mathbf{x} \succeq_K \mathbf{y}` means :math:`\mathbf{x}-\mathbf{y}\in K`
for a cone :math:`K`. For the sake of compactness, we allow generalized
inequalities over the trivial cone :math:`K_{eq} = \{\mathbf{0}\}`, such that
:math:`A \mathbf{x} \succeq_{K_{eq}} \mathbf{b}` represents an equality
constraint :math:`A \mathbf{x} = \mathbf{b}`.


The dual conic problem can be written as follows:

.. math::
   :nowrap:

   \begin{center}
   $\begin{array}{cll}
   \mbox{maximize}   & \sum_{i\in I} \mathbf{b}_i^T \mathbf{y}_i + \gamma \\
   \mbox{subject to} & \sum_{i\in I} A_i^*(\mathbf{y}_i) = \mathbf{c}, \\
                     & \mathbf{y}_i \succeq_{K_i^*} 0,\ \forall i \in I,
   \end{array}$
   \end{center}

where :math:`A^*` denotes the adjoint operator of :math:`A` and :math:`K^*`
denotes the the dual cone of :math:`K` (see the note below for a list of cones
that are supported in PICOS, together with their dual).

After an optimization problem has been solved, we can query the optimal dual
variable :math:`y_i \in K_i^*` of a conic constraint ``con`` over the cone
:math:`K_i` with its :attr:`~.constraint.Constraint.dual` attribute, i.e.,
``con.dual``.

When an optimization problem ``P`` can be reformulated to a conic program ``C``
of the above form by PICOS, we can use its :meth:`~.problem.Problem.dual`
attribute to return a :class:`~.problem.Problem` object ``D=P.dual`` which
contains the dual conic program of ``C``. It is also possible to solve ``P`` via
its dual by using the :ref:`dualize <option_dualize>` option: This passes
problem ``D`` to the solver, and the optimal primal and dual variables of ``P``
will be retrieved from the optimal solution of ``D``.


Supported cones and their dual
------------------------------

PICOS can provide dual information for problems involving the following cones:


.. rubric:: Trivial cone

The trivial cone :math:`K_{eq} = \{\mathbf{0}\}\subset \mathbb{R}^n`,
whose dual cone is the entire space :math:`K_{eq}^* = \mathbb{R}^n`.
This means that the dual variable :math:`\mathbf{y}` for an equality constraint
is unconstrained.

.. rubric:: Nonnegative Orthant

The nonnegative orthant :math:`\mathbb{R}_+^n` is self dual:
:math:`(\mathbb{R}_+^n)^* = \mathbb{R}_+^n`. Therefore the dual variable for a
set of linear inequalities is a vector :math:`\mathbf{y}\geq\mathbf{0}`.

.. rubric:: Lorentz Cone

The :ref:`Lorentz cone <lorentz>` :math:`\mathcal{Q}^n=` :math:`\{(t,\mathbf{x})
\in \mathbb{R}\times \mathbb{R}^{n-1}: \|\mathbf{x}\| \leq t \}`, which is used
to model second-order cone inequalities, is self-dual: :math:`(\mathcal{Q}^n)^*
= \mathcal{Q}^n`. This means that the dual variable for a second order cone
inequality of the form

.. math::
    :nowrap:

    $\| A \mathbf{x} - \mathbf{b} \| \leq \mathbf{h}^T \mathbf{x} - g
    \iff
    \left[
        \begin{array}{c} \mathbf{h}^T\\ A \end{array}
    \right] \mathbf{x}
    \succeq_{\mathcal{Q}^n}
    \left[
        \begin{array}{c} g\\ \mathbf{b} \end{array}
    \right]$

is a vector of the form :math:`[\lambda, \mathbf{z}^T]^T` such that
:math:`\|\mathbf{z}\| \leq \lambda`.

.. rubric:: Rotated Second-order Cone

The (widened or narrowed) :ref:`rotated second order cone <rotatedcone>` is

.. math::
    :nowrap:

    $\mathcal{R}_p^n =\{(u,v,\mathbf{x})\in\mathbb{R}\times\mathbb{R}\times\mathbb{R}^{n-2}:
    \|\mathbf{x}\|^2 \leq p\cdot u \cdot v,\ u,v\geq 0 \}$

for some :math:`p>0`, and its dual cone is :math:`(\mathcal{R}_{p}^n)^* =
\mathcal{R}_{4/p}^n`. In particular, :math:`\mathcal{R}_p^n` is self-dual for
:math:`p=2`. For example, the dual variable for the constraint :math:`\| A
\mathbf{x} - \mathbf{b} \|^2 \leq (\mathbf{h}^T \mathbf{x} - g)(\mathbf{e}^T
\mathbf{x} - f)` with :math:`(\mathbf{h}^T \mathbf{x} - g)\geq 0` and
:math:`(\mathbf{e}^T \mathbf{x} - f)\geq 0`, i.e.,

.. math::
    :nowrap:

    $
    \left[
        \begin{array}{c} \mathbf{h}^T\\ \mathbf{e}^T\\ A \end{array}
    \right] \mathbf{x}
    \succeq_{\mathcal{R}_1^n}
    \left[
        \begin{array}{c} g\\ f\\ \mathbf{b} \end{array}
    \right]$

is a vector of the form :math:`[\alpha, \beta, \mathbf{z}^T]^T` such that
:math:`\|\mathbf{z}\|^2 \leq 4 \alpha \beta`;

.. rubric:: Positive Semi-definite Cone

The positive semidefinite cone :math:`\mathbb{S}_+^n` is self dual:
:math:`(\mathbb{S}_+^n)^* = \mathbb{S}_+^n`. This means that the dual variable
for a linear matrix inequality :math:`\sum_i x_i M_i \succeq M_0` is a positive
semidefinite matrix :math:`Y \succeq 0`;

.. rubric:: Exponential Cone

PICOS can also reformulate several constraints using the *exponential cone*
:class:`~picos.expressions.ExponentialCone`, as it is the case for example for
:class:`~picos.constraints.KullbackLeiblerConstraint`. PICOS provides dual
values for :class:`~picos.constraints.ExpConeConstraint`, as computed by the
solver, but dualization of those constraints is not yet supported.
.. _examples:

Examples
========

The :ref:`quick examples <quick_examples>` are all self-contained and can be
copy-pasted to a source file or Python console to reproduce them. Most of the
remaining examples have a tutorial character and are presented in multiple
dependent code sections.

.. toctree::
   :maxdepth: 2

   quick
   graphs
   complex
   optdes
   constraints
.. _userguide:

Usage Notes
===========

.. toctree::
   :maxdepth: 2

   numscipy.rst
   slicing.rst
   duals.rst
   tolerances.rst
   cheatsheet.rst
.. _quick_examples:

Quick examples
==============

The short examples below are all self-contained and can be copied to a Python
source file or pasted into a Python console.


Projection onto a convex hull
-----------------------------

We solve the problem

.. math::

  \underset{x \in \mathbb{R}^n}{\text{minimize}}\quad&\lVert Ax - b \rVert \\
  \text{subject to}\quad&\sum_{i=1}^n x_i = 1, \\
  &x \succeq 0,

which asks for the projection :math:`Ax` of the point :math:`b \in \mathbb{R}^m`
onto the convex hull of the columns of :math:`A \in \mathbb{R}^{m \times n}`:

.. plot::
  :include-source:

  #!/usr/bin/env python3

  import numpy as np
  import picos as pc
  from matplotlib import pyplot
  from scipy import spatial

  # Make the result reproducible.
  np.random.seed(12)

  # Define the data.
  n = 20
  A = np.random.rand(2, n)
  b = np.array([1, 0])

  # Define the decision variable.
  x = pc.RealVariable("x", n)

  # Define and solve the problem.
  P = pc.Problem()
  P.minimize = abs(A*x - b)
  P += pc.sum(x) == 1, x >= 0
  P.solve(solver="cvxopt")

  # Obtain the projection point.
  p = (A*x).np

  # Plot the results.
  V = spatial.ConvexHull(A.T).vertices
  figure = pyplot.figure(figsize=(8.7, 4))
  figure.gca().set_aspect("equal")
  pyplot.axis("off")
  pyplot.fill(A.T[V, 0], A.T[V, 1], "lightgray")
  pyplot.plot(A.T[:, 0], A.T[:, 1], "k.")
  pyplot.plot(*zip(b, p), "k.--")
  pyplot.annotate("$\mathrm{conv} \{a_1, \ldots, a_n\}$", [0.25, 0.5])
  pyplot.annotate("$b$", b + 1/100)
  pyplot.annotate("$Ax$", p + 1/100)
  pyplot.tight_layout()
  pyplot.show()

.. rubric:: Example notes

- The Python builtin function :func:`abs` (absolute value) is understood as the
  default norm. For real vectors, this is the Euclidean norm.
- The attribute :attr:`~picos.valuable.Valuable.np` returns the value of a PICOS
  expression as a NumPy type.
- The choice of the CVXOPT solver is optional. Explicit solver choice is made
  throughout the documentation to make its automatic validation more reliable.


Worst-case projection
---------------------

We solve the same problem as before but now we assume that the point :math:`b`
to be projected is only known to live inside an ellipsoid around its original
location. In this case we cannot hope to obtain an exact projection but we may
compute a point :math:`p` on the convex hull of the columns of :math:`A` that
minimizes the worst-case distance to :math:`b`. This approach is known as
`robust optimization <https://en.wikipedia.org/wiki/Robust_optimization>`_.
Formally, we solve the min-max problem

.. math::

  \underset{x \in \mathbb{R}^n}{\text{minimize}}\quad&\max_{\theta \in
    \Theta}~\lVert Ax - (b + \theta) \rVert \\
  \text{subject to}\quad&\sum_{i=1}^n x_i = 1, \\
  &x \succeq 0, \\

where :math:`\Theta = \{\theta \mid L\theta \leq 1\}` is an ellipsoidal
*perturbation set* (for some invertible matrix :math:`L`):

.. plot::
  :include-source:

  #!/usr/bin/env python3

  import numpy as np
  import picos as pc
  from matplotlib import pyplot
  from matplotlib.patches import Ellipse
  from scipy import spatial

  # Make the result reproducible.
  np.random.seed(12)

  # Define the data.
  n = 20
  A = np.random.rand(2, n)
  b = np.array([1, 0])

  # Define an ellipsoidal uncertainty set Θ and a perturbation parameter θ.
  # The perturbation is later added to the data, rendering it uncertain.
  Theta = pc.uncertain.ConicPerturbationSet("θ", 2)
  Theta.bound(  # Let ‖Lθ‖ ≤ 1.
    abs([[ 5,  0],
         [ 0, 10]] * Theta.element) <= 1
  )
  theta = Theta.compile()

  # Define the decision variable.
  x = pc.RealVariable("x", n)

  # Define and solve the problem.
  P = pc.Problem()
  P.minimize = abs(A*x - (b + theta))
  P += pc.sum(x) == 1, x >= 0
  P.solve(solver="cvxopt")

  # Obtain the projection point.
  p = (A*x).np

  # Plot the results.
  V = spatial.ConvexHull(A.T).vertices
  E = Ellipse(b, 0.4, 0.2, color="lightgray", ec="k", ls="--")
  figure = pyplot.figure(figsize=(8.7, 4))
  axis = figure.gca()
  axis.add_artist(E)
  axis.set_aspect("equal")
  axis.set_xlim(0.5, 1.21)
  axis.set_ylim(-0.11, 0.5)
  pyplot.axis("off")
  pyplot.fill(A.T[V, 0], A.T[V, 1], "lightgray")
  pyplot.plot(A.T[:, 0], A.T[:, 1], "k.")
  pyplot.plot(*zip(b, p), "k.")
  pyplot.annotate("$\mathrm{conv} \{a_1, \ldots, a_n\}$", [0.25, 0.5])
  pyplot.annotate("$b$", b + 1/200)
  pyplot.annotate("$Ax$", p + 1/200)
  pyplot.tight_layout()
  pyplot.show()

.. rubric:: Example notes

- One can also scale and shift the parameter obtained from a
  :class:`~picos.uncertain.UnitBallPerturbationSet` to obtain ellipsoidal
  uncertainty. Its parent class :class:`~picos.uncertain.ConicPerturbationSet`
  that we showcased is more versatile and can represent any conically bounded
  perturbation set through repeated use of its
  :meth:`~picos.expressions.uncertain.pert_conic.ConicPerturbationSet.bound`
  method.
- A report of the robust and distributionally robust optimization models
  supported by PICOS and their mathematical background is found in
  :ref:`[1] <quick_refs>`.


Optimal Minecraft mob farm
--------------------------

Minecraft is a popular sandbox video game in which some players aim to build
efficient automated factories, referred to as *farms*. One type of farm waits
for hostile creatures (*mobs*) to appear on a platform, then pushes them off the
platform with a water dispenser in the center to collect any valuables that they
might carry. Such a farm is threatened by the possibility of Spiders to appear,
which are too large for the collection mechanism to handle. Fortunately, the
Spider requires a :math:`3 \times 3` area to spawn on while the other mobs
require just a single free :math:`1 \times 1` cell, so Spider spawns can be
prevented by blocking off some of the platform's cells.

In the following we compute an optimal platform that maximizes the number of
cells that mobs can spawn on while admitting no :math:`3 \times 3` spawnable
region for Spiders. We further compute an optimal highly symmetric (w.r.t. both
axes and diagonals) solution for those who value looks over efficiency:

.. plot::
  :include-source:

  #!/usr/bin/env python3

  import picos as pc
  from matplotlib import colors, pyplot

  # Represent the spawning platform by a 15×15 binary matrix variable S where a
  # one represents a spawnable field and a zero one that is not spawnable.
  S = pc.BinaryVariable("S", (15, 15))

  # Maximize the number of spawnable blocks.
  P = pc.Problem("Optimal Mob Farm")
  P.maximize = pc.sum(S)

  # The actual platform is shaped like a diamond of cells with taxicab distance
  # of at most seven from the center block. Mark all other cells not spawnable.
  P += [
      S[x, y] == 0
      for x in range(S.shape[0])
      for y in range(S.shape[1])
      if abs(x - 7) + abs(y - 7) > 7
  ]

  # The center block is not spawnable due to the water dispenser.
  P += S[7, 7] == 0

  # Additionally, we require that there is no 3x3 spawnable area.
  P += [
      sum([
          S[a, b]
          for a in range(x - 1, x + 2)
          for b in range(y - 1, y + 2)
      ]) <= 8
      for x in range(1, S.shape[0] - 1)
      for y in range(1, S.shape[1] - 1)
  ]

  # Solve the problem and store the optimal platform.
  P.solve(solver="glpk")
  S_opt = S.np

  # Now modify the problem to require a more symmetric solution.
  P += [S[x, :] == S[14 - x, :] for x in range(S.shape[0] // 2)]  # Vertical.
  P += S == S.T  # Diagonal.

  # Re-solve the updated problem.
  P.solve()
  S_sym = S.np

  # Display both solutions.
  figure, axes = pyplot.subplots(ncols=2, figsize=(8.7, 5))
  titles = ["An optimal platform", "An optimal symmetric platform"]
  cmap = colors.ListedColormap(["#1c1c1c", "#78ae00", "#d35e1a"])

  for axis, title, solution in zip(axes, titles, [S_opt, S_sym]):
    solution[7, 7] = 2  # Mark the center.
    axis.axis("off")
    axis.set_title(title)
    axis.pcolormesh(solution, edgecolor="#2f2f2f", linewidth=0.5, cmap=cmap)

  pyplot.tight_layout()
  pyplot.show()

.. rubric:: Example notes

- Excluding the center, the platform has 112 cells. The solutions show that an
  optimal platform has 9 obstacles and 103 free cells (92.0%) while an optimal
  symmetric platform has 12 obstacles and thus only 100 free cells (89.3%).
- The two symmetry conditions require symmetry along one axis and one main
  diagonal, respectively. Symmetry along the remaining axis and diagonal is
  obtained implicitly. With an adjustment it can be seen that only requiring
  axial symmetry does not increase efficiency.


.. _quick_refs:

References
----------

  1. "`Robust conic optimization in Python
     <https://www.static.tu.berlin/fileadmin/www/10005693/Publications/Stahlberg20.pdf>`_",
     M. Stahlberg, Master's thesis, 2020.
.. _welcome:

.. include:: badges.rst

A Python interface to conic optimization solvers
================================================

|gitlab| • |pypi| |anaconda| |aur| • |license| |cov|

Welcome to the documentation of PICOS, a powerful and user friendly Python API
for convex and mixed integer optimization that dispatches your problem to the
best fit solver that is available at runtime. A `PDF version <picos.pdf>`_ of
this documentation is available for offline use. Here's a quick example:

>>> import picos as pc
>>> x = pc.RealVariable("x", 5)
>>> a = pc.Constant("a", range(5))
>>> P = pc.Problem()
>>> P.minimize = abs(x - a)                            # abs() - Euclidean norm
>>> P += pc.sum(x) == 1                                # Add a constraint
>>> opt = P.solve(solver="cvxopt")                     # Optional: Solver choice
>>> print(x.T)                                         # .T - Transpose
[-1.80e+00 -8.00e-01  2.00e-01  1.20e+00  2.20e+00]
>>> round(P.value, 3)
4.025

.. _quickstart:

Quickstart guide
----------------

- If you are **new to PICOS**, head to the :ref:`introduction
  <introduction>`, the :ref:`tutorial <tutorial>`, or see our :ref:`examples
  <examples>`.
- As an **experienced user**, check out the :ref:`changelog <changelog>` or dive
  into the :ref:`API documentation <api>`.
- If you want to report a **bug** or **contribute to PICOS**, the
  :ref:`contribution guide <contributing>` has you covered.
- If you still have a **question**, we're happy to receive
  `your mail <incoming+picos-api/picos@incoming.gitlab.com>`_!


.. _contents:

Documentation outline
---------------------

.. toctree::
   :maxdepth: 1

   Introduction <introduction>
   tutorial
   examples
   notes
   api
   changelog
   contributing
.. _introduction:

.. include:: ../README.rst
.. warning::

    This part of the documentation has not been touched for a while. It might
    be incomplete, reference deprecated functions or make a claim that does not
    apply to the latest version of PICOS any more. On the bright side, code
    listings are validated and still work. Watch your step!


.. _cheatsheet:

Cheat Sheet
===========

Manipulate expressions
----------------------

+--------------+--------------------------------+
| **Operator** | **Interpretation**             |
+==============+================================+
|    ``+``     | addition                       |
+--------------+--------------------------------+
|    ``+=``    | inplace addition               |
+--------------+--------------------------------+
|    ``-``     | substraction                   |
+--------------+--------------------------------+
|    ``*``     | multiplication                 |
+--------------+--------------------------------+
|    ``^``     | Hadamard (elementwise) product |
+--------------+--------------------------------+
|    ``@``     | Kronecker product              |
+--------------+--------------------------------+
|    ``|``     | scalar product                 |
+--------------+--------------------------------+
|    ``/``     | division                       |
+--------------+--------------------------------+
|    ``**``    | exponentiation                 |
+--------------+--------------------------------+
|    ``abs()`` | Euclidean (or Frobenius) norm  |
+--------------+--------------------------------+
|    ``[]``    | slicing                        |
+--------------+--------------------------------+
|    ``&``     | horizontal concatenation       |
+--------------+--------------------------------+
|    ``//``    | vertical concatenation         |
+--------------+--------------------------------+
|    ``.T``    | transposition                  |
+--------------+--------------------------------+
|    ``.H``    | Hermitian transposition        |
+--------------+--------------------------------+
|    ``.Tx``   | partial transposition          |
+--------------+--------------------------------+
|    ``.conj`` | complex conjugate              |
+--------------+--------------------------------+
|    ``.real`` | real part                      |
+--------------+--------------------------------+
|    ``.imag`` | imaginary part                 |
+--------------+--------------------------------+

Create constraints
------------------

+-----------------+-----------------------------------+
| **Operator**    | **Interpretation**                |
+=================+===================================+
| ``<`` or ``<=`` | less or equal                     |
+-----------------+-----------------------------------+
| ``>`` or ``>=`` | larger or equal                   |
+-----------------+-----------------------------------+
| ``==``          | equal                             |
+-----------------+-----------------------------------+
| ``<<``          | Löwner ordering  :math:`\preceq`, |
|                 | or set membership  :math:`\in`    |
+-----------------+-----------------------------------+
| ``>>``          | Löwner ordering  :math:`\succeq`, |
|                 | or set membership   :math:`\ni`   |
+-----------------+-----------------------------------+

Create affine expressions
-------------------------

+---------------------------------+-------------------------------------------+
| **function**                    | **short doc**                             |
+=================================+===========================================+
|:func:`~picos.sum`               | sums a list of affine expressions         |
+---------------------------------+-------------------------------------------+
|:func:`~picos.diag`              | diagonal matrix defined by its diagonal   |
+---------------------------------+-------------------------------------------+
|:func:`~picos.diag_vect`         | vector of diagonal elements of a matrix   |
+---------------------------------+-------------------------------------------+
|:func:`~picos.new_param`         | constant affine expression                |
+---------------------------------+-------------------------------------------+
|:func:`~picos.trace`             | trace of a square affine expression       |
+---------------------------------+-------------------------------------------+
|:func:`~picos.partial_transpose` | partial transposition                     |
+---------------------------------+-------------------------------------------+
|:func:`~picos.partial_trace`     | partial trace                             |
+---------------------------------+-------------------------------------------+

Create convex expressions
-------------------------

+-------------------------------------+-----------------------------------------+
| **function**                        | **short doc**                           |
+=====================================+=========================================+
|:func:`~picos.geomean`               | geometric mean                          |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.norm`                  | (generalized) :math:`L_p-` norm         |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.tracepow`              | trace of a *p*-th matrix power          |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.detrootn`              | *n*-th root of determinant              |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.sum_k_largest`         | sum of k largest elements               |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.sum_k_smallest`        | sum of k smallest elements              |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.sum_k_largest_lambda`  | sum of k largest eigenvalues            |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.sum_k_smallest_lambda` | sum of k smallest eigenvalues           |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.lambda_max`            | largest eigenvalue                      |
+-------------------------------------+-----------------------------------------+
|:func:`~picos.lambda_min`            | smallest eigenvalue                     |
+-------------------------------------+-----------------------------------------+

Create sets
-----------

+-------------------------------------------------------+--------------------------------------------------------------------+
| **function**                                          | **short doc**                                                      |
+=======================================================+====================================================================+
|:func:`ball(r,p) <picos.ball>`                         | a :math:`L_p`- ball of radius ``r``                                |
+-------------------------------------------------------+--------------------------------------------------------------------+
|:func:`simplex(a) <picos.simplex>`                     | a standard simplex                                                 |
|                                                       | :math:`\{x\geq 0: \Vert x \Vert_1 \leq a \}`                       |
+-------------------------------------------------------+--------------------------------------------------------------------+
|:func:`truncated_simplex(a) <picos.truncated_simplex>` | a set of the form                                                  |
|                                                       | :math:`\{ 0\leq x\leq 1: \Vert x \Vert_1 \leq a\}`, or             |
|                                                       | :math:`\{x: \Vert x \Vert_\infty \leq 1; \Vert x \Vert_1 \leq a\}` |
+-------------------------------------------------------+--------------------------------------------------------------------+

Get information on a problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+--------------------------------------------------------------------------------------------+-------------------------------------------+
| **function**                                                                               | **short doc**                             |
+============================================================================================+===========================================+
|:meth:`get_variable(name) <.problem.Problem.get_variable>`                                  | gets the variable object ``name``         |
+--------------------------------------------------------------------------------------------+-------------------------------------------+
|:meth:`get_valued_variable(name) <.problem.Problem.get_valued_variable>`                    | gets the value of the variable ``name``   |
+--------------------------------------------------------------------------------------------+-------------------------------------------+
|:meth:`check_current_value_feasibility() <.problem.Problem.check_current_value_feasibility>`| are the current variable value feasible?  |
+--------------------------------------------------------------------------------------------+-------------------------------------------+
|:meth:`obj_value() <.problem.Problem.obj_value>`                                            | objective for the current variable values |
+--------------------------------------------------------------------------------------------+-------------------------------------------+
|:attr:`.type <.problem.Problem.type>`                                                       | returns problem's type                    |
+--------------------------------------------------------------------------------------------+-------------------------------------------+

Miscellaneous
~~~~~~~~~~~~~

+---------------------------------------------------------+-------------------------------------------+
| **function**                                            | **short doc**                             |
+=========================================================+===========================================+
|:func:`available_solvers() <picos.available_solvers>`    | lists installed solvers                   |
+---------------------------------------------------------+-------------------------------------------+
|:func:`import_cbf() <picos.import_cbf>`                  | imports data from a .cbf file             |
+---------------------------------------------------------+-------------------------------------------+
|:meth:`write_to_file() <.problem.Problem.write_to_file>` | writes problem to a file                  |
+---------------------------------------------------------+-------------------------------------------+
.. |_| unicode:: 0xA0
   :trim:

.. _tolerances:

Numeric Tolerances
==================

PICOS allows you to fine-tune how accurate your solution needs to be.
Tolerances fall in three categories:

- **Feasibility** tolerances, abbreviated ``fsb``, control the magnitude of
  constraint violation that is tolerated. The :ref:`integrality tolerance
  <option_integrality_tol>` also falls into this category.
- **Optimality** tolerances, abbreviated ``opt``, control the maximum allowed
  deviation from the mathematically exact optimum solution and serve as a
  termination criterion. An exception is the the Simplex algorithm that uses the
  :ref:`dual feasibility <option_abs_dual_fsb_tol>` as its stopping criterion.
- The remaining tolerances are used at intermediate steps of specific
  algorithms, such as the :ref:`Markowitz threshold <option_markowitz_tol>` used
  in a pivoting strategy of the Simplex algoritm.

Solvers differ in how they measure deviations from the ideal values. Some bound
**absolute** values while others consider the deviation **in relation** to the
magnitude of the numbers that occur in the problem.
PICOS abbreviates the former measurement with ``abs`` and the latter with
``rel``.
If both measurements are supported by a solver, then the standard approach is to
allow values if they are sufficiently accurate according to either one.

If solvers use a single value for **primal** and **dual** feasibility but PICOS
is configured to use differing accuracies, supplied in the options with the
``prim`` and ``dual`` abbreviations respectively, it will supply the smaller of
both values to such solvers.

By default, PICOS overrides the solver's default accuracies with common values,
so that the choice of solver becomes transparent to you.
Given that ``P`` is your problem instance, you can make PICOS respect the
solvers' individual choices as follows:

>>> import picos
>>> P = picos.Problem()
>>> P.options["*_tol"] = None

Comparison Table
----------------

The table shows what tolerance :class:`options <picos.Options>` are supported by
PICOS and each solver, and what their respective default value is.

.. list-table::
  :header-rows: 1

  * - Option |_| |_| |_| |_| |_| |_| |_| |_| |_| |_| |_| |_| |_| |_| |_| |_| |_|
    - PICOS |_| |_| |_| |_| |_|
    - CPLEX |_| |_| |_| |_| |_|
    - CVXOPT |_|
    - ECOS |_| |_| |_| |_| |_| |_| |_|
    - GLPK |_| |_| |_| |_| |_| |_| |_|
    - Gurobi |_| |_| |_| |_| |_| |_|
    - MOSEK |_| |_| |_| |_|
    - SCIP |_| |_| |_| |_| |_| |_| |_| |_|
    - SMCP |_| |_| |_| |_| |_| |_|

  * - :ref:`abs_prim_fsb_tol <option_abs_prim_fsb_tol>`
    - :math:`10^{-8}`
    - :mathlink:`\text{SX:}~10^{-6} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpRHS.html>`
    - unused
    - unused ?
    - :mathlink:`\text{SX:}~10^{-7}~? <https://fossies.org/linux/glpk/doc/glpk.pdf>`
    - :mathlink:`10^{-6}~? <https://www.gurobi.com/documentation/8.1/refman/feasibilitytol.html#parameter:FeasibilityTol>`
    - :mathlink:`\text{SX:}~10^{-6} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.basis_tol_x>`
      :mathlink:`\text{LP:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_tol_pfeas>`
      :mathlink:`\text{CP:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_co_tol_pfeas>`
      :mathlink:`\text{QP:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_qo_tol_pfeas>`
      :mathlink:`\text{NL:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_nl_tol_pfeas>`
      :mathlink:`\text{IP:}~10^{-6}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.mio_tol_feas>`
    - :mathlink:`\text{SX:}~10^{-6} <https://scip.zib.de/doc/html/PARAMETERS.php>`
    - unused

  * - :ref:`rel_prim_fsb_tol <option_rel_prim_fsb_tol>`
    - :math:`10^{-8}`
    - :mathlink:`\text{LQ:}~10^{-8} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarEpComp.html>`
      :mathlink:`\text{QC:}~10^{-8} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarQCPEpComp.html>`
    - :mathlink:`\text{CP:}~10^{-7} <https://cvxopt.org/userguide/coneprog.html#algorithm-parameters>`
      :mathlink:`\text{NL:}~10^{-7} <https://cvxopt.org/userguide/solvers.html#algorithm-parameters>`
    - :mathlink:`10^{-8}~? <https://github.com/embotech/ecos/blob/develop/include/ecos.h>`
    - unused ?
    - unused ?
    - .. MOSEK :mathlink:` <>`
    - :mathlink:`10^{-6} <https://scip.zib.de/doc/html/FAQ.php#feasibilitycomparison>`
    - :mathlink:`10^{-8} <https://smcp.readthedocs.io/en/latest/documentation/#smcp.solvers.chordalsolver_esd>`

  * - :ref:`abs_dual_fsb_tol <option_abs_dual_fsb_tol>`
    - :math:`10^{-8}`
    - :mathlink:`\text{SX:}~10^{-6} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpOpt.html>`
    - unused
    - unused ?
    - :mathlink:`\text{SX:}~10^{-7}~? <https://fossies.org/linux/glpk/doc/glpk.pdf>`
    - :mathlink:`10^{-6} <https://www.gurobi.com/documentation/8.1/refman/optimalitytol.html#parameter:OptimalityTol>`
    - :mathlink:`\text{SX:}~10^{-6} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.basis_tol_s>`
      :mathlink:`\text{LP:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_tol_dfeas>`
      :mathlink:`\text{CP:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_co_tol_dfeas>`
      :mathlink:`\text{QP:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_qo_tol_dfeas>`
      :mathlink:`\text{NL:}~10^{-8}~? <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_nl_tol_dfeas>`
    - :mathlink:`\text{SX:}~10^{-7} <https://scip.zib.de/doc/html/PARAMETERS.php>`
    - unused

  * - :ref:`rel_dual_fsb_tol <option_rel_dual_fsb_tol>`
    - :math:`10^{-8}`
    - :mathlink:`\text{LQ:}~10^{-8} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarEpComp.html>`
      :mathlink:`\text{QC:}~10^{-8} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarQCPEpComp.html>`
    - :mathlink:`\text{CP:}~10^{-7} <https://cvxopt.org/userguide/coneprog.html#algorithm-parameters>`
      :mathlink:`\text{NL:}~10^{-7} <https://cvxopt.org/userguide/solvers.html#algorithm-parameters>`
    - :mathlink:`10^{-8}~? <https://github.com/embotech/ecos/blob/develop/include/ecos.h>`
    - unused ?
    - unused
    - :mathlink:`\text{SX:}~10^{-12} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.basis_rel_tol_s>`
    - :mathlink:`10^{-6} <https://scip.zib.de/doc/html/FAQ.php#feasibilitycomparison>`
    - :mathlink:`10^{-8} <https://smcp.readthedocs.io/en/latest/documentation/#smcp.solvers.chordalsolver_esd>`

  * - :ref:`abs_ipm_opt_tol <option_abs_ipm_opt_tol>`
    - :math:`10^{-8}`
    - unused
    - :mathlink:`\text{CP:}~10^{-7} <https://cvxopt.org/userguide/coneprog.html#algorithm-parameters>`
      :mathlink:`\text{NL:}~10^{-7} <https://cvxopt.org/userguide/solvers.html#algorithm-parameters>`
    - :mathlink:`10^{-8} <https://github.com/embotech/ecos/blob/develop/include/ecos.h>`
    - unused
    - unused
    - unused
    - :mathlink:`0 <https://scip.zib.de/doc/html/PARAMETERS.php>`
    - :mathlink:`10^{-6} <https://smcp.readthedocs.io/en/latest/documentation/#smcp.solvers.chordalsolver_esd>`

  * - :ref:`rel_ipm_opt_tol <option_rel_ipm_opt_tol>`
    - :math:`10^{-8}`
    - :mathlink:`\text{LQ:}~10^{-8} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarEpComp.html>`
      :mathlink:`\text{QC:}~10^{-8} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarQCPEpComp.html>`
    - :mathlink:`\text{CP:}~10^{-6} <https://cvxopt.org/userguide/coneprog.html#algorithm-parameters>`
      :mathlink:`\text{NL:}~10^{-6} <https://cvxopt.org/userguide/solvers.html#algorithm-parameters>`
    - :mathlink:`10^{-8} <https://github.com/embotech/ecos/blob/develop/include/ecos.h>`
    - unused
    - :mathlink:`\text{CO:}~10^{-8} <https://www.gurobi.com/documentation/8.1/refman/barconvtol.html#parameter:BarConvTol>`
      :mathlink:`\text{QC:}~10^{-6} <https://www.gurobi.com/documentation/8.1/refman/barqcpconvtol.html#parameter:BarQCPConvTol>`
    - :mathlink:`\text{LP:}~10^{-8} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_tol_rel_gap>`
      :mathlink:`\text{CP:}~10^{-7} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_co_tol_rel_gap>`
      :mathlink:`\text{QP:}~10^{-7} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_qo_tol_rel_gap>`
      :mathlink:`\text{NL:}~10^{-6} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.intpnt_nl_tol_rel_gap>`
    - :mathlink:`0 <https://scip.zib.de/doc/html/PARAMETERS.php>`
    - :mathlink:`10^{-6} <https://smcp.readthedocs.io/en/latest/documentation/#smcp.solvers.chordalsolver_esd>`

  * - :ref:`abs_bnb_opt_tol <option_abs_bnb_opt_tol>`
    - :math:`10^{-6}`
    - :mathlink:`10^{-6} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpAGap.html>`
    - no IP
    - :mathlink:`10^{-6} <https://github.com/embotech/ecos/blob/develop/include/ecos_bb.h#L37>`
    - unused
    - :mathlink:`10^{-10} <https://www.gurobi.com/documentation/8.1/refman/mipgapabs.html#parameter:MIPGapAbs>`
    - :mathlink:`0 <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.mio_tol_abs_gap>`
    - :mathlink:`0 <https://scip.zib.de/doc/html/PARAMETERS.php>`
    - no IP

  * - :ref:`rel_bnb_opt_tol <option_rel_bnb_opt_tol>`
    - :math:`10^{-4}`
    - :mathlink:`10^{-4} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html>`
    - no IP
    - :mathlink:`10^{-3} <https://github.com/embotech/ecos/blob/develop/include/ecos_bb.h#L38>`
    - :mathlink:`0 <https://fossies.org/linux/glpk/doc/glpk.pdf>`
    - :mathlink:`10^{-4} <https://www.gurobi.com/documentation/8.1/refman/mipgap2.html#parameter:MIPGap>`
    - :mathlink:`10^{-4} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.mio_tol_rel_gap>`
    - :mathlink:`0 <https://scip.zib.de/doc/html/PARAMETERS.php>`
    - no IP

  * - :ref:`integrality_tol <option_integrality_tol>`
    - :math:`10^{-5}`
    - :mathlink:`10^{-5} <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpInt.html>`
    - no IP
    - :mathlink:`10^{-4} <https://github.com/embotech/ecos/blob/develop/include/ecos_bb.h#L40>`
    - :mathlink:`10^{-5} <https://fossies.org/linux/glpk/doc/glpk.pdf>`
    - :mathlink:`10^{-5} <https://www.gurobi.com/documentation/8.1/refman/intfeastol.html#parameter:IntFeasTol>`
    - :mathlink:`10^{-5} <https://docs.mosek.com/8.1/pythonapi/parameters.html#mosek.dparam.mio_tol_abs_relax_int>`
    - unused
    - no IP

  * - :ref:`markowitz_tol <option_markowitz_tol>`
    - ``None``
    - :mathlink:`0.01 <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpMrk.html>`
    - no SX
    - no SX
    - :mathlink:`0.1 <https://fossies.org/linux/glpk/doc/glpk.pdf>`
    - :mathlink:`2^{-7} <https://www.gurobi.com/documentation/8.1/refman/markowitztol.html#parameter:MarkowitzTol>`
    - unused ?
    - unused
    - no SX

.. rubric:: Pooled options

- ECOS, CVXOPT, SCIP and SMCP merge :ref:`rel_prim_fsb_tol
  <option_rel_prim_fsb_tol>` and :ref:`rel_dual_fsb_tol
  <option_rel_dual_fsb_tol>`.
- CPLEX merges :ref:`rel_prim_fsb_tol <option_rel_prim_fsb_tol>`,
  :ref:`rel_dual_fsb_tol <option_rel_dual_fsb_tol>` and
  :ref:`rel_ipm_opt_tol <option_rel_ipm_opt_tol>`.
- SCIP appears to merge :ref:`abs_ipm_opt_tol <option_abs_ipm_opt_tol>` with
  :ref:`abs_bnb_opt_tol <option_abs_bnb_opt_tol>` and :ref:`rel_ipm_opt_tol
  <option_rel_ipm_opt_tol>` with :ref:`rel_bnb_opt_tol
  <option_rel_bnb_opt_tol>` with its ``limits/absgap`` and ``limits/gap``
  options, respectively.

.. rubric:: Legend

.. list-table::
  :widths: auto

  * - ?
    - It is unclear whether an absolute or relative measure is used,
      or if an option is not available.
  * - SX
    - Linear Programs via Simplex
  * - LP
    - Linear Programs via Interior-Point Method
  * - CP
    - Conic Programs
  * - LQ
    - Linear and Quadratic Programs
  * - QP
    - Quadratic Programs
  * - QC
    - Quadratically Constrained (Quadratic) Programs
  * - NL
    - Nonlinear Programs
  * - IP
    - (Mixed) Integer Programs
.. TODO: Replace all testcode/testoutput blocks with interactive listings so
..       that test.py can validate the examples.


.. _complex:

Complex Semidefinite Programming
================================

PICOS supports complex semidefinite programming as of version 1.0.1. It was
overhauled in version 2.0 to provide some of the features showcased below.
This extension of semidefinite programming to the complex domain was introduced
by Goemans and Williamson :ref:`[1] <complex_refs>` in order to pose relaxations
of combinatorial optimization problems.
Applications include quantum information theory :ref:`[2] <complex_refs>` and
the phase recovery problem in signal processing :ref:`[3] <complex_refs>`.

Complex problems can be defined in PICOS using the complex-valued variable types
:class:`~picos.ComplexVariable` and :class:`~picos.HermitianVariable`:

>>> from picos import ComplexVariable, HermitianVariable
>>> z = ComplexVariable("z", 4)
>>> H = HermitianVariable("H", 4)
>>> z
<4×1 Complex Variable: z>
>>> H
<4×4 Hermitian Variable: H>
>>> z.real
<4×1 Real Linear Expression: Re(z)>
>>> z.imag
<4×1 Real Linear Expression: Im(z)>

Their value can be set and retrieved as in the real case but may contain an
imaginary part:

>>> z.value = [1, 2+2j, 3+3j, 4j]
>>> z.value  # Note the CVXOPT typecode of 'z'.
<4x1 matrix, tc='z'>
>>> print(z)
[ 1.00e+00-j0.00e+00]
[ 2.00e+00+j2.00e+00]
[ 3.00e+00+j3.00e+00]
[ 0.00e+00+j4.00e+00]
>>> z.real.value
<4x1 matrix, tc='d'>
>>> print(z.real)
[ 1.00e+00]
[ 2.00e+00]
[ 3.00e+00]
[ 0.00e+00]
>>> print(z.imag)
[ 0.00e+00]
[ 2.00e+00]
[ 3.00e+00]
[ 4.00e+00]

Just like real variables are the simplest form of an
:class:`~picos.expressions.AffineExpression`, complex variables are represented
to you as instances of :class:`~picos.expressions.ComplexAffineExpression`.
Most notably this gives access to complex conjugation and hermitian
transposition:

>>> z.conj
<4×1 Complex Linear Expression: conj(z)>
>>> z.H
<1×4 Complex Linear Expression: zᴴ>

Internally complex variables are represented as real variable vectors:

>>> z.dim  # Twice its dimension on the complex field.
8
>>> H.dim  # The same dimension as an arbitrary real matrix of same shape.
16

Note that in the hermitian case, we get away with just :math:`4 \cdot 4 = 16`
*real* scalar variables due to the vectorization used. This leads to a smaller
footprint when the problem is passed to a solver.

Unlike real-valued variables, :class:`~picos.ComplexVariable` and
:class:`~picos.HermitianVariable` do not accept variable bounds at creation, and
any properly complex expression formed from them cannot appear on either side of
an affine inequality constraint or as an objective function. However, PICOS
detects when you supply a real-valued expression in any of these places even if
it was created from complex expressions:

>>> A = ~z*~z.H  # Use the current value of z to create a constant 4×4 matrix.
>>> A
<4×4 Complex Constant: [z]·[zᴴ]>
>>> A.hermitian  # By construction this matrix is hermitian.
True
>>> (H|A)  # Create a complex expression involving H.
<1×1 Complex Linear Expression: ⟨H, [z]·[zᴴ]⟩>
>>> (H|A).isreal  # On closer inspection, it is always real-valued.
True
>>> (H|A).refined  # This means it can be "refined" to a real expression.
<1×1 Real Linear Expression: ⟨H, [z]·[zᴴ]⟩>
>>> (H|A) >= 0  # Refinement happens automatically wherever necessary.
<1×1 Affine Constraint: ⟨H, [z]·[zᴴ]⟩ ≥ 0>
>>> H == A  # Equalities involving complex expressions can be posed as normal.
<4×4 Complex Equality Constraint: H = [z]·[zᴴ]>

Complex linear matrix inequalities are created just as in the real case with the
overloaded ``<<`` and ``>>`` operators representing the Loewner order:

>>> H >> 0
<4×4 Complex LMI Constraint: H ≽ 0>

Since solvers at this time generally do not support complex optimization, PICOS
transforms such a constraint to an equivalent real LMI during solution search.
Only to demonstrate this behavior, we do it manually:

>>> from picos import Options
>>> from picos.constraints import ComplexLMIConstraint
>>> P = ComplexLMIConstraint.RealConversion.convert(H >> 0, Options())
>>> P.get_constraint(0)
<8×8 LMI Constraint: [Re(H), -Im(H); Im(H), Re(H)] ≽ 0>


Fidelity in Quantum Information Theory
--------------------------------------

The material of this section is inspired by a lecture of John Watrous
:ref:`[4] <complex_refs>`.

The fidelity between two (hermitian) positive semidefinite operators :math:`P`
and :math:`Q` is defined as

.. math::
    F(P,Q)
    = \left\Vert P^{\frac{1}{2}} Q^{\frac{1}{2}} \right\Vert_{\text{tr}}
    = \max_U \left|
        \operatorname{trace}\left(P^{\frac{1}{2}} U Q^{\frac{1}{2}}\right)
    \right|,

where the trace norm :math:`\Vert \cdot \Vert_{\text{tr}}` is the sum of the
singular values, and the maximization goes over the set of all unitary matrices
:math:`U`.
This quantity can be expressed as the optimal value of the following
complex-valued SDP:

.. math::
    :nowrap:

    \begin{eqnarray*}
        &\underset{Z \in \mathbb{C}^{n \times n}}{\mbox{maximize}}
        &\frac{1}{2}\operatorname{trace}(Z + Z^*)
    \\
        &\mbox{subject to}
        &\left(\begin{array}{cc}
            P & Z \\
            Z^* & Q
        \end{array}\right) \succeq 0
    \end{eqnarray*}

This model can be implemented in PICOS as follows:

.. testcode::

    import numpy
    import picos

    # Create a positive semidefinite constant P.
    _P = picos.Constant([
        [ 1  -1j,  2  +2j,  1     ],
        [     3j,     -2j, -1  -1j],
        [ 1  +2j, -0.5+1j,  1.5   ]])
    P = (_P*_P.H).renamed("P")

    # Create a positive semidefinite constant Q.
    _Q = picos.Constant([
        [-1  -2j,      2j,  1.5   ],
        [ 1  +2j,     -2j,  2.0-3j],
        [ 1  +2j, -1  +1j,  1  +4j]])
    Q = (_Q*_Q.H).renamed("Q")

    # Define the problem.
    F = picos.Problem()
    Z = picos.ComplexVariable("Z", P.shape)
    F.set_objective("max", 0.5*picos.trace(Z + Z.H))
    F.add_constraint(((P & Z) // (Z.H & Q)) >> 0)

    print(F)

    # Solve the problem.
    F.solve(solver = "cvxopt")

    print("\nOptimal value:", round(F, 4))
    print("Optimal Z:", Z.value, sep="\n")

    # Also compute the fidelity via NumPy for comparison.
    PP  = numpy.matrix(P.value)
    QQ  = numpy.matrix(Q.value)
    S,U = numpy.linalg.eig(PP)
    sqP = U * numpy.diag([s**0.5 for s in S]) * U.H  # Square root of P.
    S,U = numpy.linalg.eig(QQ)
    sqQ = U * numpy.diag([s**0.5 for s in S]) * U.H  # Square root of Q.
    Fnp = sum(numpy.linalg.svd(sqP * sqQ)[1])  # Trace-norm of sqrt(P)·sqrt(Q).

    print("Fidelity F(P,Q) computed by NumPy:", round(Fnp, 4))

.. testoutput::

    Complex Semidefinite Program
      maximize 0.5·tr(Z + Zᴴ)
      over
        3×3 complex variable Z
      subject to
        [P, Z; Zᴴ, Q] ≽ 0

    Optimal value: 39.8938
    Optimal Z:
    [ 1.06e+01+j2.04e+00 -7.21e+00+j5.77e+00  3.58e+00-j8.10e+00]
    [-8.26e+00-j2.13e+00  1.65e+01+j3.61e-01  8.59e-02-j2.29e+00]
    [-1.38e+00+j6.42e+00 -5.65e-01+j1.55e+00  1.28e+01-j2.40e+00]

    Fidelity F(P,Q) computed by NumPy: 39.8938


Phase Recovery in Signal Processing
-----------------------------------

This section is inspired by :ref:`[3] <complex_refs>`.

The goal of the phase recovery problem is to reconstruct the complex phase of a
vector given only the magnitudes of some linear measurements.
This problem can be formulated as a non-convex optimization problem, and the
authors of :ref:`[3] <complex_refs>` have proposed a complex semidefinite
relaxation similar to the well known relaxation of the **Max-Cut Problem**:
Given a linear operator :math:`A` and a vector :math:`b` of measured amplitudes,
define the positive semidefinite hermitian matrix

.. math::
    M = \operatorname{Diag}(b) (I - AA^\dagger) \operatorname{Diag}(b).

The **Phase-Cut Problem** is:

.. math::
    :nowrap:

    \begin{eqnarray*}
        &\underset{U \in \mathbb{H}_n}{\mbox{minimize}}
        &\langle U, M \rangle
    \\
        &\mbox{subject to}
        &\operatorname{diag}(U) = 1
    \\
        &&U \succeq 0
    \end{eqnarray*}

Note that :math:`U` must be hermitian (:math:`U \in \mathbb{H}_n` ).
We obtain an exact solution :math:`u` to the phase recovery problem if
:math:`U = uu^*` has rank one.
Otherwise, the leading singular vector of :math:`U` is used as an approximation.

This problem can be implemented as follows using PICOS:

.. testcode::

    import cvxopt
    import numpy
    import picos

    # Make the output reproducible.
    cvxopt.setseed(1)

    # Generate an arbitrary rank-deficient hermitian matrix M.
    n, rank = 5, 4
    m = cvxopt.normal(n, rank) + 1j*cvxopt.normal(n, rank)
    M = picos.Constant("M", m*m.H)

    # Define the problem.
    P = picos.Problem()
    U = picos.HermitianVariable("U", n)
    P.set_objective("min", (U | M))
    P.add_constraint(picos.maindiag(U) == 1)
    P.add_constraint(U >> 0)

    print(P)

    # Solve the problem.
    P.solve(solver="cvxopt")

    print("\nOptimal U:", U, sep="\n")

    # Determine the rank of U.
    S, V = numpy.linalg.eig(U.value)
    Urnk = len([s for s in S if abs(s) > 1e-6])

    print("\nrank(U) =", Urnk)

.. testoutput::

    Complex Semidefinite Program
      minimize ⟨U, M⟩
      over
        5×5 hermitian variable U
      subject to
        maindiag(U) = [1]
        U ≽ 0

    Optimal U:
    [ 1.00e+00-j0.00e+00  6.31e-01-j7.76e-01 -8.84e-01+j4.68e-01  6.23e-01-j7.82e-01  7.52e-01+j6.59e-01]
    [ 6.31e-01+j7.76e-01  1.00e+00-j0.00e+00 -9.20e-01-j3.91e-01  1.00e+00-j9.69e-03 -3.75e-02+j9.99e-01]
    [-8.84e-01-j4.68e-01 -9.20e-01+j3.91e-01  1.00e+00-j0.00e+00 -9.17e-01+j4.00e-01 -3.56e-01-j9.34e-01]
    [ 6.23e-01+j7.82e-01  1.00e+00+j9.69e-03 -9.17e-01-j4.00e-01  1.00e+00-j0.00e+00 -4.72e-02+j9.99e-01]
    [ 7.52e-01-j6.59e-01 -3.75e-02-j9.99e-01 -3.56e-01+j9.34e-01 -4.72e-02-j9.99e-01  1.00e+00-j0.00e+00]

    rank(U) = 1


.. _complex_refs:

References
----------

    1. "Approximation algorithms for MAX-3-CUT and other problems via complex
       semidefinite programming",
       M.X. Goemans and D. Williamson.
       In Proceedings of the thirty-third annual
       *ACM symposium on Theory of computing*,
       pp. 443-452. ACM, 2001.

    2. "Semidefinite programs for completely bounded norms",
       J. Watrous,
       arXiv preprint 0901.4709, 2009.

    3. "Phase recovery, maxcut and complex semidefinite programming",
       I. Waldspurger, A. d'Aspremont, and S. Mallat.
       *Mathematical Programming*, pp. 1-35, 2012.

    4. "Semidefinite programs for fidelity and optimal measurements",
       J. Watrous.
       In the script of a course on Theory of Quantum Information,
       https://cs.uwaterloo.ca/~watrous/LectureNotes/CS766.Fall2011/08.pdf.
.. TODO: Once #161 is resolved, document how to assess the solution status.
.. TODO: Bring back the commented-out section on writing to file once it works.


.. _tutorial:

Tutorial
========

First of all, let us import PICOS:

>>> import picos


.. rubric:: Output settings

PICOS makes heavy use of unicode symbols to generate pretty output.
If you find that some of these symbols are not available on your terminal, you
can call :func:`~picos.ascii` or :func:`~picos.latin1` to restrict the charset
used:

>>> X = picos.SymmetricVariable("X", 4)       # Create a dummy variable.
>>> print(X >> 0)                             # Default representation of X ≽ 0.
X ≽ 0
>>> picos.latin1()                            # Limit to ISO 8859-1 (Latin-1).
>>> print(X >> 0)
X » 0
>>> picos.ascii()                             # Limit to pure ASCII.
>>> print(X >> 0)
X >> 0

For the sake of this tutorial, we return to the full beauty of unicode:

>>> picos.default_charset()  # The same as picos.unicode().


Variables
---------

Every optimization endeavor starts with variables. As of PICOS 2.0, the
preferred way to create variables is to create an instance of the desired
variable class:

>>> from picos import RealVariable, BinaryVariable
>>> t = RealVariable("t")                     # A scalar.
>>> x = RealVariable("x", 4)                  # A column vector with 4 elements.
>>> Y = RealVariable("Y", (2, 4))             # A 2×4 matrix.
>>> Z = RealVariable("Z", (4, 2))             # A 4×2 matrix.
>>> w = BinaryVariable("w")                   # A binary scalar.

Now, let's inspect these variables:

>>> w
<1×1 Binary Variable: w>
>>> Y
<2×4 Real Variable: Y>
>>> x.shape
(4, 1)
>>> Z.name
'Z'


.. rubric:: Assigning a value

Assigning values to variables is usually the solver's job, but we can do it
manually:

>>> t.value = 2
>>> # In the case of a binary variable, we can only assign a (near) 0 or 1:
>>> w.value = 0.5  # doctest: +NORMALIZE_WHITESPACE
Traceback (most recent call last):
  ...
ValueError: Failed to assign a value to mutable w: Data is not near-binary with
    absolute tolerance 1.0e-04: Largest difference is 5.0e-01.
>>> print(Z)
Z
>>> Z.value = range(8)
>>> print(Z)  # If a variable is valued, prints the value instead.
[ 0.00e+00  4.00e+00]
[ 1.00e+00  5.00e+00]
[ 2.00e+00  6.00e+00]
[ 3.00e+00  7.00e+00]

As you can see from the last example, PICOS uses column-major order when
loading one-dimensional data such as a Python :class:`range` into a matrix.
The documentation of :func:`~picos.expressions.data.load_data` explains PICOS'
data loading concept in greater detail.


Affine expressions
------------------

The fundamental building blocks of optimization models are affine (matrix)
expressions. Each entry of such an expression is simply a linear combination of
any number of scalar variables plus a constant offset. The variable objects that
we have defined above are special cases of affine expression that refer to
themselves via an identity transformation.

We can now use our variables to create more advanced affine expressions, which
are stored as instances of :class:`~picos.expressions.ComplexAffineExpression`
or of its subclass :class:`~picos.expressions.AffineExpression`. For instance,
we may transpose a matrix variable using the suffix ``.T``:

>>> Y
<2×4 Real Variable: Y>
>>> Y.T
<4×2 Real Linear Expression: Yᵀ>

PICOS expression types overload the standard Python math operators so that you
can denote, for instance, the sum of two expressions as follows:

>>> Z + Y.T
<4×2 Real Linear Expression: Z + Yᵀ>

The overloaded operators will convert arbitrary data on the fly:

>>> t + 1
<1×1 Real Affine Expression: t + 1>
>>> x + 1  # The 1 is broadcasted to a 4×1 vector of all ones.
<4×1 Real Affine Expression: x + [1]>


.. rubric:: Constants

Constants are simply affine expressions with no linear part and are more
commonly referred to as *data*. By default, PICOS uses a short dummy string to
represent multidimensional constants, and reshapes them as needed:

>>> Y + [1, -2, 3, -4, 5, -6, 7, 8]           # Load list as a 2×4 matrix.
<2×4 Real Affine Expression: Y + [2×4]>

If you want to give your constant data a meaningful name and fix its shape for
more type safety, you can do this using :func:`~picos.expressions.Constant`:

>>> from picos import Constant
>>> alpha = Constant("α", 23)                 # Load 23 under the name α.
>>> b = Constant("b", range(4))               # Load as a column vector.
>>> C = Constant("C", [1, -2, 3, -4, 5, -6, 7, 8], (2, 4)); C
<2×4 Real Constant: C>
>>> Y + C
<2×4 Real Affine Expression: Y + C>

The data loading semantics of :func:`~picos.expressions.Constant` or when
loading data on the fly are the same as when valuing variables
(:func:`~picos.expressions.data.load_data`). In particular, you can seamlessly
input CVXOPT or NumPy matrices:

>>> import numpy
>>> Y + numpy.array([[1, 2, 3, 4], [5, 6, 7, 8]])
<2×4 Real Affine Expression: Y + [2×4]>


.. _overloads:

Overloaded operators
--------------------

Now that we have some variables (:math:`t`, :math:`x`, :math:`w`, :math:`Y` and
:math:`Z`) and a couple of constant parameters (:math:`\alpha`, :math:`b`,
:math:`C`), let us create some more affine expressions with them:

>>> C.shape, Z.shape                          # Recall the shapes.
((2, 4), (4, 2))
>>> C*Z                                       # Left multiplication.
<2×2 Real Linear Expression: C·Z>
>>> Z*C                                       # Right multiplication.
<4×4 Real Linear Expression: Z·C>
>>> C*Z*C                                     # Left and right multiplication.
<2×4 Real Linear Expression: C·Z·C>
>>> alpha*Y                                   # Scalar multiplication.
<2×4 Real Linear Expression: α·Y>
>>> t/alpha - alpha/2                         # Division and subtraction.
<1×1 Real Affine Expression: t/α - α/2>
>>> (b | x)                                   # Dot product.
<1×1 Real Linear Expression: ⟨b, x⟩>
>>> # Generalized dot product for matrices: ⟨A, B⟩ = tr(A·Bᴴ):
>>> (C | Y)
<1×1 Real Linear Expression: ⟨C, Y⟩>
>>> b^x                                       # Hadamard (element-wise) product.
<4×1 Real Linear Expression: b⊙x>
>>> C@Z                                       # Kronecker product.
<8×8 Real Linear Expression: C⊗Z>


.. rubric:: Slicing

Python slicing notation can be used to extract single elements or submatrices:

>>> Y[0, 1]                                   # Element in 1st row, 2nd column.
<1×1 Real Linear Expression: Y[0,1]>
>>> x[1:3]                                    # 2nd and 3rd element of x.
<2×1 Real Linear Expression: x[1:3]>
>>> x[-1]                                     # Last element of x.
<1×1 Real Linear Expression: x[-1]>
>>> Y[1,:]                                    # 2nd row of Y.
<1×4 Real Linear Expression: Y[1,:]>
>>> C[:, 1:3]*Y[:, -2::-2]                    # Extended slicing with step size.
<2×2 Real Linear Expression: C[:,1:3]·Y[:,-2::-2]>

In the last example, we select only the second and third column of :math:`C` as
well as the columns of :math:`Y` with an even index considered in reverse order.
The full power and notation of slicing is explained in :ref:`slicing`.


.. rubric:: Concatenation

We can also create larger affine expressions by concatenating them vertically
with ``&`` or horizontally with ``//``:

>>> (b & 2*b & x & C.T*C*x) // x.T
<5×4 Real Affine Expression: [b, 2·b, x, Cᵀ·C·x; xᵀ]>

You have to be a little careful when it comes to operator precedence, because
Python has the binding strength of ``&`` and ``//`` built into its grammar with
logical disjunction and integral division in mind. When in doubt, use
parenthesis around your blocks.


.. rubric:: Broadcasting and reshaping

To recall an example we've seen earlier with variables, scalars are broadcasted
to the necessary shape to allow an addition or subtraction to take place:

>>> 5*x - alpha
<4×1 Real Affine Expression: 5·x - [α]>

Note, however, that apart from this simple broadcasting rule, the shape of a
PICOS constant (loaded via :func:`~picos.Constant`) is already fixed. You can't
add a :math:`8 \times 1` vector to a :math:`4 \times 2` matrix:

>>> Z + (x // b)  # doctest: +NORMALIZE_WHITESPACE
Traceback (most recent call last):
  ...
TypeError: Invalid operation BiaffineExpression.__add__(Z, [x; b]):
    The operand shapes of 4×2 and 8×1 do not match.

The reason is simply that PICOS does not know *which* side to reshape. You can
make the example work by being more explicit:

>>> Z + (x // b).reshaped((4, 2))
<4×2 Real Affine Expression: Z + reshaped([x; b], 4×2)>


.. rubric:: Summing multiple expressions

Since affine expressions overload ``+``, you could use Python's :func:`sum` to
add a bunch of them. However, the string representation can become rather long:

>>> # Create a sequence of matrix constants with sensible names:
>>> A = [Constant("A[{}]".format(i), range(i, i + 8), (2, 4)) for i in range(5)]
>>> A[0]
<2×4 Real Constant: A[0]>
>>> sum([A[i]*Z for i in range(5)])
<2×2 Real Linear Expression: A[0]·Z + A[1]·Z + A[2]·Z + A[3]·Z + A[4]·Z>

To obtain a shorter representation, use :func:`picos.sum` instead:

>>> picos.sum([A[i]*Z for i in range(5)])
<2×2 Real Linear Expression: ∑(A[i]·Z : i ∈ [0…4])>

This works for all kinds of expressions and will look hard to find some pattern
in the summands' string descriptions.


Norms and quadratics
--------------------

.. rubric:: Norms

The norm of an affine expression can be expressed using Python's built-in
:func:`abs` function. If :math:`x` is an affine vector, ``abs(x)`` denotes its
Euclidean norm :math:`\sqrt{x^T x}`:

>>> abs(x)
<Euclidean Norm: ‖x‖>

If the affine expression is a matrix, :func:`abs` returns its Frobenius norm
:math:`\Vert M \Vert_F = \sqrt{\operatorname{trace} (M^T M)}`:

>>> abs(Z - 2*C.T)
<Frobenius Norm: ‖Z - 2·Cᵀ‖>

The absolute value of a scalar is expressed in the same way:

>>> abs(t)
<Absolute Value: |t|>

As is the modulus of a complex expression:

>>> t + 1j
<1×1 Complex Affine Expression: t + 1j>
>>> abs(t + 1j)
<Complex Modulus: |t + 1j|>

Additional norms are available through the :class:`~picos.Norm` class.

.. rubric:: Quadratics

Quadratic expressions can be formed in several ways:

>>> abs(x)**2                                 # Squared norm.
<Squared Norm: ‖x‖²>
>>> t**2 - x[1]*x[2] + 2*t - alpha            # Sum involving quadratic terms.
<Quadratic Expression: t² - x[1]·x[2] + 2·t - α>
>>> (x[1] - 2) * (t + 4)                      # Product of affine expressions.
<Quadratic Expression: (x[1] - 2)·(t + 4)>
>>> Y[0,:]*x                                  # Row vector times column vector.
<Quadratic Expression: Y[0,:]·x>
>>> (x + 2 | Z[:,1])                          # Scalar product.
<Quadratic Expression: ⟨x + [2], Z[:,1]⟩>
>>> (t & alpha) * C * x                       # Quadratic form.
<Quadratic Expression: [t, α]·C·x>

Note that there is no natural way to define a vector or matrix of quadratic
expressions. In PICOS, only affine expressions can be multidimensional.


Defining a problem
------------------

Now that we know how to construct affine and quadratic expressions and norms, it
is time to use them as part of an optimization problem:

>>> from picos import Problem
>>> P = Problem()
>>> P.set_objective("min", (t - 5)**2 + 2)
>>> print(P)
Quadratic Program
  minimize (t - 5)² + 2
  over
    1×1 real variable t

Next we'll search a solution for this problem, but first we configure that only
the solver `CVXOPT <https://cvxopt.org/>`_ may be used so that the documentation
examples are reproducible. We can do this by assigning to the problem's
:attr:`~.problem.Problem.options` attribute:

>>> P.options.solver = "cvxopt"

We can now obtain a solution by calling :meth:`~.problem.Problem.solve`:

>>> solution = P.solve()
>>> solution
<feasible primal solution (claimed optimal) from cvxopt>
>>> solution.primals# doctest: +SKIP
{<1×1 Real Variable: t>: [4.999997568104307]}

Unless disabled by passing ``apply_solution=False`` to
:meth:`~.problem.Problem.solve`, the solution is automatically applied to the
variables involved in the problem definition, so that the entire Problem is now
valued:

>>> round(t, 5)
5.0
>>> round(P, 5)
2.0

The Python functions :func:`round`, :class:`int`, :class:`float` and
:class:`complex` are automatically applied to the ``value`` attribute of
variables, expressions and problems.


Setting options
---------------


We've already seen the ``solver`` option used which allows you to take control
over which of the available solvers should be used. You can display all
available options and their default values by printing the
:attr:`~.problem.Problem.options` instance (we've cut some from the output):

>>> print(P.options)# doctest: +ELLIPSIS
Modified solver options:
  solver              = cvxopt (default: None)
<BLANKLINE>
Default solver options:
  ...
  apply_solution      = True
  ...
  verbosity           = 0
  ...

If you want to change an option only for a single solution attempt, you can also
pass it to :meth:`~.problem.Problem.solve` as a keyword argument:

>>> # Solve again but don't apply the result.
>>> solution = P.solve(apply_solution=False)


Constraints
-----------

Constrained optimization is only half the fun without the constraints. PICOS
again provides overloaded operators to define them:

>>> t <= 5
<1×1 Affine Constraint: t ≤ 5>
>>> x[0] == x[-1]
<1×1 Affine Constraint: x[0] = x[-1]>
>>> abs(x)**2 <= t
<Squared Norm Constraint: ‖x‖² ≤ t>
>>> abs(x)**2 >= t
<Nonconvex Quadratic Constraint: ‖x‖² ≥ t>

Unless there are solvers or reformulation strategies that can deal with a
certain nonconvex constraint type, as is the case for the
:math:`\lVert x \rVert^2 \geq t` constranint above, PICOS will raise a
:exc:`TypeError` to let you know that such a constraint is not supported:

>>> abs(x) <= t
<5×1 SOC Constraint: ‖x‖ ≤ t>
>>> abs(x) >= t
Traceback (most recent call last):
  ...
TypeError: Cannot lower-bound a nonconcave norm.

When working with multidimensional affine expressions, the inequality operators
``>=`` and ``<=`` are understood element-wise (or to put it more mathy, they
represent conic inequality with respect to the nonnegative orthant):

>>> Y >= C
<2×4 Affine Constraint: Y ≥ C>

It is possible to define linear matrix inequalities for use in semidefinite
programming with the operators ``>>`` and ``<<`` denoting the Loewner order:

>>> from picos import SymmetricVariable
>>> S = SymmetricVariable("S", 4)
>>> S >> C.T*C
<4×4 LMI Constraint: S ≽ Cᵀ·C>

Other conic inequalities do not have a Python operator of their own, but you can
denote set membership of an affine expression in a cone. To make this possible,
the operator ``<<`` is also overloaded to denote "is element of":

>>> abs(x) <= t            # Recall that this is a second order cone inequality.
<5×1 SOC Constraint: ‖x‖ ≤ t>
>>> t // x << picos.soc()  # We can also write it like this.
<5×1 SOC Constraint: ‖[t; x][1:]‖ ≤ [t; x][0]>

Here :func:`~picos.soc` is a shorthand for :class:`~picos.SecondOrderCone`,
defined as the convex set

.. math::

    \mathcal{Q}^n = \left\{
        x \in \mathbb{R}^n
    ~\middle|~
        x_1 \geq \sqrt{\sum_{i = 2}^n x_i^2}
    \right\}.

Similarly, we can constrain an expression to be in the rotated second order cone

.. math::

    \mathcal{R}_p^n = \left\{
        x \in \mathbb{R}^n
    ~\middle|~
        p x_1 x_2 \geq \sum_{i = 2}^n x_i^2 \land x_1, x_2 \geq 0
    \right\}

parameterized by :math:`p`:

>>> picos.rsoc(p=1) >> x
<4×1 RSOC Constraint: ‖x[2:]‖² ≤ x[0]·x[1] ∧ x[0], x[1] ≥ 0>

Other sets you can use like this include :class:`~picos.Ball`,
:class:`~picos.Simplex` and the :class:`~picos.ExponentialCone`.


Constrained optimization
------------------------

Let's get back to our quadratic program :math:`P`, which we have already solved
to optimality with :math:`t = 5`:

>>> print(P)
Quadratic Program
  minimize (t - 5)² + 2
  over
    1×1 real variable t


.. rubric:: Adding constraints

We can now add the constraints that :math:`t` must be the sum over all elements
of :math:`x` and that every element of :math:`x` may be at most :math:`1`:

>>> Csum = P.add_constraint(t == x.sum)
>>> Cone = P.add_constraint(x <= 1)
>>> print(P)
Quadratic Program
  minimize (t - 5)² + 2
  over
    1×1 real variable t
    4×1 real variable x
  subject to
    t = ∑(x)
    x ≤ [1]

Now let's solve the problem again and see what we get:

>>> P.solve()
<primal feasible solution pair (claimed optimal) from cvxopt>
>>> round(P, 5)
3.0
>>> round(t, 5)
4.0
>>> x.value
<4x1 matrix, tc='d'>
>>> print(x.value)
[ 1.00e+00]
[ 1.00e+00]
[ 1.00e+00]
[ 1.00e+00]
<BLANKLINE>

Note that multidimensional values such as that of :math:`x` are returned as
`CVXOPT matrix types <https://cvxopt.org/userguide/matrices.html>`_.


.. rubric:: Slack and duals

Since our problem has constraints, we now have slack values and a dual solution
as well:

>>> Csum.slack# doctest: +SKIP
-0.0
>>> Csum.dual# doctest: +SKIP
2.000004393989704
>>> print(Cone.slack)# doctest: +SKIP
[ 9.31e-12]
[ 9.31e-12]
[ 9.31e-12]
[ 9.31e-12]
<BLANKLINE>
>>> print(Cone.dual)# doctest: +SKIP
[ 2.00e+00]
[ 2.00e+00]
[ 2.00e+00]
[ 2.00e+00]
<BLANKLINE>

We did not round the values this time, to showcase that solvers don't always
produce exact solutions even if the problem is "easy". The variable :math:`t` is
also not exactly :math:`4`:

>>> t.value# doctest: +SKIP
3.999999999962744

To learn more about dual values, see :ref:`duals`. For controlling the numeric
precision requirements, see :ref:`tolerances`.


.. rubric:: Removing constraints

Let's say we are not happy with our upper bound on :math:`x` and we'd rather
constrain it to be inside a unit simplex. We can remove the former constraint as
follows:

>>> P.remove_constraint(Cone)

Instead of the constraint itself, we could also have supplied its index in the
problem, as constraints remain in the order in which you add them. Now let's add
the new constraint:

>>> Csimplex = P.add_constraint(x << picos.Simplex())
>>> print(P)
Quadratic Program
  minimize (t - 5)² + 2
  over
    1×1 real variable t
    4×1 real variable x
  subject to
    t = ∑(x)
    x ∈ {x ≥ 0 : ∑(x) ≤ 1}

If we solve again we expect :math:`t` to be :math:`1`:

>>> solution = P.solve()
>>> round(t, 5)
1.0

If the selected solver supports this, changes to a problem's constraints and
objective are passed in the form of updates to the solver's internal state which
can make successive solution searches much faster. Unfortunately, CVXOPT is
stateless so we don't get an advantage here.


.. rubric:: Grouping constraints

You can also add and remove constraints as a group. Let's compute four real
numbers between :math:`0` and :math:`1`, represented by :math:`x_1` to
:math:`x_4` (``x[0]`` to ``x[3]``), such that their minimum distance is
maximized:

>>> from pprint import pprint
>>> P.reset()                                 # Reset the problem, keep options.
>>> d = RealVariable("d", 3)                  # A vector of distances.
>>> P.set_objective("max", picos.min(d))      # Maximize the minimum distance.
>>> C1 = P.add_constraint(x[0] >= 0)          # Numbers start at 0.
>>> C2 = P.add_constraint(x[3] <= 1)          # And end at 1.
>>> # Use constraint groups to order the x[i] and map their distance to y:
>>> G1 = P.add_list_of_constraints([x[i - 1] <= x[i] for i in range(4)])
>>> G2 = P.add_list_of_constraints([d[i] == x[i+1] - x[i] for i in range(3)])
>>> pprint(G1)                                # Show the constraints added.
[<1×1 Affine Constraint: x[-1] ≤ x[0]>,
 <1×1 Affine Constraint: x[0] ≤ x[1]>,
 <1×1 Affine Constraint: x[1] ≤ x[2]>,
 <1×1 Affine Constraint: x[2] ≤ x[3]>]
>>> pprint(G2)
[<1×1 Affine Constraint: d[0] = x[1] - x[0]>,
 <1×1 Affine Constraint: d[1] = x[2] - x[1]>,
 <1×1 Affine Constraint: d[2] = x[3] - x[2]>]
>>> print(P)
Optimization Problem
  maximize min(d)
  over
    3×1 real variable d
    4×1 real variable x
  subject to
    x[0] ≥ 0
    x[3] ≤ 1
    x[i-1] ≤ x[i] ∀ i ∈ [0…3]
    d[i] = x[i+1] - x[i] ∀ i ∈ [0…2]

This looks promising and the constraint groups are nicely formatted, let's solve
the problem and see what we get:

>>> P.solve()
<primal feasible solution pair (claimed optimal) from cvxopt>
>>> print(x)# doctest: +SKIP
[ 5.00e-01]
[ 5.00e-01]
[ 5.00e-01]
[ 5.00e-01]
>>> print(d)# doctest: +SKIP
[ 1.88e-11]
[ 1.88e-11]
[ 1.88e-11]

Apparently there is an error! Revisiting our problem definition, it seems the
first constraint in :math:`G_1`, that is ``x[-1] <= x[0]``, was unnecessary and
forces all :math:`x_i` to take the same value. Luckily, we can remove it from
the group by first specifying the group to access (counting single constraints
as groups of size one) and then the constraint to remove from it:

>>> P.remove_constraint((2, 0))          # Remove 1st constraint from 3rd group.
>>> pprint(P.get_constraint((2,)))       # Show the modified 3rd group.
[<1×1 Affine Constraint: x[0] ≤ x[1]>,
 <1×1 Affine Constraint: x[1] ≤ x[2]>,
 <1×1 Affine Constraint: x[2] ≤ x[3]>]

Now it should work:

>>> print(P)
Optimization Problem
  maximize min(d)
  over
    3×1 real variable d
    4×1 real variable x
  subject to
    x[0] ≥ 0
    x[3] ≤ 1
    x[i] ≤ x[i+1] ∀ i ∈ [0…2]
    d[i] = x[i+1] - x[i] ∀ i ∈ [0…2]
>>> _ = P.solve()  # Don't show or save the solution object.
>>> print(x)#  doctest: +ELLIPSIS
[ ...]
[ 3.33e-01]
[ 6.67e-01]
[ 1.00e+00]
>>> print(d)
[ 3.33e-01]
[ 3.33e-01]
[ 3.33e-01]

(If you see an ellipsis `...` in an example that means we've cut out a near-zero
to allow the other values to be validated automatically.)


.. Problem Export
.. --------------
..
.. Lastly, we show how you can export a problem to a file, in this case in the
.. ``.lp`` format:
..
.. >>> P.reset()
.. >>> P.set_objective("min", t)
.. >>> P.add_constraint(x[0] >= 1.5)
.. >>> P.add_constraint(t - x[0] >= 0.7)
.. >>> print(P)
.. -----------------------
.. Linear Program
..   minimize t
..   over
..     1×1 real variable t
..     4×1 real variable x
..   subject to
..     x[0] ≥ 1.5
..     t - x[0] ≥ 0.7
.. -----------------------
.. >>> P.write_to_file(".helloworld.lp")
.. >>> with open(".helloworld.lp", "r") as fp:
.. ...     print(fp.read())
.. ???
.. >>> import os
.. >>> os.unlink(".helloworld.lp")
.. _numscipy:

NumPy and SciPy
===============

As a lightweight computer algebra system, PICOS sits one level above numerics
libraries such as NumPy and SciPy and acts in concert with them. Let's define a
variable and some data:

>>> import picos, numpy, scipy.sparse
>>> x = picos.RealVariable("x", 4)
>>> N = numpy.reshape(range(16), (4, 4))
>>> type(N)
<class 'numpy.ndarray'>
>>> S = scipy.sparse.spdiags(range(4), 0, 4, 4)
>>> type(S)
<class 'scipy.sparse._dia.dia_matrix'>

.. rubric:: Taking input from NumPy or SciPy

PICOS also allows loading of NumPy and SciPy data on the fly, with one caveat to
watch out for:

>>> x.T*N
<1×4 Real Linear Expression: xᵀ·[4×4]>
>>> N*x
<4×1 Real Linear Expression: [4×4]·x>
>>> x.T*S
<1×4 Real Linear Expression: xᵀ·[4×4]>
>>> S*x
Traceback (most recent call last):
    [...]
picos.valuable.NotValued: Mutable x is not valued.

The last command fails as SciPy sparse matrices `do not currently respect the
__array_priority__ attribute <https://github.com/scipy/scipy/issues/4819>`__, so
that SciPy tries to load ``x`` as an array as opposed to conceding the operation
to PICOS like NumPy does. You can fix this behavior as follows:

>>> picos.patch_scipy_array_priority()
>>> S*x
<4×1 Real Linear Expression: [4×4]·x>

Note that this `monkey-patches <https://en.wikipedia.org/wiki/Monkey_patch>`__
SciPy, so that applications importing your code calling
:func:`~picos.valuable.patch_scipy_array_priority` will also see a patched
version of SciPy.

.. rubric:: Returning NumPy or SciPy data as output

PICOS uses CVXOPT as a numerics backend and thus outputs numeric values as
CVXOPT (sparse) matrices or Python scalar types by default:

>>> x.value = range(4)
>>> x.value
<4x1 matrix, tc='d'>
>>> type(x.value)
<class 'cvxopt.base.matrix'>

However, all objects that can be valued, in particular expressions and problem
instances, also offer properties to query that value as a NumPy type, namely
:attr:`~picos.valuable.Valuable.np` and :attr:`~picos.valuable.Valuable.np2d`:

>>> x.np  # Returns a NumPy scalar, 1D, or 2D array.
array([0., 1., 2., 3.])
>>> type(x.np)
<class 'numpy.ndarray'>
>>> x.np.shape
(4,)
>>> x.np2d  # Always returns a 2D array.
array([[0.],
       [1.],
       [2.],
       [3.]])
>>> x.np2d.shape
(4, 1)

For SciPy, the :attr:`~picos.valuable.Valuable.sp` property returns a sparse
matrix whenever the data stored by PICOS internally is sparse and a NumPy 2D
array otherwise:

>>> I = picos.I(3)
>>> print(I)
[ 1.00e+00     0         0    ]
[    0      1.00e+00     0    ]
[    0         0      1.00e+00]
>>> type(I.sp)
<class 'scipy.sparse._csc.csc_matrix'>
>>> J = picos.J(3, 3)
>>> print(J)
[ 1.00e+00  1.00e+00  1.00e+00]
[ 1.00e+00  1.00e+00  1.00e+00]
[ 1.00e+00  1.00e+00  1.00e+00]
>>> type(J.sp)
<class 'numpy.ndarray'>

A full list of methods for returning values in different formats can be found in
the documentation of the :class:`~picos.valuable.Valuable` base class... _contributing:

.. include:: ../CONTRIBUTING.rst
.. _graphs:

Graph flow and cut problems
===========================

The code below initializes the graph used in all the examples of this page.
It should be run prior to any of the codes presented in this page.
The packages `networkx <http://networkx.lanl.gov/index.html>`_
and `matplotlib <http://matplotlib.sourceforge.net>`_ are required.
We use a graph generated by the LCF generator of the networkx package. The graph
and the edge capacities are deterministic, so that you can compare your results.

.. plot::
  :context:
  :include-source:

  import picos as pc
  import networkx as nx
  import pylab
  import random

  # Use a fixed RNG seed so the result is reproducible.
  random.seed(1)

  # Number of nodes.
  N=20

  # Generate a graph using LCF notation.
  G=nx.LCF_graph(N,[1,3,14],5)
  G=nx.DiGraph(G) #edges are bidirected

  # Generate edge capacities.
  c={}
  for e in sorted(G.edges(data=True)):
    capacity = random.randint(1, 20)
    e[2]['capacity'] = capacity
    c[(e[0], e[1])]  = capacity

  # Convert the capacities to a PICOS expression.
  cc=pc.new_param('c',c)

  # Manually set a layout for which the graph is planar.
  pos={
    0:  (0.07, 0.70), 1:  (0.18, 0.78), 2:  (0.26, 0.45), 3:  (0.27, 0.66),
    4:  (0.42, 0.79), 5:  (0.56, 0.95), 6:  (0.60, 0.80), 7:  (0.64, 0.65),
    8:  (0.55, 0.37), 9:  (0.65, 0.30), 10: (0.77, 0.46), 11: (0.83, 0.66),
    12: (0.90, 0.41), 13: (0.70, 0.10), 14: (0.56, 0.16), 15: (0.40, 0.17),
    16: (0.28, 0.05), 17: (0.03, 0.38), 18: (0.01, 0.66), 19: (0.00, 0.95)
  }

  # Set source and sink nodes for flow computation.
  s=16
  t=10

  # Set node colors.
  node_colors=['lightgrey']*N
  node_colors[s]='lightgreen' # Source is green.
  node_colors[t]='lightblue'  # Sink is blue.

  # Define a plotting helper that closes the old and opens a new figure.
  def new_figure():
    try:
      global fig
      pylab.close(fig)
    except NameError:
      pass
    fig=pylab.figure(figsize=(11,8))
    fig.gca().axes.get_xaxis().set_ticks([])
    fig.gca().axes.get_yaxis().set_ticks([])

  # Plot the graph with the edge capacities.
  new_figure()
  nx.draw_networkx(G, pos, node_color=node_colors)
  labels={
    e: '{} | {}'.format(c[(e[0], e[1])], c[(e[1], e[0])])
    for e in G.edges if e[0] < e[1]}
  nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
  pylab.show()

The first number on an edge label denotes the capacity from the node with the
smaller number to the node with the larger number; the second number denotes the
capacity for the other direction. Source and sink that we will use for flow
computations are drawn in green and blue, respectively.

Max-flow (LP)
-------------

Given a directed graph :math:`G(V,E)`, with a capacity :math:`c(e)` on each edge
:math:`e \in E`, a source node :math:`s` and a sink node :math:`t`, the
**max-flow** problem is to find a flow from :math:`s` to :math:`t` of maximum
value. Recall that a flow :math:`s` to :math:`t` is a mapping from :math:`E` to
:math:`\mathbb{R}^+` such that

- the capacity of each edge is respected,
  :math:`\forall e \in E,\ f(e) \leq c(e)`, and

- the flow is conserved at each non-terminal node,
  :math:`\forall n \in V \setminus \{s,t\},\ \sum_{(i,n)\in E} f((i,n)) = \sum_{(n,j)\in E} f((n,j))`.

Its value is defined as the volume passing from :math:`s` to :math:`t`:

.. math::

  \mathrm{value} (f) = \sum_{(s,j)\in E} f((s,j)) - \sum_{(i,s)\in E} f((i,s)) = \sum_{(i,t)\in E} f((i,t)) - \sum_{(t,j)\in E} f((t,j)).

This problem has a linear programming formulation, which we solve below for
``s=16`` and ``t=10``:

.. plot::
  :context:
  :nofigs:
  :include-source:

  maxflow=pc.Problem()

  # Add the flow variables.
  f={}
  for e in G.edges():
    f[e]=maxflow.add_variable('f[{0}]'.format(e))

  # Add another variable for the total flow.
  F=maxflow.add_variable('F')

  # Enforce edge capacities.
  maxflow.add_list_of_constraints([f[e] <= cc[e] for e in G.edges()])

  # Enforce flow conservation.
  maxflow.add_list_of_constraints([
      pc.sum([f[p,i] for p in G.predecessors(i)])
      == pc.sum([f[i,j] for j in G.successors(i)])
      for i in G.nodes() if i not in (s,t)])

  # Set source flow at s.
  maxflow.add_constraint(
    pc.sum([f[p,s] for p in G.predecessors(s)]) + F
    == pc.sum([f[s,j] for j in G.successors(s)]))

  # Set sink flow at t.
  maxflow.add_constraint(
    pc.sum([f[p,t] for p in G.predecessors(t)])
    == pc.sum([f[t,j] for j in G.successors(t)]) + F)

  # Enforce flow nonnegativity.
  maxflow.add_list_of_constraints([f[e] >= 0 for e in G.edges()])

  # Set the objective.
  maxflow.set_objective('max', F)

  # Solve the problem.
  maxflow.solve(solver='glpk')

.. _newversion:

An equivalent and faster way to define this problem is to use the class
:func:`~picos.flow_Constraint`:

.. plot::
  :context:
  :nofigs:
  :include-source:

  maxflow2=pc.Problem()

  # Add the flow variables.
  f={}
  for e in G.edges():
    f[e]=maxflow2.add_variable('f[{0}]'.format(e))

  # Add another variable for the total flow.
  F=maxflow2.add_variable('F')

  # Enforce all flow constraints at once.
  maxflow2.add_constraint(pc.flow_Constraint(
    G, f, source=16, sink=10, capacity='capacity', flow_value=F, graphName='G'))

  # Set the objective.
  maxflow2.set_objective('max', F)

  # Solve the problem.
  maxflow2.solve(solver='glpk')

Let us now draw the maximum flow computed with the second approach:

.. plot::
  :context:
  :include-source:

  # Close the old figure and open a new one.
  new_figure()

  # Determine which edges carry flow.
  flow_edges=[e for e in G.edges() if f[e].value > 1e-4]

  # Draw the nodes and the edges that don't carry flow.
  nx.draw_networkx(G, pos, edge_color='lightgrey', node_color=node_colors,
    edgelist=[e for e in G.edges
      if e not in flow_edges and (e[1], e[0]) not in flow_edges])

  # Draw the edges that carry flow.
  nx.draw_networkx_edges(G, pos, edgelist=flow_edges)

  # Show flow values and capacities on these edges.
  labels={e: '{0}/{1}'.format(f[e], c[e]) for e in flow_edges}
  nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)

  # Show the maximum flow value.
  fig.suptitle("Maximum flow value: {}".format(F), fontsize=16, y=0.95)

  # Show the figure.
  pylab.show()

The graph shows the source in blue, the sink in green, and the value of the flow
together with the capacity on each edge that carries flow.

Min-cut (LP)
------------

Given a directed graph :math:`G(V,E)`, with a capacity :math:`c(e)` on each edge
:math:`e \in E`, a source node :math:`s` and a sink node :math:`t`, the
**min-cut** problem is to find a partition of the nodes in two sets
:math:`(S,T)`, such that :math:`s\in S`, :math:`t \in T`, and the total capacity
of the cut,
:math:`\mathrm{capacity}(S,T)=\sum_{(i,j)\in E \cap S \times T} c((i,j)),` is
minimized.

It can be seen that binary solutions :math:`d\in\{0,1\}^E,\ p\in\{0,1\}^V`
of the following linear program yield a minimum cut:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\substack{d \in \mathbb{R}^E\\
                             p \in \mathbb{R}^V}}
                {\mbox{minimize}}
                      & \sum_{e \in E} c(e) d(e)\\
   &\mbox{subject to} & \forall (i,j) \in E,\ d((i,j)) \geq p(i)-p(j)\\
   &                  & p(s) = 1\\
   &                  & p(t) = 0\\
   &                  & \forall n \in V,\ p(n) \geq 0\\
   &                  & \forall e \in E,\ d(e) \geq 0
   \end{eqnarray*}
   \end{center}

Remarkably, this LP is the dual of the max-flow LP, and the max-flow-min-cut
theorem (also known as Ford-Fulkerson theorem :ref:`[1] <graph_refs>`) states
that the capacity of the minimum cut is equal to the value of the maximum flow.
This means that the above LP always has an optimal solution in which :math:`d`
is binary. In fact, the matrix defining this LP is *totally unimodular*, from
which we know that every extreme point of the polyhedron defining the feasible
region is integral, and hence the simplex algorithm will return a minimum cut.

We solve the min-cut problem below, again for ``s=16`` and ``t=10``:

.. plot::
  :context:
  :nofigs:
  :include-source:

  mincut=pc.Problem()

  # Add cut indicator variables.
  d={}
  for e in G.edges():
    d[e]=mincut.add_variable('d[{0}]'.format(e))

  # Add variables for the potentials.
  p=mincut.add_variable('p', N)

  # State the potential inequalities.
  mincut.add_list_of_constraints([d[i,j] >= p[i]-p[j] for (i,j) in G.edges()])

  # Set the source potential to one.
  mincut.add_constraint(p[s] == 1)

  # Set the sink potential to zero.
  mincut.add_constraint(p[t] == 0)

  # Enforce nonnegativity.
  mincut.add_constraint(p >= 0)
  mincut.add_list_of_constraints([d[e] >= 0 for e in G.edges()])

  # Set the objective.
  mincut.set_objective('min', pc.sum([cc[e]*d[e] for e in G.edges()]))

  mincut.solve(solver='glpk')

  # Determine the cut edges and node sets.
  # Rounding is done because solvers might return near-optimal solutions due to
  # numerical precision issues.
  cut=[e for e in G.edges() if abs(d[e].value-1) < 1e-6]
  S  =[n for n in G.nodes() if abs(p[n].value-1) < 1e-6]
  T  =[n for n in G.nodes() if abs(p[n].value  ) < 1e-6]

Let us now draw the minimum cut:

.. plot::
  :context:
  :include-source:

  # Close the old figure and open a new one.
  new_figure()

  # Draw the nodes and the edges that are not in the cut.
  nx.draw_networkx(G, pos, node_color=node_colors,
    edgelist=[e for e in G.edges() if e not in cut and (e[1], e[0]) not in cut])

  # Draw edges that are in the cut.
  nx.draw_networkx_edges(G, pos, edgelist=cut, edge_color='r')

  # Show capacities for cut edges.
  labels={e: '{}'.format(c[e]) for e in cut}
  nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_color='r')

  # Show the minimum cut value and the partition.
  fig.suptitle("Minimum cut value: {}\nS: {}, T: {}".format(
    mincut.obj_value(), S, T), fontsize=16, y=0.97)

  # Show the figure.
  pylab.show()

Note that the minimum-cut can also be derived from the dual variables of the
max-flow LP:

.. plot::
  :context:
  :nofigs:
  :include-source:

  # capacited flow constraint
  capaflow = maxflow.get_constraint((0,))
  dualcut = [
    e for i, e in enumerate(G.edges()) if abs(capaflow[i].dual - 1) < 1e-6]

  # flow conservation constraint
  consflow = maxflow.get_constraint((1,))

  Sdual = [s] + [
    n for i, n in enumerate([n for n in G.nodes() if n not in (s,t)])
    if abs(consflow[i].dual - 1) < 1e-6]

  Tdual = [t] + [
    n for i, n in enumerate([n for n in G.nodes() if n not in (s,t)])
    if abs(consflow[i].dual) < 1e-6]

Let's see how this dual-derived cut looks like:

.. plot::
  :context:
  :include-source:

  # Close the old figure and open a new one.
  new_figure()

  # Draw the nodes and the edges that are not in the dual cut.
  nx.draw_networkx(G, pos, node_color=node_colors, edgelist=[
    e for e in G.edges() if e not in dualcut and (e[1], e[0]) not in dualcut])

  # Draw edges that are in the dual cut.
  nx.draw_networkx_edges(G, pos, edgelist=dualcut, edge_color='b')

  # Show capacities for dual cut edges.
  labels={e: '{}'.format(c[e]) for e in dualcut}
  nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_color='b')

  # Show the dual cut value and the partition.
  fig.suptitle("Minimum cut value: {}\nS: {}, T: {}".format(
    sum(cc[e] for e in dualcut), Sdual, Tdual), fontsize=16, y=0.97)

  # Show the figure.
  pylab.show()

The graph shows the source in blue, the sink in green, and the edges defining
the cut in red, with their capacities.

Multicut (MIP)
--------------

Multicut is a generalization of the min-cut problem, in which several pairs of
nodes must be disconnected. The goal is to find a cut of minimal capacity, such
that for all pairs :math:`(s,t) \in\mathcal{P}=\{(s_1,t_1),\ldots,(s_k,t_k))\}`,
there is no path from :math:`s` to :math:`t` in the graph obtained by removing
the cut edges.

We can obtain a MIP formulation of the multicut problem via a small modification
of the min-cut LP. The idea is to introduce a different potential for every
node that is the source of a pair in :math:`\mathcal{P}`, that is

.. math::

   \forall s \in \mathcal{S}=\{s\in V: \exists t \in V\ (s,t)\in\mathcal{P}\},
   p_s \in \mathbb{R}^V,

and to constrain the cut indicator variables to be binary.

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{\substack{y \in \{0,1\}^E\\
                        \forall s \in \mathcal{S},\ p_s \in \mathbb{R}^V}}
                {\mbox{minimize}}
                      & \sum_{e \in E} c(e) y(e)\\
   &\mbox{subject to} & \forall (i,j),s \in E\times\mathcal{S},\ y((i,j)) \geq p_s(i)-p_s(j)\\
   &                  & \forall s \in \mathcal{S},\ p_s(s) = 1\\
   &                  & \forall (s,t) \in \mathcal{P},\ p_s(t) = 0\\
   &                  & \forall (s,n) \in \mathcal{S} \times V,\ p_s(n) \geq 0
   \end{eqnarray*}
   \end{center}

Unlike the min-cut problem, the LP obtained by relaxing the integer constraint
:math:`y \in \{0,1\}^E` is not guaranteed to have an integral solution
(see e.g. :ref:`[2] <graph_refs>`).

We solve the multicut problem below, for the terminal pairs
:math:`\mathcal{P}=\{(0,12),(1,5),(1,19),(2,11),(3,4),(3,9),(3,18),(6,15),(10,14)\}`.

.. plot::
  :context:
  :nofigs:
  :include-source:

  multicut=pc.Problem()

  # Define the pairs to be separated.
  pairs=[(0,12),(1,5),(1,19),(2,11),(3,4),(3,9),(3,18),(6,15),(10,14)]

  # Extract the sources and sinks.
  sources=set([p[0] for p in pairs])
  sinks=set([p[1] for p in pairs])

  # Define the cut indicator variables.
  y={}
  for e in G.edges():
    y[e]=multicut.add_variable('y[{0}]'.format(e), vtype='binary')

  # Define one potential for each source.
  p={}
  for s in sources:
    p[s]=multicut.add_variable('p[{0}]'.format(s), N)

  # State the potential inequalities.
  multicut.add_list_of_constraints(
    [y[i,j] >= p[s][i]-p[s][j] for s in sources for (i,j) in G.edges()])

  # Set the source potentials to one.
  multicut.add_list_of_constraints([p[s][s] == 1 for s in sources])

  # Set the sink potentials to zero.
  multicut.add_list_of_constraints([p[s][t] == 0 for (s,t) in pairs])

  # Enforce nonnegativity.
  multicut.add_list_of_constraints([p[s] >= 0 for s in sources])

  # Set the objective.
  multicut.set_objective('min', pc.sum([cc[e]*y[e] for e in G.edges()]))

  # Solve the problem.
  multicut.solve(solver='glpk')

  # Extract the cut.
  cut=[e for e in G.edges() if round(y[e]) == 1]

Let us now draw the multicut:

.. plot::
  :context:
  :include-source:

  # Close the old figure and open a new one.
  new_figure()

  # Define matching colors for the pairs.
  colors=[
    ('#4CF3CE','#0FDDAF'), # turquoise
    ('#FF4D4D','#FF0000'), # red
    ('#FFA64D','#FF8000'), # orange
    ('#3ABEFE','#0198E1'), # topaz
    ('#FFDB58','#FFCC11'), # mustard
    ('#BCBC8F','#9F9F5F')  # khaki
  ]

  # Assign the colors.
  node_colors=['lightgrey']*N
  for i,s in enumerate(sources):
    node_colors[s]=colors[i][0]
    for t in [t for (s0,t) in pairs if s0==s]:
      node_colors[t]=colors[i][1]

  # Draw the nodes and the edges that are not in the cut.
  nx.draw_networkx(G, pos, node_color=node_colors,
    edgelist=[e for e in G.edges() if e not in cut and (e[1], e[0]) not in cut])

  # Draw the edges that are in the cut.
  nx.draw_networkx_edges(G, pos, edgelist=cut, edge_color='r')

  # Show capacities for cut edges.
  labels={e: '{}'.format(c[e]) for e in cut}
  nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_color='r')

  # Show the cut capacity.
  fig.suptitle("Multicut value: {}"
    .format(multicut.obj_value()), fontsize=16, y=0.95)

  # Show the figure.
  pylab.show()

The graph shows terminal nodes with matching hue. Sources are a tad lighter than
sinks to make them distinguishable. The edges defining the cut are drawn in red
and show their capacities. The colors for the source nodes are, in order:
Turquoise, red, orange, topaz, mustard and khaki.

Maxcut relaxation (SDP)
-----------------------

The goal of the **maxcut** problem is to find a partition (S,T) of the nodes of
an *undirected* graph :math:`G(V,E)`, such that the capacity of the cut,
:math:`\mathrm{capacity}(S,T)=\sum_{\{i,j\} \in E \cap (S \Delta T)} c((i,j))`,
is maximized.

Goemans and Williamson have designed a famous 0.878-approximation algorithm
:ref:`[3] <graph_refs>` for this NP-hard problem based on semidefinite
programming. The idea is to introduce a variable :math:`x \in \{-1,1\}^V` where
:math:`x(n)` takes the value :math:`+1` or :math:`-1` depending on whether
:math:`n \in S` or :math:`n \in T`. Then, it can be seen that the value of the
cut is equal to :math:`\frac{1}{4} x^T L x`, where :math:`L` is the Laplacian
of the graph. If we define the matrix :math:`X=xx^T`, which is positive
semidefinite and of rank 1, we obtain an SDP by relaxing the rank-one constraint
on :math:`X`:

.. math::
   :nowrap:

   \begin{center}
   \begin{eqnarray*}
   &\underset{X \in \mathbb{S}_{|V|}}
                {\mbox{maximize}}
                      & \frac{1}{4} \langle L, X \rangle \\
   &\mbox{subject to} & \mbox{diag}(X) = \mathbf{1}\\
   &                  & X \succeq 0
   \end{eqnarray*}
   \end{center}

Then, Goemans and Williamson have shown that if we project the solution
:math:`X` onto a random hyperplan, we obtain a cut whose expected capacity is at
least 0.878 times the optimum. We give a simple implementation of their
algorithm. First, let us define and solve the SDP relaxation:

.. plot::
  :context:
  :nofigs:
  :include-source:

  import cvxopt as cvx
  import cvxopt.lapack
  import numpy as np

  # Make G undirected.
  G=nx.Graph(G)

  # Allocate weights to the edges.
  for (i,j) in G.edges():
    G[i][j]['weight']=c[i,j]+c[j,i]

  maxcut = pc.Problem()

  # Add the symmetric matrix variable.
  X=maxcut.add_variable('X', (N,N), 'symmetric')

  # Retrieve the Laplacian of the graph.
  LL = 1/4.*nx.laplacian_matrix(G).todense()
  L=pc.new_param('L', LL)

  # Constrain X to have ones on the diagonal.
  maxcut.add_constraint(pc.diag_vect(X) == 1)

  # Constrain X to be positive semidefinite.
  maxcut.add_constraint(X >> 0)

  # Set the objective.
  maxcut.set_objective('max', L|X)

  #print(maxcut)

  # Solve the problem.
  maxcut.solve(solver='cvxopt')

  #print('bound from the SDP relaxation: {0}'.format(maxcut.obj_value()))

Then, we perform the random projection:

.. plot::
  :context:
  :nofigs:
  :include-source:

  # Use a fixed RNG seed so the result is reproducable.
  cvx.setseed(1)

  # Perform a Cholesky factorization.
  V=X.value
  cvxopt.lapack.potrf(V)
  for i in range(N):
    for j in range(i+1,N):
      V[i,j]=0

  # Do up to 100 projections. Stop if we are within a factor 0.878 of the SDP
  # optimal value.
  count=0
  obj_sdp=maxcut.obj_value()
  obj=0
  while (count < 100 or obj < 0.878*obj_sdp):
    r=cvx.normal(20,1)
    x=cvx.matrix(np.sign(V*r))
    o=(x.T*L*x).value
    if o > obj:
      x_cut=x
      obj=o
    count+=1
  x=x_cut

  # Extract the cut and the seperated node sets.
  S1=[n for n in range(N) if x[n]<0]
  S2=[n for n in range(N) if x[n]>0]
  cut = [(i,j) for (i,j) in G.edges() if x[i]*x[j]<0]
  leave = [e for e in G.edges if e not in cut]

Let us now draw this cut:

.. plot::
  :context:
  :include-source:

  # Close the old figure and open a new one.
  new_figure()

  # Assign colors based on set membership.
  node_colors=[('lightgreen' if n in S1 else 'lightblue') for n in range(N)]

  # Draw the nodes and the edges that are not in the cut.
  nx.draw_networkx(G, pos, node_color=node_colors, edgelist=leave)
  labels={e: '{}'.format(G[e[0]][e[1]]['weight']) for e in leave}
  nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)

  # Draw the edges that are in the cut.
  nx.draw_networkx_edges(G, pos, edgelist=cut, edge_color='r')
  labels={e: '{}'.format(G[e[0]][e[1]]['weight']) for e in cut}
  nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_color='r')

  # Show the relaxation optimum value and the cut capacity.
  rval = maxcut.obj_value()
  sval = sum(G[e[0]][e[1]]['weight'] for e in cut)
  fig.suptitle(
    'SDP relaxation value: {0:.1f}\nCut value: {1:.1f} = {2:.3f}×{0:.1f}'
    .format(rval, sval, sval/rval), fontsize=16, y=0.97)

  # Show the figure.
  pylab.show()

The graph shows the edges defining the cut in red. The nodes are colored blue or
green depending on the partition that they belong to.

.. _graph_refs:

References
----------

  1. "Maximal Flow through a Network", LR Ford Jr and DR Fulkerson,
     *Canadian journal of mathematics*, 1956.

  2. "Analysis of LP relaxations for multiway and multicut problems",
     D.Bertsimas, C.P. Teo and R. Vohra,
     *Networks*, 34(2), p. *102-114*, 1999.

  3. "Improved approximation algorithms for maximum cut and satisfiability
     problems using semidefinite programming",
     M.X. Goemans and D.P. Williamson,
     *Journal of the ACM*, 42(6), p. *1115-1145*, 1995.
.. warning::

    This part of the documentation has not been touched for a while. It might
    be incomplete, reference deprecated functions or make a claim that does not
    apply to the latest version of PICOS any more. On the bright side, code
    listings are validated and still work. Watch your step!


.. _constraints:

Additional constraints
======================

This section introduces additional expression and constraint types that didn't
fit into the tutorial. Again, let us import PICOS:

>>> import picos as pc

We replicate some of the variables and data used in the tutorial:

>>> from picos import Constant, RealVariable
>>> t = RealVariable("t")
>>> x = RealVariable("x", 4)
>>> Y = RealVariable("Y", (2, 4))
>>> alpha = Constant("α", 23)
>>> A = [Constant("A[{}]".format(i), range(i, i + 8), (2, 4)) for i in range(5)]

This time, :math:`Z` and :math:`b` are lists:

>>> Z = [RealVariable("Z[{0}]".format(i), (4,2)) for i in range(5)]
>>> b = ([0, 2, 0, 3], [1, 1, 0, 5], [-1, 0, 2, 4], [0, 0, -2, -1],
...      [1, 1, 0, 0])
>>> b = [Constant("b[{}]".format(i), b[i]) for i in range(len(b))]


.. _flowcons:

Graph flow constraints
----------------------

Flow constraints in graphs are entered using a Networkx_ graph. The following
example finds a (trivial) maximal flow from ``'S'`` to ``'T'`` in ``G``.

.. _Networkx: https://networkx.github.io/

>>> import networkx as nx
>>> G = nx.DiGraph()
>>> G.add_edge('S','A', capacity=1)
>>> G.add_edge('A','B', capacity=1)
>>> G.add_edge('B','T', capacity=1)
>>> pb = pc.Problem()
>>> # Adding the flow variables
>>> f={}
>>> for e in G.edges():
...     f[e]=pb.add_variable('f[{},{}]'.format(e[0], e[1]), 1)
>>> # A variable for the value of the flow
>>> F = pb.add_variable('F',1)
>>> # Creating the flow constraint
>>> flowCons = pb.add_constraint(pc.flow_Constraint(
...     G, f, 'S', 'T', F, capacity='capacity', graphName='G'))
>>> pb.set_objective('max',F)
>>> sol = pb.solve(solver='cvxopt')
>>> flow = {key: var.value for key, var in f.items()}

Picos allows you to define single source - multiple sinks problems.
You can use the same syntax as for a single source - single sink problems.
Just add a list of sinks and a list of flows instead.

.. TODO: Get rid of these testcode/testoutput so we can test via unittest.

.. testcode::

    import picos as pc
    import networkx as nx

    G=nx.DiGraph()
    G.add_edge('S','A', capacity=2); G.add_edge('S','B', capacity=2)
    G.add_edge('A','T1', capacity=2); G.add_edge('B','T2', capacity=2)

    pbMultipleSinks=pc.Problem()
    # Flow variable
    f={}
    for e in G.edges():
        f[e]=pbMultipleSinks.add_variable('f[{},{}]'.format(e[0], e[1]), 1)

    # Flow value
    F1=pbMultipleSinks.add_variable('F1',1)
    F2=pbMultipleSinks.add_variable('F2',1)

    flowCons = pc.flow_Constraint(
        G, f, source='S', sink=['T1','T2'], capacity='capacity',
        flow_value=[F1, F2], graphName='G')

    pbMultipleSinks.add_constraint(flowCons)
    pbMultipleSinks.set_objective('max',F1+F2)

    # Solve the problem
    pbMultipleSinks.solve(solver='cvxopt')

    print(pbMultipleSinks)
    print()
    print('The optimal flow F1 has value {:.1f}'.format(F1.value))
    print('The optimal flow F2 has value {:.1f}'.format(F2.value))

.. testoutput::
    :options: +NORMALIZE_WHITESPACE

    Linear Program
      maximize F1 + F2
      over
        1×1 real variables F1, F2, f[A,T1], f[B,T2], f[S,A], f[S,B]
      subject to
        Feasible S-(T1,T2)-flow in G has values F1, F2.

    The optimal flow F1 has value 2.0
    The optimal flow F2 has value 2.0

A similar syntax can be used for multiple sources-single sink flows.


Second Order Cone constraints
-----------------------------

.. TODO: Remove the warning and document the difference between a conic
..       quadratic constraint and an explicit SOC/RSOC constraint. Needs #157.

.. warning::

    This section in particular is outdated: The only direct way to create
    rotated second order cone constraints is now via the :func:`~picos.rsoc`
    set-generating function. If you input such a constraint as below, then you
    will receive either a convex or a conic quadratic constraint. The former is
    handled depending on the solver used. The latter will be transformed either
    to a rotated conic constraint (which implicitly adds additional constraints
    that are part of the rotated second order cone definition) or it will remain
    nonconvex quadratic, depending on an option.

There are two types of second order cone constraints supported in PICOS.

    * The constraints of the type :math:`\Vert x \Vert \leq t`, where :math:`t`
      is a scalar affine expression and :math:`x` is
      a multidimensional affine expression (possibly a matrix, in which case the
      norm is Frobenius). This inequality forces
      the vector :math:`[t; x]` to belong to a Lorrentz-Cone (also called
      *ice-cream cone*).
    * The constraints of the type :math:`\Vert x \Vert^2 \leq t u,\ t \geq 0`,
      where :math:`t` and :math:`u` are scalar affine expressions and
      :math:`x` is a multidimensional affine expression, which constrain
      the vector :math:`[t; u; x]` inside a rotated version of the Lorretz cone.

A few examples:

>>> # A simple ice-cream cone constraint
>>> abs(x) < (2|x-1)
<5×1 SOC Constraint: ‖x‖ ≤ ⟨[2], x - [1]⟩>
>>> # SOC constraint with Frobenius norm
>>> abs(Y+Z[0].T) < t+alpha
<9×1 SOC Constraint: ‖Y + Z[0]ᵀ‖ ≤ t + α>
>>> # Rotated SOC constraint
>>> abs(Z[1][:,0])**2 < (2*t-alpha)*(x[2]-x[-1])
<Conic Quadratic Constraint: ‖Z[1][:,0]‖² ≤ (2·t - α)·(x[2] - x[-1])>
>>> # t**2 is internally represented as the squared norm of [t]
>>> t**2 < alpha + t
<Squared Scalar Constraint: t² ≤ α + t>
>>> # 1 is understood as the squared norm of [1]
>>> 1 < (t-1)*(x[2]+x[3])
<Conic Quadratic Constraint: (t - 1)·(x[2] + x[3]) ≥ 1>


Semidefinite constraints
------------------------

Linear matrix inequalities (LMI) can be entered thanks to an overload of the
operators ``<<`` and ``>>``. For example, the LMI

.. math::
    :nowrap:

    \begin{equation*}
        \sum_{i=0}^3 x_i b_i b_i^T \succeq b_4 b_4^T,
    \end{equation*}

where :math:`\succeq` is used to denote the Löwner ordering, is passed to PICOS
by writing:

>>> pc.sum([x[i]*b[i]*b[i].T for i in range(4)]) >> b[4]*b[4].T
<4×4 LMI Constraint: ∑(x[i]·b[i]·b[i]ᵀ : i ∈ [0…3]) ≽ b[4]·b[4]ᵀ>

Note the difference with

>>> pc.sum([x[i]*b[i]*b[i].T for i in range(4)]) > b[4]*b[4].T
<4×4 Affine Constraint: ∑(x[i]·b[i]·b[i]ᵀ : i ∈ [0…3]) ≥ b[4]·b[4]ᵀ>

which yields an elementwise inequality.

For convenience, it is possible to add a symmetric matrix variable ``X``,
by specifying the option ``vtype=symmetric``. This has the effect to
store all the affine expressions which depend on ``X`` as a function
of its lower triangular elements only.

>>> sdp = pc.Problem()
>>> X = sdp.add_variable('X',(4,4),vtype='symmetric')
>>> C = sdp.add_constraint(X >> 0)
>>> print(sdp)
Feasibility Problem
  find an assignment
  for
    4×4 symmetric variable X
  subject to
    X ≽ 0

In this example, you see indeed that the problem has 10=(4*5)/2 variables,
which correspond to the lower triangular elements of ``X``.

.. warning::

     When a constraint of the form ``A >> B`` is passed to PICOS, it is not
     enforced that :math:`A - B` is symmetric. How the constraint is passed then
     depends on the solver, for instance it could be that the lower or upper
     triangular part is ignored. You can add a constraint of the form
     ``A - B == (A - B).T`` to enforce symmetry.


Inequalities involving geometric means
--------------------------------------

It is possible to enter an inequality of the form

.. math::
    t \leq \prod_{i=1}^n x_i^{1/n}

in PICOS, where :math:`t` is a scalar affine expression and :math:`x` is an
affine expression of dimension :math:`n` (possibly a matrix, in which case
:math:`x_i` is counted in column major order). This inequality is internally
converted to an equivalent set of second order cone inequalities, by using
standard techniques (cf. e.g. :ref:`[1] <tuto_refs>`).

Many convex constraints can be formulated using inequalities that involve
a geometric mean. For example, :math:`t \leq x_1^{2/3}` is equivalent
to :math:`t \leq t^{1/4} x_1^{1/4} x_1^{1/4}`, which can be entered in PICOS
thanks to the function :func:`~picos.geomean`:

  >>> t < pc.geomean(t //x[1] //x[1] //1)
  <Geometric Mean Constraint: geomean([t; x[1]; x[1]; 1]) ≥ t>

Note that the latter example can also be passed to PICOS in a more simple way,
thanks to an overloading of the ``**`` exponentiation operator:

  >>> t < x[1]**(2./3)
  <Power Constraint: x[1]^(2/3) ≥ t>

Such a power constraint will be reformulated as a geometric mean inequality when
the problem is solved, which in turn will be translated to conic inequalities.


Inequalities involving real powers or trace of matrix powers
------------------------------------------------------------

As mentionned above, the ``**`` exponentiation operator has been overloaded
to support real exponents. A rational approximation of the exponent is used,
and the inequality are internally reformulated as a set of equivalent SOC
inequalities. Note that only inequalities defining a convex regions can be
passed:

>>> t**0.6666 > x[0]
<Power Constraint: t^(2/3) ≥ x[0]>
>>> t**-0.5 < x[0]
<Power Constraint: t^(-1/2) ≤ x[0]>
>>> t**-0.5 > x[0]
Traceback (most recent call last):
  ...
TypeError: Cannot lower-bound a nonconcave (trace of) power.

More generally, inequalities involving trace of matrix powers can be passed to
PICOS, by using the :func:`~picos.tracepow` function. The following example
creates the constraint

.. math::

    \operatorname{trace}\ \big(x_0 A_0 A_0^T + x_2 A_2 A_2^T\big)^{2.5} \leq 3.

>>> pc.tracepow(x[0] * A[0]*A[0].T + x[2] * A[2]*A[2].T, 2.5) <= 3
<Trace of Power Constraint: tr((x[0]·A[0]·A[0]ᵀ + x[2]·A[2]·A[2]ᵀ)^(5/2)) ≤ 3>

.. Warning::

    when a power expression :math:`x^p` (resp. the trace of matrix power
    :math:`\operatorname{trace}\ X^p` ) is used, the base :math:`x` is forced
    to be nonnegative (resp. the base :math:`X` is forced to be positive
    semidefinite) by picos.

When the exponent is :math:`0<p<1`,
it is also possible to represent constraints of the form
:math:`\operatorname{trace}(M X^p) \geq t`
with SDPs, where :math:`M\succeq 0`, see :ref:`[2] <tuto_refs>`.

>>> pc.tracepow(X, 0.6666, coef = A[0].T*A[0]+"I") >= t
<Trace of Scaled Power Constraint: tr((A[0]ᵀ·A[0] + I)·X^(2/3)) ≥ t>

As for geometric means, inequalities involving real powers yield their internal
representation via the ``constraints`` and ``variables`` attributes.


.. _pnorms:

Inequalities involving generalized p-norm
-----------------------------------------

Inequalities of the form :math:`\Vert x \Vert_p \leq t` can be entered by using the
function :func:`~picos.norm`. This function is also defined for :math:`p < 1`
by the usual formula :math:`\Vert x \Vert_p :=  \Big(\sum_i |x_i|^p \Big)^{1/p}`.
The norm function is convex over :math:`\mathbb{R}^n` for all :math:`p\geq 1`, and
concave over the set of vectors with nonnegative coordinates for :math:`p \leq 1`.

>>> pc.norm(x,3) < t
<Vector p-Norm Constraint: ‖x‖_3 ≤ t>
>>> pc.norm(x,'inf') < 2
<Maximum Norm Constraint: ‖x‖_max ≤ 2>
>>> pc.norm(x,0.5) > x[0]-x[1]
<Generalized p-Norm Constraint: ‖x‖_(1/2) ≥ x[0] - x[1] ∧ x ≥ 0>

.. Warning::

    Note that when a constraint of the form ``norm(x,p) >= t`` is entered (with
    :math:`p \leq 1` ), PICOS forces the vector ``x`` to be nonnegative
    (componentwise).

Inequalities involving the generalized :math:`L_{p,q}` norm of
a matrix can also be handled with picos, cf. the documentation of
:func:`~picos.norm` .

As for geometric means, inequalities involving p-norms yield their internal
representation via the ``constraints`` and ``variables`` attributes.


Inequalities involving the nth root of a determinant
----------------------------------------------------

The function :func:`~picos.detrootn`
can be used to enter the :math:`n`-th root of the determinant of a
:math:`(n \times n)`-symmetric positive semidefinite matrix:

>>> M = sdp.add_variable('M',(5,5),'symmetric')
>>> t < pc.detrootn(M)
<n-th Root of a Determinant Constraint: det(M)^(1/5) ≥ t>

.. warning::

    Note that when a constraint of the form ``t < pc.detrootn(M)`` is entered
    (with :math:`p \leq 1`), PICOS forces the matrix ``M`` to be positive
    semidefinite.

As for geometric means, inequalities involving the nth root of a determinant
yield their internal representation via the ``constraints`` and ``variables``
attributes.


Set membership
--------------

Since Picos 1.0.2, there is a :class:`Set <picos.expressions.Set>` class that
can be used to pass constraints as membership of an affine expression to a set.

Following sets are currently supported:

    * :math:`L_p-` balls representing the set
      :math:`\{x: \Vert x \Vert_p \leq r\}` can be constructed with the function
      :func:`~picos.ball`
    * The standard simplex (scaled by a factor :math:`\gamma`)
      :math:`\{x \geq 0: \sum_i x_i \leq r \}` can be constructed with the
      function :func:`~picos.simplex`
    * Truncated simplexes :math:`\{0 \leq x \leq 1: \sum_i x_i \leq r \}`
      and symmetrized Truncated simplexes
      :math:`\{x: \Vert x \Vert_\infty \leq 1, \Vert x \Vert_1\leq r \}`
      can be constructed with the function :func:`~picos.truncated_simplex`

Membership of an affine expression to a set can be expressed with the overloaded
operator ``<<``. This returns a temporary object that can be passed to a picos
problem with the function :meth:`~.problem.Problem.add_constraint`.

>>> x << pc.simplex(1)
<Unit Simplex Constraint: x ∈ {x ≥ 0 : ∑(x) ≤ 1}>
>>> x << pc.truncated_simplex(2)
<Box-Truncated Simplex Constraint: x ∈ {0 ≤ x ≤ 1 : ∑(x) ≤ 2}>
>>> x << pc.truncated_simplex(2,sym=True)
<Box-Truncated 1-norm Ball Constraint: x ∈ {-1 ≤ x ≤ 1 : ∑(|x|) ≤ 2}>
>>> x << pc.ball(3)
<5×1 SOC Constraint: ‖x‖ ≤ 3>
>>> pc.ball(2,'inf') >> x
<Maximum Norm Constraint: ‖x‖_max ≤ 2>
>>> x << pc.ball(4,1.5)
<Vector p-Norm Constraint: ‖x‖_(3/2) ≤ 4>


.. _tuto_refs:

References
----------

1. "`Applications of second-order cone programming`",
   M.S. Lobo, L. Vandenberghe, S. Boyd and H. Lebret,
   *Linear Algebra and its Applications*,
   284, p. *193-228*, 1998.

2. "`On the semidefinite representations of real functions applied to symmetric
   matrices <http://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/1751>`_"
   , G. Sagnol,
   *Linear Algebra and its Applications*,
   439(10), p. *2829-2843*, 2013.
.. |docs| image:: https://img.shields.io/website.svg?down_color=red&down_message=offline&label=docs&logo=read-the-docs&logoColor=white&up_color=brightgreen&up_message=online&url=https%3A%2F%2Fpicos-api.gitlab.io%2Fpicos
  :alt: Documentation
  :target: https://picos-api.gitlab.io/picos/

.. |gitlab| image:: https://img.shields.io/badge/GitLab-source-blue?logo=gitlab&logoColor=white
  :alt: GitLab
  :target: https://gitlab.com/picos-api/picos

.. |pypi| image:: https://img.shields.io/pypi/v/picos?label=PyPI&logo=pypi&logoColor=white
  :alt: PyPI
  :target: https://pypi.org/project/PICOS/

.. |anaconda| image:: https://img.shields.io/conda/vn/picos/picos?label=Anaconda&logo=anaconda&logoColor=white
  :alt: Anaconda
  :target: https://anaconda.org/picos/picos

.. |aur| image:: https://img.shields.io/aur/version/python-picos?label=AUR&logo=arch-linux&logoColor=white
  :alt: AUR
  :target: https://aur.archlinux.org/packages/?SeB=n&K=python-picos

.. |freenode| image:: https://img.shields.io/badge/freenode-%23picos-green.svg
  :alt: Freenode
  :target: https://webchat.freenode.net/?nick=picosuser&channels=%23picos

.. |cov| image:: https://img.shields.io/gitlab/coverage/picos-api/picos/master
  :alt: Coverage
  :target: coverage/

.. |license| image:: https://img.shields.io/badge/copying-GPL--3.0-brightgreen
  :alt: License
  :target: https://choosealicense.com/licenses/gpl-3.0/.. _slicing:

Matrix Slicing
==============

Affine matrix expressions form the core of PICOS' modeling toolbox: All
:class:`constant <picos.Constant>` and :class:`variable
<picos.expressions.variables>` expressions that you enter, including integral
variables, and any linear combination of these objects, are stored as instances
of the multidimensional :class:`~picos.expressions.ComplexAffineExpression` or
its real subclass :class:`~picos.expressions.AffineExpression`. Their common
base class :class:`~picos.expressions.BiaffineExpression` implements plenty of
algebraic operations to combine and modify your initial expressions to yield the
desired statements. One of these operations is slicing, denoted by ``A[·]`` for
an affine expression ``A``.

.. rubric:: Preliminaries

Unlike in NumPy, all multidimensional expressions are strictly matrices. In
particular, there are no flat arrays but only row and column vectors, and any
scalar expression is also a :math:`1 \times 1` matrix. PICOS does not support
tensors or higher order expressions, but it does support the `Kronecker product
<https://en.wikipedia.org/wiki/Kronecker_product>`_ as well as :meth:`partial
trace <.exp_biaffine.BiaffineExpression.partial_trace>` and :meth:`partial
transposition <.exp_biaffine.BiaffineExpression.partial_transpose>` operations
to enable some of the optimization problems naturally defined on tensors. If you
enter data in the form of a flat array (e.g. a Python :class:`list` or a NumPy
:class:`~numpy:numpy.ndarray` with one axis), it will be read as a column
vector.

In PICOS, all indices start from zero.

.. rubric:: Slicing modes

PICOS has two modes for slicing: :ref:`arbitrary_access` and
:ref:`proper_slicing`.

Arbitrary Access lets you select individual elements from a vector or matrix
expression and put them in a column vector in the desired order.
:meth:`Transposition <.exp_biaffine.BiaffineExpression.T>`, :meth:`reshaping
<.exp_biaffine.BiaffineExpression.reshaped>` and :meth:`broadcasting
<.exp_biaffine.BiaffineExpression.broadcasted>` can then be used to put the
selection into the desired shape. Arbitrary Access has the form ``A[·]`` where
``·`` stands for an integer, a Python :class:`slice`, a flat collection of
integers such as a :class:`list`, or a dictionary storing sparse index pairs.

Proper Slicing refers to selecting certain rows and columns of a matrix, and
forming a new matrix where all elements that are not selected are removed.
It has the form ``A[·,·]`` where each ``·`` stands for an integer, a
:class:`slice`, or a flat collection of integers.

To demonstrate the different possibilities, we use a constant :math:`5 \times 5`
expression:

>>> from picos import Constant
>>> A = Constant("A", range(25), (5,5))
>>> A
<5×5 Real Constant: A>
>>> print(A)
[ 0.00e+00  5.00e+00  1.00e+01  1.50e+01  2.00e+01]
[ 1.00e+00  6.00e+00  1.10e+01  1.60e+01  2.10e+01]
[ 2.00e+00  7.00e+00  1.20e+01  1.70e+01  2.20e+01]
[ 3.00e+00  8.00e+00  1.30e+01  1.80e+01  2.30e+01]
[ 4.00e+00  9.00e+00  1.40e+01  1.90e+01  2.40e+01]

.. _arbitrary_access:

Arbitrary Access
----------------

.. rubric:: By integer

If a single integer or a single flat collection of integers is given, then these
indices refer to the column-major vectorization of the matrix, represented by
the order of the numbers in the demonstration matrix ``A``.

The most common case is selecting a single element via an integer index:

>>> A[0]  # Select the first element as a scalar expression.
<1×1 Real Constant: A[0]>
>>> print(A[0])  # Print its value.
0.0
>>> print(A[7])  # The value of the eighth element.
7.0
>>> # Negative indices are counted from the rear; -1 refers to the last element:
>>> print(A[-1])
24.0

.. rubric:: By slice

Python slices allow you to compactly specify a structured sequence of elements
to extract.
A Python slice has the form ``a:b`` or ``a:b:s`` with :math:`a` the inclusive
start index, :math:`b` the exclusive stop index and :math:`s` a step size.
Negative :math:`a` and :math:`b`, as in the integer index case, are counted from
the rear, while a negative step size reverses the order.
All of :math:`a`, :math:`b` and :math:`s` may be omitted. Then, the defaults are

.. math::

  s &= 1, \\
  a &= \begin{cases}
    0,~&\text{if}~s > 0, \\
    \text{len}(A) - 1,~&\text{if}~s < 0,
  \end{cases} \\
  b &= \begin{cases}
    \text{len}(A),~&\text{if}~s > 0, \\
    \textbf{None},~&\text{if}~s < 0.
  \end{cases}

Note the :obj:`None` in the statement above: When going backwards, this special
token is the only way to stop at the first element with index :math:`0` as
:math:`-1` refers to the last element. For example, the first two elements in
reverse order are selected via the slice ``1:None:-1`` or just ``1::-1``.

>>> A[:2]  # The first two elements as a column vector.
<2×1 Real Constant: A[:2]>
>>> print(A[:2])
[ 0.00e+00]
[ 1.00e+00]
>>> print(A[1::-1])  # The first two elements reversed (indices 1 and 0).
[ 1.00e+00]
[ 0.00e+00]
>>> print(A[-2:])  # The last two elements.
[ 2.30e+01]
[ 2.40e+01]
>>> print(A[2:7].T)  # The third to seventh element (transposed).
[ 2.00e+00  3.00e+00  4.00e+00  5.00e+00  6.00e+00]
>>> print(A[2:7:2].T)  # As before, but with a step size of 2.
[ 2.00e+00  4.00e+00  6.00e+00]

You could use this to vectorize :math:`A` in column-major order, but ``A.vec``
is both individually faster and has its result cached:

>>> A[:].equals(A.vec)
True
>>> A.vec is A.vec  # Cached.
True
>>> A[:] is A[:]  # Computed again as new expression.
False

.. rubric:: By integer sequence

By providing a :class:`list` or a similar vector of integers, you can select
arbitrary elements in any order, including duplicates:

>>> print(A[[0,1,0,1,-1]])
[ 0.00e+00]
[ 1.00e+00]
[ 0.00e+00]
[ 1.00e+00]
[ 2.40e+01]

Note that you cannot provide a :class:`tuple` instead of a list, as ``A[(·,·)]``
is understood by Python as ``A[·,·]`` (see :ref:`proper_slicing`).
Any other object that the function :func:`~picos.expressions.data.load_data`
with ``typecode="i"`` loads as an integer row or column vector works, including
integral NumPy arrays.

.. _sparse_index_dict:
.. rubric:: By sparse index pair dictionary

If you provide a dictionary with exactly two keys that can be compared via
``<`` and whose values are integer sequences of same length (anything recognized
by :func:`~picos.expressions.data.load_data` as an integer vector), PICOS
interprets the sequence corresponding to the smaller key as row indices and the
sequence corresponding to the greater key as the corresponding column indices:

>>> print(A[{"x": range(3), "y": [1]*3}])  # Select (0,1), (1,1) and (2,1).
[ 5.00e+00]
[ 6.00e+00]
[ 7.00e+00]
>>> print(A[{"y": range(3), "x": [1]*3}])  # Transposed selection, as "x" < "y".
[ 1.00e+00]
[ 6.00e+00]
[ 1.10e+01]

You could use this to extract the main diagonal of :math:`A`, but ``A.maindiag``
is both individually faster and has its result cached:

>>> indices = dict(enumerate([range(min(A.shape))]*2))
>>> indices
{0: range(0, 5), 1: range(0, 5)}
>>> A[indices].equals(A.maindiag)
True
>>> A.maindiag is A.maindiag  # Cached.
True
>>> A[indices] is A[indices]  # Computed again as new expression.
False

.. _proper_slicing:

Proper Slicing
--------------

If you provide not one but two integers, slices, or integer sequences separated
by a comma or given as a :obj:`tuple`, then they are understood as row and
column indices, respectively.
Unlike when providing a sparse index pair by dictionary, these indices select
*entire* rows and columns and PICOS returns the matrix of all elements that are
selected twice (both by row and by column):

>>> print(A[1,2])  # The single element at (1,2) (second row, third column).
11.0
>>> print(A[0,:])  # The first row of the matrix.
[ 0.00e+00  5.00e+00  1.00e+01  1.50e+01  2.00e+01]
>>> print(A[range(3),-1])  # The first three elements of the last column.
[ 2.00e+01]
[ 2.10e+01]
[ 2.20e+01]
>>> print(A[[0,1],[0,1]])  # The first second-order principal submatrix.
[ 0.00e+00  5.00e+00]
[ 1.00e+00  6.00e+00]
>>> print(A[1:-1,1:-1])  # Cut away the outermost pixels of an image.
[ 6.00e+00  1.10e+01  1.60e+01]
[ 7.00e+00  1.20e+01  1.70e+01]
[ 8.00e+00  1.30e+01  1.80e+01]
>>> print(A[::2,::2])  # Sample every second element.
[ 0.00e+00  1.00e+01  2.00e+01]
[ 2.00e+00  1.20e+01  2.20e+01]
[ 4.00e+00  1.40e+01  2.40e+01]

You can even select the entire matrix to effectively create a copy of it, though
this is discouraged as expressions are supposed to be immutable so that reusing
an expression in multiple places is always safe.

>>> A[:,:].equals(A)
True
>>> A[:,:] is A
False

We refer to this as proper slicing because you cut out the rows that you want,
throwing away the rest, then cut the desired columns out from the remainder.
It's like cutting a square cake except that you can also duplicate the pieces!

.. note::
  In NumPy, ``A[[0,1],[0,1]]`` would create a flat array with the elements
  ``A[0,0]`` and ``A[1,1]`` while PICOS creates a submatrix from the first two
  rows and columns as in the example above. If you want to mirror NumPy's
  behavior in PICOS, see :ref:`sparse_index_dict`.
.. _api:

API Reference
=============

PICOS is organized in a number of submodules and subpackages, most of which you
do not need to access directly when solving optimization problems. It is usually
sufficient to ``import picos`` and use the functions and classes provided in the
:mod:`picos` namespace. A notable exception are the tools for handling uncertain
data that are found in the :mod:`picos.uncertain` namespace.

.. rubric:: Modules

.. toctree::
  :includehidden:
  :maxdepth: 2

  api/picos
{# THEMODULE #}
{{ node.name }}
{{ '=' * node.name|length }}

.. automodule:: {{ node.name }}
  :no-special-members:

{# SET CURRENT MODULE #}
.. currentmodule:: {{ node.name }}
{# EXCEPTIONS #}
{%- block exceptions -%}
{%- if node.exceptions %}

.. rubric:: Exceptions

{% for item, obj in node.exceptions.items() %}

{% if node.name == obj|module %}
.. autoexception:: {{ item }}
  :no-special-members:
{% else %}
.. exception:: {{ item }}

  See :py:exc:`{{ obj|module }}.{{ obj|qualname }}`.
{% endif %}
{% endfor -%}
{%- endif -%}
{%- endblock -%}
{# CLASSES #}
{%- block classes -%}
{%- if node.classes %}

.. rubric:: Classes

{% for item, obj in node.classes.items() %}

{% if node.name == obj|module %}
.. autoclass:: {{ item }}
  :members:
  :undoc-members:
{% else %}
.. class:: {{ item }}

  See :py:class:`{{ obj|module }}.{{ obj|qualname }}`.
{% endif %}
{% endfor -%}
{%- endif -%}
{%- endblock -%}
{# FUNCTIONS #}
{%- block functions -%}
{%- if node.functions %}

.. rubric:: Functions

{% for item, obj in node.functions.items() %}

{% if node.name == obj|module %}
.. autofunction:: {{ item }}
{% else %}
.. function:: {{ item }}

  See :py:func:`{{ obj|module }}.{{ obj|qualname }}`.
{% endif %}
{% endfor -%}
{%- endif -%}
{%- endblock -%}
{# OBJECTS #}
{%- block objects -%}
{%- if node.variables %}

.. rubric:: Objects

{% for item, obj in node.variables.items() %}

.. autodata:: {{ item }}
  :annotation:

  :Default value:

    .. code-block:: text

      {{ obj|pprint|truncate|indent(6) }}
{% endfor -%}
{%- endif -%}
{%- endblock -%}
{# TOCTREE #}
{%- block toctree -%}
{%- if subnodes %}

.. toctree::
  :hidden:

{% for subnode in subnodes %}
  {{ subnode.name }}
{% endfor -%}
{%- endif -%}
{%- endblock -%}

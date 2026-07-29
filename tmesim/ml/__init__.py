"""ML-based panel-gap imputation.

Additive package: nothing here is imported, or run, by any existing
reconstruction path unless a caller explicitly opts in (see
``panel_imputation.impute_gene_score`` and ``state_mode='ml_imputed'`` in
``reconstruct.real_pipeline``). No file under ``tmesim/core``, ``tmesim/physics``
or ``tmesim/sim`` is touched by this package.
"""

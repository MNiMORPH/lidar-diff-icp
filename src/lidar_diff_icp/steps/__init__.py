"""Pipeline steps: the producers the workflow graph declares, as package modules.

These were scripts under ``analysis/`` until 2026-09-06. Three of them --
``trace_ridgelines``, ``convexity_dod_landcover``, ``curvature_diffusion`` -- lived in
``analysis/ridgelines/``, which reads as a concluded investigation: of its 61 scripts only
those 3 were reachable from the pipeline, so the folder could not be retired without
losing half the step chain. A producer the graph declares is part of the product, and it
now lives with the product.

Each module is run as ``python -m lidar_diff_icp.steps.<name>``; ``workflow.py`` holds the
argument list. They execute on import by design, being scripts -- do not import them for
their functions.
"""

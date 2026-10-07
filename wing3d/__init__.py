"""wing3d: full-3D panel + viscous-coupled wing simulator.

Submodules
----------
geometry : NACA wing mesh generation, STL import, mesh decimation
panel    : Morino Dirichlet source/doublet panel solver (Hess & Smith family)
ibl      : surface-streamline integral boundary layer (Thwaites/Michel/Head)
coupled  : viscous-inviscid interaction via transpiration (blowing) velocity
forces   : pressure integration, Trefftz induced drag, Squire-Young drag
openfoam : OpenFOAM case generation + Docker result parsing/comparison
viz      : visualizer set (surface Cp, cuts, polars, convergence, VTK export)
"""

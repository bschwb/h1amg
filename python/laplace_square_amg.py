"""Solve simple laplace equation on a square."""
# pylint: disable=no-member

from ctypes import CDLL

import ngsolve as ngs
from ngsolve import grad
from netgen.geom2d import unit_square

CDLL('libh1amg.so')

with ngs.TaskManager():
    mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=0.2))

    fes = ngs.H1(mesh, dirichlet=[1, 2, 3], order=1)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # rhs
    f = ngs.LinearForm(fes)
    f += ngs.SymbolicLFI(v)

    # lhs
    a = ngs.BilinearForm(fes, symmetric=True)
    a += ngs.SymbolicBFI(grad(u) * grad(v))

    c = ngs.Preconditioner(a, 'h1amg', flags={'test': True, 'semi_smoothed': True})

    gfu = ngs.GridFunction(fes)

    bvp = ngs.BVP(bf=a, lf=f, gf=gfu, pre=c)

    while True:
        fes.Update()
        gfu.Update()
        a.Assemble()
        f.Assemble()

        bvp.Do()
        ngs.Draw(gfu, mesh, 'solution')

        input('Hit enter for refinement')
        mesh.Refine()

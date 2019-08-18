"""Solve complex equation on a square."""
# pylint: disable=no-member

from ctypes import CDLL

from netgen.geom2d import unit_square

import ngsolve as ngs
from ngsolve import grad

CDLL('libh1amg.so')


with ngs.TaskManager():
    mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=0.15))

    fes = ngs.H1(mesh, dirichlet=[1, 2, 3], order=1, complex=True)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # rhs
    f = ngs.LinearForm(fes)
    f += ngs.SymbolicLFI(v)

    # lhs
    a = ngs.BilinearForm(fes, symmetric=True)
    a += ngs.SymbolicBFI(grad(u) * grad(v) + 1j * u * v)

    c = ngs.Preconditioner(a, 'h1amg2')
    # c = ngs.Preconditioner(a, 'direct')

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

"""Solve simple elliptic equation on a cube."""
# pylint: disable=no-member

from ctypes import CDLL

import ngsolve as ngs
from ngsolve import x, grad
from netgen.libngpy._meshing import NgException

from cube_geo import cube_geo

CDLL('libh1amg.so')
# ngs.ngsglobals.pajetrace = 100000000

with ngs.TaskManager():
    try:
        mesh = ngs.Mesh('cube.vol')
    except NgException:
        ngmesh = cube_geo().GenerateMesh(maxh=0.1)
        ngmesh.Save('cube.vol')
        mesh = ngs.Mesh('cube.vol')

    print(mesh.GetMaterials())
    print(mesh.GetBoundaries())

    fes = ngs.H1(mesh, dirichlet='dirichlet', order=1)
    print('Dofs:', fes.ndof)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # rhs
    f = ngs.LinearForm(fes)
    f += ngs.SymbolicLFI(x*x*x*x * v)
    f.Assemble()

    # lhs
    a = ngs.BilinearForm(fes, symmetric=False)
    a += ngs.SymbolicBFI(grad(u) * grad(v) + u * v)

    c = ngs.Preconditioner(a, 'h1amg')

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

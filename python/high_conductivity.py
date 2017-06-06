"""Solve simple laplace equation on a square."""
# pylint: disable=invalid-name, no-member

from ctypes import CDLL

import netgen.geom2d as geom2d

import ngsolve as ngs
from ngsolve import grad

CDLL('libh1amg.so')


def square_conductivity_geo():
    """Return geometry for a square with high conductivity area."""
    geo = geom2d.SplineGeometry()  # pylint: disable=no-member

    p1, p2, p3, p4 = [geo.AppendPoint(x, y) for x, y in
                      [(0, 0), (1, 0), (1, 1), (0, 1)]]

    geo.Append(['line', p1, p2], leftdomain=1, rightdomain=0, bc='dirichlet')
    geo.Append(['line', p2, p3], leftdomain=1, rightdomain=0, bc='dirichlet')
    geo.Append(['line', p3, p4], leftdomain=1, rightdomain=0, bc='dirichlet')
    geo.Append(['line', p4, p1], leftdomain=1, rightdomain=0)

    geo.AddRectangle((0.2, 0.6), (0.8, 0.7), leftdomain=2, rightdomain=1)
    return geo


with ngs.TaskManager():

    square_conductivity = square_conductivity_geo()
    mesh = ngs.Mesh(square_conductivity.GenerateMesh(maxh=0.05))

    fes = ngs.H1(mesh, dirichlet='dirichlet', order=1)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # rhs
    f = ngs.LinearForm(fes)
    f += ngs.SymbolicLFI(v)

    # lhs
    lam = ngs.DomainConstantCF([1, 1000])
    a = ngs.BilinearForm(fes, symmetric=True)
    a += ngs.SymbolicBFI(lam * grad(u) * grad(v))
    c = ngs.Preconditioner(a, 'h1amg', flags={'test': True, 'levels': 10})

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

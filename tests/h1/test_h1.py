"""Test solving of h1 equations by amg on a square."""
# pylint: disable=no-member, invalid-name

from ctypes import CDLL

from nose.tools import assert_less_equal, assert_greater

import ngsolve as ngs
from ngsolve import grad
from netgen.geom2d import unit_square


def setup_module():
    """Load library"""
    CDLL('libh1amg.so')


def test_h1_real():
    """Test h1 amg for real example."""
    with ngs.TaskManager():
        mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=0.2))

        fes = ngs.H1(mesh, dirichlet=[1, 2, 3], order=1)

        u = fes.TrialFunction()
        v = fes.TestFunction()

        # rhs
        f = ngs.LinearForm(fes)
        f += ngs.SymbolicLFI(v)
        f.Assemble()

        # lhs
        a = ngs.BilinearForm(fes, symmetric=True)
        a += ngs.SymbolicBFI(grad(u) * grad(v))

        c = ngs.Preconditioner(a, 'h1amg2')
        a.Assemble()

        solver = ngs.CGSolver(mat=a.mat, pre=c.mat)

        gfu = ngs.GridFunction(fes)
        gfu.vec.data = solver * f.vec

    assert_greater(solver.GetSteps(), 0)
    assert_less_equal(solver.GetSteps(), 4)

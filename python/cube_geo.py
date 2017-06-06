"""Generate a cube mesh for amg hcurl testing."""
# pylint: disable=no-member

import netgen.csg as csg
import ngsolve as ngs

with ngs.TaskManager():
    origin = csg.Pnt(0, 0, 0)

    side = 1
    box = csg.OrthoBrick(origin, csg.Pnt(side, 2*side, side)).bc('cube_outer')

    normal_vec = csg.Vec3d(0, 1, 0)
    topplane = csg.Plane(csg.Pnt(0, 1, 0), normal_vec).bc('dirichlet')

    cube = (box * topplane).mat('cube_mat')

    cube_geom = csg.CSGeometry()
    cube_geom.Add(cube)

    mesh = cube_geom.GenerateMesh(maxh=0.02)
    mesh.Save('cube.vol')

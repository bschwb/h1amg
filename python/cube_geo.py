"""Generate a cube mesh for amg hcurl testing."""
# pylint: disable=no-member

import netgen.csg as csg
import ngsolve as ngs


def cube_geo():
    origin = csg.Pnt(0, 0, 0)

    side = 1
    box = csg.OrthoBrick(origin, csg.Pnt(side, 2*side, side)).bc('cube_outer')

    normal_vec = csg.Vec(0, 1, 0)
    topplane = csg.Plane(csg.Pnt(0, 1, 0), normal_vec).bc('dirichlet')

    cube = (box * topplane).mat('cube_mat')

    cube_geom = csg.CSGeometry()
    cube_geom.Add(cube)
    return cube_geom

if __name__ == '__main__':
    with ngs.TaskManager():
        mesh = cube_geo().GenerateMesh(maxh=0.1)
        mesh.Save('cube.vol')

from netgen.geom2d import SplineGeometry

n_layers = 5
max_h = 0.1

geo = SplineGeometry()

regs = ['even', 'odd']
l = 1
c = 1/n_layers
d = 1/n_layers
geo.AddRectangle((-c, -c), (c, c), leftdomain=1+l, rightdomain=2-l,
    bcs=('int_bnd', 'int_bnd', 'int_bnd', 'int_bnd'))
l = 1-l

for layer in range(n_layers-2):
    c = c+d
    geo.AddRectangle((-c, -c), (c, c), leftdomain=1+l, rightdomain=2-l,
        bcs=('int_bnd', 'int_bnd', 'int_bnd', 'int_bnd'))
    l = 1-l

c = c+d
geo.AddRectangle((-c, -c), (c, c), leftdomain=1+l, rightdomain=0,
        bcs=('ext_bnd', 'ext_bnd', 'ext_bnd', 'ext_bnd'))

mesh = geo.GenerateMesh(maxh=max_h)

mesh.Save('2d_matrioshka_l{}_h0_1.vol'.format(n_layers))

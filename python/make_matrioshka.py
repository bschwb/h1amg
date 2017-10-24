from netgen.csg import CSGeometry, OrthoBrick, Point3d


n_layers = 5
max_h = 0.1

geo = CSGeometry()

regs = ['even', 'odd']
l = 1
c = 1/n_layers
d = 1/n_layers
brick0 = OrthoBrick(Point3d(-c, -c, -c), Point3d(c, c, c)).bc('inner_bnd')
geo.Add(brick0.mat(regs[l]))
l = 1-l

last_brick = brick0
for layer in range(n_layers-1):
    c = c+d
    if layer < n_layers-2:
        brick = OrthoBrick(Point3d(-c, -c, -c), Point3d(c, c, c)).bc('inner_bnd')
    else:
        brick = OrthoBrick(Point3d(-c, -c, -c), Point3d(c, c, c)).bc('ext_bnd')
    geo.Add((brick-last_brick).mat(regs[l]))
    l = 1-l
    last_brick = brick

mesh = geo.GenerateMesh(maxh=max_h)

mesh.Save('3d_matrioshka_l'+str(n_layers)+'_h0_1.vol')

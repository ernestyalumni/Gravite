# Gravite
Notes and symbolic and numerical computations and implementations on Gravity

Includes: 
* Gravity_Notes_grande.tex
* /pdfs/Gravity_Notes_grande.pdf
* `Rn.sage` 

## Rn.sage - Euclidean spaces as manifolds using [sagemanifolds](http://sagemanifolds.obspm.fr)
*Features*
* R^2,R^3,R^n as a manifold, with a chart atlas

```
sage: load(‘‘Rn.sage’’)  
sage: R2eg = R2() 
sage: R3eg = R3() 
sage: R4 = Rn(4) 
sage: R2eg.M.atlas()
[Chart (R2, (x, y)), Chart (U, (x, y)), Chart (U, (r, ph))]
sage: R3eg.M.atlas()
[Chart (R3, (x, y, z)),
 Chart (U, (x, y, z)),
 Chart (U, (rh, th, ph)),
 Chart (U, (r, phi, zc))]
sage: R4.M.atlas()
[Chart (R4, (x1, x2, x3, x4)),
 Chart (U, (x1, x2, x3, x4)),
 Chart (U, (rh, th1, th2, ph)),
 Chart (U, (r, the1, phi, z))]
```

* (carefully) define a spherical coordinate and cylindrical coordinate chart on Euclidean spaces, e.g.

```
sage: R2eg.transit_sph_to_cart.display() 
x = r*cos(ph)
y = r*sin(ph)
sage: R3eg.transit_sph_to_cart.display() 
x = rh*cos(ph)*sin(th)
y = rh*sin(ph)*sin(th)
z = rh*cos(th)
sage: R3eg.transit_cyl_to_cart.display() 
x = r*cos(phi)
y = r*sin(phi)
z = zc
sage: R4.transit_sph_to_cart.display() 
x1 = rh*cos(ph)*sin(th1)*sin(th2)
x2 = rh*sin(ph)*sin(th1)*sin(th2)
x3 = rh*cos(th2)*sin(th1)
x4 = rh*cos(th1)
```

* calculate the **Jacobian!**

```
sage: to_orthonormal2 , e2, Jacobians2 = R2eg.make_orthon_frames(R2eg.sph_ch) 
sage: Jacobians2[0].inverse()[:,R2eg.sph_ch]
[ cos(ph) -r*sin(ph)]
[ sin(ph) r*cos(ph)]
sage: to_orthonormal3sph, e3sph, Jacobians3sph = R3eg.make_orthon_frames(R3eg.sph_ch) 
sage: to_orthonormal3cyl, e3cyl, Jacobians3cyl = R3eg.make_orthon_frames(R3eg.cyl_ch) 
sage: Jacobians3sph[0].inverse()[:,R3eg.sph_ch]
[ cos(ph)*sin(th) rh*cos(ph)*cos(th) -rh*sin(ph)*sin(th)]
[ sin(ph)*sin(th) rh*cos(th)*sin(ph) rh*cos(ph)*sin(th) ]
[ cos(th)         -rh*sin(th)        0                  ]
sage: Jacobians3cyl[0].inverse()[:,R3eg.cyl_ch]
[ cos(phi) -r*sin(phi) 0]
[ sin(phi) r*cos(phi) 0] 
[ 0 0 1]
```

* equip the Euclidean space manifold with a *metric* *g* and calculate the metric automatically:

```
sage: R2eg.equip_metric() 
sage: R3eg.equip_metric() 
sage: R4.equip_metric()

sage: R2eg.g.display(R2eg.sph_ch.frame(),R2eg.sph_ch)
g = dr*dr + r^2 dph*dph
sage: R3eg.g.display(R3eg.sph_ch.frame(),R3eg.sph_ch)
g = drh*drh + rh^2 dth*dth + rh^2*sin(th)^2 dph*dph
sage: R3eg.g.display(R3eg.cyl_ch.frame(),R3eg.cyl_ch)
g = dr*dr + r^2 dphi*dphi + dzc*dzc
sage: R4.g.display(R4.sph_ch.frame(),R4.sph_ch)
g = drh*drh + rh^2 dth1*dth1 + rh^2*sin(th1)^2 dth2*dth2 + rh^2*sin(th1)^2*sin(th2)^2 dph*dph
sage: R4.g.display(R4.cyl_ch.frame(),R4.cyl_ch)
g = dr*dr + r^2 dthe1*dthe1 + r^2*sin(the1)^2 dphi*dphi + dz*dz
```
* Calculate the so-called orthonormal non-coordinate basis vectors in terms of the (local) coordinate basis vectors, showing clearly and distinctively the difference between the two (concepts)

```
sage: e2[1].display( R2eg.sph_ch.frame(), R2eg.sph_ch)
e_1 = d/dr
sage: e2[2].display( R2eg.sph_ch.frame(), R2eg.sph_ch)
e_2 = 1/r d/dph
sage: for i in range(1,3+1):                                                         
    e3sph[i].display( R3eg.sph_ch.frame(), R3eg.sph_ch )
....:     
e_1 = d/drh
e_2 = 1/rh d/dth
e_3 = 1/(rh*sin(th)) d/dph
sage: for i in range(1,3+1):
    e3cyl[i].display( R3eg.cyl_ch.frame(), R3eg.cyl_ch )
....:     
e_1 = d/dr
e_2 = 1/r d/dphi
e_3 = d/dzc
```

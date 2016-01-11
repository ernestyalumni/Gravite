# Gravite
Notes and symbolic and numerical computations and implementations on Gravity

Includes: 
* Gravity_Notes_grande.tex
* /pdfs/Gravity_Notes_grande.pdf

## Rn.sage Euclidean spaces as manifolds using [sagemanifolds](http://sagemanifolds.obspm.fr)
*Features*
* R^2,R^3,R^n as a manifold, with a chart atlas

```
sage: load(‘‘Rn.sage’’)  
sage: R2eg = R2() 
sage: R3eg = R3() 
sage: R4 = Rn(4) 
```

* (carefully) define a spherical coordinate and cylindrical coordinate chart on Euclidean spaces, e.g.

```
sage: R2eg.transit_sph_to_cart.display() x = r*cos(ph)
y = r*sin(ph)
sage: R3eg.transit_sph_to_cart.display() x = rh*cos(ph)*sin(th)
y = rh*sin(ph)*sin(th)
z = rh*cos(th)
sage: R3eg.transit_cyl_to_cart.display() x = r*cos(phi)
y = r*sin(phi)
z = zc
sage: R4.transit_sph_to_cart.display() x1 = rh*cos(ph)*sin(th1)*sin(th2)
x2 = rh*sin(ph)*sin(th1)*sin(th2)
x3 = rh*cos(th2)*sin(th1)
x4 = rh*cos(th1)
```

* calculate the **Jacobian!**

```
sage: to_orthonormal2 , e2, Jacobians2 = R2eg.make_orthon_frames(R2eg.sph_ch) sage: Jacobians2[0].inverse()[:,R2eg.sph_ch]
[ cos(ph) -r*sin(ph)]
[ sin(ph) r*cos(ph)]
sage: to_orthonormal3sph, e3sph, Jacobians3sph = R3eg.make_orthon_frames(R3eg.sph_ch) sage: to_orthonormal3cyl, e3cyl, Jacobians3cyl = R3eg.make_orthon_frames(R3eg.cyl_ch) sage: Jacobians3sph[0].inverse()[:,R3eg.sph_ch]
[ cos(ph)*sin(th) rh*cos(ph)*cos(th) -rh*sin(ph)*sin(th)]
[ sin(ph)*sin(th) rh*cos(th)*sin(ph) rh*cos(ph)*sin(th) ]
[ cos(th)         -rh*sin(th)        0                  ]
sage: Jacobians3cyl[0].inverse()[:,R3eg.cyl_ch]
[ cos(phi) -r*sin(phi) 0]
[ sin(phi) r*cos(phi) 0] 
[ 0 0 1]
```

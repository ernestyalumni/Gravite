# Gravite - Gravity, Gravité
Notes and symbolic and numerical computations and implementations on Gravity

Includes: 
* Gravity_Notes_grande.tex
* /pdfs/Gravity_Notes_grande.pdf
* `Rn.sage`
* Installation of `sagemanifolds` into Sage Math

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

### Installation of `sagemanifolds` into Sage Math

#### EY : 20160611 note:

I'm on Fedora Linux, Fedora 23 Workstation, 64-bit (`x86_64`).  I have successfully installed Sage Math 7.2 developer's version from its github source.  Then I successfully installed `sagemanifolds` following the `sagemanifolds`'s webpage's [instructions](http://sagemanifolds.obspm.fr/download.html) for installation on Linux, that is *not* Debian Live.

However, I am unsuccess at installing `sagemanifolds` onto Sage Math if Sage Math is a prebuilt Linux version.  I had followed the [instructions](http://doc.sagemath.org/html/en/installation/binary.html) and downloaded the appropriate version for my "distro" (Fedora) and architecture (64-bit), i.e. `sage-7.2-Fedora_23-x86_64.tar.bz2` from [here](http://files.sagemath.org/linux/64bit/index.html).  I unpacked it, with Archive Manager, extracting it into my user account (I'm not on admin).  At the first time, I ran `./sage` and it built some `.so` files.  It runs fine at this point, but without `sagemanifolds`.  Then I follow the same instructions as I had eluded to above for `sagemanifolds` installation.

However, as `sagemanifolds` was building, and it is a very lengthy build, seeming like at least 440 builds must be done, by 129, I am receiving a fatal error and by 130, it stops; here it is:

```
[129/440] creating build/temp.linux-x86_64-2.7/home/topolo/Public/SageMath/src/build/cythonized/sage/libs/fplll
gcc -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wno-unused -fPIC -I/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/cysignals -I/home/topolo/Public/SageMath/local/include -I/home/topolo/Public/SageMath/local/include/python2.7 -I/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/numpy-1.11.0-py2.7-linux-x86_64.egg/numpy/core/include -I/home/topolo/Public/SageMath/src -I/home/topolo/Public/SageMath/src/sage/ext -I/home/topolo/Public/SageMath/src/build/cythonized -I/home/topolo/Public/SageMath/src/build/cythonized/sage/ext -I/home/topolo/Public/SageMath/local/include/python2.7 -c /home/topolo/Public/SageMath/src/build/cythonized/sage/libs/fplll/fplll.cpp -o build/temp.linux-x86_64-2.7/home/topolo/Public/SageMath/src/build/cythonized/sage/libs/fplll/fplll.o -DFPLLL_V3_COMPAT -fno-strict-aliasing
In file included from /home/topolo/Public/SageMath/local/include/fplll/nr.h:8:0,
                 from /home/topolo/Public/SageMath/src/build/cythonized/sage/libs/fplll/fplll.cpp:346:
/home/topolo/Public/SageMath/local/include/fplll/defs.h:29:0: warning: "FPLLL_V3_COMPAT" redefined
 #define FPLLL_V3_COMPAT
 ^
<command-line>:0:0: note: this is the location of the previous definition
In file included from /home/topolo/Public/SageMath/local/include/fplll/nr.h:26:0,
                 from /home/topolo/Public/SageMath/src/build/cythonized/sage/libs/fplll/fplll.cpp:346:
/home/topolo/Public/SageMath/local/include/fplll/nr_FP_dd.inl:9:24: fatal error: qd/dd_real.h: No such file or directory
compilation terminated.
[130/440] creating build/temp.linux-x86_64-2.7/home/topolo/Public/SageMath/src/build/cythonized/sage/libs/glpk
gcc -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wno-unused -fPIC -I/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/cysignals -I/home/topolo/Public/SageMath/local/include -I/home/topolo/Public/SageMath/local/include/python2.7 -I/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/numpy-1.11.0-py2.7-linux-x86_64.egg/numpy/core/include -I/home/topolo/Public/SageMath/src -I/home/topolo/Public/SageMath/src/sage/ext -I/home/topolo/Public/SageMath/src/build/cythonized -I/home/topolo/Public/SageMath/src/build/cythonized/sage/ext -I/home/topolo/Public/SageMath/local/include/python2.7 -c /home/topolo/Public/SageMath/src/build/cythonized/sage/libs/glpk/error.c -o build/temp.linux-x86_64-2.7/home/topolo/Public/SageMath/src/build/cythonized/sage/libs/glpk/error.o -fno-strict-aliasing
error: command 'gcc' failed with exit status 1
Makefile:6: recipe for target 'sage' failed
make: *** [sage] Error 1

Installation of SageManifolds 0.9 completed!

```

I try to run Sage Math with `./sage` in Sage Math's "root" directory.  Non-sagemanifolds function work fine.  However, I obtain this error when I type into sage, "Manifold":

```
sage: Manifold
---------------------------------------------------------------------------
ImportError                               Traceback (most recent call last)
/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/core/formatters.pyc in __call__(self, obj)
    697                 type_pprinters=self.type_printers,
    698                 deferred_pprinters=self.deferred_printers)
--> 699             printer.pretty(obj)
    700             printer.flush()
    701             return stream.getvalue()

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/lib/pretty.pyc in pretty(self, obj)
    381                             if callable(meth):
    382                                 return meth(obj, self, cycle)
--> 383             return _default_pprint(obj, self, cycle)
    384         finally:
    385             self.end_group()

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/lib/pretty.pyc in _default_pprint(obj, p, cycle)
    501     if _safe_getattr(klass, '__repr__', None) not in _baseclass_reprs:
    502         # A user-provided repr. Find newlines and replace them with p.break_()
--> 503         _repr_pprint(obj, p, cycle)
    504         return
    505     p.begin_group(1, '<')

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/lib/pretty.pyc in _repr_pprint(obj, p, cycle)
    692     """A pprint that just redirects to the normal repr function."""
    693     # Find newlines and replace them with p.break_()
--> 694     output = repr(obj)
    695     for idx,output_line in enumerate(output.splitlines()):
    696         if idx:

/home/topolo/Public/SageMath/src/sage/misc/lazy_import.pyx in sage.misc.lazy_import.LazyImport.__repr__ (/home/topolo/Public/SageMath/src/build/cythonized/sage/misc/lazy_import.c:3694)()
    398             'Integer Ring'
    399         """
--> 400         return repr(self._get_object())
    401 
    402     def __str__(self):

/home/topolo/Public/SageMath/src/sage/misc/lazy_import.pyx in sage.misc.lazy_import.LazyImport._get_object (/home/topolo/Public/SageMath/src/build/cythonized/sage/misc/lazy_import.c:2232)()
    244         elif self._at_startup and not startup_guard:
    245             print('Option ``at_startup=True`` for lazy import {0} not needed anymore'.format(self._name))
--> 246         self._object = getattr(__import__(self._module, {}, {}, [self._name]), self._name)
    247         alias = self._as_name or self._name
    248         if self._deprecation is not None:

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/sage/manifolds/manifold.py in <module>()
    297 from sage.rings.integer import Integer
    298 from sage.manifolds.subset import ManifoldSubset
--> 299 from sage.manifolds.structure import TopologicalStructure, \
    300                                      RealTopologicalStructure
    301 

ImportError: No module named structure
<repr(<sage.misc.lazy_import.LazyImport at 0x7f75cb26cb40>) failed: ImportError: No module named structure>
sage: Manifold
---------------------------------------------------------------------------
ImportError                               Traceback (most recent call last)
/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/core/formatters.pyc in __call__(self, obj)
    697                 type_pprinters=self.type_printers,
    698                 deferred_pprinters=self.deferred_printers)
--> 699             printer.pretty(obj)
    700             printer.flush()
    701             return stream.getvalue()

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/lib/pretty.pyc in pretty(self, obj)
    381                             if callable(meth):
    382                                 return meth(obj, self, cycle)
--> 383             return _default_pprint(obj, self, cycle)
    384         finally:
    385             self.end_group()

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/lib/pretty.pyc in _default_pprint(obj, p, cycle)
    501     if _safe_getattr(klass, '__repr__', None) not in _baseclass_reprs:
    502         # A user-provided repr. Find newlines and replace them with p.break_()
--> 503         _repr_pprint(obj, p, cycle)
    504         return
    505     p.begin_group(1, '<')

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/IPython/lib/pretty.pyc in _repr_pprint(obj, p, cycle)
    692     """A pprint that just redirects to the normal repr function."""
    693     # Find newlines and replace them with p.break_()
--> 694     output = repr(obj)
    695     for idx,output_line in enumerate(output.splitlines()):
    696         if idx:

/home/topolo/Public/SageMath/src/sage/misc/lazy_import.pyx in sage.misc.lazy_import.LazyImport.__repr__ (/home/topolo/Public/SageMath/src/build/cythonized/sage/misc/lazy_import.c:3694)()
    398             'Integer Ring'
    399         """
--> 400         return repr(self._get_object())
    401 
    402     def __str__(self):

/home/topolo/Public/SageMath/src/sage/misc/lazy_import.pyx in sage.misc.lazy_import.LazyImport._get_object (/home/topolo/Public/SageMath/src/build/cythonized/sage/misc/lazy_import.c:2232)()
    244         elif self._at_startup and not startup_guard:
    245             print('Option ``at_startup=True`` for lazy import {0} not needed anymore'.format(self._name))
--> 246         self._object = getattr(__import__(self._module, {}, {}, [self._name]), self._name)
    247         alias = self._as_name or self._name
    248         if self._deprecation is not None:

/home/topolo/Public/SageMath/local/lib/python2.7/site-packages/sage/manifolds/manifold.py in <module>()
    297 from sage.rings.integer import Integer
    298 from sage.manifolds.subset import ManifoldSubset
--> 299 from sage.manifolds.structure import TopologicalStructure, \
    300                                      RealTopologicalStructure
    301 

ImportError: No module named structure
<repr(<sage.misc.lazy_import.LazyImport at 0x7f75cb26cb40>) failed: ImportError: No module named structure>
```  
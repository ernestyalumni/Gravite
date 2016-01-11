## Rn.sage
## Sage Math Sagemanifolds implementation of Euclidean R^n as classes
## namespace or names follow closely the Tutorial pdf on sagemanifolds webpage: 
## http://sagemanifolds.obspm.fr/examples/pdf/SM_tutorial.pdf
############################################################################ 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
## 20160109
##                                                                               
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, or   
## (at your option) any later version.                                        
##     
## This program is distributed in the hope that it will be useful,               
## but WITHOUT ANY WARRANTY; without even the implied warranty of              
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 
## GNU General Public License for more details.                             
##                                                                          
## linkedin     : ernestyalumni                                                    
## wordpress    : ernestyalumni                                                    
############################################################################
t = var('t')
assume(t,"real")

class R1(object):
      def __init__(self):
      	  self.M = Manifold(2,'R1',r'\mathbb{R}',start_index=1)
	  self.cart_ch = self.M.chart('x')
	  

class R2(object):
      def __init__(self):
      	  self.M = Manifold(2,'R2',r'\mathbb{R}^2',start_index=1)
	  self.cart_ch = self.M.chart('x y')
	  self.U       = self.M.open_subset('U',
				coord_def={self.cart_ch: 
					(self.cart_ch[1]<0, self.cart_ch[2]!=0)})	  # cf. http://sagemanifolds.obspm.fr/examples/pdf/SM_tutorial.pdf "Introducing a second chart on the manifold" the condition AND written with [] instead of ()

	  self.sph_ch  = self.U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
	  self.cart_ch_U = self.cart_ch.restrict(self.U)
	  self.transit_sph_to_cart = self.sph_ch.transition_map(self.cart_ch_U, 
	  			   [self.sph_ch[1]*cos(self.sph_ch[2]), 
				    self.sph_ch[1]*sin(self.sph_ch[2])])
	  Eucnorm = sqrt( sum([self.cart_ch_U[i[0]]**2 for i in self.M.index_generator(1)]) )
	  self.transit_sph_to_cart.set_inverse( Eucnorm, 
	  					atan2(self.cart_ch_U[2],self.cart_ch_U[1]) )
      def equip_metric(self):
      	  self.g = self.M.riemannian_metric('g')
	  for i in self.M.index_generator(1):
	      self.g[i[0],i[0]] = 1  

      def make_orthon_frames(self,ch):
      	  """
	  make_orthon_frames=make_orthon_frames(self,ch)
	  This method creates a change-of-basis matrix for orthonormal coordinates as 
	       to_orthonormal
	  and a new orthonormal frame from the input of a (spherical coordinates) chart, ch

	  PARAMETERS (INPUTS):
	  ch - <a sagemanifolds chart>

	  OUTPUT
	  to_orthonormal - a change-of-basis matrix
	  eo - new orthonormal frame 
	  
	  EXAMPLES of USAGE:
	  R2eg = R2()	  
	  R2.equip_metric()
	  R2.make_orthon_frames(R2.sph_ch)	  

	  """
      	  try:
	      to_orthonormal = ch.domain().automorphism_field()
	      for i,j in self.M.index_generator(2):
	      	  if self.g[ch.frame(),i,j,ch]!= 0:
		     g_ij = self.g[ch.frame(),i,j,ch]
		     to_orthonormal[ch.frame(),i,j,ch] = 1/sqrt(g_ij)
	      esph = ch.frame().new_frame(to_orthonormal,'e') 		     

	      # cf. https://github.com/sagemanifolds/SageManifolds/blob/master/Worksheets/v0.9/SM_Cartesian_spherical-3D.ipynb for explanation on the change of frame
	      Jacobian_sph_to_cart = ch.domain().change_of_frame(ch.frame(), ch.domain().default_chart().frame() )
	      Jacobian_och_to_sph  = ch.domain().change_of_frame(esph, ch.frame() )

	      ch.domain().set_change_of_frame(ch.domain().default_frame(), esph, 
	      					Jacobian_och_to_sph.inverse()*Jacobian_sph_to_cart.inverse())

	      ch.domain().set_change_of_frame( esph, ch.domain().default_frame(),
	      				       Jacobian_sph_to_cart*Jacobian_och_to_sph )

	      return to_orthonormal, esph, [Jacobian_sph_to_cart, Jacobian_och_to_sph]

	  except AttributeError:
	      print "Equip the manifold with a metric g by doing the method .equip_metric() first!"


class R3(object):
      def __init__(self):
      	  self.M = Manifold(3,'R3',r'\mathbb{R}^3',start_index=1)
	  self.cart_ch = self.M.chart('x y z')
	  
	  self.U = self.M.open_subset('U',coord_def={self.cart_ch: (self.cart_ch[1]<0, self.cart_ch[2]!=0)})
	  self.cart_ch_U = self.cart_ch.restrict(self.U)
	  
	  self.sph_ch  = self.U.chart(r'rh:(0,+oo):\rho th:(0,pi):\theta ph:(0,2*pi):\phi')
	  rh,th,ph = [self.sph_ch[i[0]] for i in self.M.index_generator(1)] 
	  self.transit_sph_to_cart = self.sph_ch.transition_map(self.cart_ch_U,
					[rh*sin(th)*cos(ph),rh*sin(th)*sin(ph),rh*cos(th)])
	  Sphnorm = sqrt(sum([self.cart_ch_U[i[0]]**2 for i in self.M.index_generator(1)]))
	  self.transit_sph_to_cart.set_inverse( Sphnorm,atan2( sqrt( sum([ self.cart_ch_U[i]**2 for i in range(1,3)])), self.cart_ch_U[3]), atan2( self.cart_ch_U[2],self.cart_ch_U[1]) )

	  
	  self.cyl_ch = self.U.chart(r'r:(0,+oo) phi:(0,2*pi):\phi zc')
	  r,phi,zc = [self.cyl_ch[i[0]] for i in self.M.index_generator(1)]
	  self.transit_cyl_to_cart = self.cyl_ch.transition_map(self.cart_ch_U,
	  			        [r*cos(phi),r*sin(phi),zc])
	  self.transit_cyl_to_cart.set_inverse( sqrt( self.cart_ch_U[1]**2 + self.cart_ch_U[2]**2 ) , atan2( self.cart_ch_U[2],self.cart_ch_U[1]), self.cart_ch_U[3] )
	  
      def equip_metric(self):
      	  self.g = self.M.riemannian_metric('g')
	  for i in self.M.index_generator(1):
	      self.g[i[0],i[0]] = 1

      def make_orthon_frames(self,ch):
      	  """
	  make_orthon_frames=make_orthon_frames(self,ch)
	  This method creates a change-of-basis matrix for orthonormal coordinates as 
	       to_orthonormal
	  and a new orthonormal frame from the input of a (spherical coordinates) chart, ch

	  PARAMETERS (INPUTS):
	  ch - <a sagemanifolds chart>

	  OUTPUT
	  to_orthonormal - a change-of-basis matrix
	  eo - new orthonormal frame 
	  
	  EXAMPLES of USAGE:
	  R3eg = R3()	  
	  R3.equip_metric()
	  R3.make_orthon_frames(R3.sph_ch)	  

	  """
      	  try:
	      to_orthonormal = ch.domain().automorphism_field()
	      for i,j in self.M.index_generator(2):
	      	  if self.g[ch.frame(),i,j,ch]!= 0:
		     g_ij = self.g[ch.frame(),i,j,ch]
		     to_orthonormal[ch.frame(),i,j,ch] = 1/sqrt(g_ij)
	      eoch = ch.frame().new_frame(to_orthonormal,'e') 		     
	      
	      # cf. https://github.com/sagemanifolds/SageManifolds/blob/master/Worksheets/v0.9/SM_Cartesian_spherical-3D.ipynb for explanation on the change of frame
	      Jacobian_ch_to_cart = ch.domain().change_of_frame(ch.frame(), ch.domain().default_chart().frame() )
	      Jacobian_och_to_ch  = ch.domain().change_of_frame(eoch, ch.frame() )

	      ch.domain().set_change_of_frame(ch.domain().default_frame(), eoch, 
	      					Jacobian_och_to_ch.inverse()*Jacobian_ch_to_cart.inverse())

	      ch.domain().set_change_of_frame( eoch, ch.domain().default_frame(),
	      				       Jacobian_ch_to_cart*Jacobian_och_to_ch )


	      return to_orthonormal, eoch, [Jacobian_ch_to_cart, Jacobian_och_to_ch]

	  except AttributeError:
	      print "Equip the manifold with a metric g by doing the method .equip_metric() first!"


class Rn(object):
      def __init__(self,n):
      	  assert n>0
      	  if n == 2:
	     print "Use the class R2"
	  elif n == 3:
	     print "Use the class R3"
	  else:
	     self.M = Manifold(n,'R'+str(n),r'\mathbb{R}^'+str(n),start_index=1)
	     self.cart_ch = self.M.chart(r" ".join([r"x"+str(i) for i in range(1,n+1)]))
	     xis = [self.cart_ch[i[0]] for i in self.M.index_generator(1)]

	     self.U = self.M.open_subset('U',coord_def={self.cart_ch:(xis[0]<0,xis[1]!=0)})
	     self.cart_ch_U = self.cart_ch.restrict(self.U)

	     # spherical coordinates	     
	     self.sph_ch = self.U.chart(r'rh:(0,+oo):\rho '+r" ".join([r"th"+str(i)+r":(0,pi)" for i in range(1,n+1-2)])+r' ph:(0,2*pi):\phi')
	     sphs = [self.sph_ch[i[0]] for i in self.M.index_generator(1)]
	     self.transit_sph_to_cart = self.sph_ch.transition_map(self.cart_ch_U,
	     			      [sphs[0]*prod([sin(sphs[i]) for i in range(1,n-1)])*cos(sphs[-1]), sphs[0]*prod([sin(sphs[i]) for i in range(1,n-1)])*sin(sphs[-1])]+
				       [sphs[0]*prod([sin(sphs[i]) for i in range(1,j)])*cos(sphs[j]) for j in range(n-2,1,-1)]+[sphs[0]*cos(sphs[1]),])
	     gen_transit_list_sph = [ sqrt(sum([ xis[i]**2 for i in range(len(xis))])),] + [atan2( sqrt( sum([xis[i]**2 for i in range(j)])),xis[j]) for j in range(n-1,1,-1)]+[atan2(xis[1],xis[0]),] 
	     self.transit_sph_to_cart.set_inverse(*gen_transit_list_sph)

	     # cylindrical coordinates
	     self.cyl_ch = self.U.chart(r'r:(0,+oo) '+r" ".join([r"the"+str(i)+r":(0,pi)" for i in range(1,n+1-3)])+r' phi:(0,2*pi):\varphi z')
	     cyls = [self.cyl_ch[i[0]] for i in self.M.index_generator(1)]
	     self.transit_cyl_to_cart = self.cyl_ch.transition_map(self.cart_ch_U, [cyls[0]*prod([sin(cyls[i]) for i in range(1,n-2)])*cos(cyls[-2]), cyls[0]*prod([sin(cyls[i]) for i in range(1,n-2)])*sin(cyls[-2])]+ [cyls[0]*prod([sin(cyls[i]) for i in range(1,j)])*cos(cyls[j]) for j in range(n-3,1,-1)]+[cyls[0]*cos(cyls[1]),cyls[-1]] ) 
	     gen_transit_list_cyl = [ sqrt(sum([ xis[i]**2 for i in range(len(xis)-1)])),]+[atan2(sqrt( sum([xis[i]**2 for i in range(j)])),xis[j]) for j in range(n-2,1,-1)]+[atan2(xis[1],xis[0]),xis[n-1]]
	     self.transit_cyl_to_cart.set_inverse(*gen_transit_list_cyl)

      def equip_metric(self):
      	  self.g=self.M.riemannian_metric('g')
	  for i in self.M.index_generator(1):
	      self.g[i[0],i[0]]=1	     

      def make_orthon_frames(self,ch):
      	  """
	  make_orthon_frames=make_orthon_frames(self,ch)
	  This method creates a change-of-basis matrix for orthonormal coordinates as 
	       to_orthonormal
	  and a new orthonormal frame from the input of a (spherical coordinates) chart, ch

	  PARAMETERS (INPUTS):
	  ch - <a sagemanifolds chart>

	  OUTPUT
	  to_orthonormal - a change-of-basis matrix
	  eoch - new orthonormal frame 
	  
	  EXAMPLES of USAGE:
	  R4 = Rn(4)	  
	  R4.equip_metric()
	  R4.make_orthon_frames(R4.sph_ch)	  

	  """
      	  try:
	      to_orthonormal = ch.domain().automorphism_field()
	      for i,j in self.M.index_generator(2):
	      	  if self.g[ch.frame(),i,j,ch]!= 0:
		     g_ij = self.g[ch.frame(),i,j,ch]
		     to_orthonormal[ch.frame(),i,j,ch] = 1/sqrt(g_ij)
	      eoch = ch.frame().new_frame(to_orthonormal,'e') 		     
	      
	      # cf. https://github.com/sagemanifolds/SageManifolds/blob/master/Worksheets/v0.9/SM_Cartesian_spherical-3D.ipynb for explanation on the change of frame
	      Jacobian_ch_to_cart = ch.domain().change_of_frame(ch.frame(), ch.domain().default_chart().frame() )
	      Jacobian_och_to_ch  = ch.domain().change_of_frame(eoch, ch.frame() )

	      ch.domain().set_change_of_frame(ch.domain().default_frame(), eoch, 
	      					Jacobian_och_to_ch.inverse()*Jacobian_ch_to_cart.inverse())

	      ch.domain().set_change_of_frame( eoch, ch.domain().default_frame(),
	      				       Jacobian_ch_to_cart*Jacobian_och_to_ch )


	      return to_orthonormal, eoch, [Jacobian_ch_to_cart,Jacobian_och_to_ch]

	  except AttributeError:
	      print "Equip the manifold with a metric g by doing the method .equip_metric() first!"
      
	  
def make_pt(ch):
    """
    make_pt = make_pt(ch)
    INPUT
    ch = sagemanifold chart

    EXAMPLES of USAGE
    p = make_pt(R3.cart_ch)
    """
    coords = ch[:]
    farglst = ['p',]+list(coords)
    p = ch.scalar_field( function(*farglst) )
    return p


def make_u(ch):
    """
    make_u = make_u(ch)
    make_u creates a time-INDEPENDENT velocity vector field
    
    INPUT
    ch = sage manifold chart

    EXAMPLEs of USAGE:
    R2 = Rd(2)
    u2 = make_u(R2.X_U)
   
    R3 = Rd(3)
    u3 = make_u(R3.X_U)

    u3[1].expr().diff(t) # 0 ; this demonstrates that this velocity vector is time-INDEPENDENT
    """				  
    n_0 = ch.domain().manifold().dim()
    # ucomplst components list of u  				  
    ucomplst = []
    for i in ch.domain().manifold().index_generator(1):
    	farglst = ['u'+str(i[0]),] + list(ch[:]) 
	ucomplst.append( function( *farglst ) )
    u = ch.domain().vector_field()
    u[ch.frame(),:,ch] = ucomplst
    return u


def make_ut(ch):
    """
    make_ut = make_ut(ch)
    
    INPUT
    ch = sage manifold chart

    EXAMPLEs of USAGE:
    R2 = Rd(2)
    ut2 = make_ut(R2.X_U)
   
    R3 = Rd(3)
    ut3 = make_ut(R3.X_U)
    """				  
    n_0 = ch.domain().manifold().dim()
    # ucomplst components list of u  				  
    ucomplst = []
    for i in ch.domain().manifold().index_generator(1):
    	farglst = ['u'+str(i[0]),] + [t,] + list(ch[:]) 
	ucomplst.append( function( *farglst ) )
    u = ch.domain().vector_field()
    u[ch.frame(),:,ch] = ucomplst
    return u


def make_material_der(u, ch):
    """
    make_material_der = make_material_der(u,ch)

    EXAMPLES of USAGE:
    R3 = Rd(3)
    u3t = make_ut(R3.X_U)
    udu = make_material_der(u3t, R3.X_U)
    """
    uedcomp = []
    for ui in u[ch.frame(),:,ch]:
    	uidict = dict( [(ch,ui),])
	uedcomp.append( u( ch.domain().scalar_field( uidict ) ))
    X = sum( [ uedcomp[i[0]-1]*ch.frame()[i[0]] for i in ch.domain().manifold().index_generator(1) ] )
    return X


def div(u,g):
    """	
    div = div(u,g)
    Return the divergence of vector field u \in \mathfrak{X}(M), given the metric g for the manifold M
    """
    uflat = g['_ij']*u['^j']
    return xder( uflat.hodge_star(g) )

def grad(p,g):
    """
    grad = grad(p,g)

    EXAMPLE of USAGE
    R3 = Rd(3)
    p = make_pt(R3.M)
    grad(p,R3.g)
    """
    dp = xder( p )
    gradp = g.inverse()['^ij']*dp['_j']
    return gradp

def curl(u,g):
    """
    curl = curl(u,g)
    Return the curl of vector field u \in \mathfrak{X}(M), given the metric g for the manifold M
    """
    uflat = g['_ij']*u['^j']
    duflat = xder( uflat )
    return duflat.hodge_star(g)

def buildrho(ch):
    """
    buildrho = buildrho(ch)
    build a time-dependent $\rho$ the mass density, as a scalar function on a chart of a manifold
    
    EXAMPLE of USAGE:
    R2=Rd(2)
    rho2=buildrho(R2.X_U)
    """
    n_0 = ch.domain().manifold().dim()
    variables = [t,]+[ch[i] for i in range(1,n_0+1)]
    rho = ch.domain().scalar_field(function('rho',*variables),name='rho',latex_name=r'\rho' )
    return rho



##############################
## Usage Examples
##############################
"""
R2eg = R2()
R2eg.transit_sph_to_cart.display()
R2eg.equip_metric()
R2eg.g.display(R2eg.sph_ch.frame(),R2eg.sph_ch)

to_orthonormal2, e2, Jacobians2 = R2eg.make_orthon_frames(R2eg.sph_ch)
to_orthonormal2.display(R2eg.sph_ch.frame(),R2eg.sph_ch)
e2[1].display( R2eg.sph_ch.frame(), R2eg.sph_ch)
e2[2].display( R2eg.sph_ch.frame(), R2eg.sph_ch)
Jacobians2[0].inverse()[:,R2eg.sph_ch]
Jacobians2[1].inverse()[R2eg.sph_ch.frame(),:,R2eg.sph_ch]

R3eg = R3()
R3eg.transit_sph_to_cart.display()
R3eg.transit_cyl_to_cart.display()
R3eg.equip_metric()
R3eg.g.display(R3eg.sph_ch.frame(),R3eg.sph_ch)
R3eg.g.display(R3eg.cyl_ch.frame(),R3eg.cyl_ch)

to_orthonormal3sph, e3sph, Jacobians3sph = R3eg.make_orthon_frames(R3eg.sph_ch)
to_orthonormal3cyl, e3cyl, Jacobians3cyl = R3eg.make_orthon_frames(R3eg.cyl_ch)
to_orthonormal3sph.display(R3eg.sph_ch.frame(),R3eg.sph_ch)
to_orthonormal3cyl.display(R3eg.cyl_ch.frame(),R3eg.cyl_ch)

for i in range(1,3+1):	
    e3sph[i].display( R3eg.sph_ch.frame(), R3eg.sph_ch )
for i in range(1,3+1):
    e3cyl[i].display( R3eg.cyl_ch.frame(), R3eg.cyl_ch )
Jacobians3sph[0].inverse()[:,R3eg.sph_ch]
Jacobians3cyl[0].inverse()[:,R3eg.cyl_ch]

R4 = Rn(4)
R4.transit_sph_to_cart.display()
R4.transit_cyl_to_cart.display()
R4.equip_metric()
R4.g.display(R4.sph_ch.frame(),R4.sph_ch)
R4.g.display(R4.cyl_ch.frame(),R4.cyl_ch)

to_orthonormal4sph, e4sph, Jacobians4sph = R4.make_orthon_frames(R4.sph_ch)
to_orthonormal4cyl, e4cyl, Jacobians4cyl = R4.make_orthon_frames(R4.cyl_ch)
to_orthonormal4sph.display(R4.sph_ch.frame(),R4.sph_ch)
to_orthonormal4cyl.display(R4.cyl_ch.frame(),R4.cyl_ch)

for i in range(1,4+1):	
    e4sph[i].display( R4.sph_ch.frame(), R4.sph_ch )
for i in range(1,4+1):
    e4cyl[i].display( R4.cyl_ch.frame(), R4.cyl_ch )
Jacobians4sph[0].inverse()[:,R4.sph_ch]
Jacobians4cyl[0].inverse()[:,R4.cyl_ch]

"""
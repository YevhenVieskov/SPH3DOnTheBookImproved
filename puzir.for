	subroutine puzir(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     bubble rizing
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =2   water
c     h-- smoothing lengths of particles                           [out]
c     ntotal-- total particle number                               [out]
      use interf
      use ex
      implicit none     
      include 'param.inc'
      
      integer itype(:), ntotal
      double precision x(:, :), vx(:, :), mass(:),
     &     rho(:), p(:), u(:), hsml(:)  !, c(:)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy
      double precision radius,cr,xr,yr,r2

c     Giving mass and smoothing length as well as other data.
	radius=1.0
      m =40+1
      n =80+1
      mp = m-1
      np = n-1
      ntotal = mp * np
      xl =6.0*radius    
      yl =10.0*radius   
      dx = xl/mp
      dy = yl/np
      dx = xl/mp
      dy = yl/np

      do i = 1, mp
	do j = 1, np
	  k = j + (i-1)*np
	  x(1, k) = (i-1)*dx + dx/2.
	  x(2, k) = (j-1)*dy + dy/2.
        enddo
      enddo
	
	xr=3.0*radius
	yr=2.0*radius
	r2=radius*radius

      do i = 1, mp*np
	vx(1, i) = 0.
	vx(2, i) = 0.
	cr=(x(1,i)-xr)*(x(1,i)-xr)+(x(2,i)-yr)*(x(2,i)-yr)
	
	if(cr>r2)then     
        rho (i) = 1000. 
	  itype(i) = 2
        u(i)=357.1
      else
        rho (i) = 1.
	  itype(i) = 1 
	  u(i)=357.1 
      end if
        mass(i) = dx*dy*rho(i)  
        p(i)= 9.81d0*rho(i)*(yl-x(2,i))   !+1.0d5   
        
        
        hsml(i) = dx
      enddo
	 
c	initial particle distribution for tensile correction	
	 dpt=dx

      end
	
	
	

	subroutine rt(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     Raley-Teylor instability
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
     &     rho(:), p(:), u(:), hsml(:)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy !,c(:)
      double precision yr

c     Giving mass and smoothing length as well as other data.
	
      m =20+1
      n =40+1
      mp = m-1
      np = n-1
      ntotal = mp * np
      xl =1.0    
      yl =2.0   
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
	
	

      do i = 1, mp*np
	vx(1, i) = 0.
	vx(2, i) = 0.
	
	yr=1.0-sin(2.0*pi*x(1,i))
	if(x(2,i)>yr)then     
        rho (i) = 1.8 
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
	
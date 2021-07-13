      subroutine ext_force(ntotal,mass,x,niac,pair_i,pair_j,
     &           itype,hsml,dvxdt)

c--------------------------------------------------------------------------
c     Subroutine to calculate the external forces, e.g. gravitational forces.      
c     The forces from the interactions with boundary virtual particles 
c     are also calculated here as external forces.

c     here as the external force. 
c     ntotal  : Number of particles                                 [in]
c     mass    : Particle masses                                     [in]
c     x       : Coordinates of all particles                        [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     itype   : type of particles                                   [in]
c     hsml   : Smoothing Length                                     [in]
c     dvxdt   : Acceleration with respect to x, y and z            [out] 
      use interf
      implicit none
      include 'param.inc'
      
      integer ntotal, itype(:), niac,
     &        pair_i(:), pair_j(:)
      double precision mass(:), x(:,:), hsml(:),          
     &       dvxdt(:,:)   !,c(:)
      integer i, j, k, d
      double precision  rr, f, rr0, dd, p1, p2,q,am,gam     
      double precision,allocatable::dx(:),c(:)
      allocate(dx(dim),c(maxn))
      do i = 1, ntotal
        do d = 1, dim
          dvxdt(d, i) = 0.
	enddo
      enddo
        
c     Consider self-gravity or not ?

      if (self_gravity) then
        do i = 1, ntotal
          dvxdt(dim, i) = -9.8
        enddo
      endif 

c     Boundary particle force and penalty anti-penetration force. 
	if(apf.eq.1) then
      rr0 = 1.25e-5
      dd = 1.e-2
      p1 = 12
      p2 = 4
      
      do  k=1,niac
        i = pair_i(k)
        j = pair_j(k)  
        if(itype(i).gt.0.and.itype(j).lt.0) then  
          rr = 0.      
          do d=1,dim
            dx(d) =  x(d,i) -  x(d,j)
            rr = rr + dx(d)*dx(d)
          enddo  
          rr = sqrt(rr)
          if(rr.lt.rr0) then
            f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
            do d = 1, dim
              dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
            enddo
          endif
        endif        
      enddo 

	end if
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!antipenetration force Becker, Teshner
	if(apf.eq.2) then
	do  k=1,niac
        i = pair_i(k)
        j = pair_j(k)  
        if(itype(i).gt.0.and.itype(j).lt.0) then
	    do d=1,dim
            dx(d) =  x(d,i) -  x(d,j)
            rr = rr + dx(d)*dx(d)
		  rr=sqrt(rr)		  
          enddo
		q=rr/hsml(i)
		if(q.gt.0d0.and.q.lt.2d0/3d0) then 
		  am=2d0/3d0 
	    else if(q.gt.2d0/3d0.and.q.lt.1d0) then
	      am=2d0*q-3d0/2d0*q*q
		else if(q.gt.1d0.and.q.lt.2d0) then
	      am=0.5d0*(2d0-q)*(2d0-q)
		else
	      am=0d0
	    end if
		gam=0.02d0*c(i)*c(i)*am/rr
		f=mass(j)/(mass(i)+mass(j))*gam/rr
		do d = 1, dim
            dvxdt(d, i) = dvxdt(d, i)+(x(d,i) -  x(d,j))*f
          enddo
	  end if
	end do  
      
	end if	
      deallocate(dx,c)
      end         
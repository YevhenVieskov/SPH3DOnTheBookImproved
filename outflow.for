	subroutine outflow(Xleft,Xright,Ytop,dx,dy,vel,ntotal,nbound, 
     &							nout,x,vx, mass, rho, p, u,itype, hsml)
c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2 d shear driven cavity probem with Re = 1
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
      
      integer itype(maxn), ntotal,nbound,nout
      double precision x(:, :), vx(:, :), mass(:),
     &     rho(:), p(:), u(:), hsml(:),vel
      double precision Xleft,Xright,Ytop,Ybottom,dx,dy
      integer i, j, d, m, n, mp, np, k,scale_k
      double precision xl, yl,c(maxn) !!!!
	double precision outlength,treshold 
      double precision,allocatable:: xtemp(:, :)
      allocate(xtemp(dim, maxn))
      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
        scale_k = 3 
      endif
	treshold=scale_k*maxval(hsml) 
	Ybottom=Ytop-treshold
	n=ntotal+nbound
	mp=nint((Xright-Xleft)/dx)
	np=nint((Ytop-Ybottom)/dy)      
      nout = mp * np
c      N=nint((XPright-XPleft)/dx2)
c	L=nint((YPtop-YPbottom)/dy2)
c	plate bottom wall
c	N=nint((XPright-XPleft)/dx2)
c	inner layer
c	do i=1,N-1
c	  nghost=nghost+1
c	  x(1,ntotal+nghost)=XPleft+dx+(i-1)*dx2
c	  x(2,ntotal+nghost)=YPbottom+dy2	  
c	end do
c	outer layer
c	do i=1,N-1
c	  nghost=nghost+1
c	  x(1,ntotal+nghost)=XPleft+dx+(i-1)*dx2
c	  x(2,ntotal+nghost)=YPbottom

c	XMright=1.0
      do i = 1, mp
	  do j = 1, np
	    k = j + (i-1)*np	  	
c	  xtemp(1, k) =(Xleft)+(i-1)*dx+dx/2.
c	  xtemp(2, k) =(Ytop)-(j-1)*dy+dy/2.
          
	      x(1, n+k) =(Xleft)+(i-1)*dx !+dx/2.
	      x(2, n+k) =(Ytop)-(j-1)*dy-dy  !/2.
		 
        enddo
      enddo

c	do i=1,nout
c	 if(xtemp(1, i)<(1.0-dx/2.)) then
c	  x(1, n+i)=xtemp(1, i)
c	  x(2, n+i)=xtemp(2, i)
c	 end if
c	end do

	do i=n+1,n+nout
	 vx(1,i)=0.0
       vx(2,i)=vel
	 p(i)=0.0
	 u(i)=0.0
       rho(i)=1000.
	 mass(i)=rho(i)*dx*dy
	 itype(i)=-6
      end do
      deallocate(xtemp)
	end


	subroutine sort_outflow(Xleft,Xright,Ytop,dx,dy,vel,ntotal, 
     &							x,vx, mass, rho, p, u,itype, hsml)
c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2 d shear driven cavity probem with Re = 1
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
      
      integer itype(:), ntotal,nbound,nout
      double precision x(:, :), vx(:, :), mass(:),
     &     rho(:), p(:), u(:), hsml(:),vel
      double precision Xleft,Xright,Ytop,Ybottom,dx,dy
      integer i, j, d, m, n, mp, np, k,scale_k
      double precision xl, yl,c(maxn)!!!!!!
	double precision outlength,treshold
      double precision,allocatable::xtemp(:, :), vxtemp(:, :), 
     &    masstemp(:), rhotemp(:), ptemp(:), utemp(:), hsmltemp(:)
      integer,allocatable::mark(:),itypetemp(:)
	integer ntotalm1
	integer,parameter::ifluid=1,ibound=2,iout=3
      double precision Ymax
      allocate(xtemp(dim, maxn), vxtemp(dim, maxn), masstemp(maxn),
     &     rhotemp(maxn), ptemp(maxn), utemp(maxn), hsmltemp(maxn))
      allocate(mark(maxn),itypetemp(maxn))
      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
        scale_k = 3 
      endif
	treshold=scale_k*maxval(hsml)
	ntotalm1=ntotal 
	n=ntotal+nbound+nout
      ntotal=0
	nbound=0
	nout=0
	do i=1,maxn
	  xtemp(1,i)=0.0
	  xtemp(2,i)=0.0
	end do
	
	if(Ytop<=0.0) then
	  Ymax=Ytop-treshold/2.0
	else
        Ymax=Ytop+treshold/2.0
	end if
	
	do i=1,ntotalm1
	 if(((x(2,i)>=Ymax).and.(itype(i).eq.2))) then
	   ntotal=ntotal+1
	   xtemp(1,ntotal)=x(1,i)
         xtemp(2,ntotal)=x(2,i)
         vxtemp(1,ntotal)=vx(1,i)
         vxtemp(2,ntotal)=vx(2,i)
	   ptemp(ntotal)=p(i)
	   utemp(ntotal)=u(i)
         rhotemp(ntotal)=rho(i)
	   masstemp(ntotal)= mass(i)
	   hsmltemp(ntotal)=  hsml(i)
	   itypetemp(ntotal)=itype(i)
	 end if       
	end do

	do i=1,maxn
	  x(1,i)=0.0
         x(2,i)=0.0
         vx(1,i)=0.0
        vx(2,i)=0.0
	   p(i)=0.0
	   u(i)=0.0
         rho(i)=0.0
	    mass(i)=0.0
	    hsml(i)=0.0
	    itype(i)=0.0
      end do

      do i=1,maxn
	  x(1,i)=xtemp(1,i)
        x(2,i)=xtemp(2,i)
        vx(1,i)= vxtemp(1,i)
        vx(2,i)=vxtemp(2,i)
	  p(i)=ptemp(i)
	  u(i)=utemp(i)
        rho(i)= rhotemp(i)
	  mass(i)=masstemp(i)
	  hsml(i)=hsmltemp(i)
	  itype(i)=itypetemp(i)
      end do
      deallocate(xtemp, vxtemp, masstemp,
     &     rhotemp, ptemp, utemp, hsmltemp)
       deallocate(mark,itypetemp)
	end
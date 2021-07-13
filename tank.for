	subroutine tank(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

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
      use ex
      implicit none     
      include 'param.inc'
      
      integer itype(:), ntotal
      double precision x(:, :), vx(:, :), mass(:),
     &     rho(:), p(:), u(:), hsml(:) !,c(:)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy
      double precision XMleft,XMright,YMbottom,YMtop,XBleft,XBright
      double precision YBbottom,YBtop,XPleft,XPright,YPbottom,YPtop 
	double precision  dx2,dy2
	logical cplate,cmleft,cmright,cnull
      double precision,allocatable:: xtemp(:, :)
      allocate(xtemp(dim, maxn))
c     Giving mass and smoothing length as well as other data.
	
	XMleft=-1.0
	XMright=1.0
	YMbottom=-5.0
	YMtop=0.0

	XBleft=-5.0
	XBright=5.0
	YBbottom=0.0
	YBtop=5.0

	XPleft=-2.0
	XPright=2.0
	YPbottom=1.0
	YPtop=1.5
	 
	dx=0.25
	dy=0.25
	dx2=dx/2.0
	dy2=dy/2.0

	x(1:2,1:maxn)=0.0

	mp=nint((XBright-XBleft)/dx)
	np=nint((YBtop-YMbottom)/dy)
	
	
	
      
      
      
      do i = 1, mp
	do j = 1, np
	  k = j + (i-1)*np	  	
	  xtemp(1, k) =(XBleft)+(i-1)*dx+dx/2.
	  xtemp(2, k) =(YMbottom)+(j-1)*dy   !+dy/2.
        enddo
      enddo


      ntotal=0
     
	do k = 1, mp*np

	 cplate=(xtemp(1, k)>XPleft.and.xtemp(1, k)<XPright)
     & .and.(xtemp(2, k)>YPbottom.and.xtemp(2, k)<YPtop)

	 cmleft=(xtemp(1, k)>=XBleft.and.xtemp(1, k)<=XMleft)
     & .and.(xtemp(2, k)>=YMbottom.and.xtemp(2, k)<=YBbottom)

	 cmright=(xtemp(1, k)>=XMright.and.xtemp(1, k)<=XBright)
     & .and.(xtemp(2, k)>=YMbottom.and.xtemp(2, k)<=YBbottom)
	
c	cnull=(xtemp(1, k)>=-dx2).and.(xtemp(1, k)<=dx2).and.
c     &(xtemp(2, k)>=-dx2).and.(xtemp(2, k)<=dx2)

c	cnull=(xtemp(1, k).eq.0.0).and.(xtemp(1, k).eq.0.0).and.
c     &(xtemp(2, k).eq.0.0).and.(xtemp(2, k).eq.0.0)

c	cnull=(abs(xtemp(1, k)-0.0)<=0.01).and.
c     &(abs(xtemp(1, k)-0.0)<=0.01).and.
c     &(abs(xtemp(2, k)-0.0)<=0.01).and.(abs(xtemp(2, k)-0.0)<=0.01)

	 
       if((.not.(cplate.or.cmleft.or.cmright))) then
	   ntotal=ntotal+1
	   x(1, ntotal) = xtemp(1, k)
	   x(2, ntotal) = xtemp(2, k)
       end if
      enddo

c	do i=1,maxn
c	  xtemp(1, i)=0.0
c	  xtemp(2, i)=0.0
c	end do
      
c	k=0
c	do i=1,ntotal
c	if((x(1, i).eq.0.0.and.x(2, i).eq.0.0)) then
c	  k=k+1
c	  xtemp(1, k)=x(1, k)
c	  xtemp(2, k)=x(2, k)
c	end if 	
c	end do
      

c	do i=1,ntotal
c        x(1, i) = 0.0
c	  x(2, i) = 0.0
c      end do

c	do i=1,ntotal
c        x(1, i) = xtemp(1, i)
c	  x(2, i) = xtemp(2, i)
c      end do
	
	
      do i = 1, ntotal
	   vx(1, i) = 0.
	   vx(2, i) = 0.      
        rho (i) = 1000.   
        mass(i) = dx*dy*rho(i)  
        p(i)= 9.81d0*rho(i)*(Ybtop-x(2,i))   !+1.0d5   
        u(i)=357.1
        itype(i) = 2
        hsml(i) = dx
      enddo
	 
c	initial particle distribution for tensile correction	
	 dpt=dx
c	call vtk(x, vx, mass, rho, p, u, c, itype, hsml, ntotal,
c     &		1111)
      deallocate(xtemp) 
      end	
	
	
	
	
	
	
	
	
	
	
	
	subroutine dbp_part_tank(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	 use ex
       use interf
	implicit none
      include 'param.inc'
      integer itimestep, ntotal, nghost, itype(:)
      double precision hsml(:),mass(:),x(:,:),vx(:,:),
     &                 rho(:), u(:), p(:)
      integer i, j, d, im, mp, scale_k,mp2,M,N,L
      double precision xleft,xright,ytop,ybottom
	double precision dr,dx,dy,xl,v_inf,dr1,dr2,xl2  !, dxiac(maxn) 
	double precision XMleft,XMright,YMbottom,YMtop,XBleft,XBright
      double precision YBbottom,YBtop,XPleft,XPright,YPbottom,YPtop 
	double precision  dx2,dy2
	nghost = 0
      XMleft=-1.0
	XMright=1.0
	YMbottom=-5.0
	YMtop=0.0

	XBleft=-5.0
	XBright=5.0
	YBbottom=0.0
	YBtop=5.0

	XPleft=-2.0
	XPright=2.0
	YPbottom=1.0
	YPtop=1.5
	 
	dx=0.25
	dy=0.25
	dx2=dx/2.0
	dy2=dy/2.0
c	magistral 
      N=nint((XMright-XMleft)/dx2)
	L=nint((YMtop-YMbottom)/dy2)

c	magistral left wall
c	inner layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XMleft
	  x(2,ntotal+nghost)=YMbottom+(i-1)*dy2	  
	end do
c	outer layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XMleft-dx2
	  x(2,ntotal+nghost)=YMbottom+(i-1)*dy2	  
	end do
c	magistral right wall
c	inner layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XMright
	  x(2,ntotal+nghost)=YMbottom+(i-1)*dy2	  
	end do
c	outer layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XMright-dx2
	  x(2,ntotal+nghost)=YMbottom+(i-1)*dy2	  
	end do
c	tank
	N=nint((XBright-XBleft)/dx2)
	L=nint((YBtop-YBbottom)/dy2)
c	tank left wall
c	inner layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XBleft
	  x(2,ntotal+nghost)=YBbottom-dy2+(i-1)*dy2	  
	end do
c	outer layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XBleft-dx2
	  x(2,ntotal+nghost)=YBbottom-dy2+(i-1)*dy2	  
	end do
c	tank right wall
c	inner layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XBright-dx2
	  x(2,ntotal+nghost)=YBbottom-dy2+(i-1)*dy2	  
	end do
c	outer layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XBright
	  x(2,ntotal+nghost)=YBbottom-dy2+(i-1)*dy2	  
	end do
c	bottom wall
c	left side
	N=nint((XMleft-XBleft)/dx2)
c	inner layer
	do i=1,N-2
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XBleft+dx2+(i-1)*dx2
	  x(2,ntotal+nghost)=YBbottom	  
	end do
c	outer layer
	do i=1,N-2
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XBleft+dx2+(i-1)*dx2
	  x(2,ntotal+nghost)=YBbottom-dy2	 
	end do

c	right side
	N=nint((XBright-XMright)/dx2)
c	inner layer
	do i=1,N-2
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XMright+dx2+(i-1)*dx2
	  x(2,ntotal+nghost)=YBbottom	  
	end do
c	outer layer
	do i=1,N-2
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XMright+dx2+(i-1)*dx2
	  x(2,ntotal+nghost)=YBbottom-dy2	 
	end do

c	plate
	N=nint((XPright-XPleft)/dx2)
	L=nint((YPtop-YPbottom)/dy2)

c	plate left wall
c	inner layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPleft+dx2
	  x(2,ntotal+nghost)=YPbottom+(i-1)*dy2	  
	end do
c	outer layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPleft
	  x(2,ntotal+nghost)=YPbottom+(i-1)*dy2	  
	end do
c	plate right wall
c	inner layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPright-dx2
	  x(2,ntotal+nghost)=YPbottom+(i-1)*dy2	 
	end do
c	outer layer
      do i=1,L+1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPright
	  x(2,ntotal+nghost)=YPbottom+(i-1)*dy2	  
	end do
c	plate bottom wall
	N=nint((XPright-XPleft)/dx2)
c	inner layer
	do i=1,N-1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPleft+dx+(i-1)*dx2
	  x(2,ntotal+nghost)=YPbottom+dy2	  
	end do
c	outer layer
	do i=1,N-1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPleft+dx+(i-1)*dx2
	  x(2,ntotal+nghost)=YPbottom	  
	end do
c	plate top wall
c	inner layer
	do i=1,N-1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPleft+dx+(i-1)*dx2
	  x(2,ntotal+nghost)=YPtop-dy2	  
	end do
c	outer layer
	do i=1,N-1
	  nghost=nghost+1
	  x(1,ntotal+nghost)=XPleft+dx+(i-1)*dx2
	  x(2,ntotal+nghost)=YPtop	  
	end do

	do i = 1, nghost
	  rho (ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
	  p(ntotal + i) = 0.
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -3
	  hsml(ntotal + i) = dx
	  vx(1,ntotal + i)=0.0
	  vx(2,ntotal + i)=0.0
      enddo
      rho0=1000.0  
     

      if (mod(itimestep,save_step).eq.0) then
        open(1,file="../data/xv_gp.dat")
        open(2,file="../data/state_gp.dat")
        open(3,file="../data/other_gp.dat")            
        write(1,*) nghost
        do i = ntotal + 1, ntotal + nghost         
          write(1,1001) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)              
          write(2,1002) i, mass(i), rho(i), p(i), u(i)
          write(3,1003) i, itype(i), hsml(i)                               
        enddo       
1001    format(1x, I6, 6(2x, e14.8))
1002    format(1x, I6, 7(2x, e14.8)) 
1003    format(1x, I6, 2x, I4, 2x, e14.8)
        close(1)
        close(2) 
        close(3) 
      endif 

      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
         print *,' >> Statistics: Ghost boundary particles:'
         print *,'          Number of ghost particles:',nghost
        endif     
      endif


	end
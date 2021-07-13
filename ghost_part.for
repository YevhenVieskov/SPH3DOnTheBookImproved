	     
      subroutine g_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	use interf
      implicit none
      include 'param.inc'
      integer itimestep, ntotal, nghost, itype(:)
      double precision hsml(:),mass(:),x(:,:),vx(:,:),
     &                 rho(:), u(:), p(:)
      integer i, j, d, im, mp, scale_k
      double precision xleft,xright,ytop,ybottom
	double precision dr,dx,xl,v_inf,dr1,dr2,dr3,dr4,skh      
	logical cond
      double precision,allocatable::dxiac(:)
      allocate(dxiac(maxn))
	xleft=0d0
	xright=1d0
	ytop=1d0
	ybottom=0d0

	

	mp =2*40
	xl = 0.1
	dx = xl / mp
	v_inf = 2.



	if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
        scale_k = 3 
      endif 

	do i=1,ntotal
	  if(itype(i)>0) then
	  skh=8.0*hsml(i)
c	left wall ghost particles
c		horisontal ghost part
          
	    dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.le.skh.and.dr2.gt.skh.and.dr3.gt.skh.and.dr4.gt.skh
		if(cond) then
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xleft-x(1,i) 
            x(2, ntotal + nghost) = 2.0*x(2,i)-x(2,i)  
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

	    end if

c	right wall ghost particles
c		horisontal ghost part
	    dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.gt.skh.and.dr2.le.skh.and.dr3.gt.skh.and.dr4.gt.skh
		if(cond) then
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xright-x(1,i) 
            x(2, ntotal + nghost) = 2.0*x(2,i)-x(2,i)  
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx


	    end if

c	bottom wall ghost particles
c		vertical ghost part
	    dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.gt.skh.and.dr2.gt.skh.and.dr3.le.skh.and.dr4.gt.skh
		if(cond) then
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*x(1,i)-x(1,i) 
            x(2, ntotal + nghost) = 2.0*ybottom-x(2,i)  
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx



	    end if

c	top wall ghost particles
c		vertical ghost part
	    dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.gt.skh.and.dr2.gt.skh.and.dr3.gt.skh.and.dr4.le.skh
		if(cond) then
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*x(1,i)-x(1,i) 
            x(2, ntotal + nghost) = 2.0*ytop-x(2,i)  
            vx(1, ntotal + nghost) = 2.0*v_inf-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

	    end if
c       corner particles
c       left lower corner
	    dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.le.skh.and.dr2.gt.skh.and.dr3.le.skh.and.dr4.gt.skh
		if(cond) then
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xleft-x(1,i) 
            x(2, ntotal + nghost) =2.0*ybottom-x(2,i)  
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx
cccccccccccccccccccccccccccccccccccccccccccccc
            nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xleft-x(1,i) 
            x(2, ntotal + nghost) = 2.0*x(2,i)-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*x(1,i)-x(1,i) 
            x(2, ntotal + nghost) = 2.0*ybottom-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

	    end if
c       right lower corner
	    dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.gt.skh.and.dr2.le.skh.and.dr3.le.skh.and.dr4.gt.skh
		if(cond) then
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xright-x(1,i)  
            x(2, ntotal + nghost) = 2.0*ybottom-x(2,i) 
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx
cccccccccccccccccccccccccccccccccccccccccccccc
            nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xright-x(1,i) 
            x(2, ntotal + nghost) = 2.0*x(2,i)-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*x(1,i)-x(1,i) 
            x(2, ntotal + nghost) = 2.0*ybottom-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*0.0-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

	   end if
c       left top corner
	    dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.le.skh.and.dr2.gt.skh.and.dr3.gt.skh.and.dr4.le.skh
		if(cond) then
		nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xleft-x(1,i)  
            x(2, ntotal + nghost) = 2.0*ytop-x(2,i)  
            vx(1, ntotal + nghost) = 2.0*v_inf-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx
cccccccccccccccccccccccccccccccccccccccccccccc
            nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xleft-x(1,i) 
            x(2, ntotal + nghost) = 2.0*x(2,i)-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*v_inf-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*x(1,i)-x(1,i) 
            x(2, ntotal + nghost) = 2.0*ytop-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*v_inf-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

	end if
c       right top corner
	  dr1=abs(x(1,i)-xleft)
	    dr2=abs(x(1,i)-xright)
	    dr3=abs(x(2,i)-ybottom)
	    dr4=abs(x(2,i)-ytop)
	    cond=dr1.gt.skh.and.dr2.le.skh.and.dr3.gt.skh.and.dr4.le.skh
		if(cond) then
		nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xright-x(1,i)
            x(2, ntotal + nghost) = 2.0*ytop-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*v_inf-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx
cccccccccccccccccccccccccccccccccccccccccccccc
            nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*xright-x(1,i) 
            x(2, ntotal + nghost) = 2.0*x(2,i)-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*v_inf-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
	      nghost=nghost+1
		  x(1, ntotal + nghost) = 2.0*x(1,i)-x(1,i) 
            x(2, ntotal + nghost) = 2.0*ytop-x(2,i)   
            vx(1, ntotal + nghost) = 2.0*v_inf-vx(1,i) 
	      vx(2, ntotal + nghost) = 2.0*0.0-vx(2,i)
		  rho (ntotal + nghost) = 1000.
	      mass(ntotal + nghost) = rho (ntotal + nghost) * dx * dx
	      p(ntotal + nghost) = p(i)
	      u(ntotal + nghost) = u(i)
	      itype(ntotal + nghost) = -4
	      hsml(ntotal + nghost) = dx

	end if
	  end if
	end do
	
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
        deallocate(dxiac)
	end
      
      


	subroutine dbp_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	 use ex
	implicit none
      include 'param.inc'
      integer itimestep, ntotal, nghost, itype(:)
      double precision hsml(:),mass(:),x(:,:),vx(:,:),
     &                 rho(:), u(:), p(:)
      integer i, j, d, im, mp, scale_k,mp2
      double precision xleft,xright,ytop,ybottom
	double precision dr,dx,xl,v_inf,dr1,dr2,xl2
      double precision,allocatable:: dxiac(:)
      allocate(dxiac(maxn))
	nghost = 0
        mp = 16
	xl = 0.6
	dx = xl / mp
	v_inf = 1.
	xl2=xl+dx
	mp2=xl2/dx

c     Monaghan type virtual particle on the Upper side

        do i = 1, 2*mp+1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = (i-1)*dx/2.0 
          x(2, ntotal + nghost) = xl  
          vx(1, ntotal + nghost) =0.0
	  vx(2, ntotal + nghost) = 0.
        enddo

	  do i = 1, 2*mp+1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = (i-1)*dx/2.0 
          x(2, ntotal + nghost) = xl+dx/2  
          vx(1, ntotal + nghost) =0.0
	  vx(2, ntotal + nghost) = 0.
        enddo

c     Monaghan type virtual particle on the Lower side

        do i = 1, 2*mp+1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = (i-1)*dx/2.0 
          x(2, ntotal + nghost) = 0.  
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo

	  do i = 1, 2*mp+1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = (i-1)*dx/2.0 
          x(2, ntotal + nghost) =-dx/2.0  
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo

c     Monaghan type virtual particle on the Left side

        do i = 1, 2*mp-1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = 0. 
          x(2, ntotal + nghost) = i*dx/2.0
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo
        
	 do i = 1, 2*mp-1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) =-dx/2.0 
         x(2, ntotal + nghost) = i*dx/2.0
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo
	 
	
c     Monaghan type virtual particle on the Right side

        do i = 1, 2*mp-1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = xl 
          x(2, ntotal + nghost) = i*dx/2.  
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo

	  do i = 1, 2*mp-1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = xl+dx/2.0 
          x(2, ntotal + nghost) = i*dx/2.  
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo



	do i = 1, nghost
	  rho (ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
	  p(ntotal + i) = 0.
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -3
	  hsml(ntotal + i) = dx
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

       deallocate(dxiac)
	end

	subroutine db_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	 use ex
	implicit none
      include 'param.inc'
      integer itimestep, ntotal, nghost, itype(maxn)
      double precision hsml(maxn),mass(maxn),x(dim,maxn),vx(dim,maxn),
     &                 rho(maxn), u(maxn), p(maxn)
      integer i, j, d, im, mp, scale_k,mp2
      double precision xleft,xright,ytop,ybottom
	double precision dxiac(maxn),dr,dx,xl,v_inf,dr1,dr2,xl2
	nghost = 0
        mp = 80
	xl = 0.6
	dx = xl / mp
	v_inf = 1.
	xl2=xl+dx
	mp2=xl2/dx

c     Monaghan type virtual particle on the Upper side

        do i = 1, 2*mp+1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = (i-1)*dx/2.0 
          x(2, ntotal + nghost) = xl  
          vx(1, ntotal + nghost) =0.0
	  vx(2, ntotal + nghost) = 0.
        enddo

	  

c     Monaghan type virtual particle on the Lower side

        do i = 1, 2*mp+1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = (i-1)*dx/2.0 
          x(2, ntotal + nghost) = 0.  
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo

	  

c     Monaghan type virtual particle on the Left side

        do i = 1, 2*mp-1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = 0. 
          x(2, ntotal + nghost) = i*dx/2.0
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo
        
	 
	 
	
c     Monaghan type virtual particle on the Right side

        do i = 1, 2*mp-1
   	  nghost = nghost + 1
	  x(1, ntotal + nghost) = xl 
          x(2, ntotal + nghost) = i*dx/2.  
          vx(1, ntotal + nghost) = 0.
	  vx(2, ntotal + nghost) = 0.
        enddo

	  



	do i = 1, nghost
	  rho (ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
	  p(ntotal + i) = 0.
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -3
	  hsml(ntotal + i) = dx
        enddo
      rho0=1000.0  
     

      if (mod(itimestep,save_step).eq.0) then
        open(1,file="../data/xv_db.dat")
        open(2,file="../data/state_db.dat")
        open(3,file="../data/other_db.dat")            
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
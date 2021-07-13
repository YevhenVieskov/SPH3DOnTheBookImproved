	subroutine gmod_part(itimestep, ntotal,ngm,hsml,mass,x,vx,
     &           rho,u,p,itype) 

c----------------------------------------------------------------------
c   Subroutine to determine the information of virtual particles
c   Here only the Monaghan type virtual particles for the 2D shear
c   cavity driven problem are generated.
c     itimestep : Current time step                                 [in]
c     ntotal : Number of particles                                  [in]
c     ndrp  : Number of virtual particles                         [out]
c     hsml   : Smoothing Length                                 [in|out]
c     mass   : Particle masses                                  [in|out]
c     x      : Coordinates of all particles                     [in|out]
c     vx     : Velocities of all particles                      [in|out]
c     rho    : Density                                          [in|out]
c     u      : internal energy                                  [in|out]
c     itype   : type of particles                               [in|out]
      use interf
      implicit none
      include 'param.inc'
      integer itimestep, ntotal, ndrp, itype(:),ngm
      double precision hsml(:),mass(:),x(:,:),vx(:,:),
     &                 rho(:), u(:), p(:)
      integer i, j, d, im, mp
      double precision xl, dx, v_inf

	double precision xmax,xmin,ymax,ymin

      if (vp_input) then          
                        
        open(1,file="../data/xv_dp.dat")
        open(2,file="../data/state_dp.dat")
        open(3,file="../data/other_dp.dat")            
        read(1,*) ndrp
        do j = 1, ndrp   
          i = ntotal + j      
          read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
          read(2,*)im, mass(i), rho(i), p(i), u(i)        
          read(3,*)im, itype(i), hsml(i)                            
        enddo  
        close(1)
        close(2) 
        close(3) 
      
	else 
       
	ngm = 0
        mp = 40
	xl = 1.0
	dx = xl / mp
	xl=xl+dx/2.0
	v_inf = 1.
	xmin=0.0
	xmax=xl
	ymin=xmin
	ymax=xmax


c     ndrp=0

c     Monaghan type virtual particle on the Upper side

        do i = 1, 2*mp+1
   	  ndrp = ndrp + 1
	  x(1, ntotal + ndrp) = (i-1)*dx/2 
          x(2, ntotal + ndrp) = xl  
          vx(1, ntotal + ndrp) = v_inf
	  vx(2, ntotal + ndrp) = 0.
        enddo



        do i=1,ntotal
       if(itype(i)==-2) then
	   if(x(1,i)>=xmax) then
			ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)+dx
	        x(2,ntotal+ndrp)=x(2,i)
			ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)+2d0*dx
	        x(2,ntotal+ndrp)=x(2,i)
	  
        
	   end if 

	if(x(1,i)<=xmin) then
			ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)-dx
	        x(2,ntotal+ndrp)=x(2,i)
			ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)-2d0*dx
	        x(2,ntotal+ndrp)=x(2,i)
	  
	   end if 

 

	if(x(2,i)<=ymin) then
			ndrp=ndrp+1
	        x(2,ntotal+ndrp)=x(2,i)-dx
	        x(1,ntotal+ndrp)=x(1,i)
			ndrp=ndrp+1
	        x(2,ntotal+ndrp)=x(2,i)-2d0*dx
	        x(1,ntotal+ndrp)=x(1,i)
	  
	   end if
	   if(x(1,i)==0d0.and.x(2,i)==0d0) then
	         ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)-dx
	        x(2,ntotal+ndrp)=x(2,i)-dx
              ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)-2.d0*dx
	        x(2,ntotal+ndrp)=x(2,i)-dx
	        ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)-dx
	        x(2,ntotal+ndrp)=x(2,i)-2.d0*dx
	        ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)-2.d0*dx
	        x(2,ntotal+ndrp)=x(2,i)-2.d0*dx
	   end if 
	   if(x(1,i)==xmax.and.x(2,i)==0d0) then
	         ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)+dx
	        x(2,ntotal+ndrp)=x(2,i)-dx
	        ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)+2.d0*dx
	        x(2,ntotal+ndrp)=x(2,i)-dx
	        ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)+dx
	        x(2,ntotal+ndrp)=x(2,i)-2.d0*dx
	        ndrp=ndrp+1
	        x(1,ntotal+ndrp)=x(1,i)+2.d0*dx
	        x(2,ntotal+ndrp)=x(2,i)-2.d0*dx
	   end if 
	end if
	end do
 

	

	do i = 1, ndrp
	  rho (ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
	  p(ntotal + i) = 0.
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -3
	  hsml(ntotal + i) = dx
        enddo
        
      endif   

      if (mod(itimestep,save_step).eq.0) then
        open(1,file="../data/xv_dp.dat")
        open(2,file="../data/state_dp.dat")
        open(3,file="../data/other_dp.dat")            
        write(1,*) ndrp
        do i = ntotal + 1, ntotal + ndrp         
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
         print *,' >> Statistics: Dalrimple boundary particles:'
         print *,'          Number of Dalrimple particles:',ndrp
        endif     
      endif

      end
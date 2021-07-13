		subroutine sym_euler(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
     
c----------------------------------------------------------------------
c      x-- coordinates of particles                       [input/output]
c      vx-- velocities of particles                       [input/output]
c      mass-- mass of particles                                  [input]
c      rho-- dnesities of particles                       [input/output]
c      p-- pressure  of particles                         [input/output]
c      u-- internal energy of particles                   [input/output]
c      c-- sound velocity of particles                          [output]
c      s-- entropy of particles, not used here                  [output]
c      e-- total energy of particles                            [output]
c      itype-- types of particles                               [input]
c           =1   ideal gas
c           =2   water
c           =3   tnt
c      hsml-- smoothing lengths of particles              [input/output]
c      ntotal-- total particle number                            [input]  
c      maxtimestep-- maximum timesteps                           [input]
c      dt-- timestep                                             [input]
      use interf
	use ex
      implicit none     
      include 'param.inc'
      
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt
      integer i, j, k, itimestep, d, current_ts, nstart        
             
      double precision  time, temp_rho, temp_u
       
		
	double precision dt2
	
      double precision,allocatable::x_min(:, :), v_min(:, :), u_min(:),
     &       rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:)
      double precision,allocatable::xn(:, :), vxn(:, :),  
     &       rhon(:), pn(:), un(:)
      double precision,allocatable::xn12(:, :), vxn12(:, :),  
     &       rhon12(:), pn12(:), un12(:)
      
      allocate(x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),
     &       rho_min(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn),  
     &       drho(maxn),  av(dim, maxn), ds(maxn), 
     &       t(maxn), tdsdt(maxn),
     &       xn(dim, maxn), vxn(dim, maxn),  
     &       rhon(maxn), pn(maxn), un(maxn),
     &       xn12(dim, maxn), vxn12(dim, maxn),  
     &       rhon12(maxn), pn12(maxn), un12(maxn) )
      
      do i = 1, ntotal
        do d = 1, dim
          av(d, i) = 0.
        enddo
      enddo  
     
      x_min=0.;  v_min=0.;  
     &       u_min=0.; rho_min=0.;  dx=0.; dvx=0.;  du=0.;   
     &       drho=0.;   av=0.;  ds=0.;  
     &       t=0.;  tdsdt=0.  
     
      nstart=0;   current_ts=0; time=0.

      
      
	dt2=dt/2.0

      do itimestep = nstart+1, nstart+maxtimestep   
	   
        current_ts=current_ts+1
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif      
      
	
	do i=1,ntotal
	  rhon(i)=rho(i)
	  un(i)=u(i)       
	  pn(i)=p(i)
	  
        do d=1,dim
          vxn(d,i)=vx(d,i)	
		xn(d,i)=x(d,i)
	  end do
	end do


c     first step sym_euler
c	v=r+f(t,r)*(h/2)
c	r=r+v*h

	call single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u, s, 
     &       rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)
      
	 

	 
	do i=1,ntotal
	  
	  
        do d=1,dim
          vx(d,i)=vxn(d,i)+dt2*dvx(d,i)	
		x(d,i)=xn(d,i)+dt*vx(d,i)
	  end do
		u(i)=un(i)+dt*du(i)
	  if(.not.summation_density) then       
	   rho(i)=rhon(i)+dt*drho(i) 
	  end if
	end do


c	second step
c	r=v+f(t,r)*h/2
      
	call single_step(itimestep+1, dt, ntotal, hsml, mass, x, vx, u, s, 
     &       rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)

	do i=1,ntotal	  
	 
        do d=1,dim
          vx(d,i)=vx(d,i)+dt2*dvx(d,i)	
		
	  end do
	end do




	time = time + dt  
        

	if (mod(itimestep,save_step).eq.0) then
          call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
		call vtk(x, vx, mass, rho, p, u, c, itype, hsml, ntotal,
     &		itimestep) 

	endif 

        if (mod(itimestep,print_step).eq.0) then
          write(*,*)
          write(*,101)'x','velocity', 'dvx'    
          write(*,100)x(1,moni_particle), vx(1,moni_particle), 
     &                dvx(1,moni_particle)    
        endif
        
101     format(1x,3(2x,a12))	 
100     format(1x,3(2x,e12.6))
	 
      enddo

      nstart=current_ts
      deallocate(x_min, v_min, u_min,
     &       rho_min, dx, dvx, du,  
     &       drho,  av, ds, 
     &       t, tdsdt,
     &       xn, vxn,  
     &       rhon, pn, un,
     &       xn12, vxn12,  
     &       rhon12, pn12, un12 )
      
      end

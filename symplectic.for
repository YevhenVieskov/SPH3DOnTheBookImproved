      subroutine symplectic(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt )
     
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
   
	use ex
      use interf
      implicit none     
      include 'param.inc'
      
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt,time
      integer i, j, k, itimestep, d, current_ts, nstart        
           
      double precision,allocatable:: x_min(:, :), v_min(:, :), 
     &    u_min(:), rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:),
     &       xn(:, :), vxn(:, :), 
     &       rhon(:), pn(:), un(:), 
     &         xn1(:, :), vxn1(:, :), 
     &       rhon1(:), pn1(:), un1(:)
             
      	double precision dt2, epsilon_rdot
          
          allocate  (x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),
     &       rho_min(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn),  
     &       drho(maxn),  av(dim, maxn), ds(maxn), 
     &       t(maxn), tdsdt(maxn),
     &       xn(dim, maxn), vxn(dim, maxn), 
     &       rhon(maxn), pn(maxn), un(maxn), 
     &        xn1(dim, maxn), vxn1(dim, maxn), 
     &       rhon1(maxn), pn1(maxn), un1(maxn) )  
	   

	
	dt2=0.5*dt   
	         
      do i = 1, ntotal
        do d = 1, dim
          av(d, i) = 0.
        enddo
      enddo  


	
     
      do itimestep = nstart+1, nstart+maxtimestep   
	   
        current_ts=current_ts+1
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif      
      
c     If not first time step, then update thermal energy, density and 
c     velocity half a time step  
	

	
	do i=1,ntotal
	  
        
	  rhon(i)=rho(i)
	  un(i)=u(i)       
	  pn(i)=p(i)

	  
        do d=1,dim          

	    vxn(d,i)=vx(d,i)	
		xn(d,i)=x(d,i)
	  end do
	end do
	
      call single_step(itimestep, dt2, ntotal, hsml,mass,x,vx,u,s, 
     &       rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)

	do i=1,ntotal
	   if(.not.summation_density) then
	    rho(i)=rhon(i)+dt2*drho(i)
	   end if
		u(i)=un(i)+dt2*du(i)
	   

          do d=1,dim
	      vx(d,i)=vxn(d,i)+dt2*dvx(d,i)
            x(d,i)=xn(d,i)+dt2*vxn(d,i)
	    end do
	  end do

c	!Hughes & Grahams 2009 correction 
	
	if(ihg_correct) then
	  if(mod(itimestep,ihg_step).eq.0)then
	    do i=1,m_nb
		  rho(i)=rhon(i)+dt2*drho(i) 
		  if(rho(i).lt.rho0) rho(i)=rho0
		  if ((itype(i)).eq.-30) then
            call p_gas(rho(i), u(i), p(i),c(i))  
	      else if (abs(itype(i)).eq.-3) then	     
	        call p_art_water(rho(i), p(i), c(i))
		  end if
		end do
	  end if
	end if


c	calculating pressure
	

      call single_step(itimestep, dt2, ntotal, hsml,mass,x,vx,u,s, 
     &       rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)

	do i=1,ntotal
	   if(.not.summation_density) then
	    epsilon_rdot=-drho(i)/rho(i)
	    rho(i)=rhon(i)*(2.0-epsilon_rdot*dt)/(2.0+epsilon_rdot*dt)
	   end if
		u(i)=un(i)+dt*du(i)
	   

          do d=1,dim
	      vx(d,i)=vxn(d,i)+dt*dvx(d,i)
            x(d,i)=xn(d,i)+dt2*(vxn(d,i)+vx(d,i))
	    end do
	  end do

      do i=1,ntotal
	if (abs(itype(i)).eq.1) then
          call p_gas(rho(i), u(i), p(i),c(i))  
	else if (abs(itype(i)).eq.2) then	     
	  call p_art_water(rho(i), p(i), c(i))
      endif 
      end do 

c	!Hughes & Grahams 2009 correction 
	
	if(ihg_correct) then
	  if(mod(itimestep,ihg_step).eq.0)then
	    do i=1,m_nb
             epsilon_rdot=-drho(i)/rho(i)
	     rho(i)=rhon(i)*(2.0-epsilon_rdot*dt)/(2.0+epsilon_rdot*dt)
		  if(rho(i).lt.rho0) rho(i)=rho0
		  if ((itype(i)).eq.-30) then
            call p_gas(rho(i), u(i), p(i),c(i))  
	      else if (abs(itype(i)).eq.-3) then	     
	        call p_art_water(rho(i), p(i), c(i))
		  end if
		end do
	  end if
	end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
  
     
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
      deallocate  (x_min, v_min, u_min,
     &       rho_min, dx, dvx, du,  
     &       drho,  av, ds, 
     &       t, tdsdt,
     &       xn, vxn, 
     &       rhon, pn, un, 
     &        xn1, vxn1, 
     &       rhon1, pn1, un1 )  
      end
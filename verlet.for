      subroutine verlet(x,vx, mass, rho, p, u, c, s, e, itype, 
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
   
      use interf
      implicit none     
      include 'param.inc'
      
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt
      integer i, j, k, itimestep, d, current_ts, nstart        
            
      double precision  time, temp_rho, temp_u
         
	  
	   
      double precision,allocatable:: x_min(:, :), v_min(:, :),
     &        u_min(:),rho_min(:), dx(:,:), dvx(:, :), 
     &        du(:),drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:) , 
     &       xn(:, :), vxn(:, :), 
     &       rhon(:), pn(:), un(:), rhon1(:) ,un1(:),pn1(:),
     &       vxn1(:,:),xn1(:,:) 
      allocate (x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),
     &       rho_min(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn),  
     &       drho(maxn),  av(dim, maxn), ds(maxn), 
     &       t(maxn), tdsdt(maxn) , 
     &       xn(dim, maxn), vxn(dim, maxn), 
     &       rhon(maxn), pn(maxn), un(maxn),
     &       rhon1(maxn),un1(maxn),pn1(maxn),
     &       vxn1(dim,maxn),xn1(dim,maxn) )
	         
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
	  rhon1(i)=rhon(i)
	  un1(i)=un(i)       
	  pn1(i)=pn(i)
        
	  rhon(i)=rho(i)
	  un(i)=u(i)       
	  pn(i)=p(i)

	  
        do d=1,dim
          vxn1(d,i)=vxn(d,i)	
		xn1(d,i)=xn(d,i)

	    vxn(d,i)=vx(d,i)	
		xn(d,i)=x(d,i)
	  end do
	end do
	

      if(itimestep.eq.1.or.(.not.mod(itimestep,50)))then
	
	  call single_step(itimestep, dt, ntotal, hsml,mass,x,vx,u,s, 
     &       rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)

	  do i=1,ntotal
	   if(.not.summation_density) then
	    rho(i)=rhon(i)+dt*drho(i)
	   end if
		u(i)=un(i)+dt*du(i)
	   

          do d=1,dim
	      vx(d,i)=vxn(d,i)+dt*dvx(d,i)
            x(d,i)=x(d,i)+dt*vxn(d,i)+0.5*dt*dt*dvx(d,i)
	    end do
	  end do

	else

	call single_step(itimestep, dt, ntotal, hsml, mass, x,vx,u,s, 
     &       rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)
	
	do i=1,ntotal
	   if(.not.summation_density) then
	    rho(i)=rhon1(i)+2.0*dt*drho(i)
	   end if

		u(i)=un1(i)+2.0*dt*du(i)
	  
          do d=1,dim
	      vx(d,i)=vxn1(d,i)+2.0*dt*dvx(d,i)
            x(d,i)=x(d,i)+dt*vxn(d,i)+0.5*dt*dt*dvx(d,i)
	    end do
	  end do

	
	 
	end if 
	 
  
     
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
      deallocate (x_min, v_min, u_min,
     &       rho_min, dx, dvx, du,  
     &       drho,  av, ds, 
     &       t, tdsdt , 
     &       xn, vxn, 
     &       rhon, pn, un,rhon1,un1,pn1,
     &       vxn1,xn1)
      end
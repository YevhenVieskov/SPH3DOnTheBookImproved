	subroutine runge(x,vx, mass, rho, p, u, c, s, e, 
     &	itype,hsml, ntotal,maxtimestep , dt )

		
     
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
   

      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal, maxtimestep,num
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn), 
     &       rho(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), 
     &       hsml(maxn), dt
      integer i, j, k, itimestep, d, current_ts, nstart        
      double precision  xtemp(dim, maxn), vxtemp(dim, maxn),utemp(maxn),
     &       rhotemp(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn),  
     &       drho(maxn),  av(dim, maxn), ds(maxn), 
     &       t(maxn), tdsdt(maxn) 
      double precision xk1(dim,maxn), vxk1(dim, maxn), uk1(maxn), 
     &	rhok1(maxn)  
      double precision xk2(dim,maxn), vxk2(dim, maxn), uk2(maxn), 
     &	rhok2(maxn) 
      double precision xk3(dim,maxn), vxk3(dim, maxn), uk3(maxn), 
     &	rhok3(maxn) 
      double precision xk4(dim,maxn), vxk4(dim, maxn), uk4(maxn), 
     &	rhok4(maxn)       
      double precision  time, temp_rho, temp_u,tmtemp
c      double precision  up(ntotal), rhop(ntotal),vxp(ntotal), xp(ntotal)
 
      double precision tm,tlast,dt2
c      open(15,file="../data/xk1.dat")
	 dt2=dt/2.0
c    *********************************************************************************	 
	do itimestep = nstart+1, nstart+maxtimestep   
	   
        current_ts=current_ts+1
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif      
     
	  
	  call single_step(itimestep, dt, ntotal, hsml, mass, x, vx,  
     &           u, s,rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)
        do i=1,ntotal
	     uk1(i)=du(i)*dt
		 rhok1(i)=drho(i)*dt
		 do d = 1, dim
	        vxk1(d, i) = dvx(d, i)*dt 
              xk1(d, i) = vx(d, i)*dt
		
            enddo 
	  end do 
c	do i=1,ntotal
c        write(15,*) xk1(1, i),xk1(2,i)
c         write(15,*) dvx(1, i) !ne mozet bit raven nulu
c      end do
	  
	  do i=1,ntotal
	    utemp(i)=u(i)+0.5*uk1(i)
          rhotemp(i)=rho(i)+0.5*rhok1(i)
          do d = 1, dim
	        vxtemp(d, i) =vx(d, i)+0.5*vxk1(d, i)  
              xtemp(d, i) = x(d, i)+0.5*xk1(d, i)
              
            enddo 
 	  end do


	call single_step(itimestep, dt2, ntotal, hsml, mass, x, vx,  
     &           u, s,rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)

        
	   do i=1,ntotal
	     uk2(i)=du(i)*dt
		 rhok2(i)=drho(i)*dt
		 do d = 1, dim
	         vxk2(d, i) = dvx(d, i)*dt
              xk2(d, i) = vxtemp(d, i)*dt
               
            enddo 
	    end do 
      do i=1,ntotal
	    utemp(i)=u(i)+0.5*uk2(i)
          rhotemp(i)=rho(i)+0.5*rhok2(i)
          do d = 1, dim
	        vxtemp(d, i) =vx(d, i)+0.5*vxk2(d, i)  
              xtemp(d, i) = x(d, i)+0.5*xk2(d, i)
              
            enddo 
 	  end do

      call single_step(itimestep, dt2, ntotal, hsml,mass,xtemp,vxtemp,  
     &   utemp, s,rhotemp, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)

	   do i=1,ntotal
	     uk3(i)=du(i)*dt
		 rhok3(i)=drho(i)*dt
		 do d = 1, dim
	         vxk3(d, i) = dvx(d, i)*dt
              xk3(d, i) = vxtemp(d, i)*dt
              
            enddo 
	    end do 

	do i=1,ntotal
	    utemp(i)=u(i)+uk3(i)
          rhotemp(i)=rho(i)+rhok3(i)
          do d = 1, dim
	        vxtemp(d, i) =vx(d, i)+vxk3(d, i)  
              xtemp(d, i) = x(d, i)+xk3(d, i)
              
            enddo 
 	  end do

      call single_step(itimestep, dt, ntotal, hsml, mass, xtemp,vxtemp,  
     &   utemp, s,rhotemp, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)

	   do i=1,ntotal
	     uk4(i)=du(i)*dt
		 rhok4(i)=drho(i)*dt
		 do d = 1, dim
	         vxk4(d, i) = dvx(d, i)*dt
              xk4(d, i) = vxtemp(d, i)*dt
               
            enddo 
	    end do 

     
	  do i=1,ntotal
	    u(i)=u(i)+1.0/6.0*(uk1(i)+2.0*(uk2(i)+uk3(i))+uk4(i))
          rho(i)=rho(i)+
     &		1.0/6.0*(rhok1(i)+2.0*(rhok2(i)+rhok3(i))+rhok4(i))
          do d = 1, dim
              vx(d, i) =vx(d, i)+
     &   1.0/6.0*(vxk1(d, i)+2.0*(vxk2(d, i)+vxk3(d, i))+vxk4(d, i))
               x(d, i) = x(d, i)+
     &	1.0/6.0*(xk1(d, i)+2.0*(xk2(d, i)+xk3(d, i))+xk4(d, i))  
            enddo 
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



      end

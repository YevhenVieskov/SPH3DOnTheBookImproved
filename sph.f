      program SPH

c----------------------------------------------------------------------
c     This is a three dimensional SPH code. the followings are the 
c     basic parameters needed in this codeor calculated by this code

c     mass-- mass of particles                                      [in]
c     ntotal-- total particle number ues                            [in]
c     dt--- Time step used in the time integration                  [in]
c     itype-- types of particles                                    [in]
c     x-- coordinates of particles                              [in/out]
c     vx-- velocities of particles                              [in/out]
c     rho-- dnesities of particles                              [in/out]
c     p-- pressure  of particles                                [in/out]
c     u-- internal energy of particles                          [in/out]
c     hsml-- smoothing lengths of particles                     [in/out]
c     c-- sound velocity of particles                              [out]
c     s-- entropy of particles                                     [out]
c     e-- total energy of particles                                [out]
      use interf
      implicit none 
       
      include 'param.inc'

      integer ntotal, maxtimestep, d, m, i, yesorno      
      double precision  dt
      double precision s1, s2
	
      
      integer,allocatable::itype(:)[:]
      double precision,allocatable:: x(:,:)[:],vx(:,:)[:],mass(:)[:],
     & rho(:)[:],p(:)[:],u(:)[:],c(:)[:],s(:)[:],e(:)[:],hsml(:)[:]
      double precision,allocatable::pm1(:)[:], vxm1(:, :)[:], 
     &rhom1(:)[:] 
      allocate(itype(maxn)[*])
      allocate(x(dim,maxn)[*], vx(dim, maxn)[*], mass(maxn)[*],
     &rho(maxn)[*], p(maxn)[*], u(maxn)[*], c(maxn)[*], s(maxn)[*], 
     &e(maxn)[*], hsml(maxn)[*], pm1(maxn) [*], vxm1(dim, maxn)[*],
     &rhom1(maxn)[*])
      sync all
      x=0.; vx=0.; mass=0.; rho=0.
      p=0.;u=0.;c=0.;s=0.;e=0.;hsml=0.
      if(this_image()==1) then
        call time_print
        call time_elapsed(s1)      

        if (shocktube)   dt = 0.005
        if (shearcavity) dt = 5.e-5
	  if(breakdam) dt = 5.e-5
        if(drop) dt = 5.e-5
        if(itank) dt = 5.e-5
	end if
      call input(x, vx, mass, rho, p, u, itype, hsml, ntotal) 
      if(this_image()==1) then
 1    write(*,*)'  ***************************************************' 
      write(*,*)'          Please input the maximal time steps '
      write(*,*)'  ***************************************************'
      read(*,*) maxtimestep
      end if
	if(ileapfrog) then     
      call time_integration(x, vx, mass, rho, p, u, c, s, e, itype, 
     &     hsml, ntotal, maxtimestep, dt )
	end if
	if(iverlet) then  
	call verlet(x, vx, mass, rho, p, u, c, s, e, itype, 
     &     hsml, ntotal, maxtimestep, dt )
	end if
	if(iprcorr) then  
      call prcorr(x, vx, mass, rho, p, u, c, s, e, itype, 
     &     hsml, ntotal, maxtimestep, dt )
	end if
	if(ibeeman) then  
       call beeman(x, vx, mass, rho, p, u, c, s, e, itype, 
     &     hsml, ntotal, maxtimestep, dt )
	end if
	if(isymplectic) then 
	  call symplectic(x, vx, mass, rho, p, u, c, s, e, itype, 
     &     hsml, ntotal, maxtimestep, dt )
	end if
      if(isym_euler) then
	call sym_euler(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
	end if
	if(isym_verlet) then
	call sym_verlet(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
	end if

      if(isym_verlet2) then
	call sym_verlet2(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
	end if
      call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal) 
      if(this_image()==1) then
      write(*,*)'  ***************************************************'
      write(*,*) 'Are you going to run more time steps ? (0=No, 1=yes)'
      write(*,*)'  ***************************************************'
      read (*,*) yesorno     
      if (yesorno.ne.0) go to 1
      call time_print
      call time_elapsed(s2)      
      write (*,*)'        Elapsed CPU time = ', s2-s1
      end if
      deallocate(itype)
      deallocate(x, vx, mass, rho, p, u, c, s, 
     &e, hsml, pm1, vxm1, rhom1)                     
      end

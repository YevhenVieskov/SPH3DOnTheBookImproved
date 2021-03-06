      subroutine single_step(itimestep, dt, ntotal, hsml, mass, x, vx,  
     &           u, s,rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av) 

c----------------------------------------------------------------------
c   Subroutine to determine the right hand side of a differential 
c   equation in a single step for performing time integration 

c   In this routine and its subroutines the SPH algorithms are performed.
c     itimestep: Current timestep number                            [in]
c     dt       : Timestep                                           [in]
c     ntotal   :  Number of particles                               [in]
c     hsml     :  Smoothing Length                                  [in]
c     mass     :  Particle masses                                   [in]
c     x        :  Particle position                                 [in]
c     vx       :  Particle velocity                                 [in]
c     u        :  Particle internal energy                          [in]
c     s        :  Particle entropy (not used here)                  [in]
c     rho      :  Density                                       [in/out]
c     p        :  Pressure                                         [out]
c     t        :  Temperature                                   [in/out]
c     tdsdt    :  Production of viscous entropy t*ds/dt            [out]
c     dx       :  dx = vx = dx/dt                                  [out]
c     dvx      :  dvx = dvx/dt, force per unit mass                [out]
c     du       :  du  = du/dt                                      [out]
c     ds       :  ds  = ds/dt                                      [out]     
c     drho     :  drho =  drho/dt                                  [out]
c     itype    :  Type of particle                                 [in]
c     av       :  Monaghan average velocity                        [out]
      use interf
	use ex
      implicit none
      include 'param.inc'

      integer itimestep, ntotal, itype(:)    
      double precision dt, hsml(:), mass(:), x(:,:),     
     &       vx(:,:), u(:), s(:), rho(:), p(:) !, eta(:) !, c(:)
      double precision tdsdt(:), t(:), dx(:,:), dvx(:,:), 
     &       du(:), ds(:), drho(:), av(:, :) 
      integer i, d, nvirt, niac
      
             
     
      integer nghost,ndbp,nout                        
      integer,allocatable:: pair_i(:), pair_j(:), ns(:)
      
      double precision,allocatable:: w(:), dwdx(:,:),  
     &       indvxdt(:,:),exdvxdt(:,:),ardvxdt(:,:),  
     &       avdudt(:), ahdudt(:),c(:),eta(:)
      
      allocate(pair_i(max_interaction),
     &        pair_j(max_interaction), ns(maxn))
      
      allocate( w(max_interaction), dwdx(dim,max_interaction),  
     &       indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),  
     &       avdudt(maxn), ahdudt(maxn),c(maxn),eta(maxn))
      
      
      do  i=1,ntotal
        avdudt(i) = 0.
        ahdudt(i) = 0.
        do  d=1,dim
          indvxdt(d,i) = 0.
          ardvxdt(d,i) = 0.
          exdvxdt(d,i) = 0.
        enddo
      enddo  
 
c---  Positions of virtual (boundary) particles: 

      nvirt = 0
      if (virtual_part) then 
c        call virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
c     &       rho,u,p,itype)
      endif 
	nghost=0

	if (idbp_part) then
	if(.not.itank) then
	call dbp_part(itimestep, ntotal,ndbp,hsml,mass,x,vx,
     &           rho,u,p,itype)
	end if
      if(itank) then
	  call dbp_part_tank(itimestep, ntotal,ndbp,hsml,mass,x,vx,
     &           rho,u,p,itype)
	end if
		
	nvirt=ndbp
	m_nb=ndbp
	end if
      
      if (ghost_part) then
	call g_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	nvirt=nghost
	m_nb=nghost
      end if

	if (igdb_part) then
	call g_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	call db_part(itimestep, ntotal+nghost,nvirt,hsml,mass,x,vx,
     &           rho,u,p,itype)
	nvirt=nvirt+nghost
	m_nb=nvirt
      end if
	
	

	if(itank) then
	call outflow(-1.0d0,1.0d0,-5.d0,0.25d0,0.25d0,-2.0d0,ntotal,nvirt, 
     &						nout,x,vx, mass, rho, p, u,itype, hsml)
	m_nb=m_nb+nout
	end if

      

c---  Interaction parameters, calculating neighboring particles
c     and optimzing smoothing length

      if (nnps.eq.1) then 
        call direct_find(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      else if (nnps.eq.2) then
        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      else if (nnps.eq.3) then 
        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
     &       pair_j,w,dwdx,ns)

		
      endif         
                        
c---  Density approximation or change rate

	if(ishepard) then
	if(mod(itimestep,istepshepard).eq.0) then
      call shepard(ntotal+nvirt,hsml,mass,niac,pair_i,pair_j,w,
     &       itype,rho) 
      end if  
	end if
      
      if (summation_density) then      
        call sum_density(ntotal+nvirt,hsml,mass,niac,pair_i,pair_j,w,
     &       itype,rho)          
      else             
        call con_density(ntotal+nvirt,mass,niac,pair_i,pair_j,
     &       dwdx,vx, itype,x,rho, drho)         
      endif

c---  Dynamic viscosity:

      if (visc.or.visc_morris.or.visc_combine.or.visc_cleary ) then
	    call viscosity(ntotal+nvirt,itype,x,rho,eta)
      end if 
c---  Internal forces:
 
      call int_force(itimestep,dt,ntotal+nvirt,hsml,mass,vx,niac,rho,
     &     eta, pair_i,pair_j,w,dwdx,u,itype,x,t,c,p,indvxdt,tdsdt,du) 
                  
c---  Artificial viscosity:

      if (visc_artificial) call art_visc(ntotal+nvirt,hsml,
     &      mass,x,vx,niac,rho,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt)

      if(visc_combine) call combine_visc(ntotal+nvirt,hsml,
     &   mass,x,vx,niac,rho,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt,eta)
      
c---  External forces:

      if (ex_force) call ext_force(ntotal+nvirt,mass,x,niac,
     &                   pair_i,pair_j,itype, hsml, exdvxdt)


c     Calculating the neighboring particles and undating HSML
      
      if (sle.ne.0) call h_upgrade(dt,ntotal, mass, vx, rho, niac, 
     &                   pair_i, pair_j, dwdx, hsml)

      if (heat_artificial) call art_heat(ntotal+nvirt,hsml,
     &         mass,x,vx,niac,rho,u, c,pair_i,pair_j,w,dwdx,ahdudt)
     
c     Calculating average velocity of each partile for avoiding penetration

      if (average_velocity) call av_vel(ntotal,mass,niac,pair_i,
     &                           pair_j, w, vx, rho, av) 

c---  Convert velocity, force, and energy to f and dfdt  

      do i=1,ntotal
        do d=1,dim
          dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i)
        enddo
        du(i) = du(i) + avdudt(i) + ahdudt(i)
      enddo

      if (mod(itimestep,print_step).eq.0) then      
          write(*,*)
          write(*,*) '**** Information for particle ****', 
     &        		moni_particle         
          write(*,101)'internal a ','artifical a=',
     &         		'external a ','total a '   
          write(*,100)indvxdt(1,moni_particle),ardvxdt(1,moni_particle),
     &                exdvxdt(1,moni_particle),dvx(1,moni_particle)          
      endif
101   format(1x,4(2x,a12))      
100   format(1x,4(2x,e12.6))      
      
      deallocate(pair_i, pair_j, ns)
      
      deallocate( w, dwdx, indvxdt, exdvxdt, ardvxdt, avdudt,
     &ahdudt)
      end


      subroutine single_step2(itimestep, dt, ntotal, hsml, mass, x, vx,  
     &           u, s,rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av) 

c----------------------------------------------------------------------
c   Subroutine to determine the right hand side of a differential 
c   equation in a single step for performing time integration 

c   In this routine and its subroutines the SPH algorithms are performed.
c     itimestep: Current timestep number                            [in]
c     dt       : Timestep                                           [in]
c     ntotal   :  Number of particles                               [in]
c     hsml     :  Smoothing Length                                  [in]
c     mass     :  Particle masses                                   [in]
c     x        :  Particle position                                 [in]
c     vx       :  Particle velocity                                 [in]
c     u        :  Particle internal energy                          [in]
c     s        :  Particle entropy (not used here)                  [in]
c     rho      :  Density                                       [in/out]
c     p        :  Pressure                                         [out]
c     t        :  Temperature                                   [in/out]
c     tdsdt    :  Production of viscous entropy t*ds/dt            [out]
c     dx       :  dx = vx = dx/dt                                  [out]
c     dvx      :  dvx = dvx/dt, force per unit mass                [out]
c     du       :  du  = du/dt                                      [out]
c     ds       :  ds  = ds/dt                                      [out]     
c     drho     :  drho =  drho/dt                                  [out]
c     itype    :  Type of particle                                 [in]
c     av       :  Monaghan average velocity                        [out]
	use ex
      use interf
      implicit none
      include 'param.inc'

      integer itimestep, ntotal, itype(maxn)    
      double precision dt, hsml(maxn), mass(maxn), x(dim,maxn),
     &       vx(dim,maxn), u(maxn), s(maxn), rho(maxn), p(maxn), 
     &       t(maxn), tdsdt(maxn), dx(dim,maxn), dvx(dim,maxn), 
     &       du(maxn), ds(maxn), drho(maxn), av(dim, maxn)              
      integer i, d, nvirt, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), ns(maxn)
      double precision w(max_interaction), dwdx(dim,max_interaction),  
     &       indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),  
     &       avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn)  
     
      integer nghost,ndbp                         

      do  i=1,ntotal
        avdudt(i) = 0.
        ahdudt(i) = 0.
        do  d=1,dim
          indvxdt(d,i) = 0.
          ardvxdt(d,i) = 0.
          exdvxdt(d,i) = 0.
        enddo
      enddo  
 
c---  Positions of virtual (boundary) particles: 

      nvirt = 0
      if (virtual_part) then 
c        call virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
c     &       rho,u,p,itype)
      endif 
	nghost=0

	if (idbp_part) then
	call dbp_part(itimestep, ntotal,ndbp,hsml,mass,x,vx,
     &           rho,u,p,itype)	
	nvirt=ndbp
	m_nb=ndbp
	end if
      
      if (ghost_part) then
	call g_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	nvirt=nghost
	m_nb=nghost
      end if

c---  Interaction parameters, calculating neighboring particles
c     and optimzing smoothing length

      

c---  Convert velocity, force, and energy to f and dfdt  

      do i=1,ntotal
        do d=1,dim
          dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i)
        enddo
        du(i) = du(i) + avdudt(i) + ahdudt(i)
      enddo

      if (mod(itimestep,print_step).eq.0) then      
          write(*,*)
          write(*,*) '**** Information for particle ****', 
     &        		moni_particle         
          write(*,101)'internal a ','artifical a=',
     &         		'external a ','total a '   
          write(*,100)indvxdt(1,moni_particle),ardvxdt(1,moni_particle),
     &                exdvxdt(1,moni_particle),dvx(1,moni_particle)          
      endif
101   format(1x,4(2x,a12))      
100   format(1x,4(2x,e12.6))      

      end

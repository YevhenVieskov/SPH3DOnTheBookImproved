c---------------------------------------------------------
c     Including file for parameters and constants used 
c     in the entire SPH software packages.
c---------------------------------------------------------

c     dim : Dimension of the problem (1, 2 or 3)
      integer dim
      parameter ( dim = 2)

c     maxn    : Maximum number of particles
c     max_interation : Maximum number of interaction pairs
      integer maxn,max_interaction
      parameter ( maxn    = 15000    ,
     &            max_interaction = 100 * maxn )

c     Parameters for the computational geometry,  
c     x_maxgeom : Upper limit of allowed x-regime 
c     x_mingeom : Lower limit of allowed x-regime 
c     y_maxgeom : Upper limit of allowed y-regime 
c     y_mingeom : Lower limit of allowed y-regime 
c     z_maxgeom : Upper limit of allowed z-regime 
c     z_mingeom : Lower limit of allowed z-regime 
      double precision x_maxgeom,x_mingeom,y_maxgeom,
     &                 y_mingeom,z_maxgeom,z_mingeom
      parameter ( x_maxgeom =  10.e0     ,
     &            x_mingeom = -10.e0     ,
     &            y_maxgeom =  10.e0     ,
     &            y_mingeom = -10.e0     ,
     &            z_maxgeom =  10.e0     ,
     &            z_mingeom = -10.e0     )
    
c     SPH algorithm for particle approximation (pa_sph)
c     pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
c              2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
      integer pa_sph 
      parameter(pa_sph = 2)

c     Nearest neighbor particle searching (nnps) method
c     nnps = 1 : Simplest and direct searching
c            2 : Sorting grid linked list
c            3 : Tree algorithm
      integer nnps 
      parameter(nnps = 3 )

c     Smoothing length evolution (sle) algorithm
c     sle = 0 : Keep unchanged,
c           1 : h = fac * (m/rho)^(1/dim)
c           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
c           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) ) 

      integer sle 
      parameter(sle = 0)

c     Smoothing kernel function 
c     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c         = 2, Gauss kernel   (Gingold and Monaghan 1981) 
c         = 3, Quintic kernel (Morris 1997)
      integer skf 
      parameter(skf = 1)

c     Switches for different senarios

c     summation_density = .TRUE. : Use density summation model in the code, 
c                        .FALSE.: Use continuiity equation
c     average_velocity = .TRUE. : Monaghan treatment on average velocity,
c                       .FALSE.: No average treatment.
c     config_input = .TRUE. : Load initial configuration data,
c                   .FALSE.: Generate initial configuration.
c     virtual_part = .TRUE. : Use vritual particle,
c                   .FALSE.: No use of vritual particle.
c     vp_input = .TRUE. : Load virtual particle information,
c               .FALSE.: Generate virtual particle information.
c     visc = .true. : Consider viscosity,
c           .false.: No viscosity.
c     ex_force =.true. : Consider external force,
c               .false.: No external force.
c     visc_artificial = .true. : Consider artificial viscosity,
c                      .false.: No considering of artificial viscosity.
c     heat_artificial = .true. : Consider artificial heating,
c                      .false.: No considering of artificial heating.
c     self_gravity = .true. : Considering self_gravity,
c                    .false.: No considering of self_gravity
c     nor_density =  .true. : Density normalization by using CSPM,
c                    .false.: No normalization.
      
      logical summation_density, average_velocity, config_input,
     &        virtual_part, vp_input, visc, ex_force, heat_artificial,
     & visc_artificial, self_gravity, nor_density,ghost_part,idbp_part
      parameter ( summation_density  = .false. )
      parameter ( average_velocity  = .false. )
      parameter ( config_input  = .false. )
      parameter ( virtual_part  = .false. )
	  parameter ( ghost_part  = .true. )
	  parameter ( idbp_part  = .false. )
      parameter ( vp_input  = .false.  )
      parameter ( visc  = .false.  )
      parameter ( ex_force  = .true.)
      parameter ( visc_artificial  = .true. )
      parameter ( heat_artificial  = .false. )
      parameter ( self_gravity  = .true. )      
      parameter ( nor_density  = .false. )      

c     Symmetry of the problem
c     nsym = 0 : no symmetry,
c          = 1 : axis symmetry,
c          = 2 : center symmetry.     
      integer    nsym
      parameter ( nsym = 0)

c     Control parameters for output 
c     int_stat = .true. : Print statistics about SPH particle interactions.
c                        including virtual particle information.
c     print_step: Print Timestep (On Screen)
c     save_step : Save Timestep    (To Disk File)
c     moni_particle: The particle number for information monitoring.
      logical int_stat
      parameter ( int_stat = .true. )
      integer print_step, save_step, moni_particle
      parameter ( print_step = 100 ,
     &            save_step = 100,
     &            moni_particle = 1600   )
           
      double precision pi
      parameter ( pi = 3.14159265358979323846 )


c	  Shepard kernel correction - ishepard, step correction -istepshepard
		logical ishepard
		integer istepshepard
        parameter ( ishepard  = .true.)
		parameter ( istepshepard =30)

c      integrator parameter
       logical ileapfrog,iverlet,iprcorr,ibeeman,isymplectic
       parameter ( ileapfrog  = .true.)
	   parameter ( iverlet  = .false.)
	   parameter ( iprcorr  = .false.)
	   parameter ( ibeeman  = .false.)
	   parameter ( isymplectic  = .false.)

c      symplectic integrator for SPH
	   logical isym_euler,isym_verlet,isym_verlet2
	   parameter ( isym_euler  = .false.)
	   parameter ( isym_verlet  = .false.)
	   parameter ( isym_verlet2  = .false.)
c    Hughes & Grahams 2009 correction integrator, only for dynamic boundary particles
		logical ihg_correct
		parameter(ihg_correct=.false.)
		integer ihg_step
		parameter(ihg_step=30)

c      Morris viscous model, for morris model visc=.false., visc_morris=.true.
       
	   logical visc_morris
       parameter ( visc_morris  = .false.)

c   viscous model combined morris and artificial viscosity
       logical visc_combine
       parameter ( visc_combine  = .false.)

c      Cleary viscous model, for Cleary model visc=.false., visc_cleary=.true.
       
	   logical visc_cleary
       parameter ( visc_cleary  = .false.)

c	   Monaghan tensile instability correction
	   logical tinst
	   parameter(tinst=.false.)

c     Simulation cases
c     shocktube = .true. : carry out shock tube simulation
c     shearcavity = .true. : carry out shear cavity simulation
      logical shocktube,shearcavity,breakdam,drop
      parameter ( shocktube  = .false.)
      parameter ( shearcavity  = .false.)
	  parameter ( breakdam  = .true. )
	  parameter (drop=.false.)

c	  antipenetration force select apf=1- original, apf=2-new
	  integer apf
	  parameter(apf=0)
c    particles type : itype(i) 1-liquid, 2-gas, -2- virt part liquid, -3-dbp part, -4-ghost_part,
c    for gas -20,-30,-40 respectively
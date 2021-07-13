      subroutine p_gas(rho, u, p, c)
      
c----------------------------------------------------------------------
c   Gamma law EOS: subroutine to calculate the pressure and sound  
 
c     rho    : Density                                              [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
          
      implicit none
      double precision rho, u, p, c   
      double precision gamma 
          
c      For air (idea gas)

      gamma=1.4
      p = (gamma-1) * rho * u     
      c = sqrt((gamma-1) * u) 
     
      end         
      
      subroutine p_art_water(rho, p, c)
      
c----------------------------------------------------------------------
c   Artificial equation of state for the artificial compressibility 

c     rho    : Density                                              [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
c     Equation of state for artificial compressibility   

      implicit none
      double precision rho, u, p, c,b,c0
      double precision gamma, rho0,h_SWL
	integer coef

c     Artificial EOS, Form 1 (Monaghan, 1994)
      gamma=7.
      rho0=1000.       
c      b = 1.013e5
c	b=c0*c0*rho0/gamma
c	press_EoS = B * ( (rho_EoS/rho0)**i_gamma - 1.)
c     cs_EoS = cs0*((rho_EoS/rho0)**3)
	h_SWL=0.15
	coef=10       
	b=coef*coef*rho0*h_SWL/gamma 
	c0 =10.0      !sqrt(gamma*B/rho0)
c      b=200.0*9.81*h_SWL/(rho0*gamma)
c	c0=sqrt(200.0*9.81*h_SWL)
c	c0=10.0*sqrt(2.0*9.81*h_SWL)
c	b=rho0*c0*c0/gamma
      p = b*((rho/rho0)**gamma-1)      
      c = c0*((rho/rho0)**3)

c     Artificial EOS, Form 2 (Morris, 1997)
c      c = 0.01
c      p = c**2 * rho      
      
      end 
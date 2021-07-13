      double precision function tensile(wk,p_i,p_j,rho_i,rho_j) 

c----------------------------------------------------------------------
c   Subroutine to calculate the internal forces on the right hand side 
c   of the Navier-Stokes equations, i.e. the pressure gradient and the
c   gradient of the viscous stress tensor, used by the time integration. 
c   Moreover the entropy production due to viscous dissipation, tds/dt, 
c   and the change of internal energy per mass, de/dt, are calculated. 
 
c     wk         : Kernel                                               [in]
c     pi,pj      : pressure i,j paryicle                                [in]
c     rhoi,rhoj  : pressure i,j paryicle                                [in]     

      use interf
      use ex
      implicit none
      include 'param.inc'
      
      
      double precision  wk,p_i,p_j,rho_i,rho_j
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision wdp,od_wdp,fab,Ra,Rb,R
      double precision,allocatable:: dwdxp(:),dx(:)
      allocate(dwdxp(dim),dx(dim))
      dwdxp=0.
      dx=1.
      tensile=0.
c     Smoothing kernel function 
c     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c         = 2, Gauss kernel   (Gingold and Monaghan 1981) 
c         = 3, Quintic kernel (Morris 1997)

      if(pa_sph.eq.2.and.skf.le.2)  then 
        if(skf.eq.1) then
            call kernel(1d0,dx,1.3d0,wdp,dwdxp)
          else if(skf.eq.2) then
            call kernel(sqrt(2d0),dx,2d0,wdp,dwdxp)
          else
             wdp=1. 
          end if
          od_wdp=1./wdp
          
        
          fab=wk*od_wdp
          fab=fab*fab
          fab=fab*fab
                    
          if(p_i.gt.0d0) then
            Ra=0.01*p_i/rho_i**2
          else
            Ra=0.2*abs(p_i/rho_i**2)
          end if
          
          if(p_j.gt.0d0) then
            Rb=0.01*p_j/rho_j**2
          else
            Rb=0.2*abs(p_j/rho_j**2)
          end if
          R=Ra+Rb
          tensile=R*fab        
      else 
          tensile=0d0
      end if 
          
      deallocate(dwdxp,dx)
      !return
      end
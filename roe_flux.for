      subroutine roe_flux(i,j,vol,vx,u,rf) 
       
c----------------------------------------------------------------------------------------------------- 
c       Subroutine to calculate Roe's numerical flux between  
c         each pair of interacting particles. It is a vector with 
c                   3 components.       
c----------------------------------------------------------------------------------------------------- 
 
      implicit none 
      include 'param.inc' 
      integer i, j, x, y, isign 
      double precision vol(:), vx(:,:), u(:) 
      double precision voll, volr, vl, vr, el, er, pl, pr 
      double precision cmean, rf(3),af(3), a(3), b(3) 
      double precision pmean, vmean, pdrho, pde, r(3,3)  
      double precision lflux(3), rflux(3), c 
       
c     The following variables are calculated by imposing particle i as 
c     the left state and particle j as the right state 
       
      voll=vol(i) 
      volr=vol(j) 
      vl=vx(1,i) 
      vr=vx(1,j) 
      el=u(i) 
      er=u(j)  
             
c     The pressure at both states is calculated using the Tait 
c     equation evaluated at each state 
 
c     Water      
c     call p_art_water(rhol,pl,c) 
c     call p_art_water(rhor,pr,c) 
 
c     Perfect gas     
      call p_gas(1./voll,el,pl,c) 
      call p_gas(1./volr,er,pr,c) 
       
c     The eigenvalues, eigenvectors, wave strength coefficients and  
c     fluxes at the left and right states of the interface are  
c     calculated 
       
      call reigenvalues(voll,volr,vl,vr,el,er,pl,pr,a) 
      call reigenvectors(voll,volr,vl,vr,el,er,pl,pr,r) 
      call ws_coeff(voll,volr,vl,vr,el,er,pl,pr,b) 
      call lrflux(vl,vr,pl,pr,lflux,rflux) 
       
       
      do x=1,3 
         af(x)=0.0 
         do  y=1,3 
            af(x)=af(x)+(dabs(a(y))*b(y)*r(y,x)) 
         enddo 
   
         rf(x)=0.5*(lflux(x)+rflux(x))-0.5*af(x) 
      enddo       
      end
            
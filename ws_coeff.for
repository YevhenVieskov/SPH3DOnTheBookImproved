      subroutine ws_coeff(voll,volr,vl,vr,el,er,pl,pr,b) 
       
c----------------------------------------------------------------------------------------------------- 
c       Subroutine to calculate the wave strength coefficients       
c----------------------------------------------------------------------------------------------------- 
 
 

 
      implicit none 
      include 'param.inc'   
      double precision pmean, vmean, pdrho, pde, ter, tel 
      double precision voll, volr, vl, vr, x1, x2, x3, x4, b(3) 
      double precision el, er, pl, pr, cmean, rhomean, emean 
             
      call sound_speed(voll,volr,vl,vr,el,er,pl,pr,vmean, 
     &           pmean,emean,pde,pdrho,cmean) 
       
c  Calculation of specific total energy at left and right states  
c  based on internal energy and kinetic energy on both states  
 
      tel=el+0.5*(vl**2) 
      ter=er+0.5*(vr**2) 
       
      x1=-pde/(2*cmean**2) 
      x2=ter-tel 
      x3=(pmean-vmean*cmean-pdrho/pde)*(vr-vl)/cmean 
      x4=(volr-voll)*pdrho/pde 
       
      b(3)=x1*(x2+x3+x4) 
      b(2)=volr-voll-(vr-vl)/cmean-2.*b(3) 
      b(1)=(vr-vl)/cmean+b(3)      
 
      end 

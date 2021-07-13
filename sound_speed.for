      subroutine sound_speed(voll,volr,vl,vr,el,er,pl,pr,vmean, 
     &           pmean,emean,pde,pdrho,cmean) 
       
c----------------------------------------------------------------------------------------------------- 
c       Subroutine to calculate the Lagrangian sound velocity       
c----------------------------------------------------------------------------------------------------- 
 
 
      implicit none 
      include 'param.inc' 
       
      double precision voll, volr, vl, vr, el, er, pl, pr, volmean 
      double precision vmean, emean, cmean, cmeans, pmean, pde, pdrho 
      double precision gamma 
             
       
      call pressure_pd(voll,volr,vl,vr,el,er,pl,pr,vmean, 
     &           pmean,emean,pde,pdrho) 
       
c      cmeans=pmean*pde-pdrho 
      gamma = 1.4 
      volmean = 0.5*(voll + volr) 
      cmeans = gamma * pmean / volmean 
C sound speed Cl in Lagragian coordinate 
      cmean=sqrt(cmeans) 
         
      end 

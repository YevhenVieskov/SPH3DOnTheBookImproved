	subroutine reigenvectors(voll,volr,vl,vr,el,er,pl,pr,r)

C-----------------------------------------------------------------------------------------------

C	Subroutine to calculate the corresponding right
C	eigenvectors of the Roe matrix

C-----------------------------------------------------------------------------------------------

	implicit none
	include 'param.inc'
      
      double precision cmean,voll,volr,vl,vr,el,er,pl,pr
      double precision pmean,rhomean,vmean,emean,pdrho,pde,r(3,3)
      
	call sound_speed(voll,volr,vl,vr,el,er,pl,pr,vmean,
     &                 pmean,emean,pde,pdrho,cmean)

c	r1(3)=(1;cmean;vmean*cmean-pmean)
c     r2(3)=(1;0;rhomean**2pdrho/pde) 
c     r3(3)=(1;-cmean;-vmean*cmean-pmean)

      r(1,1)=1.
	r(1,2)=cmean
	r(1,3)=vmean*cmean-pmean

      r(2,1)=1.
	r(2,2)=0.
	r(2,3)=-pdrho/pde

	r(3,1)=1.
	r(3,2)=-cmean
	r(3,3)=-vmean*cmean-pmean
	end
	
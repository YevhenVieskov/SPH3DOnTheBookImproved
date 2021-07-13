	subroutine reigenvalues(voll,volr,vl,vr,el,er,pl,pr,a)

C-----------------------------------------------------------------------------------------------

C	Subroutine to calculate the eigenvalues of the Roe matrix

C-----------------------------------------------------------------------------------------------

	implicit none
	include 'param.inc'
	
	double precision cmean,voll,volr,vl,vr,el,er,pl,pr
	double precision rhomean,vmean,pmean,emean,pde,pdrho
	double precision a(3)
	call sound_speed(voll,volr,vl,vr,el,er,pl,pr,vmean,
     &                 pmean,emean,pde,pdrho,cmean)
	
	a(1)=-cmean
	a(2)=0.
	a(3)=cmean
	end
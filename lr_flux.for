	subroutine lrflux(vl,vr,pl,pr,lflux,rflux)

c-----------------------------------------------------------------
c	Subroutine to calculate the flux at the left and right
c	states corresponding to the pair of interacting particles
c	It then has to be corrected to find the Roe flux
	implicit none
	include 'param.inc'
	double precision vl, vr, pl, pr    
	double precision lflux(3), rflux(3) 

c     The left and right fluxes are defined

	lflux(1)=-vl
	lflux(2)=pl
	lflux(3)=vl*pl

	rflux(1)=-vr
	rflux(2)=pr
	rflux(3)=vr*pr
		
	end
	subroutine pressure_pd(voll,volr,vl,vr,el,er,pl,pr,vmean,
     &	                        pmean,emean,pde,pdrho)
C-----------------------------------------------------------------------------------------------

C	Subroutine to calculate the mean partial derivatives  of the

C     pressure at the intermediate state
C-----------------------------------------------------------------------------------------------
      use interf
	implicit none
	include 'param.inc'
	double precision rho0,gamma,a0,c0,deltap,deltarho,deltae
	double precision voll,volr,vl,vr,el,er,pl,pr,pde,pdrho
	double precision vmean,emean,pmean,x1,x2,b,volmean

	call mv_interf(voll,volr,vl,vr,el,er,pl,pr,vmean,
     &                pmean,emean)
	rho0=1000.

c	Water
c	gamma=7
      gamma=1.4
	b=1.013e5
c	a0=10.*c0
	deltap=pr-pl
	deltarho=volr-voll
	deltae=er-el
	if(deltae.le.0.) deltae=1e-6
	if(deltap.le.0.) deltap=1e-6
	if(deltarho.le.0.) deltarho=1e-6
c	***********WATER
cc    pdrho=rho0*(a0**2)*((rhomean/rho0)**(gamma-1))

cc    pdrho=b*((rhomean/rho0)**(gamma-1))
c     pdrho=0.01**2
c	pde=(deltap-pdrho*deltarho)/deltae

c**************PERFECT GAS

	volmean=0.5*(voll+volr)
	pdrho=-pmean/volmean
	pde=(gamma-1)/volmean
c     pde=(deltap-pdrho*deltarho)/deltae
      end
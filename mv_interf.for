	subroutine mv_interf(voll,volr,vl,vr,el,er,pl,pr,vmean,
     &	                     pmean,emean)
C-----------------------------------------------------------------------------------------------

C	Subroutine to calculate the mean values of the

C     variables at the interface state
C-----------------------------------------------------------------------------------------------

	implicit none
	include 'param.inc'
	double precision cmean,pmean,rhomean,vmean,emean,c
	double precision voll,volr,sl,sr,vl,vr,el,er,pl,pr
C	Now I should call the equation of state to calculate the
C	mean pressure at the interface state
C	Water
      call p_art_water(rhomean,pmean,c)
C	perfect gas 
C	call p_gas(rhomean,emean,pmean,c)
C	sl=sqrt(rhol)
C	sr=sqrt(rhor)
C	rhomean=sl*sr
C	vmean=(sl*vl+sr*vr)/(sl+sr)
C	emean=(sl*el+sr*er)/(sl+sr)
C	pmean=pmean
	pmean=0.5*(pl+pr)
	vmean=0.5*(vl+vr)
	end
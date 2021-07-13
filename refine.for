	subroutine refine(x, vx, mass, rho, p, u, c,itype,hsml,ilevel, 
     &ntotal,nbound,vcrit)
c	particle refinement by Lopez and Roose
      use interf
	use ex
      implicit none
      include 'param.inc'
c      !inout- x, vx, mass, rho, p, u, c, itype, hsml, ntotal
c	!in- vcrit
      integer  ntotal,nbound, itype(:)    
      double precision  hsml(:), mass(:), x(:,:),
     &       vx(:,:), u(:), rho(:), p(:), c(:)
     &                    
c      double precision w(:), dwdx(:,:)
c      integer i, d, niac, pair_i(:),
c     &        pair_j(:), ns(:)
      integer ilevel(:), i,d
      
	double precision vcrit 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	logical crit_type !crit_type=0-velocity, crit_type=1-velocity
c	temporary array for store initial particle distribution
	integer   ntotalt  
      double precision sp,sr,dref,size,xc,yc,size4
	
      integer,allocatable::itypet(:), iref(:),ilevelt(:),ireft(:)
      double precision,allocatable:: hsmlt(:), masst(:), xt(:,:),
     &  vxt(:,:), ut(:), rhot(:), pt(:), ct(:)
c	ilevel-level of refinement,iref=0- if not refine,iref=1- if  refine
c	sp-separation parameter,sr-smoothed ratio,dref-initial particle distribution
      allocate(itypet(maxn))
      allocate(hsmlt(maxn), masst(maxn), xt(dim,maxn),
     &  vxt(dim,maxn), ut(maxn), rhot(maxn), pt(maxn), ct(maxn))
	sp=0.5
	sr=0.5
	dref=dpt
	crit_type=.false.
	do i=1,ntotal
	  ilevel(i)=0
	  iref(i)=0
	end do
c	storing parameters
      do i=1,ntotal
	  itypet(i)=itype(i)
        masst(i)=mass(i)
	  ut(i)=u(i) 
	  rhot(i)=rho(i) 
	  pt(i)=p(i) 
	  ct(i)=c(i)
	  hsmlt(i)=hsml(i)
        ilevelt(i)=ilevel(i)
c	  ireft(i)=iref(i)
	  do d=1,dim
	    vxt(d,i)=vx(d,i)
	    vxt(d,i)=vx(d,i)
	  end do
	end do	    
	
	ntotalt=ntotal
	

	if(crit_type) then
	else
c       velocity criterion
	  do i=1,ntotal
	    if((vx(1,i)>=vcrit.or.vx(2,i)>=vcrit).and.ilevel(i)<1) then 
	      iref(i)=1
	    end if
	  end do
	end if

      do i=1,ntotal
	  ireft(i)=iref(i)
	end do


c	delete particles with iref(i)=1 and reforming particles arrays
      ntotal=0
	do i=1,ntotalt
	  if(iref(i).eq.1) then
	    itype(i)=0
          mass(i)=0.
		u(i)=0.
		rho(i)=0. 
		p(i)=0.
		c(i)=0.
	    hsml(i)=0.0
		ilevel(i)=0
		iref(i)=0
		do d=1,dim
	      vx(d,i)=0
	      vx(d,i)=0
		end do
	  else 
	    ntotal=ntotal+1
	    itype(ntotal)=itypet(i)
          mass(ntotal)=masst(i)
		u(ntotal)=ut(i) 
		rho(ntotal)=rhot(i) 
		p(ntotal)=pt(i) 
		c(ntotal)=ct(i)
	    hsml(ntotal)=hsmlt(i)
		ilevel(ntotal)=ilevelt(i)
		iref(ntotal)=0
		do d=1,dim
	      vx(d,ntotal)=vxt(d,i)
	      vx(d,ntotal)=vxt(d,i)
		end do  
	  end if
	end do
c	refine particles	
	do i=1,ntotalt
	  if(ireft(i).eq.1) then
	    size=sp*dref
	    size4=size/4.0
          xc=xt(1,i)
	    yc=xt(2,i)
c		lower left corner
	    ntotal=ntotal+1
          x(1,ntotal)=xc-size4
	    x(2,ntotal)=yc-size4
          itype(ntotal)=itypet(i)
          mass(ntotal)=masst(i)/4.0
		u(ntotal)=ut(i) 
		rho(ntotal)=rhot(i) 
		p(ntotal)=pt(i) 
		c(ntotal)=ct(i)
	    hsml(ntotal)=sr*hsmlt(ntotal)
		ilevel(ntotal)=ilevelt(i)+1
		iref(ntotal)=0
c		top left corner
	    ntotal=ntotal+1
          x(1,ntotal)=xc-size4
	    x(2,ntotal)=yc-size4
          itype(ntotal)=itypet(i)
          mass(ntotal)=masst(i)/4.0
		u(ntotal)=ut(i) 
		rho(ntotal)=rhot(i) 
		p(ntotal)=pt(i) 
		c(ntotal)=ct(i)
	    hsml(ntotal)=sr*hsmlt(ntotal)
		ilevel(ntotal)=ilevelt(i)+1
		iref(ntotal)=0
c		top right corner
	    ntotal=ntotal+1
          x(1,ntotal)=xc+size4
	    x(2,ntotal)=yc+size4
          itype(ntotal)=itypet(i)
          mass(ntotal)=masst(i)/4.0
		u(ntotal)=ut(i) 
		rho(ntotal)=rhot(i) 
		p(ntotal)=pt(i) 
		c(ntotal)=ct(i)
	    hsml(ntotal)=sr*hsmlt(ntotal)
		ilevel(ntotal)=ilevelt(i)+1
		iref(ntotal)=0
c		lower right corner
	    ntotal=ntotal+1
          x(1,ntotal)=xc+size4
	    x(2,ntotal)=yc-size4
          itype(ntotal)=itypet(i)
          mass(ntotal)=masst(i)/4.0
		u(ntotal)=ut(i) 
		rho(ntotal)=rhot(i) 
		p(ntotal)=pt(i) 
		c(ntotal)=ct(i)
	    hsml(ntotal)=sr*hsmlt(ntotal)
		ilevel(ntotal)=ilevelt(i)+1
		iref(ntotal)=0
	  end if
	end do
	deallocate(itypet)
      deallocate(hsmlt, masst, xt,
     &  vxt, ut, rhot, pt, ct)	 
	end subroutine refine
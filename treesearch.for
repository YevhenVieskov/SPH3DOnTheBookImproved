      subroutine tree_search(itimestep, ntotal,hsml,x,niac,pair_i,
     &           pair_j,w,dwdx,countiac)
	use treedefs
      use interf
c	use param
	INCLUDE 'param.inc'
      LOGICAL INTERSECT
	INTEGER SPTR,Q,P,NIAC
	REAL(8) PSIZE(MXBODY),DRIAC,MHSML
	REAL(8) DX,DY,DZ,DPC
	INTEGER PAIR_I(:), PAIR_J(:)
	integer i, j, d,  sumiac, maxiac, miniac, noiac,
     &        maxp, minp, scale_k ,countiac(:)
	double precision r
	double precision  x(:,:), w(:),
     &       dwdx(:,:),HSML(:)
		real(4) startt, endt
       double precision, allocatable:: dxiac(:), STACK(:), tdwdx(:)
       allocate(dxiac(dim), STACK(max_interaction), tdwdx(dim))
	call cpu_time(startt)
      	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
        scale_k = 3 
      endif 
      NBODY=ntotal
      do i=1,ntotal
        countiac(i) = 0
      enddo
C      
C      ÂÛÏÎËÍÈÒÜ ÏÅÐÅÂÎÄ ÌÀÑÑÈÂÀ X Â ÌÀÑÑÈÂ    POS   
      POS(1:NBODY,1)=x(1,1:NBODY)
      POS(1:NBODY,2)=x(2,1:NBODY)
	IF(dim==3) THEN
	  POS(1:NBODY,3)=x(3,1:NBODY)
      ELSE
        POS(1:NBODY,3)=0.0
      END IF
      RSIZE=4.0
C     ********************************************8
	CALL EXPBOX(NBODY)
C	--------------------------
C       Load bodies into the tree.
C	--------------------------
	CALL LDTREE(NBODY)
C     ---------------------
C	create interaction list
C	---------------------
	NIAC=0
C	PRINT*,NBODY
	DO P=1, NBODY    
	SPTR=1
	STACK(SPTR)=ROOT
   10 IF(SPTR.GT.0) THEN
	Q=STACK(SPTR)
	SPTR=SPTR-1
	IF(Q .LT. INCELL) THEN

C     ---------------------
C	process body-body interaction
C	---------------------
      IF(Q.NE.P) THEN
	IF(Q.GT.P) THEN
C	
       DX=POS(P,1)-POS(Q,1)
       DY=POS(P,2)-POS(Q,2)
	 DZ=POS(P,3)-POS(Q,3)
	 DRIAC=DX*DX+DY*DY+DZ*DZ
	 MHSML=(HSML(P)+HSML(Q))/2.0
c	PRINT*,HSML(P)
	 IF(SQRT(DRIAC).LT.scale_k*MHSML)THEN
		NIAC=NIAC+1
	    PAIR_I(NIAC)=P
	    PAIR_J(NIAC)=Q
	    countiac(P) = countiac(Q) + 1
          countiac(Q) = countiac(Q) + 1
	    
C     COMPUTE KERNEL
          r=sqrt(DRIAC)
	    countiac(P) = countiac(P) + 1
          countiac(Q) = countiac(Q) + 1
          dxiac(1) = x(1,P) - x(1,Q)
		if(dim>=2) dxiac(2) = x(2,P) - x(2,Q)
		if(dim==3) dxiac(3) = x(3,P) - x(3,Q)
		call kernel(r,dxiac,mhsml,w(niac),tdwdx)
              do d=1,dim
                dwdx(d,niac) = tdwdx(d)
              enddo   

	 END IF
	 END IF
	END IF	
	ELSE
		DX=ABS(POS(P,1)-MID(Q,1))
		DY=ABS(POS(P,2)-MID(Q,2))
		DZ=ABS(POS(P,3)-MID(Q,3))
C	PRINT*,HSML
		PSIZE(P)=scale_k*HSML(P)*2.0  !!!!!!!!!!!!!!!!!!!!!
	    DPC=(PSIZE(P)+CLSIZE(Q))/2.0
		INTERSECT=DX.LE.DPC.AND.DY.LE.DPC.AND.DZ.LE.DPC
		IF(INTERSECT) THEN
C      ÏÐÎÏÓÑÊÀÅÌ ÂÑÒÀÂÊÓ Â ÑÒÅÊ ÅÑËÈ ÍÅ ÏÅÐÅÑÅÊËÎÑÜ
          DO 20 K=1,NSUBC
	    IF(SUBP(Q,K).NE.NULL) THEN
			SPTR=SPTR+1
			STACK(SPTR)=SUBP(Q,K)
		END IF
   20 CONTINUE
	END IF
      ENDIF
	GOTO 10
	ENDIF
	END DO
C     ÑÒÀÒÈÑÒÈÊÀ
      sumiac = 0
      maxiac = 0
      miniac = 1000
      noiac  = 0
      do i=1,ntotal
        sumiac = sumiac + countiac(i)
        if (countiac(i).gt.maxiac) then
	  maxiac = countiac(i)
	  maxp = i
	endif
	if (countiac(i).lt.miniac) then 
	  miniac = countiac(i)
          minp = i
	endif
        if (countiac(i).eq.0)      noiac  = noiac + 1
      enddo
 
c      if (mod(itimestep,print_step).eq.0) then
c        if (int_stat) then
c          print *,' >> Statistics: interactions per particle:'
c          print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
c          print *,'**** Particle:',minp, ' minimal interactions:',miniac
c          print *,'**** Average :',real(sumiac)/real(ntotal)
c          print *,'**** Total pairs : ',niac
c          print *,'**** Particles with no interactions:',noiac
c        endif     
c      endif  
      deallocate(dxiac,STACK,tdwdx)
    	end

C ================================================================
C TREELOAD: routines to handle tree construction and verification.
C ================================================================
 

 
C -----------------------------------------------------
C EXPBOX: enlarge system cube to include all particles.
C -----------------------------------------------------
 
      SUBROUTINE EXPBOX(NBODY)
 
c        INCLUDE 'treedefs.fi'
      use treedefs
	REAL*8 XYZMAX
        INTEGER J, I, NBODY
 
C       ------------------------------
C       Find maximum coordinate value.
C       ------------------------------
	XYZMAX = 0
        DO 10 J = 1, NDIM
	  DO 20 I = 1, NBODY
	    XYZMAX = MAX(XYZMAX, ABS(POS(I,J)))
  20	  CONTINUE
  10	CONTINUE
C	-----------------------------------------
C	Expand box by factors of 2 until it fits.
C	-----------------------------------------
  30	IF (XYZMAX .GE. RSIZE/2.0) THEN
	  RSIZE = 2.0 * RSIZE
	  GOTO 30
	ENDIF
      END

C ================================================================
C TREELOAD: routines to handle tree construction and verification.
C ================================================================
 
C ----------------------------------------------------------------
C MKTREE: initialize the tree structure for the force calculation.
C ----------------------------------------------------------------
 
      
C -----------------------------------------------------
C EXPBOX: enlarge system cube to include all particles.
C -----------------------------------------------------
 
      


C ---------------------------------------------------------------------
C LDTREE: construct tree body-by-body.  This phase initializes the SUBP
C array and loads the geometric midpoint into the MID of each cell.
C ---------------------------------------------------------------------

      SUBROUTINE LDTREE(NBODY)

	use treedefs
        INTEGER MKCELL, K, P

C       ---------------------------------------
C       Deallocate current tree, begin new one.
C       ---------------------------------------
        NCELL = 0
        ROOT = MKCELL()
C	------------------------------------------
C	Initialize midpoint and size of root cell.
C	------------------------------------------
	DO 10 K = 1, NDIM
	  MID(ROOT,K) = 0.0
  10	CONTINUE
	CLSIZE(ROOT) = RSIZE
C       ---------------------------------------------
C       Load bodies into the new tree, one at a time.
C       ---------------------------------------------
        DO 20 P = 1, NBODY
          CALL LDBODY(P)
  20    CONTINUE


c      do i=INCELL,MXNODE
c	  write(1,*) (subp(i,j),j=1,4)
c      end do


c	do i=INCELL,MXNODE
c	  write(2,*) (subp(i,j),j=5,8)
c      end do
      
      END
 
C ----------------------------------------
C LDBODY: load body P into tree structure.
C ----------------------------------------
 
      SUBROUTINE LDBODY(P)
	use treedefs
      INTEGER P
             
      INTEGER Q, QIND, SBINDX, MKCELL, C, K, P0

C       ---------------------------------------------
C       Start Q,QIND pair in correct subcell of root.
C       ---------------------------------------------
        Q = ROOT
        QIND = SBINDX(P, Q)
C       -----------------------------------------------------
C       Loop descending tree until an empty subcell is found.
C       -----------------------------------------------------
  20    IF (SUBP(Q, QIND) .NE. NULL) THEN
C         --------------------------------------
C         On reaching another body, extend tree.
C         --------------------------------------
          IF (SUBP(Q, QIND) .LT. INCELL) THEN
C           -------------------------------------------
C           Allocate an empty cell to hold both bodies.
C           -------------------------------------------
            C = MKCELL()
C           ------------------------------------------------------
C           Locate midpoint of new cell wrt. parent, and set size.
C           ------------------------------------------------------
            DO 30 K = 1, NDIM
              IF (POS(P,K) .GE. MID(Q,K)) THEN
                MID(C,K) = MID(Q,K) + CLSIZE(Q)/4.0
              ELSE
                MID(C,K) = MID(Q,K) - CLSIZE(Q)/4.0
              ENDIF
  30        CONTINUE
	    CLSIZE(C) = CLSIZE(Q) / 2.0
C           ------------------------------------------------------
C           Store old body in appropriate subcell within new cell.
C           ------------------------------------------------------
            P0 = SUBP(Q, QIND)
            SUBP(C, SBINDX(P0, C)) = P0
C           ---------------------------------------------
C           Link new cell into tree in place of old body.
C           ---------------------------------------------
            SUBP(Q, QIND) = C
          ENDIF
C         --------------------------------------------------------
C         At this point, the node indexed by Q,QIND is known to be
C	  a cell, so advance to the next level of tree, and loop.
C         --------------------------------------------------------
          Q = SUBP(Q, QIND)
          QIND = SBINDX(P, Q)
          GOTO 20
        ENDIF
C       ---------------------------------------------
C       Found place in tree for P, so store it there.
C       ---------------------------------------------
        SUBP(Q, QIND) = P
      END
 
C -------------------------------------------------------
C SBINDX: compute subcell index for node P within cell Q.
C -------------------------------------------------------
 
      INTEGER FUNCTION SBINDX(P, Q)
	use treedefs
      INTEGER P, Q
 
        
        INTEGER K
 
C       ---------------------------------------------------
C       Initialize subindex to point to lower left subcell.
C       ---------------------------------------------------
        SBINDX = 1
C       ---------------------------------
C       Loop over all spatial dimensions.
C       ---------------------------------
        DO 10 K = 1, NDIM
          IF (POS(P,K) .GE. MID(Q,K))
     &      SBINDX = SBINDX + 2 ** (3 - K)
 10     CONTINUE
      END
 
C ---------------------------------------------------------
C MKCELL: function to allocate a cell, returning its index.
C ---------------------------------------------------------
 
      INTEGER FUNCTION MKCELL()
 
        use treedefs
        INTEGER I
 
C       ----------------------------------------------------------
C       Terminate simulation if no remaining space for a new cell.
C       ----------------------------------------------------------
        IF (NCELL .GE. MXCELL)
     &    CALL TERROR(' MKCELL: NO MORE MEMORY')
C       ----------------------------------------------------
C       Increment cell counter, initialize new cell pointer.
C       ----------------------------------------------------
        NCELL = NCELL + 1
        MKCELL = NCELL + MXBODY
C       --------------------------------------
C       Zero pointers to subcells of new cell.
C       --------------------------------------
        DO 10 I = 1, NSUBC
          SUBP(MKCELL,I) = NULL
 10     CONTINUE
      END
 



      SUBROUTINE TERROR(MSG)
	use treedefs
      CHARACTER*(*) MSG

	

C       ------------------------------------
C       Write error message to the log file.
C       ------------------------------------
        CALL OUTERR(MSG)
C       ---------------------------------------------------
C       Output timing data, close files, terminate the run.
C       ---------------------------------------------------
C        CALL OUTCPU
C        CALL STPOUT
	STOP
      END

	SUBROUTINE OUTERR(MSG)
	use treedefs
      CHARACTER*(*) MSG

	
 
C       ------------------------------------------------------
C       Write the message, surrounded by stars for visibility.
C       ------------------------------------------------------
C        WRITE (ULOG, '(/,1X,72(''*''))')
C        WRITE (ULOG, '(/,A)') MSG
C        WRITE (ULOG, '(/,1X,72(''*''))')
C       ---------------------------
C       Repeat message to terminal.
C	---------------------------
        WRITE (UTERM, '(/,1X,72(''*''))')
        WRITE (UTERM, '(/,A)') MSG
        WRITE (UTERM, '(/,1X,72(''*''))')
      END



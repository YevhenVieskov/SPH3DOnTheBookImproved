module treedefs
! NULL        value denoting pointer to nothing.
! NDIM        number of spatial dimensions.
! NSUBC       number of descendents per cell.
! NQUAD       number of independent quadrupole components.
	INTEGER NULL, NDIM, NSUBC, NQUAD
        PARAMETER(NULL = 0)
        PARAMETER(NDIM = 3)
        PARAMETER(NSUBC = 2**NDIM)
        PARAMETER(NQUAD = 2*NDIM - 1)
!C ---------------------------
!C Tree structure data arrays.
!C ---------------------------
 
!C Bodies and cells are used to represent the tree.  The information for a
!C given node is distributed across a number of arrays.  Quantities defined
!C for both bodies and cells, such as mass and position, are stored in
!C arrays with indicies running from 1 to MXNODE.  Quantities defined only
!C for bodies are stored in arrays with indicies running from 1 to MXBODY,
!C while those defined only for cells are stored in arrays with indicies
!C running from INCELL to MXNODE.  With this convention, indicies can refer
!C to either bodies or cells without conflict, and the type is deduced by
!C comparing the index to the index of the first cell.  During tree
!C construction, descendent indicies are stored in the SUBP arrays:
!C
!C          +-------------------------------------------------------+
!C ROOT --> | CELL: MASS, POS, RCRIT2, QUAD, SUBP:[/,o,/,/,/,/,o,/] |
!C          +----------------------------------------|---------|----+
!C                                                   |         |
!C     +---------------------------------------------+         |
!C     |                                                       |
!C     |    +--------------------------------+                 |
!C     +--> | BODY: MASS, POS, VEL, ACC, PHI |                 |
!C          +--------------------------------+                 |
!C                                                             |
!C     +-------------------------------------------------------+
!C     |
!C     |    +-------------------------------------------------------+
!C     +--> | CELL: MASS, POS, RCRIT2, QUAD, SUBP:[o,/,/,o,/,/,o,/] |
!C          +--------------------------------------|-----|-----|----+
!C                                                etc   etc   etc
 
!C MXBODY      maximum number of bodies allowed.
!C MXCELL      maximum number of cells allowed.
!C INCELL      index of first cell in arrays.
!C MXNODE      maximum number of nodes (bodies + cells).

!C NBODY       number of bodies in the system.
!C NCELL       number of cells currently in use.
!C MASS        mass of each node.
!C POS         position of each node.
!C VEL         velocity of each body.
!C ACC         acceleration of each body
!C PHI         potential of each body.
!C RCRIT2      critical distances**2 of each cell.
!C QUAD        quad moments of each cell.
!C SUBP        descendent of each cell.
!C TPOS        current position time.
!C TVEL        current velocity time.
!C RSIZE       side-length of root cell.
!C ROOT        index of cell representing root.
 
        INTEGER MXBODY, MXCELL, INCELL, MXNODE
        PARAMETER(MXBODY = 12000)
        PARAMETER(MXCELL = MXBODY)
        PARAMETER(INCELL = MXBODY + 1)
        PARAMETER(MXNODE = MXBODY + MXCELL)
 
        INTEGER  NCELL
        REAL*8  POS(MXNODE,NDIM) !MASS(MXNODE),
!        REAL*8 VEL(MXBODY,NDIM), ACC(MXBODY,NDIM), PHI(MXBODY)
!        REAL*8 RCRIT2(INCELL:MXNODE), QUAD(INCELL:MXNODE,NQUAD)
        INTEGER SUBP(INCELL:MXNODE,NSUBC), ROOT
        REAL*8 TPOS, TVEL, RSIZE
!        COMMON /NODES/ NBODY, NCELL, MASS, POS, VEL, ACC, PHI,&
!                      RCRIT2, QUAD, SUBP, TPOS, TVEL, RSIZE, ROOT
 
!C These definitions apply only during tree construction.

!C MID         geometric center of each cell [equiv to POS].
!C CLSIZE      size of each cell [equiv to RCRIT2].

	REAL*8 MID(MXNODE,NDIM), CLSIZE(INCELL:MXNODE)
!	EQUIVALENCE (MID(1,1), POS(1,1))
!	EQUIVALENCE (CLSIZE(INCELL), RCRIT2(INCELL))


 
!C -----------------------------------------------------------------
!C Declarations for force calculation with shared interaction

!C --------------------------
!C Input/output unit numbers.
!C --------------------------
 
        INTEGER UTERM, UPARS, UBODI, UBODO, UBODF, ULOG
        PARAMETER(UTERM =  6, UPARS = 10, UBODI = 11)
        PARAMETER(UBODO = 12, UBODF = 13, ULOG  = 14)
!ccccccccccccccccccccccccccccccccc

end

 
      subroutine vtk(x, vx, mass, rho, p, u, c, itype, hsml, ntotal,num) 
      
c----------------------------------------------------------------------           
c     Subroutine for saving particle information to external disk file

c     x-- coordinates of particles                                  [in]
c     vx-- velocities of particles                                  [in]
c     mass-- mass of particles                                      [in]
c     rho-- dnesities of particles                                  [in]
c     p-- pressure  of particles                                    [in]
c     u-- internal energy of particles                              [in]
c     c-- sound velocity of particles                               [in]
c     itype-- types of particles                                    [in]
c     hsml-- smoothing lengths of particles                         [in]
c     ntotal-- total particle number                                [in]
	use interf
      use ex
      implicit none     
      include 'param.inc'
      
      integer itype(:), ntotal
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:),p(:), u(:), c(:), hsml(:)
     
      double precision,allocatable::xtemp(:, :), vxtemp(:, :)
      
      integer i, d, npart,num,IOpst     
      
      character(50)	:: szNumber
      character(41)	:: szFileName
	character(41)	:: szFileNameV
	character(41)	:: szFileNameP
	character(41)	:: szFileNameVOF
      character(200) filename1,filename2
      allocate(xtemp(3, maxn), vxtemp(3, maxn))
	IOpst=1

!      generation rezult files v zadannij moment vremeni
      if(num .ne. 0) then
	 Write(szNumber,'(i6.6)') num
	szFileName= Trim("../data/")//Trim("f_xv")//Trim(szNumber)
     &//Trim(".vtk")
c	szFileName= Trim("../data/")//Trim("f_state")//Trim(szNumber)
c     &//Trim(".vtk")
c      szFileName= Trim("../data/")//Trim("f_other")//Trim(szNumber)
c     &//Trim(".vtk")
    
      end if

      open(IOpst,file=szFileName)
      write(IOpst,'(A)')   '# vtk DataFile Version 2.0'
	write(IOpst,'(A,A,A)')   'SPH result output'
	write(IOpst,'(A)')      'ASCII'
      write(IOpst,'(A)')      'DATASET UNSTRUCTURED_GRID'
	write(IOpst,'(A,1X,I10,1X,A)') 'POINTS',ntotal+m_nb,'double'
	xtemp=0d0
      if(dim==1) then
        xtemp(1,1:maxn)=x(1,1:maxn)
	  vxtemp(1,1:maxn)=vx(1,1:maxn)
	else if(dim==2) then
	  xtemp(1,1:maxn)=x(1,1:maxn)
	  xtemp(2,1:maxn)=x(2,1:maxn)
	  vxtemp(1,1:maxn)=vx(1,1:maxn)
	  vxtemp(2,1:maxn)=vx(2,1:maxn)
      else if(dim==3) then
	  xtemp(1:dim,1:maxn)=x(1:dim,1:maxn)
	  vxtemp(1:dim,1:maxn)=vx(1:dim,1:maxn)
	end if
       
	
	 do i=1,ntotal+m_nb
	  do d=1,3
	  
	  
	   write (IOpst,*) xtemp(d,i)   
	  
	  
	 end do
	end do


      write(IOpst,'(A,1X,I10,1X,I10)') 'CELLS',ntotal+m_nb,
     &2*(ntotal+m_nb)
c
      do i=1,ntotal+m_nb
	  write(IOpst,'(I10,1X,I10)') 1,i-1
	end do

      write(IOpst,'(A,1X,I10)')  'CELL_TYPES',ntotal+m_nb

      do i=1,ntotal+m_nb
	  write(IOpst,'(I1)') 1
	end do
	
	
	
	
	 
      write(IOpst,'(A,1X,I10)')  'POINT_DATA',ntotal+m_nb

	write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'pressure',
     &'double',1
	write(IOpst,'(A)') "LOOKUP_TABLE default"
      do i=1,ntotal+m_nb
	  write (IOpst,*) p(i)
	end do

	write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'density',
     &'double',1
	write(IOpst,'(A)') "LOOKUP_TABLE default"
      do i=1,ntotal+m_nb
	  write (IOpst,*) rho(i)
	end do
      
	write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'mass',
     &'double',1
	write(IOpst,'(A)') "LOOKUP_TABLE default"
      do i=1,ntotal+m_nb
	  write (IOpst,*) mass(i)
	end do

      write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'internal_energy',
     &'double',1
	write(IOpst,'(A)') "LOOKUP_TABLE default"
      do i=1,ntotal+m_nb
	  write (IOpst,*) u(i)
	end do
      
	write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'itype',
     &'double',1
	write(IOpst,'(A)') "LOOKUP_TABLE default"
      do i=1,ntotal+m_nb
	  write (IOpst,*) itype(i)
	end do


	write(IOpst,'(A,1X,A,1X,A)') 'VECTORS', 'velocity',
     &'double'

        do i=1,ntotal+m_nb
	   
	    write (IOpst,*) vxtemp(1,i), 
     &		vxtemp(2,i),vxtemp(3,i)
	 
      end do

                          
      close(1)
      
      deallocate(xtemp, vxtemp)
      end
	
	
	
	
	

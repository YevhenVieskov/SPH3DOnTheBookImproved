	subroutine shepard(ntotal,hsml,mass,niac,pair_i,pair_j,w,
     &           itype,rho)

C----------------------------------------------------------------------
C   Subroutine to calculate the density with SPH summation algorithm.

C     ntotal : Number of particles                                  [in]
C     hsml   : Smoothing Length                                     [in]
C     mass   : Particle masses                                      [in]
C     niac   : Number of interaction pairs                          [in]
C     pair_i : List of first partner of interaction pair            [in]
C     pair_j : List of second partner of interaction pair           [in]
C     w      : Kernel for all interaction pairs                     [in]
c     itype   : type of particles                                   [in]
c     x       : Coordinates of all particles                        [in]
c     rho    : Density                                             [out]
      use interf
      implicit none
      include 'param.inc'
      
      integer ntotal, niac, pair_i(:),
     &        pair_j(:), itype(:)  
      double precision hsml(:),mass(:), w(:),
     &       rho(:) 
      integer i, j, k, d      
      double precision selfdens, r     
      double precision,allocatable::hv(:), wi(:) 
      allocate(hv(dim),  wi(maxn))
c     wi(maxn)---integration of the kernel itself
        
      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      
c     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*mass(i)/rho(i)
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        wi(i) = wi(i) + mass(j)/rho(j)*w(k)
        wi(j) = wi(j) + mass(i)/rho(i)*w(k)
      enddo

c     Secondly calculate the rho integration over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        rho(i) = selfdens*mass(i)
      enddo

c     Calculate SPH sum for rho:
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        rho(i) = rho(i) + mass(j)*w(k)
        rho(j) = rho(j) + mass(i)*w(k)
      enddo

c     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
     
      
        do i=1, ntotal
          rho(i)=rho(i)/wi(i)
        enddo
      
      deallocate (hv,  wi)
      end
      
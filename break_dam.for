      subroutine bdam(x, vx, mass, rho, p, u,
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------
c     This subroutine is used to generate initial data for the
c     2 d shear driven cavity probem with Re = 1
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =2   water
c     h-- smoothing lengths of particles                           [out]
c     ntotal-- total particle number                               [out]
c      use interf
      use ex
      implicit none
      include 'param.inc'

      integer itype(:)[*], ntotal
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*],
     &     rho(:)[*], p(:)[*], u(:)[*], hsml(:)[*]
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy !,c(:)

c     Giving mass and smoothing length as well as other data.
      if(this_image()==1) then
      m =2*20+1  !2*20+1
      n =2*40+1  !2*40+1
      mp = m-1
      np = n-1
      ntotal = mp * np
      xl =0.2    !1.e-3
      yl =0.4   !1.e-3
      dx = xl/mp
      dy = yl/np
      dx = xl/mp
      dy = yl/np

      do i = 1, mp
	do j = 1, np
	  k = j + (i-1)*np
	  x(1, k) = (i-1)*dx + dx/2.
	  x(2, k) = (j-1)*dy + dy/2.
        enddo
      enddo

      do i = 1, mp*np
	vx(1, i) = 0.
	vx(2, i) = 0.
        rho (i) = 1000.
        mass(i) = dx*dy*rho(i)
        p(i)= 9.81d0*rho(i)*(0.2-x(2,i))   !+1.0d5
        u(i)=357.1
        itype(i) = 2
        hsml(i) = dx
      enddo

c	initial particle distribution for tensile correction
	 dpt=dx
c	call vtk(x, vx, mass, rho, p, u, c, itype, hsml, ntotal,
c     &		1111)
       end if
      end

	subroutine virt_dam(tm,dt, ntotal,nvirt,hsml,mass,x,vx,
     &           rho,u,p,itype)

c----------------------------------------------------------------------
c   Subroutine to determine the information of virtual particles
c   Here only the Monaghan type virtual particles for the 2D shear
c   cavity driven problem are generated.
c     itimestep : Current time step                                 [in]
c     ntotal : Number of particles                                  [in]
c     nvirt  : Number of virtual particles                         [out]
c     hsml   : Smoothing Length                                 [in|out]
c     mass   : Particle masses                                  [in|out]
c     x      : Coordinates of all particles                     [in|out]
c     vx     : Velocities of all particles                      [in|out]
c     rho    : Density                                          [in|out]
c     u      : internal energy                                  [in|out]
c     itype   : type of particles                               [in|out]
      use interf
      implicit none
      include 'param.inc'
      integer itimestep, ntotal, nvirt, itype(:)
      double precision hsml(:),mass(:),x(:,:),vx(:,:),
     &                 rho(:), u(:), p(:)
      integer i, j, d, im, mp
      double precision xl, dx, v_inf,dt,tm,yl,dy
      if(this_image()==1) then
      if (vp_input) then

        open(1,file="../data/xv_vp.dat")
        open(2,file="../data/state_vp.dat")
        open(3,file="../data/other_vp.dat")
        read(1,*) nvirt
        do j = 1, nvirt
          i = ntotal + j
          read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)
          read(2,*)im, mass(i), rho(i), p(i), u(i)
          read(3,*)im, itype(i), hsml(i)
        enddo
        close(1)
        close(2)
        close(3)

	else

	nvirt = 0
        mp =2*40   !100      !40
	xl = 0.8d0
	yl=0.6
	dx = xl / mp
      dy=yl/mp
	v_inf = 1.e-3


c     Monaghan type virtual particle on the Lower side

        do i = 1, 2*mp+1
   	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = (i-1)*dx/2
          x(2, ntotal + nvirt) = 0.
          vx(1, ntotal + nvirt) = 0.
	  vx(2, ntotal + nvirt) = 0.
        enddo

c     Monaghan type virtual particle on the Left side

        do i = 1, 2*mp-1
   	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = 0.
          x(2, ntotal + nvirt) = i*dy/2
          vx(1, ntotal + nvirt) = 0.
	  vx(2, ntotal + nvirt) = 0.
        enddo
c     Monaghan type virtual particle on the Right side

      do i = 1, 2*mp-1
   	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = xl
          x(2, ntotal + nvirt) = i*dx/2
          vx(1, ntotal + nvirt) = 0.
	  vx(2, ntotal + nvirt) = 0.
        enddo

	do i = 1, nvirt
	  rho (ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
	  p(ntotal + i) = 1d5+x(2,i)
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -2
	  hsml(ntotal + i) = dx
        enddo

      endif

      if (mod(itimestep,save_step).eq.0) then
        open(1,file="../data/xv_vp.dat")
        open(2,file="../data/state_vp.dat")
        open(3,file="../data/other_vp.dat")
        write(1,*) nvirt
        do i = ntotal + 1, ntotal + nvirt
          write(1,1001) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)
          write(2,1002) i, mass(i), rho(i), p(i), u(i)
          write(3,1003) i, itype(i), hsml(i)
        enddo
1001    format(1x, I6, 6(2x, e14.8))
1002    format(1x, I6, 7(2x, e14.8))
1003    format(1x, I6, 2x, I4, 2x, e14.8)
        close(1)
        close(2)
        close(3)
      endif
       itimestep=int(tm/dt)
      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
         print *,' >> Statistics: Virtual boundary particles:'
         print *,'          Number of virtual particles:',NVIRT
        endif
      endif
      if(this_image()==1) then
      end

      subroutine fs_dam(ntotal,hsml,mass,niac,pair_i,pair_j,w,dwdx,
     &	vx,itype,x,rho,drhodt)

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
      integer i, j, k, d,ki
      double precision selfdens, r,
     &dwdx(:, :),
     &vx(:,:), x(:,:), drhodt(:)
      double precision    vcc
	double precision rho_ref
      double precision,allocatable:: dvx(:),rho_sum(:),drhodt_temp(:)
      double precision,allocatable:: hv(:), wi(:)
      allocate(dvx(dim),rho_sum(maxn),drhodt_temp(maxn))
      allocate(hv(maxn), wi(maxn))
      do d=1,dim
        hv(d) = 0.d0
      enddo
      do i=1,maxn
	 rho_sum(i)=rho(i)
	end do
c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      rho_ref=1000d0
c     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*mass(i)/rho_sum(i)
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        wi(i) = wi(i) + mass(j)/rho_sum(j)*w(k)
        wi(j) = wi(j) + mass(i)/rho_sum(i)*w(k)
      enddo

c     Secondly calculate the rho integration over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        rho_sum(i) = selfdens*mass(i)
      enddo

c     Calculate SPH sum for rho:
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        rho_sum(i) = rho_sum(i) + mass(j)*w(k)
        rho_sum(j) = rho_sum(j) + mass(i)*w(k)
      enddo

c     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

      if (nor_density) then
        do i=1, ntotal
          rho_sum(i)=rho_sum(i)/wi(i)
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i = 1, ntotal
        drhodt(i) = 0.
	 drhodt_temp(i) = 0.
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j)
        enddo
        vcc = dvx(1)*dwdx(1,k)
        do d=2,dim
          vcc = vcc + dvx(d)*dwdx(d,k)
        enddo
        drhodt_temp(i) = drhodt_temp(i) + mass(j)*vcc
        drhodt_temp(j) = drhodt_temp(j) + mass(i)*vcc
      enddo

      do i=1,ntotal
	if((rho_sum(i)/rho_ref<1d0).or.(rho_sum(i)/rho_ref>=1.05)) then
	 drhodt(i)=drhodt_temp(i)
      else
	 rho(i)=rho_sum(i)
	end if
      end do
      deallocate(dvx,rho_sum,drhodt_temp)
      deallocate(hv, wi)
      end


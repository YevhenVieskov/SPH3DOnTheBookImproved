 
      subroutine riemann(ntotal,hsml,x,vx,mass,vol,e,flux,  
     &     niac,pair_i,pair_j,w,dwdx,rho,p) 
       
c----------------------------------------------------------------------------------------------------- 
c       Subroutine to calculate, for each particle i, the density,  
c       velocity and internal energy Roe flux averaged over all  
c        the particles in the support domain of particle i       
c----------------------------------------------------------------------------------------------------- 
 
 
      implicit none 
      include 'param.inc' 
       
      integer ntotal, i, j, h, k,l,pp,qq 
      double precision hsml(:), x(:,:), mass(:), vol(:),rho(:), 
     &        p(:), e(:)   
      double precision w(:), dwdx(:,:) 
      integer m, scale_k, niac, pair_i(:) 
      integer pair_j(:) 
      double precision  driac, mhsml
      double precision  r(3,3),beta(3)
      double precision voll,volr,vl,vr,el,er,ul,ur,pl,pr,cl,cr 
      double precision lambda(3,3),volmean,pmean,vmean,cmean,gamma 
      double precision pdvolmean,pdumean,lf(3),rf(3),roef(3),aroe(3) 
      double precision alpharoe,vmax, Dis_R 
      external DeltaR, r_3 
      double precision DeltaR, r_3, ri(3,3),prod(3,3),mat_roe(3,3) 
      double precision selfdens,dselfdens 
      double precision,allocatable:: flux(:,:),vx(:,:)
      allocate(flux(3,maxn), vx(dim,maxn))
      gamma = 1.4 
       
      do i = 1, ntotal 
         do m = 1, 3 
            flux(m,i) = 0. 
         enddo  
      enddo 
      
      do k=1,niac   
           i = pair_i(k) 
         j = pair_j(k) 
          
C Left and right state 
         voll=vol(i) 
         volr=vol(j) 
         vl=vx(1,i) 
         vr=vx(1,j) 
         el=e(i) 
         er=e(j) 
         pl=p(i) 
         pr=p(j) 
          
          
c The eigenvalues, eigenvectors, wave strength coefficients and  
c fluxes at the left and right states of the interface are  
c calculated 
       
C Mean values 
         volmean = 0.5 * (voll + volr) 
         pmean = 0.5 * (pl + pr) 
         vmean = 0.5 * (vl + vr) 
         cmean = sqrt(gamma * pmean / volmean) 
         pdvolmean = -pmean/volmean 
         pdumean = (gamma - 1.) / volmean 
          
C Eigenvalues of Jacobian matrix 
         do pp = 1,3 
            do qq =1,3 
               lambda(pp,qq) = 0.0 
            enddo 
         enddo 
 
         lambda(1,1)=-cmean 
         if (dabs(vl).ge.dabs(vr)) then 
            vmax=dabs(vl) 
         else 
            vmax=dabs(vr) 
         endif 
         alpharoe=0.5 
         lambda(2,2)=alpharoe*vmax 
         lambda(3,3)=cmean 
          
C Eigenvectors of Jacobian matrix : R matrix 
         r(1,1)=1.0 
         r(2,1)=-lambda(1,1) 
         r(3,1)= r_3(pdvolmean,pdumean,lambda(1,1),vmean) 
          
         r(1,2)=1.0 
         r(2,2)=-lambda(2,2) 
         r(3,2)= r_3(pdvolmean,pdumean,lambda(2,2),vmean) 
          
         r(1,3)=1.0 
         r(2,3)=-lambda(3,3) 
         r(3,3)= r_3(pdvolmean,pdumean,lambda(3,3),vmean) 
 
C Discriminant of R 
 
         Dis_R = DeltaR(r(2,2),r(3,2),r(2,3),r(3,3)) - 
     &           DeltaR(r(2,1),r(3,1),r(2,3),r(3,3)) + 
     &           DeltaR(r(2,1),r(3,1),r(2,2),r(3,2)) 
 
C R-1 matrix 
         ri(1,1) = DeltaR(r(2,2),r(3,2),r(2,3),r(3,3)) / Dis_R 
         ri(1,2) = - DeltaR(r(1,2),r(3,2),r(1,3),r(3,3)) / Dis_R 
         ri(1,3) = DeltaR(r(1,2),r(2,2),r(1,3),r(2,3)) / Dis_R 
 
         ri(2,1) = - DeltaR(r(2,1),r(3,1),r(2,3),r(3,3)) / Dis_R 
         ri(2,2) = DeltaR(r(1,1),r(3,1),r(1,3),r(3,3)) / Dis_R 
         ri(2,3) = - DeltaR(r(1,1),r(2,1),r(1,3),r(2,3)) / Dis_R 
 
         ri(3,1) = DeltaR(r(2,1),r(3,1),r(2,2),r(3,2)) / Dis_R 
         ri(3,2) = - DeltaR(r(1,1),r(3,1),r(1,2),r(3,2)) / Dis_R 
         ri(3,3) = DeltaR(r(1,1),r(2,1),r(1,2),r(2,2)) / Dis_R 
 
C A = R * |Lambda| * R-1 
         lambda(1,1)=dabs(lambda(1,1)) 
         lambda(2,2)=dabs(lambda(2,2)) 
         lambda(3,3)=dabs(lambda(3,3)) 
         call Two_Matrix_Product(r,lambda,prod,3,3,3,3,3,3) 
         call Two_Matrix_Product(prod,ri,mat_roe,3,3,3,3,3,3) 
C Eigen flux   
         aroe(1) = mat_roe(1,1)*(volr - voll)+mat_roe(1,2)*(vr - vl)+ 
     %             mat_roe(1,3)*(er - el) 
         aroe(2) = mat_roe(2,1)*(volr - voll)+mat_roe(2,2)*(vr - vl)+ 
     %             mat_roe(2,3)*(er - el) 
         aroe(3) = mat_roe(3,1)*(volr - voll)+mat_roe(3,2)*(vr - vl)+ 
     %             mat_roe(3,3)*(er - el) 
 
C Left and right system flux 
         lf(1) = -vl 
         lf(2) = pl 
         lf(3) = vl*pl 
         rf(1) = -vr 
         rf(2) = pr 
         rf(3) = vr*pr 
 
C Roe flux 
         do l=1,3 
            roef(l)=0.5*(lf(l)+rf(l))-0.5*aroe(l) 
         enddo       
 
         write(10,*)rho(i),vol(i),rho(i)*vol(i) 
C Particle flux 
         do l=1,3 
            flux(l,i)=flux(l,i)+2.*mass(j)*roef(l)* 
     &           dwdx(1,k)*vol(j)*vol(i) 
            flux(l,j)=flux(l,j)-2.*mass(i)*roef(l)* 
     &           dwdx(1,k)*vol(i)*vol(j) 
         enddo 
      enddo 
      
       deallocate(flux, vx)
       
      end 
       
      double precision function DeltaR(da,db,dc,dd) 
 
C  Calculation of the discriminant 2x2 
 
      double precision da,db,dc,dd 
 
      DeltaR = da * dd - db * dc 
 
      return  
      end 
 
      double precision function r_3(dpv,dpu,dlam,dv) 
 
C  Calculation of third term of eigenvector R 
      
      double precision dpv,dpu,dlam,dv 
 
      r_3 = -(dpv + dlam**2) / dpu - dlam * dv 
      
      return 
      end 
 
      subroutine Two_Matrix_Product(mat1,mat2,dst,nbl1,nbc1,nbl2,nbc2, 
     &                              diml,dimc) 
c ---------------------------------------------------------------------------------------------------- 
c     Product of mat1 and mat2. The result goes to dst (destination)  
c     nbl1,nbc1,nbl2,nbc2,nbldst,nbcdst : nb of lines and columns of matrices 
c         diml,dimc : dimension of three matrices (nb maxi) 
c ---------------------------------------------------------------------------------------------------- 
      implicit none 
      integer ll,cc,kk,nbl1,nbc1,nbl2,nbc2,diml,dimc 
      double precision som,mat1(diml,dimc),mat2(diml,dimc) 
      double precision dst(diml,dimc) 
      if(nbc1.ne.nbl2) then 
        print *, 'Error: product of two matrices of incompatible size' 
        stop 
      endif 
      do 10 ll=1,nbl1 
        do 20 cc=1,nbc2 
          som=0 
          do 30 kk=1,nbc1 
            som=som+mat1(ll,kk)*mat2(kk,cc) 
 30      enddo 
          dst(ll,cc)=som 
 20    enddo 
 10   enddo 
 
      end 


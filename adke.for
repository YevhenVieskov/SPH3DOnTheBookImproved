      subroutine adke(ntotal,mass,itype,rho,x,hsml, w, dwdx)
         use ex
      implicit none
      include 'param.inc'

      integer itimestep, ntotal, itype(:)
      double precision  hsml(:), mass(:), x(:,:),
     &        rho(:)
      integer  nvirt, niac, pair_i(:),
     &        pair_j(:), ns(:)

      double precision w(:), dwdx(:,:)


      double precision ak,eps,gm(:),alam(:)
      integer na(maxn),k,i,j,d

      allocate(gm(maxn), alam(maxn))

         if (nnps.eq.1) then
        call direct_find(itimestep, ntotal,hsml,x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      else if (nnps.eq.2) then
        call link_list(itimestep, ntotal,hsml(1),x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      else if (nnps.eq.3) then
        call tree_search(itimestep, ntotal,hsml,x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      end if

c         call sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,
c    &           itype,rho)
       ak=1.0
       eps=0.5
       do i=1,ntotal
       na(i)=0
       gm(i)=0.0
       end do
       do k=1,niac
         i=pair_i(k)
         j=pair_j(k)
         na(i)=na(i)+1
         na(j)=na(i)+1
         gm(i)=gm(i)+log(rho(j))
         gm(j)=gm(j)+log(rho(i))
       end do
       do i=1,ntotal
         gm(i)=gm(i)/na(i)
       end do
       do i=1,ntotal

        gm(i)=exp(gm(i))
        alam(i)=ak*(rho(i)/gm(i))**(-eps)
        hsml(i)=hsml(i)*alam(i)
       end do

       deallocate(gm(maxn), alam(maxn))

      end

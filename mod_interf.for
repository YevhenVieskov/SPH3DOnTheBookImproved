      module interf
      
      interface
        subroutine art_heat(ntotal,hsml,mass,x,vx,niac,rho,u,
     &           c,pair_i,pair_j,w,dwdx,dedt)


        implicit none
      
      
        integer ntotal,niac,pair_i(:) [*],
     &        pair_j(:)[*]
        double precision hsml(:)[*], mass(:)[*], x(:,:)[*],vx(:,:)[*],
     &                 rho(:)[*], u(:)[*], c(:)[*],w(:)[*],
     &                 dwdx(:,:)[*], dedt(:)[*]
        integer i,j,k,d
        double precision dx, vr, rr, h, mc, mrho, mhsml, 
     &                 hvcc, mui, muj, muij, rdwdx, g1,g2
        double precision,allocatable::dvx(:),vcc(:)
        end
      end  interface
      
      
      interface
         subroutine art_visc(ntotal,hsml,mass,x,vx,niac,rho,c,
     &           pair_i,pair_j,w,dwdx,dvxdt,dedt) 
        implicit none
       
      
        integer ntotal, niac, pair_i(:)[*],
     &        pair_j(:)[*]            
        double precision hsml(:)[*], mass(:)[*], x(:,:)[*],vx(:,:)[*],
     &       rho(:)[*], c(:)[*], w(:)[*],
     &       dwdx(:,:)[*], dvxdt(:,:)[*], dedt(:)[*]
        integer i,j,k,d
        double precision dx, alpha, beta, etq, piv,
     &       muv, vr, rr, h, mc, mrho, mhsml
        double precision,allocatable:: dvx(:)
        end
      end  interface
      
      
      interface
        subroutine av_vel(ntotal,mass,niac,pair_i,pair_j,
     &           w, vx, rho, av)
  
      implicit none
  
      
      integer ntotal, niac, pair_i(:)[*],
     &        pair_j(:)[*]
      double precision   mass(:)[*],w(:)[*],
     &       vx(:,:)[*], rho(:)[*], av(:, :)[*]       
      integer i,j,k,d       
      double precision   vcc, epsilon
      double precision,allocatable:: dvx(:)[*]
      end
      end  interface
      
      
      interface
      subroutine beeman(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt )
     	
      implicit none     
            
      integer itype(:)[*], ntotal, maxtimestep
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*], 
     &       rho(:)[*], p(:)[*], u(:)[*], c(:)[*], s(:)[*], e(:)[*], 
     &       hsml(:)[*], dt
      integer i, j, k, itimestep, d, current_ts, nstart        
      double precision,allocatable::  x_min(:, :)[:], v_min(:, :)[:], 
     &      u_min(:)[*], rho_min(:)[:], dx(:,:)[:], dvx(:, :)[:],  
     &       du(:)[:], drho(:)[:],  av(:, :)[:], ds(:)[:], 
     &       t(:)[:], tdsdt(:)[:]         
      double precision  time, temp_rho, temp_u        
	
      double precision,allocatable:: xn(:, :)[*], vxn(:, :)[*], 
     &    rhon(:)[*], pn(:)[*], un(:)[*], xn1(:, :)[*], vxn1(:, :)[*], 
     &       rhon1(:)[*], pn1(:)[*], un1(:)[*], dvxn(:, :)[*], 
     &       drhon(:)[*], dun(:)[*],dvxn1(:, :)[*], 
     &       drhon1(:)[*], dun1(:)[*]
      end
      end  interface
      
      
      interface
      subroutine bdam(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

      implicit none     
           
      integer itype(:)[*], ntotal
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*],
     &     rho(:)[*], p(:)[*], u(:)[*], hsml(:)[*]
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy !,c(:)
      end
      end  interface
      
      
      interface
      subroutine virt_dam(tm,dt, ntotal,nvirt,hsml,mass,x,vx,
     &           rho,u,p,itype) 

      implicit none
      integer itimestep, ntotal, nvirt, itype(:)[*]
      double precision hsml(:)[*],mass(:)[*],x(:,:)[*],vx(:,:)[*],
     &                 rho(:)[*], u(:)[*], p(:)[*]
      integer i, j, d, im, mp
      double precision xl, dx, v_inf,dt,tm,yl,dy
      end
      end  interface
      
      
      interface
      subroutine fs_dam(ntotal,hsml,mass,niac,pair_i,pair_j,w,dwdx,
     &	vx,itype,x,rho,drhodt)
      
      implicit none
           
      integer ntotal, niac, pair_i(:)[*],
     &        pair_j(:)[*], itype(:)[*]  
      double precision hsml(:)[*],mass(:)[*], w(:)[*],
     &       rho(:) [*]
      integer i, j, k, d,ki      
      double precision selfdens,  r, 
     &dwdx(:, :)[*],
     &vx(:,:)[*], x(:,:)[*], drhodt(:)[*] 
      double precision    vcc
	double precision rho_ref
      double precision,allocatable:: dvx(:)[*],rho_sum(:)[*],
     &drhodt_temp(:)[*]
       double precision,allocatable:: hv(:), wi(:)
      end
      end  interface
      
      
      interface
       subroutine combine_visc(ntotal,hsml,mass,x,vx,niac,rho,c,
     &           pair_i,pair_j,w,dwdx,dvxdt,dedt,eta)


      implicit none
           
      integer ntotal, niac, pair_i(:)[*],
     &        pair_j(:) [*]           
      double precision hsml(:)[*], mass(:)[*], x(:,:)[*],vx(:,:)[*],
     &       rho(:)[*], c(:)[*], w(:)[*],
     &  dwdx(:,:)[*], dvxdt(:,:)[*],dedt(:)[*],eta(:)[*]
      integer i,j,k,d
      double precision dx, alpha, beta, etq, piv,
     &       muv, vr, rr, h, mc, mrho, mhsml
      double precision dxm,dym,dzm,drm,mer
      double precision  driac, r
      double precision,allocatable:: dvx(:),dxiac(:)
      end
      end  interface
      
      
      interface
       subroutine sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,
     &           itype,rho)


      implicit none
          
      integer ntotal, niac, pair_i(:)[*],
     &        pair_j(:)[*], itype(:) [*] 
      double precision hsml(:)[*],mass(:)[*], w(:)[*],
     &       rho(:)[*] 
      integer i, j, k, d      
      double precision selfdens,  r   
      double precision,allocatable:: hv(:)[*],wi(:)[*]
      end
      end  interface
      
      
      interface
      subroutine con_density(ntotal,mass,niac,pair_i,pair_j,
     &           dwdx,vx,itype,x,rho, drhodt)

      implicit none
         
      integer ntotal,niac,pair_i(:)[*],
     &        pair_j(:)[*], itype(:)[*]   
      double precision mass(:)[*], dwdx(:, :)[*],
     &       vx(:,:)[*], x(:,:)[*], rho(:)[*], drhodt(:)[*]
      integer i,j,k,d    
      double precision    vcc 
      double precision,allocatable:: dvx(:) 
      end
      end  interface
      
      
      interface
      subroutine direct_find(itimestep, ntotal,hsml,x,niac,pair_i,
     &           pair_j,w,dwdx,countiac)

      implicit none
      
      integer itimestep, ntotal,niac,pair_i(:)[*],
     &        pair_j(:)[*], countiac(:)[*]
      double precision hsml(:)[*], x(:,:)[*], w(:)[*],
     &       dwdx(:,:)[*]
      integer i, j, d,  sumiac, maxiac, miniac, noiac,
     &        maxp, minp, scale_k 
      double precision  driac, r, mhsml  
      double precision,allocatable:: dxiac(:),tdwdx(:)
      end
      end  interface
      
      
      interface
       subroutine ext_force(ntotal,mass,x,niac,pair_i,pair_j,
     &           itype,hsml,dvxdt)

      implicit none
     
      integer ntotal, itype(:)[*], niac,
     &        pair_i(:)[*], pair_j(:)[*]
      double precision mass(:)[*], x(:,:)[*], hsml(:)[*],          
     &       dvxdt(:,:)[*] !,c(:)
      integer i, j, k, d
      double precision  rr, f, rr0, dd, p1, p2,q,am,gam     
      double precision,allocatable::dx(:)
      end
      end  interface
      
      
      interface
      subroutine g_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	implicit none
     
      integer itimestep, ntotal, nghost, itype(:)[*]
      double precision hsml(:)[*],mass(:)[*],x(:,:)[*],vx(:,:)[*],
     &                 rho(:)[*], u(:)[*], p(:)[*]
      integer i, j, d, im, mp, scale_k
      double precision xleft,xright,ytop,ybottom
	double precision dr,dx,xl,v_inf,dr1,dr2,dr3,dr4,skh      
	logical cond
      double precision,allocatable::dxiac(:)
      end
      end  interface
      
      
      interface
      subroutine dbp_part(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	
	implicit none
   
      integer itimestep, ntotal, nghost, itype(:)[*]
      double precision hsml(:)[*],mass(:)[*],x(:,:)[*],vx(:,:)[*],
     &                 rho(:)[*], u(:)[*], p(:)[*]
      integer i, j, d, im, mp, scale_k,mp2
      double precision xleft,xright,ytop,ybottom
	double precision dr,dx,xl,v_inf,dr1,dr2,xl2
      double precision,allocatable:: dxiac(:)
      end
      end  interface
      
      
      interface
      subroutine gmod_part(itimestep, ntotal,ngm,hsml,mass,x,vx,
     &           rho,u,p,itype) 

      implicit none
      
      integer itimestep, ntotal, ndrp, itype(:),ngm
      double precision hsml(:)[*],mass(:)[*],x(:,:)[*],vx(:,:)[*],
     &                 rho(:)[*], u(:)[*], p(:)[*]
      integer i, j, d, im, mp
      double precision xl, dx, v_inf
      end
      end  interface
      
      
      interface
      subroutine grid_geom(i,x,ngridx,maxgridx,mingridx,dgeomx,xgcell)

      implicit none
      
      integer i, ngridx(:),xgcell(3)
      double precision x(:), maxgridx(:), mingridx(:), dgeomx(:)
      integer d
      end
      end  interface
      
      
      interface
      subroutine h_upgrade(dt,ntotal, mass, vx, rho, niac, 
     &           pair_i, pair_j, dwdx, hsml)

      implicit none
            
      integer ntotal, niac, pair_i(:)[*], 
     &        pair_j(:)[*]
      double precision mass(:)[*], vx(:, :)[*], rho(:)[*],
     &       dwdx(:, :)[*], hsml(:)[*]     
      integer i,j,k,d
      double precision dt, fac, hvcc   
      double precision,allocatable:: dvx(:),vcc(:), dhsml(:)
      end
      end  interface
      
      
      interface
       subroutine init_grid(ntotal,hsml,grid,ngridx,ghsmlx,
     &           maxgridx,mingridx,dgeomx)

      implicit none
     
      !integer maxngx,maxngy,maxngz      
      integer ntotal,  ngridx(:), ! grid(maxngx,maxngy,maxngz),
     &ghsmlx(:),grid(:,:,:)
      double precision hsml, maxgridx(:), mingridx(:), dgeomx(:)
   
      integer i, j, k, d,  ngrid(3)
      integer,allocatable::maxng(:)
      end
      end  interface
      
      
      interface
      subroutine input(x, vx, mass, rho, p, u, itype, hsml, ntotal)
      
      implicit none     
      integer itype(:)[*], ntotal       
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*], 
     &                 p(:)[*], u(:)[*], hsml(:)[*], rho(:)[*]
      integer i, d, im  
      end
      end  interface
      
      
      interface
      subroutine int_force(itimestep,dt,ntotal,hsml,mass,vx,niac,rho,
     &        eta,pair_i,pair_j,w,dwdx,u,itype,x,t,c,p,dvxdt,tdsdt,dedt)

      implicit none
            
      integer itimestep, ntotal,niac,pair_i(:)[*],
     &        pair_j(:)[*], itype(:) [*]
      double precision dt, hsml(:)[*], mass(:)[*], vx(:,:)[*],      
     &       rho(:)[*], eta(:)[*], dwdx(:,:)[*], u(:)[*],
     &       x(:,:)[*], t(:)[*], c(:)[*], p(:)[*], dvxdt(:,:)[*],          
     &       tdsdt(:)[*],dedt(:)[*],w(:)[*]
      integer i, j, k, d
      double precision  hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij

	double precision nu0,dx,dy,dz,dr,mer,ksi
	double precision  driac, r
      double precision Ri,Rj,Rp,nt,wp,mhsml,fp !!!!!
      double precision delm_i,delm_j,delp_ij,eps1,eps2
	integer n
      double precision,allocatable:: dvx(:), txx(:), tyy(:),
     &tzz(:), txy(:), txz(:), tyz(:), vcc(:),dxiac(:),dwdxp(:)
      end
      end  interface
      
      
      interface
      subroutine kernel(r,dx,hsml,w,dwdx)   

      implicit none            
      double precision r, dx(:), hsml, w, dwdx(:)
      integer i, j, d      
      double precision q, dw, factor
      end
      end  interface
      
      
      interface
      subroutine link_list(itimestep, ntotal,hsml,x,niac,pair_i,
     &           pair_j,w,dwdx,countiac)
     

      implicit none
      
      !integer maxngx,maxngy,maxngz        
      integer itimestep, ntotal, niac, pair_i(:)[*],
     &        pair_j(:)[*], countiac(:)
      double precision hsml, x(:,:)[*],w(:),
     &       dwdx(:,:)
      integer i, j, d, scale_k, sumiac, maxiac, noiac, miniac, maxp,minp    
c      integer grid(maxngx,maxngy,maxngz),gcell(3),
c     &     xcell,ycell,zcell,minxcell(3),maxxcell(3)
          
      double precision hsml2,dr,r
      integer,allocatable:: ngridx(:), ghsmlx(:),xgcell(:,:), 
     &celldata(:), dnxgcell(:),dpxgcell(:)
      double precision,allocatable:: dx(:), mingridx(:),maxgridx(:),
     &tdwdx(:), dgeomx(:)
      end
      end  interface
      
      
      interface
      subroutine outflow(Xleft,Xright,Ytop,dx,dy,vel,ntotal,nbound, 
     &							nout,x,vx, mass, rho, p, u,itype, hsml)

      implicit none     
      
      integer itype(:)[*], ntotal,nbound,nout
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*],
     &     rho(:)[*], p(:)[*], u(:)[*], hsml(:)[*],vel
      double precision Xleft,Xright,Ytop,Ybottom,dx,dy
      integer i, j, d, m, n, mp, np, k,scale_k
      double precision xl, yl !,c(maxn) !!!!
	double precision outlength,treshold 
      double precision,allocatable:: xtemp(:, :)[*]
      end
      end  interface
      
      
      interface
      subroutine sort_outflow(Xleft,Xright,Ytop,dx,dy,vel,ntotal, 
     &							x,vx, mass, rho, p, u,itype, hsml)

      implicit none     
      
      integer itype(:)[*], ntotal,nbound,nout
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*],
     &     rho(:)[*], p(:)[*], u(:)[*], hsml(:)[*],vel
      double precision Xleft,Xright,Ytop,Ybottom,dx,dy
      integer i, j, d, m, n, mp, np, k,scale_k
      double precision xl, yl  !,c(maxn)!!!!!!
	double precision outlength,treshold
      double precision,allocatable::xtemp(:, :), vxtemp(:, :), 
     &    masstemp(:), rhotemp(:), ptemp(:), utemp(:), hsmltemp(:)
      integer,allocatable::mark(:),itypetemp(:)
      end
      end  interface
      
      
      interface
      subroutine output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal) 
      
      implicit none     
          
      integer itype(:)[*], ntotal
      double precision x(:, :),[*] vx(:, :)[*], mass(:)[*], 
     &       rho(:)[*],p(:)[*], u(:)[*], c(:)[*], hsml(:)[*]
      integer i, d, npart     
      end
      end  interface
      
      
      interface
      subroutine prcorr(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
     
      implicit none    
      integer itype(:)[*], ntotal, maxtimestep
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*], 
     &       rho(:)[*], p(:)[*], u(:)[*], c(:)[*], s(:)[*], e(:)[*], 
     &       hsml(:)[*], dt
      integer i, j, k, itimestep, d, current_ts, nstart        
      !t,tdsdt      
      double precision  time, temp_rho, temp_u       
	 
	double precision,allocatable:: x_min(:, :), v_min(:, :), u_min(:),
     &       rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:),  xn(:, :), vxn(:, :),  
     &       rhon(:), pn(:), un(:) , xn12(:, :), vxn12(:, :),  
     &       rhon12(:), pn12(:), un12(:) 
      end
      end  interface
      
      
      interface
      subroutine puzir(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

      implicit none     
      
      integer itype(:), ntotal
      double precision x(:, :), vx(:, :), mass(:),
     &     rho(:), p(:), u(:), hsml(:)  !, c(:)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy
      double precision radius,cr,xr,yr,r2
      end
      end  interface
      
      
      interface
      subroutine rt(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

      implicit none     
           
      integer itype(:), ntotal
      double precision x(:, :), vx(:, :), mass(:),
     &     rho(:), p(:), u(:), hsml(:)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy !,c(:)
      double precision yr
      end
      end  interface
      
      
      interface
      subroutine refine(x, vx, mass, rho, p, u, c,itype,hsml,ilevel, 
     &ntotal,nbound,vcrit)

      implicit none
      
      integer  ntotal,nbound, itype(:)    
      double precision  hsml(:), mass(:), x(:,:),
     &       vx(:,:), u(:), rho(:), p(:), c(:)
     &                    
c      double precision w(:), dwdx(:,:)
c      integer i, d, niac, pair_i(:), !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     &        pair_j(:), ns(:)
      integer ilevel(:), i
      
	double precision vcrit 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	logical crit_type !crit_type=0-velocity, crit_type=1-velocity
c	temporary array for store initial particle distribution
	integer   ntotalt  
      double precision sp,sr,dref,size,xc,yc,size4
	
      integer,allocatable::itypet(:), iref(:),ilevelt(:),ireft(:)
      double precision,allocatable:: hsmlt(:), masst(:), xt(:,:),
     &  vxt(:,:), ut(:), rhot(:), pt(:), ct(:)
      end
      end  interface
                
      interface
      subroutine sym_verlet(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
     
           implicit none     
           
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt
      integer i, j, k, itimestep, d, current_ts, nstart     
      double precision  time, temp_rho, temp_u
     	double precision dt2
      double precision,allocatable::x_min(:, :), v_min(:, :), u_min(:),
     &       rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:),
     &       xn(:, :), vxn(:, :),  
     &       rhon(:), pn(:), un(:),
     &       xn12(:, :), vxn12(:, :),  
     &       rhon12(:), pn12(:), un12(:) 
      end
      end interface
      
      
      interface
      subroutine sym_verlet2(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
     

	
      implicit none     
      
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt
      integer i, j, k, itimestep, d, current_ts, nstart        
              
      double precision  time, temp_rho, temp_u
       
	      
      double precision,allocatable:: x_min(:, :), v_min(:, :), u_min(:),
     &       rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:), xn(:, :), vxn(:, :),  
     &       rhon(:), pn(:), un(:), xn12(:, :), vxn12(:, :),  
     &       rhon12(:), pn12(:), un12(:) 
      double precision dt2
      end
      end interface
     
      interface
      subroutine sym_euler(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
     

      implicit none     
   
      
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt
      integer i, j, k, itimestep, d, current_ts, nstart        
             
      double precision  time, temp_rho, temp_u
       
		
	double precision dt2
	
      double precision,allocatable::x_min(:, :), v_min(:, :), u_min(:),
     &       rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:)
      double precision,allocatable::xn(:, :), vxn(:, :),  
     &       rhon(:), pn(:), un(:)
      double precision,allocatable::xn12(:, :),  
     &vxn12(:, :),  rhon12(:), pn12(:), un12(:)
      end
      end interface
     
      interface
      subroutine tank(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

     
      implicit none     
     
      
      integer itype(:), ntotal
      double precision x(:, :), vx(:, :), mass(:),
     &     rho(:), p(:), u(:), hsml(:) !,c(:)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy
      double precision XMleft,XMright,YMbottom,YMtop,XBleft,XBright
      double precision YBbottom,YBtop,XPleft,XPright,YPbottom,YPtop 
	double precision  dx2,dy2
	logical cplate,cmleft,cmright,cnull
      double precision,allocatable:: xtemp(:, :)
      end
      end interface
     
     
c      interface
c      subroutine prcorr(x,vx, mass, rho, p, u, c, s, e, itype, 
c     &        hsml, ntotal, maxtimestep, dt)
     
c      implicit none     
          
c      integer itype(:), ntotal, maxtimestep
c      double precision x(:, :), vx(:, :), mass(:), 
c     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
c     &       hsml(:), dt
c      integer i, j, k, itimestep, d, current_ts, nstart        
      !t,tdsdt      
c      double precision  time, temp_rho, temp_u       
	 
c	double precision,allocatable:: x_min(:, :), v_min(:, :), u_min(:),
c     &       rho_min(:), dx(:,:), dvx(:, :), du(:),  
c     &       drho(:),  av(:, :), ds(:), 
c     &       t(:), tdsdt(:),  xn(:, :), vxn(:, :),  
c     &       rhon(:), pn(:), un(:) , xn12(:, :), vxn12(:, :),  
c    &       rhon12(:), pn12(:), un12(:) 
c      end
          
c      end interface
     
      interface
      subroutine single_step(itimestep, dt, ntotal, hsml, mass, x, vx,  
     &           u, s,rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av) 
      implicit none     

      integer itimestep, ntotal, itype(:)[*]    
      double precision dt, hsml(:)[*], mass(:)[*], x(:,:)[*],
     &       vx(:,:)[*], u(:)[*], s(:)[*], rho(:)[*], p(:)[*] !, eta(:) !, c(:)
      integer i, d, nvirt, niac     
                  
      integer nghost,ndbp,nout                        
      integer,save,allocatable:: pair_i(:)[*], pair_j(:)[*], ns(:)
      double precision,allocatable:: tdsdt(:), t(:), dx(:,:), dvx(:,:), 
     &       du(:), ds(:), drho(:), av(:, :) 
      double precision,allocatable:: w(:), dwdx(:,:),  
     &       indvxdt(:,:),exdvxdt(:,:),ardvxdt(:,:),  
     &       avdudt(:), ahdudt(:)
      end
      end interface
      
      interface
      
      subroutine vtk(x, vx, mass, rho, p, u, c, itype, hsml, ntotal,num) 
      implicit none   
          
      integer itype(:)[*], ntotal,num
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*], 
     &       rho(:)[*],p(:)[*], u(:)[*], c(:)[*], hsml(:)[*]
      double precision,allocatable::xtemp(:, :)[*], vxtemp(:, :)[*]
      end
      end interface
      
      interface
      subroutine dbp_part_tank(itimestep, ntotal,nghost,hsml,mass,x,vx,
     &           rho,u,p,itype)
	 
	implicit none
      integer itimestep, ntotal, nghost, itype(:)
      double precision hsml(:),mass(:),x(:,:),vx(:,:),
     &                 rho(:), u(:), p(:)
      integer i, j, d, im, mp, scale_k,mp2,M,N,L
      double precision xleft,xright,ytop,ybottom
	double precision dr,dx,dy,xl,v_inf,dr1,dr2,xl2  !, dxiac(maxn) 
	double precision XMleft,XMright,YMbottom,YMtop,XBleft,XBright
      double precision YBbottom,YBtop,XPleft,XPright,YPbottom,YPtop 
	double precision  dx2,dy2
      end
      end interface
      
      interface
      subroutine viscosity(ntotal,itype,x,rho,eta)
      implicit none
      integer ntotal,i,itype(:)[*]
      double precision x(:,:)[*],rho(:)[*],eta(:)
      end
      end interface
      
      interface
       subroutine tree_search(itimestep, ntotal,hsml,x,niac,pair_i,
     &           pair_j,w,dwdx,countiac)
	
      LOGICAL INTERSECT
	INTEGER SPTR,Q,P,NIAC
	REAL(8) DRIAC,MHSML !, PSIZE(MXBODY)
	REAL(8)DX,DY,DZ,DPC
	INTEGER PAIR_I(:), PAIR_J(:)
	integer i, j, d,  sumiac, maxiac, miniac, noiac,
     &        maxp, minp, scale_k ,countiac(:)
	double precision r
	double precision  x(:,:), w(:),
     &       dwdx(:,:),HSML(:)
		real(4) startt, endt
      double precision, allocatable:: dxiac(:), STACK(:), tdwdx(:) 
      end
      end interface
      
      interface
      subroutine shepard(ntotal,hsml,mass,niac,pair_i,pair_j,w,
     &           itype,rho)

      implicit none
      
      integer ntotal, niac, pair_i(:)[*],
     &        pair_j(:)[*], itype(:)[*]  
      double precision hsml(:)[*],mass(:)[*], w(:)[*],
     &       rho(:)[*] 
      integer i, j, k, d      
      double precision selfdens, r     
      double precision,allocatable::hv(:), wi(:) 
      end
      end interface
      
      interface
      subroutine time_integration(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt )
     
      implicit none     
            
      integer itype(:)[*], ntotal, maxtimestep
      double precision x(:, :)[*], vx(:, :)[*], mass(:)[*], 
     &       rho(:)[*], p(:)[*], u(:)[*], c(:)[*], s(:)[*], e(:)[*], 
     &       hsml(:)[*], dt
      integer i, j, k, itimestep, d, current_ts, nstart     
            
      double precision  time, temp_rho, temp_u
       double precision,allocatable::  x_min(:, :), v_min(:, :), 
     &       u_min(:),rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:) 
      end
      end interface
      
      interface
      subroutine symplectic(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt )
     

      implicit none     
            
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt,time
      integer i, j, k, itimestep, d, current_ts, nstart        
           
      double precision,allocatable:: x_min(:, :), v_min(:, :), 
     &    u_min(:), rho_min(:), dx(:,:), dvx(:, :), du(:),  
     &       drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:),
     &       xn(:, :), vxn(:, :), 
     &       rhon(:), pn(:), un(:), 
     &         xn1(:, :), vxn1(:, :), 
     &       rhon1(:), pn1(:), un1(:)
             
      	double precision dt2, epsilon_rdot
      end
      end interface
      
      interface
      subroutine verlet(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt )
     

      implicit none     
    
      
      integer itype(:), ntotal, maxtimestep
      double precision x(:, :), vx(:, :), mass(:), 
     &       rho(:), p(:), u(:), c(:), s(:), e(:), 
     &       hsml(:), dt
      integer i, j, k, itimestep, d, current_ts, nstart        
            
      double precision  time, temp_rho, temp_u
         
	  
	   
      double precision,allocatable:: x_min(:, :), v_min(:, :),
     &        u_min(:),rho_min(:), dx(:,:), dvx(:, :), 
     &        du(:),drho(:),  av(:, :), ds(:), 
     &       t(:), tdsdt(:) , 
     &       xn(:, :), vxn(:, :), 
     &       rhon(:), pn(:), un(:), rhon1(:) ,un1(:),pn1(:),
     &       vxn1(:,:),xn1(:,:) 
      end
      end interface
      
      !interface
        !double precision function tensile(wk,p_i,p_j,rho_i,rho_j)
          !implicit none     
          !double precision  wk,p_i,p_j,rho_i,rho_j
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !double precision wdp,od_wdp,fab,Ra,Rb,R
         !double precision,allocatable:: dwdxp(:),dx(:)      
        !end
      !end interface
      
      end module interf
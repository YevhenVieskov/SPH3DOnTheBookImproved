	subroutine kernel(r,dx,hsml,w,dwdx)   

c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing kernel wij and its 
c   derivatives dwdxij.
c     if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c            = 2, Gauss kernel   (Gingold and Monaghan 1981) 
c            = 3, Quintic kernel (Morris 1997)

c     r    : Distance between particles i and j                     [in]
c     dx   : x-, y- and z-distance between i and j                  [in]  
c     hsml : Smoothing length                                       [in]
c     w    : Kernel for all interaction pairs                      [out]
c     dwdx : Derivative of kernel with respect to x, y and z       [out]
      use interf
      implicit none
      include 'param.inc'
      
      double precision r, dx(:), hsml, w, dwdx(:)
      integer i, j, d      
      double precision q, dw, factor
      double precision xs,AK,BK,CK,DK,EK,FK, AN2D,AN3D,alpha,beta
      double precision PK,QK,gamma
      integer nk
	  
	  
      q = r/hsml 
      w = 0.e0
      do d=1,dim         
        dwdx(d) = 0.e0
      enddo   

      if (skf.eq.1) then     
      
        if (dim.eq.1) then
          factor = 1.e0/hsml
        elseif (dim.eq.2) then
          factor = 15.e0/(7.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 3.e0/(2.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif                                           
        if (q.ge.0.and.q.le.1.e0) then          
          w = factor * (2./3. - q*q + q**3 / 2.)
          do d = 1, dim
            dwdx(d) = factor * (-2.+3./2.*q)/hsml**2 * dx(d)       
          enddo   
        else if (q.gt.1.e0.and.q.le.2) then          
          w = factor * 1.e0/6.e0 * (2.-q)**3 
          do d = 1, dim
            dwdx(d) =-factor * 1.e0/6.e0 * 3.*(2.-q)**2/hsml * (dx(d)/r)        
          enddo              
	else
	  w=0.
          do d= 1, dim
            dwdx(d) = 0.
          enddo             
        endif     
                                    
      else if (skf.eq.2) then
      
        factor = 1.e0 / (hsml**dim * pi**(dim/2.))      
	if(q.ge.0.and.q.le.3) then
	  w = factor * exp(-q*q)
          do d = 1, dim
            dwdx(d) = w * ( -2.* dx(d)/hsml/hsml)
          enddo 
	else
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo 	   
	endif	       
	
      else if (skf.eq.3) then	
      
        if (dim.eq.1) then
          factor = 1.e0 / (120.e0*hsml)
        elseif (dim.eq.2) then
          factor = 7.e0 / (478.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif              
	if(q.ge.0.and.q.le.1) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * ( (-120 + 120*q - 50*q**2) 
     &                        / hsml**2 * dx(d) )
          enddo 
	else if(q.gt.1.and.q.le.2) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4 + 30*(2-q)**4)  
     &                       / hsml * (dx(d)/r) 
          enddo 
        else if(q.gt.2.and.q.le.3) then
          w = factor * (3-q)**5 
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4) / hsml * (dx(d)/r) 
          enddo 
        else   
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo  
        endif                      
                
      else if (skf.eq.5) then
        xs=0.2
        alpha=1.0/(xs*xs*xs-3.0*xs*xs+3.0*xs-1.0)
        AK=alpha/2.0
        BK=-alpha*(1+xs)
        CK=3.0*alpha*xs
        DK=-alpha*(-1+3.0*xs)
        EK=alpha*(2.0*xs-1.0)/2.0
        FK=AK*xs*xs*xs*xs+BK*xs*xs*xs+CK*xs*xs+DK*xs+EK+xs
        AN2D=1.0/2.0*FK*xs*xs-1.0/3.0*xs*xs*xs
        AN2D=AN2D+AK/6.0+BK/5.0+CK/4.0+DK/3.0+EK/2.0
        AN2D=AN2D-AK/6.0*xs*xs*xs*xs*xs*xs-BK/5.0*xs*xs*xs*xs*xs
        AN2D=AN2D-CK/4.0*xs*xs*xs*xs-DK/3.0*xs*xs*xs-EK/2.0*xs*xs
        AN2D=2.0*pi*AN2D 
        AN2D=1.0/AN2D      
        AN3D=1.0/3.0*FK*xs*xs*xs-1.0/4.0*xs*xs*xs*xs
        AN3D=AN3D+AK/7.0+BK/6.0+CK/5.0+DK/4.0+EK/3.0
        AN3D=AN3D-AK/7.0*xs*xs*xs*xs*xs*xs*xs
        AN3D=AN3D-BK/6.0*xs*xs*xs*xs*xs*xs
        AN3D=AN3D-CK/5.0*xs*xs*xs*xs*xs
        AN3D=AN3D-DK/4.0*xs*xs*xs*xs-EK/3.0*xs*xs*xs
        AN3D=4.0*pi*AN3D
          AN3D=1.0/AN3D
           
        !AN2D=1.798
        !AN3D=2.168
         
        if(dim.eq.2)factor=AN2D/hsml/hsml
        if(dim.eq.3)factor=AN3D/hsml/hsml/hsml
        
        if(0.0.ge.q.and.q.lt.xs) then
          w =factor*(FK-q)
        else if(xs.ge.q.and.q.lt.1)then 
           w=factor*(AK*q*q*q*q+BK*q*q*q+CK*q*q+DK*q+EK)
        else
          w=0.0
        end if
        
        if(0.0.ge.q.and.q.lt.xs) then
          do d=1,dim
            dwdx(d)=factor*(-1.0)*dx(d)/hsml/hsml/q
          end do
        else if(xs.ge.q.and.q.lt.1)then 
          do d=1,dim
            dwdx(d)=factor*(4.0*AK*q*q*q+3.0*BK*q*q+2.0*CK*q+DK)
     &              *dx(d)/hsml/hsml/q    
          end do
        else
          do d=1,dim
            dwdx(d)=0.0
          end do
        end if
        
      else if (skf.eq.6) then
      alpha=1.0/3.0
      beta=1.0+6.0*alpha*alpha-12.0*alpha*alpha*alpha
      AN2D=8.0/(pi*(6.4*alpha**5-16.0*alpha**6+1))
      factor=AN2D/(hsml*hsml*hsml)
      if(0.0.ge.q.and.q.lt.xs) then
         w=factor*((-12.0*alpha+18.0*alpha*alpha)*q+beta)
      else if (q.gt.alpha.and.q.le.0.5) then 
         w=factor*(1.0-6.0*q*q+6*q*q*q)
      else if (q.gt.0.5.and.q.le.1.0) then
        w=factor*2*(1.0-q*q*q)**3 
      else
        w=0.0
      end if
      if(q.eq.0.0) then
        do d=1,dim
          dwdx(d)=0.0
        end do
      else if(0.0.gt.q.and.q.lt.xs) then
         do d=1,dim
            dwdx(d)=factor*(-12.0*alpha+18.0*alpha*alpha)
     &              *dx(d)/hsml/hsml/q            
         end do
      else if (q.gt.alpha.and.q.le.0.5) then
         do d=1,dim
            dwdx(d)=factor*(-12.0*q+18.0*q*q)*dx(d)/hsml/hsml/q            
         end do 
      else if (q.gt.0.5.and.q.le.1.0) then
        do d=1,dim
            dwdx(d)=factor*(-6.0)*(1.0-q)**2*dx(d)/hsml/hsml/q            
         end do
      else
        do d=1,dim
            dwdx(d)=0.0           
         end do
      end if    
      else if (skf.eq.7) then
        nk=8	                        
         if(nk.eq.3) then
           AK=2.4
           BK=-9.4
           PK=-1.81
           QK=1.028
           alpha=0.317
           AN2D=3.71
         else if(nk.eq.4) then
           AK=3.2
           BK=-18.8
           PK=-2.15
           QK=0.98
           alpha=0.214
           AN2D=6.52
         else if(nk.eq.5) then
           AK=4.27
           BK=-37.6
           PK=-2.56
           QK=0.962
           alpha=0.161
           AN2D=10.4 
         else if(nk.eq.8) then 
           AK=10.1
           BK=-300.8
           PK=-3.86
           QK=0.942
           alpha=0.0927
           AN2D=30.75
         end if 
         factor=AN2D/hsml/hsml/hsml
         if(q.gt.0.and.q.le.alpha)then
           w=factor*(PK*q+QK)
         else if(q.gt.alpha.and.q.le.beta)then
           w=factor*((1.0-q)**nk+AK*(gamma-q)**nk+BK*(beta-q)**nk)
         else if(q.gt.beta.and.q.le.gamma)then
           w=factor*((1.0-q)**nk+AK*(gamma-q)**nk)
         else if(q.gt.gamma.and.q.le.1.0)then
           w=factor*(1.0-q)**nk
         else
           w=0
         end if 
         
         if(q.gt.0.and.q.le.alpha)then
           do d=1,dim
            dwdx(d)=factor*PK*dx(d)/hsml/hsml/q
           end do
         else if(q.gt.alpha.and.q.le.beta)then
           do d=1,dim
           
            dwdx(d)=factor*(-nk*(1-q)**(nk-1)-BK*nk*(beta - q)**(nk-1)
     &              -AK*nk*(gamma - q)**(nk-1))*dx(d)/hsml/hsml/q
           end do
         else if(q.gt.beta.and.q.le.gamma)then
           do d=1,dim
            dwdx(d)=factor*(-nk*(1-q)**(nk-1)-AK*nk*(gamma - q)**(nk-1))
     &                *dx(d)/hsml/hsml/q
           end do
         else if(q.gt.gamma.and.q.le.1.0)then
           do d=1,dim
            dwdx(d)=factor*(-nk*(1-q)**(nk-1))*dx(d)/hsml/hsml/q
           end do
         else
           do d=1,dim
            dwdx(d)=0.0
           end do
         end if   
      endif  
		
      end

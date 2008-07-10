!-----[--.----+----.----+----.-----------------------------------------]
      subroutine therm2d(alpha,beta,ul,xl,s,nel,nst,ndm,ndf)
!-----[--.----+----.----+----.-----------------------------------------]
      
!     Two dimensional (plane/axisymmetric) Linear Thermal Element
!             alpha del^2 u + beta u
!     Input:
!     alpha         - stiffness or conductivity 
!     beta          - mass or capacitance
!     ul(ndf,nel,*) - Element nodal solution parameters
!     xl(ndm,nel)   - Block nodal coordinate array
!     nel           - number of nodes on elem
!     nst           - size of element arrays
!     ndm           - 2: -2D
!     ndf           - 1: scalar problem
!
!     Output:
!     s(nst,nst)    - stiffness matrix
!  
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      
      logical   quad
      integer   ndf,ndm,nst, i,j, i1,j1, l,lint, tdof,nel
      real*8    xsj, a1,a2,a3,a4, hh,tinf
      real*8    ul(ndf,nel,1),xl(ndm,nel),s(nst,nst),alpha,beta
      real*8    gradt(2),dd(2,2),shp(3,9),el(4,7)
      
      save

!     Compute conductivity (stiffness) matrix
         
      if(nel.eq.3 .or. nel.eq.6 .or. nel.eq.7) then
         if( nel.eq.3 ) then
            l = 1
         else 
            l = 7
         endif
         quad = .false.
         call tint2d(l,lint,el)
      else
         stop 'nel='
!     l    =  d(5)
!     quad = .true.
!     if(l*l.ne.lint) call int2d(l,lint,sg)
      endif

!     Get global dof for thermal variable
      tdof = 1                  !  max(1,nint(d(19)))
      hh   = 0.d0                 ! d(127)
      tinf = 0.d0                 ! d(128)
      
      do l = 1,lint
         
         if(quad) then
            stop 'nel='
!     call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
!     xsj = xsj*sg(3,l)*d(14)
         else
            if( nel.eq.3 ) then
               call trishp(el(1,l),xl,ndm,1,xsj,shp)
            else 
               call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
            endif
            xsj = xsj*el(4,l)   ! *d(14) shell thickness
         endif
         
!     Compute flux
         
         call thfx2d( alpha, ul, shp, gradt, dd, ndf, nel )
         
         j1 = 1
         do j = 1,nel
            a1 = (dd(1,1)*shp(1,j) + dd(1,2)*shp(2,j))*xsj
            a2 = (dd(2,1)*shp(1,j) + dd(2,2)*shp(2,j))*xsj
            a3 = beta*shp(3,j)*xsj   ! d(4)*d(64)*  -- d(4) is 'densiity'
            a4 = 0.            !d(127)*shp(3,j)*xsj  -- d(127) is surface convection
            
!     Compute residual
            
!     p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2)
!     &              - a3*tdot
!     &              - a4*temp ! (temp - d(128))
            
!     Compute tangent
            
!            a1 = a1*ctan(1)
!            a2 = a2*ctan(1)
!            a3 = a3*ctan(3) ! 
            a4 = a3 ! a4*alpha + a3*cfac
            
!     Lumped rate terms
            
!     s(j1,j1) = s(j1,j1) + a3*lfac
            
!     Consistent rate and conductivity terms
            
            i1 = 1
            do i = 1,nel
               s(i1,j1) = s(i1,j1) + a1*shp(1,i) + a2*shp(2,i) + a4*shp(3,i)
               i1 = i1 + ndf
            end do
            j1 = j1 + ndf
         end do
      end do
      
      end
      
!-----[--.----+----.----+----.-----------------------------------------]      
      subroutine thfx2d( alpha, ul, shp, gradt, dd, ndf, nel )
!-----[--.----+----.----+----.-----------------------------------------]      
!     Compute thermal gradient
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
 
      integer   ndf,nel, i
      real*8    ul(ndf,*), shp(3,*), alpha
      real*8    psi,cs,sn,c2,s2,gradt(2),dd(2,2)
      
      save
      
!      gradt(1) = 0.0d0
!      gradt(2) = 0.0d0
!      do i = 1,nel
!         gradt(1) = gradt(1) + shp(1,i)*ul(1,i)
!         gradt(2) = gradt(2) + shp(2,i)*ul(1,i)
!      end do
      
!     Compute thermal flux
      
!      psi = 0. ! d(31)
!      cs  = 1. ! cos(psi)
!      sn  = 0. ! sin(psi)
!      c2  = cs*cs
!      s2  = sn*sn
!      cs  = cs*sn
      
      dd(1,1) = alpha
      dd(2,2) = alpha
      dd(1,2) = 0.d0 ! cs*(alpha - alpha)
      dd(2,1) = 0.d0 ! dd(1,2)

      end

!-----[--.----+----.----+----.-----------------------------------------]
      subroutine trishp(el,xl,ndm,iord, xsj,shp)
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Triangular shape function routine

!     Type:  |iord| = 1:  Linear  three-node
!     |iord| = 2:  Quadratic six-node
!     |iord| = 3:  Quadratic seven-node
!     |iord| = 4:  Quadratic + 3 bubbles (Zienkiewicz/Lefebre)
      
!     iord  > 0:  Mid-side and center node are  global  coords.
!     iord  < 0:  Mid-side and center node heirarchical coords.
      
!     Inputs:
!     el(3)     - Area coordinates for point
!     xl(ndm,*) - Nodal coordinates for element
!     ndm       - Spatial dimension of mesh
!     iord      - Order of shape functions (see above)
      
!     Outputs:
!     xsj       - Jacobian determinant at point
!     shp(3,*)  - Shape functions and derivatives at point
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      
      integer   ndm,iord, i
      real*8    one3,four9,xsj,xsjr, fel1,fel2,fel3, fel12,fel23,fel31
      real*8    x1,x2,x3,x4,x5,x6,x7, y1,y2,y3,y4,y5,y6,y7
      real*8    el(3),xl(ndm,*), shp(3,*)
      
      data one3 /0.3333333333333333d0/, four9 /0.4444444444444444d0/
      
!     Form Jacobian terms
      
      x1 = xl(1,1)
      x2 = xl(1,2)
      x3 = xl(1,3)
      
      y1 = xl(2,1)
      y2 = xl(2,2)
      y3 = xl(2,3)
      
      if(abs(iord).gt.1) then
         
         fel1 = 4.d0*el(1)
         fel2 = 4.d0*el(2)
         fel3 = 4.d0*el(3)
         
         x4   = xl(1,4)
         x5   = xl(1,5)
         x6   = xl(1,6)
         
         y4   = xl(2,4)
         y5   = xl(2,5)
         y6   = xl(2,6)
         
!     Form shape functions in total coordinates
         
         if(iord.gt.0) then
            
            x4 = x4 - 0.5d0*(x1 + x2)
            x5 = x5 - 0.5d0*(x2 + x3)
            x6 = x6 - 0.5d0*(x3 + x1)
            x7 = one3*(x1 + x2 + x3) + four9*(x4 + x5 + x6)
            
            y4 = y4 - 0.5d0*(y1 + y2)
            y5 = y5 - 0.5d0*(y2 + y3)
            y6 = y6 - 0.5d0*(y3 + y1)
            y7 = one3*(y1 + y2 + y3) + four9*(y4 + y5 + y6)
            
         endif

         x1   = x1 + x4*fel2 + x6*fel3
         x2   = x2 + x5*fel3 + x4*fel1
         x3   = x3 + x6*fel1 + x5*fel2
         
         y1   = y1 + y4*fel2 + y6*fel3
         y2   = y2 + y5*fel3 + y4*fel1
         y3   = y3 + y6*fel1 + y5*fel2
         
         fel12 = el(1)*el(2)
         fel23 = el(2)*el(3)
         fel31 = el(3)*el(1)
         
         if(iord.eq.3) then
            
            fel12 = 27.d0*fel12
            fel23 = 27.d0*fel23
            fel31 = 27.d0*fel31
            
            x7    = xl(1,7) - x7
            x1    = x1 + x7*fel23
            x2    = x2 + x7*fel31
            x3    = x3 + x7*fel12
            
            y7    = xl(2,7) - y7
            y1    = y1 + y7*fel23
            y2    = y2 + y7*fel31
            y3    = y3 + y7*fel12
            
         endif
      endif
      
      xsj  = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
      xsjr = 1.0d0
      if(xsj.ne.0.0d0) then
         xsjr = 1.0d0/xsj
      else 
         stop 'xsj==0'
      endif
      xsj  = 0.5d0*xsj
      
!     Specify shape functions and their derivatives
      
      shp(1,1) = (y2-y3)*xsjr
      shp(2,1) = (x3-x2)*xsjr
      shp(3,1) = el(1)
      
      shp(1,2) = (y3-y1)*xsjr
      shp(2,2) = (x1-x3)*xsjr
      shp(3,2) = el(2)
      
      shp(1,3) = (y1-y2)*xsjr
      shp(2,3) = (x2-x1)*xsjr
      shp(3,3) = el(3)

!     Add quadratic hierarchical functions
      
      if(abs(iord).gt.1) then
         
         shp(1,4) = shp(1,1)*fel2 + shp(1,2)*fel1
         shp(2,4) = shp(2,1)*fel2 + shp(2,2)*fel1
         shp(3,4) = el(1)*fel2
         
         shp(1,5) = shp(1,2)*fel3 + shp(1,3)*fel2
         shp(2,5) = shp(2,2)*fel3 + shp(2,3)*fel2
         shp(3,5) = el(2)*fel3
         
         shp(1,6) = shp(1,3)*fel1 + shp(1,1)*fel3
         shp(2,6) = shp(2,3)*fel1 + shp(2,1)*fel3
         shp(3,6) = el(3)*fel1
         
!     Bubble at baricenter
         
         if(abs(iord).eq.3) then
            shp(1,7) = shp(1,1)*fel23 + shp(1,2)*fel31 + shp(1,3)*fel12
            shp(2,7) = shp(2,1)*fel23 + shp(2,2)*fel31 + shp(2,3)*fel12
            shp(3,7) = fel12*el(3)
         elseif( abs(iord).eq.4 ) then
            shp(1,7) = (shp(1,1)*fel23*2.d0 + shp(1,2)*fel31 + shp(1,3)*fel12)*el(1)
            shp(2,7) = (shp(2,1)*fel23*2.d0 + shp(2,2)*fel31 + shp(2,3)*fel12)*el(1)
            shp(3,7) = fel12*fel31
            
            shp(1,8) = (shp(1,1)*fel23 + shp(1,2)*fel31*2.d0 + shp(1,3)*fel12)*el(2)
            shp(2,8) = (shp(2,1)*fel23 + shp(2,2)*fel31*2.d0 + shp(2,3)*fel12)*el(2)
            shp(3,8) = fel12*fel23
            
            shp(1,9) = (shp(1,1)*fel23 + shp(1,2)*fel31 + shp(1,3)*fel12*2.d0)*el(3)
            shp(2,9) = (shp(2,1)*fel23 + shp(2,2)*fel31 + shp(2,3)*fel12*2.d0)*el(3)
            shp(3,9) = fel31*fel23
         endif
         
!     Modify shape functions for mid-side and interior values
         
         if(iord.gt.1) then
            
!     Modify vertex and mid-side values for bubble
            
            if(iord.eq.3) then
               do i = 1,3
                  shp(i,1) = shp(i,1) -  one3*shp(i,7)
                  shp(i,2) = shp(i,2) -  one3*shp(i,7)
                  shp(i,3) = shp(i,3) -  one3*shp(i,7)
                  shp(i,4) = shp(i,4) - four9*shp(i,7)
                  shp(i,5) = shp(i,5) - four9*shp(i,7)
                  shp(i,6) = shp(i,6) - four9*shp(i,7)
               end do
            endif
            
!     Modify vertex shape functions for mid-side values
            
            do i = 1,3
               shp(i,1) = shp(i,1) - 0.5d0*(shp(i,4) + shp(i,6))
               shp(i,2) = shp(i,2) - 0.5d0*(shp(i,5) + shp(i,4))
               shp(i,3) = shp(i,3) - 0.5d0*(shp(i,6) + shp(i,5))
            end do
            
         endif
      endif
      
      end
      
!-----[--.----+----.----+----.-----------------------------------------]
      subroutine tint2d(l,lint,el)
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose: Set gauss points and weights for triangular elements
      
!     Inputs:
!     l       - Number of gauss points indicator
      
!     Outputs:
!     lint    - Total number of points
!     el(4,*) - Area coordinate points and weights for quadrature
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
            
      integer   l, lint
      real*8    el(4,*), r0,r1,r2, ww, eta
      
      data ww, eta / 0.3333333333333333d0 , 0.1666666666666667d0 /
      
!     1-point gauss integration
      
      if(l.eq.1) then
         el(1,1) = ww
         el(2,1) = ww
         el(3,1) = ww
         el(4,1) = 1.d0
         lint    = 1
         
!     3-point integration: mid-edge points
         
      elseif(l.eq.3) then
         el(1,1) = 0.d0
         el(2,1) = 0.5d0
         el(3,1) = 0.5d0
         el(4,1) = ww
         
         el(1,2) = 0.5d0
         el(2,2) = 0.d0
         el(3,2) = 0.5d0
         el(4,2) = ww
         
         el(1,3) = 0.5d0
         el(2,3) = 0.5d0
         el(3,3) = 0.d0
         el(4,3) = ww
         
         lint    = 3
         
!     3-point integration: interior points
         
      elseif(l.eq.-3) then
         
         el(1,1) = 1.0d0 - ww
         el(2,1) = eta
         el(3,1) = eta
         el(4,1) = ww
         
         el(1,2) = eta
         el(2,2) = 1.0d0 - ww
         el(3,2) = eta
         el(4,2) = ww
         
         el(1,3) = eta
         el(2,3) = eta
         el(3,3) = 1.0d0 - ww
         el(4,3) = ww
         
         lint    = 3
         
!     6-point nodal integration
         
      elseif(l.eq.6) then
         
         el(1,1) =  1.0d0
         el(2,1) =  0.0d0
         el(3,1) =  0.0d0
         el(4,1) =  eta
         
         el(1,2) =  0.0d0
         el(2,2) =  1.0d0
         el(3,2) =  0.0d0
         el(4,2) =  eta
         
         el(1,3) =  0.0d0
         el(2,3) =  0.0d0
         el(3,3) =  1.0d0
         el(4,3) =  eta
         
         el(1,4) =  0.5d0
         el(2,4) =  0.5d0
         el(3,4) =  0.0d0
         el(4,4) =  eta
         
         el(1,5) =  0.0d0
         el(2,5) =  0.5d0
         el(3,5) =  0.5d0
         el(4,5) =  eta
         
         el(1,6) =  0.5d0
         el(2,6) =  0.0d0
         el(3,6) =  0.5d0
         el(4,6) =  eta
         
         lint    =  6
         
!     6-point order 4 formula
         
      elseif(l.eq.-6) then
         
         el(1,1) = 0.816847572980459d0
         el(2,1) = 0.091576213509771d0
         el(3,1) = 0.091576213509771d0
         el(4,1) = 0.109951743655322d0
         
         el(1,2) = 0.091576213509771d0
         el(2,2) = 0.816847572980459d0
         el(3,2) = 0.091576213509771d0
         el(4,2) = 0.109951743655322d0
         
         el(2,3) = 0.091576213509771d0
         el(1,3) = 0.091576213509771d0
         el(3,3) = 0.816847572980459d0
         el(4,3) = 0.109951743655322d0
         
         el(1,4) = 0.108103018168070d0
         el(2,4) = 0.445948490915965d0
         el(3,4) = 0.445948490915965d0
         el(4,4) = 0.223381589678011d0
         
         el(1,5) = 0.445948490915965d0
         el(2,5) = 0.108103018168070d0
         el(3,5) = 0.445948490915965d0
         el(4,5) = 0.223381589678011d0
         
         el(1,6) = 0.445948490915965d0
         el(2,6) = 0.445948490915965d0
         el(3,6) = 0.108103018168070d0
         el(4,6) = 0.223381589678011d0
         
         lint    = 6
         
!     7-point gauss integration
         
      elseif(l.eq.7) then
         r0      =  sqrt(15.0d0)
         r1      =  3.d0/7.d0
         r2      =  (r0 + r0)/21.d0
         
         el(1,1) =  ww
         el(2,1) =  el(1,1)
         el(3,1) =  el(1,1)
         el(4,1) =  0.225d0
         
         el(1,2) =  r1 + r2
         el(2,2) =  0.5d0 - 0.5d0*el(1,2)
         el(3,2) =  el(2,2)
         el(4,2) =  (155.d0 - r0)/1200.d0
         
         el(1,3) =  el(2,2)
         el(2,3) =  el(1,2)
         el(3,3) =  el(2,2)
         el(4,3) =  el(4,2)
         
         el(1,4) =  el(2,2)
         el(2,4) =  el(2,2)
         el(3,4) =  el(1,2)
         el(4,4) =  el(4,2)
         
         el(1,5) =  r1 - r2
         el(2,5) =  0.5d0 - 0.5d0*el(1,5)
         el(3,5) =  el(2,5)
         el(4,5) =  (155.d0 + r0)/1200.d0
         
         el(1,6) =  el(2,5)
         el(2,6) =  el(1,5)
         el(3,6) =  el(2,5)
         el(4,6) =  el(4,5)
         
         el(1,7) =  el(2,5)
         el(2,7) =  el(2,5)
         el(3,7) =  el(1,5)
         el(4,7) =  el(4,5)
         
         lint    =  7
         
!     12-point order 6 formula
         
      elseif(l.eq.12) then
         
         el(1, 1) = 0.873821971016996d0
         el(2, 1) = 0.063089014491502d0
         el(3, 1) = 0.063089014491502d0
         el(4, 1) = 0.050844906370207d0
         
         el(1, 2) = 0.063089014491502d0
         el(2, 2) = 0.873821971016996d0
         el(3, 2) = 0.063089014491502d0
         el(4, 2) = 0.050844906370207d0
         
         el(1, 3) = 0.063089014491502d0
         el(2, 3) = 0.063089014491502d0
         el(3, 3) = 0.873821971016996d0
         el(4, 3) = 0.050844906370207d0
         
         el(1, 4) = 0.501426509658179d0
         el(2, 4) = 0.249286745170910d0
         el(3, 4) = 0.249286745170910d0
         el(4, 4) = 0.116786275726379d0
         
         el(1, 5) = 0.249286745170910d0
         el(2, 5) = 0.501426509658179d0
         el(3, 5) = 0.249286745170910d0
         el(4, 5) = 0.116786275726379d0
         
         el(1, 6) = 0.249286745170910d0
         el(2, 6) = 0.249286745170910d0
         el(3, 6) = 0.501426509658179d0
         el(4, 6) = 0.116786275726379d0
         
         el(1, 7) = 0.636502499121399d0
         el(2, 7) = 0.310352451033785d0
         el(3, 7) = 0.053145049844816d0
         el(4, 7) = 0.082851075618374d0
         
         el(1, 8) = 0.636502499121399d0
         el(2, 8) = 0.053145049844816d0
         el(3, 8) = 0.310352451033785d0
         el(4, 8) = 0.082851075618374d0
         
         el(1, 9) = 0.310352451033785d0
         el(2, 9) = 0.636502499121399d0
         el(3, 9) = 0.053145049844816d0
         el(4, 9) = 0.082851075618374d0
         
         el(1,10) = 0.053145049844816d0
         el(2,10) = 0.636502499121399d0
         el(3,10) = 0.310352451033785d0
         el(4,10) = 0.082851075618374d0
         
         el(1,11) = 0.310352451033785d0
         el(2,11) = 0.053145049844816d0
         el(3,11) = 0.636502499121399d0
         el(4,11) = 0.082851075618374d0
         
         el(1,12) = 0.053145049844816d0
         el(2,12) = 0.310352451033785d0
         el(3,12) = 0.636502499121399d0
         el(4,12) = 0.082851075618374d0
         
         lint     = 12

      else
        write(  *,2000) l
        lint    = -1
      endif
      
!     Format

 2000 format(' *ERROR* TINT2D: Wrong quadrature, l =',i3)
      
      end

!-----[--.----+----.----+----.-----------------------------------------]
      subroutine cktris(xl,shp,ndm,nel)
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Check isoparametric triangles elements for valid data

!      Inputs:
!         xl(ndm,*) - Nodal coordinates for element
!         ndm       - Spatial dimension of problem

!      Outputs:
!         None

!      Scratch:
!         shp(*)    - Storage for shape functions
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ineg, ndm, i,l, ic(2,7), nel
      real*8    xsj, shp(*),xl(ndm,*),el(3,7)

      save

      data el/0.0d0,0.0d0,1.0d0, 1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, &
              0.5d0,0.0d0,0.5d0, 0.5d0,0.5d0,0.0d0, 0.0d0,0.5d0,0.5d0, &
              0.3333333333333d0, 0.3333333333333d0, 0.3333333333333d0/

!     Check element for input errors

      ineg = 0
      do l = 1,nel
         call trishp(el(1,l),xl,ndm,1,xsj,shp)
         if(xsj.le.0.0d0) then
            ineg       = ineg + 1
            ic(1,ineg) = l
            ic(2,ineg) = 1 ! abs(ix(l))
         endif
      end do
      
      if(ineg.gt.0) then
         write(*,2001) -1,(ic(1,i),ic(2,i),i=1,ineg)
      endif
      
!     Formats
 
2001  format(' >Element',i8,' has negative jacobian at nodes:',('                    Local =',i3,' Global =',i4))
 
      end
      
!-----[--.----+----.----+----.-----------------------------------------]
!      program fe
	subroutine fe
!-----[--.----+----.----+----.-----------------------------------------]
!     Driver for therm2d
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      
      integer nel,i,j,k,ix(3,8),nst,e,ni,nj
      real*8 xl(2,9),xlt(2,3),s(3,4),alpha,beta
      real*8 arr(9,10),ul(3,4)

      data ix/1,4,2,2,4,5,3,2,5,3,5,6,5,4,7,5,7,8,6,5,8,6,8,9/
      data xl/0.d0,2.d0,1.d0,2.d0,2.d0,2.d0,0.d0,1.d0,1.d0,1.d0,2.d0,1.d0,0.d0,0.d0,1.d0,0.d0,2.d0,0.d0/

      nst = 3
      nel = 3
      alpha = 1.d0
      beta = 0.d0
      !      xl = xl * 2.d0
      !
      do j=1,10
         do i=1,9
            arr(i,j) = 0.d0
         enddo
      enddo
      !     
      do e=1,8
         !     
         do j=1,4
            do k=1,3
               s(k,j) = 0.d0
            enddo
         enddo

         do j=1,nel
            ul(j,1) = 0.d0        ! displ: u
            ul(j,4) = 0.d0        ! rate: \dot u
            xlt(1,j) = xl( 1, ix(j,e) )
            xlt(2,j) = xl( 2, ix(j,e) )
         enddo
         !     
         call therm2d(alpha,beta,ul,xlt,s,nel,nst,2,1)
         !     
         if( s(3,3) .eq. 0. .or. s(1,4) .ne. 0. ) then
            write (*,*) 'ERROR: s(3,3) =',s(3,3),'s(1,4)=',s(1,4) 
         endif
         !     
         do i=1,3
            ni = ix(i,e)
            do j=1,3
               nj = ix(j,e)
               arr(ni,nj) = arr(ni,nj) + s(i,j)
            enddo               ! j
         enddo                  ! i
         
      enddo                     ! e
      
      if( arr(9,9) .eq. 0. .or. arr(1,10) .ne. 0. ) then
         write (*,*) 'ERR: a(9,9) =',arr(9,9),'a(1,10)=',arr(1,10) 
      endif
      
      write (*,*) arr

      end



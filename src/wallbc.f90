!=============================================================================!
        subroutine wallbc(vl, rhow, pnorm) 
!=============================================================================!
        use global
        use material
        use stencil
        implicit none

        real :: vl(ny,nx,ndof), rhow(nx), Pnorm(nx)
        integer :: iz
        
!.... local wall variables

        integer :: i, j, ij, idof, jdof
        
        real :: g1v(nx,ndof), g2v(nx,ndof)
        real :: g11v(nx,ndof), g12v(nx,ndof), g22v(nx,ndof)
        real :: g1vl, g2vl, g11vl, g12vl, g22vl

        real :: g1divu(nx), g2divu(nx)
        real :: S1jj(nx), S2jj(nx)
        real :: divu(nx)
        real :: S(nx,2,2)
        real :: mu(nx),     lm(nx),     con(nx)
        real :: dmu(nx),    d2mu(nx)
        real :: dlm(nx),    d2lm(nx)
        real :: dcon(nx),   d2con(nx)
        real :: g1mu(nx),   g2mu(nx)
        real :: g1lm(nx),   g2lm(nx)
        real :: g1con(nx),  g2con(nx)
        real :: g1dmu(nx),  g2dmu(nx)
        real :: g1dlm(nx),  g2dlm(nx)
        real :: g1dcon(nx), g2dcon(nx)

!.... stuff for the wall normal pressure gradient
        
        real :: bn1(nx), bn2(nx)
!=============================================================================!

!.... Compute first derivatives of field in the mapped space

        call wgrad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                   xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... Compute second derivatives of field
        
        call wgrad2(ndof, nx, ny, vl, g2v, g11v, g12v, g22v, dxi, deta, &
                    optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall BC on the field

        if (wallt.eq.1 .or. wallt.eq.2) g2v(:,ndof) = zero

#if 0

!.... transform the gradients to physical space

        j = 1
        do idof = 1, ndof
          do i = 1, nx
            ij = j + (i-1)*ny
            g1vl  = g1v(i,idof)*m1(ij) + g2v(i,idof)*n1(ij)
            g2vl  = g1v(i,idof)*m2(ij) + g2v(i,idof)*n2(ij)
  
            g11vl = g11v(i,idof)       * m1m1(ij)       + &
                    two * g12v(i,idof) * m1n1(ij)       + &
                    g22v(i,idof)       * n1n1(ij)       + &
                    g1v(i,idof)        * m11(ij)        + &
                    g2v(i,idof)        * n11(ij)
  
            g12vl = g11v(i,idof)       * m1m2(ij)       + &
                    g12v(i,idof)       * m1n2(ij)       + &
                    g12v(i,idof)       * m2n1(ij)       + &
                    g22v(i,idof)       * n1n2(ij)       + &
                    g1v(i,idof)        * m12(ij)        + &
                    g2v(i,idof)        * n12(ij)
  
            g22vl = g11v(i,idof)       * m2m2(ij)       + &
                    two * g12v(i,idof) * m2n2(ij)       + &
                    g22v(i,idof)       * n2n2(ij)       + &
                    g1v(i,idof)        * m22(ij)        + &
                    g2v(i,idof)        * n22(ij)
  
            g1v(i,idof)  = g1vl
            g2v(i,idof)  = g2vl
            g11v(i,idof) = g11vl
            g12v(i,idof) = g12vl
            g22v(i,idof) = g22vl
          end do
        end do

!=============================================================================!
!.... setup local variables
!=============================================================================!

!.... compute the gradient of the divergence of um

        g1divu = g11v(:,2) + g12v(:,3)
        g2divu = g12v(:,2) + g22v(:,3)

!.... compute strain rate tensor for the viscous terms

        S1jj = pt5 * ( g11v(:,2) + g11v(:,2) + g22v(:,2) + &
                       g12v(:,3) )
                
        S2jj = pt5 * ( g11v(:,3) + g12v(:,2) + g22v(:,3) + &
                       g22v(:,3) )
        
!.... compute divergence
  
        divu = g1v(:,2) + g2v(:,3)

!.... compute strain rate tensor for the viscous terms

        S(:,1,1) = g1v(:,2)
        S(:,1,2) = pt5 * ( g2v(:,2) + g1v(:,3) )
        S(:,2,1) = pt5 * ( g1v(:,3) + g2v(:,2) )
        S(:,2,2) = g2v(:,3)

!.... compute the material properties

        call getmat(vl(1,:,ndof),                       &
                    mu(:),   lm(:),    con(:),          &
                    dmu(:),  d2mu(:),  dlm(:),          &
                    d2lm(:), dcon(:),  d2con(:)         )

!.... compute gradients of viscosity using chain-rule

        g1mu = dmu * g1v(:,5)
        g2mu = dmu * g2v(:,5)

        g1dmu = d2mu * g1v(:,5)
        g2dmu = d2mu * g2v(:,5)

!.... compute gradients of conductivity using chain-rule

        g1con = dcon * g1v(:,5)
        g2con = dcon * g2v(:,5)

        g1dcon = d2con * g1v(:,5)
        g2dcon = d2con * g2v(:,5)

!.... compute gradients of bulk viscosity using chain-rule

        g1lm = dlm * g1v(:,5)
        g2lm = dlm * g2v(:,5)

        g1dlm = d2lm * g1v(:,5)
        g2dlm = d2lm * g2v(:,5)

!=============================================================================!

!.... compute the boundary normal unit-vector

        j = 1
        do i = 1, nx
          ij = j + (i-1) * ny
          bn1(i)  = n1(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
          bn2(i)  = n2(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
        end do    

!.... wall normal pressure gradient from the normal momentum equation
!.... Note that Pnorm is normalized by 1/sqrt(n1^2 + n2^2)

        do i = 1, nx
          ij = j + (i-1) * ny

          Pnorm(i) = one / Re * (                                           &
            divu(i) * ( g1lm(i) * bn1(i) + g2lm(i) * bn2(i) )             + &
            lm(i)   * ( g1divu(i) * bn1(i) + g2divu(i) * bn2(i) )         + &
            two * mu(i) * ( S1jj(i) * bn1(i) + S2jj(i) * bn2(i) )         + &
            two * ( g1mu(i) * S(i,1,1) + g2mu(i) * S(i,1,2) ) * bn1(i)    + &
            two * ( g1mu(i) * S(i,2,1) + g2mu(i) * S(i,2,2) ) * bn2(i) ) /  &
            sqrt( n1(ij)**2 + n2(ij)**2 )

!         Pnorm(i) = one / Re * (                                             &
!           damp(ij) * g1lm(i) * divu(i) * bn1(i)                           + &
!           g2lm(i) * ( damp(ij) * g1v(i,2) + g2v(i,3) ) * bn2(i)           + &
!           lm(i)   * ( damp(ij) * g1divu(i) *  bn1(i)                      + &
!                        g2divu(i) * bn2(i) )                               + &
!           mu(i) * ( ( damp(ij) * two * g11v(i,2) + g22v(i,2)              + &
!                       damp(ij) * g12v(i,3) ) * bn1(i)                     + &
!                     ( damp(ij) * g11v(i,3) + damp(ij) * g12v(i,2)         + &
!                       two * g22v(i,3) ) * bn2(i) )                        + &
!           ( g1mu(i) * damp(ij) * two * g1v(i,2)                           + &
!             g2mu(i) * ( g2v(i,2) + damp(ij) * g1v(i,3) ) ) * bn1(i)       + &
!           ( g1mu(i) * ( damp(ij) * g1v(i,3) + g2v(i,2) )                  + &
!             g2mu(i) * two * g2v(i,3) ) * bn2(i) )                         / &
!            sqrt( n1(ij)**2 + n2(ij)**2 )
        end do

!.... now compute the wall density using the normal momentum equation

!.... WARNING: This assumes that T,n = 0

        j = 1
        do i = 1, nx
          ij = j + (i-1) * ny
          if (carp) then
            rhow(i) = ( deta * gamma * Ma**2 * Pnorm(i) / vl(1,i,5) -   &
                        gg2 * vl(2,i,1) -                               &
                        gg3 * vl(3,i,1) -                               &
                        gg4 * vl(4,i,1) -                               &
                        gg5 * vl(5,i,1) -                               &
                        gg6 * vl(6,i,1) ) / gg1
          else
            rhow(i) = ( deta * gamma * Ma**2 * Pnorm(i) / vl(1,i,5) -   &
                        gc2 * vl(2,i,1) -                               &
                        gc3 * vl(3,i,1) -                               &
                        gc4 * vl(4,i,1) -                               &
                        gc5 * vl(5,i,1) ) / gc1
          end if
        end do
#else
        write(*,*) "Need to fix wall.f90"
#endif
        return
        end

!=============================================================================!
        subroutine wgrad( ndof, nx, ny, v, g1v, g2v, dx, dy, optx, opty, &
                          xper, yper, lsym, rsym, bsym, tsym, carp)
!
!  Take the gradient of a 2-D field at the wall.
!
!=============================================================================!
        use stencil
        implicit none
        
        integer :: ndof, nx, ny, optx, opty
        logical :: xper, yper
        logical :: lsym, rsym, bsym, tsym, carp
        real    :: v(ny,nx,ndof), g1v(nx,ndof), g2v(nx,ndof)
        real    :: dx, dy
        
        real,parameter :: one = 1.0
        real dxinv, dyinv
        real a, b, c, w
        real gx1, gx2, gx3, gx4, gx5, gx6
        real gy1, gy2, gy3, gy4, gy5, gy6

        integer :: idof
        real :: isign
!=============================================================================!

        dxinv  = one / dx
        dyinv  = one / dy

!.... seven point stencil in x

        if (optx.eq.0) then
          c = 1.0 / 60.0
        else if (optx.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gx1 =  -c * dxinv
        gx2 =  -b * dxinv
        gx3 =  -a * dxinv
        gx4 =   a * dxinv
        gx5 =   b * dxinv
        gx6 =   c * dxinv

!.... seven point stencil in y

        if (opty.eq.0) then
          c = 1.0 / 60.0
        else if (opty.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gy1 =  -c * dyinv
        gy2 =  -b * dyinv
        gy3 =  -a * dyinv
        gy4 =   a * dyinv
        gy5 =   b * dyinv
        gy6 =   c * dyinv

!=============================================================================!
!.... compute the gradient in x
!=============================================================================!

        if (xper) then

          g1v(1,:)          = ( gx1 * v(1,nx-3,:)       + &
                                gx2 * v(1,nx-2,:)       + &
                                gx3 * v(1,nx-1,:)       + &
                                gx4 * v(1,2,:)          + &
                                gx5 * v(1,3,:)          + &
                                gx6 * v(1,4,:)  ) 
  
          g1v(2,:)          = ( gx1 * v(1,nx-2,:)       + &
                                gx2 * v(1,nx-1,:)       + &
                                gx3 * v(1,1,:)          + &
                                gx4 * v(1,3,:)          + &
                                gx5 * v(1,4,:)          + &
                                gx6 * v(1,5,:)  ) 
  
          g1v(3,:)          = ( gx1 * v(1,nx-1,:)       + &
                                gx2 * v(1,1,:)          + &
                                gx3 * v(1,2,:)          + &
                                gx4 * v(1,4,:)          + &
                                gx5 * v(1,5,:)          + &
                                gx6 * v(1,6,:)  ) 
  
          g1v(nx-2,:)       = ( gx1 * v(1,nx-5,:)       + &
                                gx2 * v(1,nx-4,:)       + &
                                gx3 * v(1,nx-3,:)       + &
                                gx4 * v(1,nx-1,:)       + &
                                gx5 * v(1,1,:)          + &
                                gx6 * v(1,2,:)  ) 
  
          g1v(nx-1,:)       = ( gx1 * v(1,nx-4,:)       + &
                                gx2 * v(1,nx-3,:)       + &
                                gx3 * v(1,nx-2,:)       + &
                                gx4 * v(1,1,:)          + &
                                gx5 * v(1,2,:)          + &
                                gx6 * v(1,3,:)  ) 
  
          g1v(nx,:)     = g1v(1,:)
          
        else
        
          g1v(1,:)      = ( gc1 * v(1,1,:)  + &
                            gc2 * v(1,2,:)  + &
                            gc3 * v(1,3,:)  + &
                            gc4 * v(1,4,:)  + &
                            gc5 * v(1,5,:)  ) * dxinv
  
          g1v(2,:)      = ( gb1 * v(1,1,:)  + &
                            gb2 * v(1,2,:)  + &
                            gb3 * v(1,3,:)  + &
                            gb4 * v(1,4,:)  + &
                            gb5 * v(1,5,:)  ) * dxinv
  
          g1v(3,:)      = ( ga1 * v(1,1,:)  + &
                            ga2 * v(1,2,:)  + &
                            ga3 * v(1,4,:)  + &
                            ga4 * v(1,5,:)  ) * dxinv
          
          g1v(nx-2,:)   = ( ga1 * v(1,nx-4,:)  + &
                            ga2 * v(1,nx-3,:)  + &
                            ga3 * v(1,nx-1,:)  + &
                            ga4 * v(1,nx  ,:)  ) * dxinv
  
          g1v(nx-1,:)  = -( gb1 * v(1,nx  ,:)  + &
                            gb2 * v(1,nx-1,:)  + &
                            gb3 * v(1,nx-2,:)  + &
                            gb4 * v(1,nx-3,:)  + &
                            gb5 * v(1,nx-4,:)  ) * dxinv
  
          g1v(nx,:)    = -( gc1 * v(1,nx  ,:)  + &
                            gc2 * v(1,nx-1,:)  + &
                            gc3 * v(1,nx-2,:)  + &
                            gc4 * v(1,nx-3,:)  + &
                            gc5 * v(1,nx-4,:)  ) * dxinv

        end if

!.... implement the symmetric conditions

        if (lsym) then
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            g1v(1,:)        = ( isign * gx1 * v(1,4,:)  + &
                                isign * gx2 * v(1,3,:)  + &
                                isign * gx3 * v(1,2,:)  + &
                                gx4 * v(1,2,:)          + &
                                gx5 * v(1,3,:)          + &
                                gx6 * v(1,4,:)          ) 
    
            g1v(2,:)        = ( isign * gx1 * v(1,3,:)  + &
                                isign * gx2 * v(1,2,:)  + &
                                gx3 * v(1,1,:)          + &
                                gx4 * v(1,3,:)          + &
                                gx5 * v(1,4,:)          + &
                                gx6 * v(1,5,:)          ) 
    
            g1v(3,:)        = ( isign * gx1 * v(1,2,:)  + &
                                gx2 * v(1,1,:)          + &
                                gx3 * v(1,2,:)          + &
                                gx4 * v(1,4,:)          + &
                                gx5 * v(1,5,:)          + &
                                gx6 * v(1,6,:)          ) 
          end do
        end if

!.... interior

        g1v(4:nx-3,:)  = ( gx1 * v(1,1:nx-6,:)  + &
                           gx2 * v(1,2:nx-5,:)  + &
                           gx3 * v(1,3:nx-4,:)  + &
                           gx4 * v(1,5:nx-2,:)  + &
                           gx5 * v(1,6:nx-1,:)  + &
                           gx6 * v(1,7:nx  ,:)  ) 

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

          g2v(:,:)     = ( gc1 * v(1,:,:)       + &
                           gc2 * v(2,:,:)       + &
                           gc3 * v(3,:,:)       + &
                           gc4 * v(4,:,:)       + &
                           gc5 * v(5,:,:)  ) * dyinv
  
          if (carp) then

          g2v(:,:)     = ( gg1 * v(1,:,:)       + &
                           gg2 * v(2,:,:)       + &
                           gg3 * v(3,:,:)       + &
                           gg4 * v(4,:,:)       + &
                           gg5 * v(5,:,:)       + &
                           gg6 * v(6,:,:)  ) * dyinv
          end if

        return
        end

!=============================================================================!
        subroutine wgrad2 (ndof, nx, ny, v, g2v, g11v, g12v, g22v, dx, dy, &
                           optx, opty, xper, yper, lsym, rsym, bsym, tsym, &
                           carp)
!
!  Take the second derivative of a 2-D field at the wall.
!  Updated to sixth order accurate differencing on the interior with
!  the option for optimized fourth order differencing.
!
!  Added g2v as input to improve efficiency
!
!=============================================================================!
        use stencil
        implicit none
        
        integer :: ndof, nx, ny, optx, opty
        logical :: xper, yper
        logical :: lsym, rsym, bsym, tsym, carp
        real    :: v(ny,nx,ndof), g2v(nx,ndof)
        real    :: g11v(nx,ndof), g12v(nx,ndof), g22v(nx,ndof)
        real    :: dx, dy

        real, parameter :: one = 1.0
        real :: dxinv, dyinv, dxsinv, dysinv
        real :: a, b, c, w
        real :: gx1, gx2, gx3, gx4, gx5, gx6
        real :: gy1, gy2, gy3, gy4, gy5, gy6
        real :: dx1, dx2, dx3, dx4, dx5, dx6, dx7
        real :: dy1, dy2, dy3, dy4, dy5, dy6, dy7

        integer :: idof, isign
!=============================================================================!

        dxinv  = one / dx
        dyinv  = one / dy
        dxsinv = one / dx**2
        dysinv = one / dy**2

!.... seven point stencil in x

        if (optx.eq.0) then
          c = 1.0 / 60.0
        else if (optx.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gx1 =  -c * dxinv
        gx2 =  -b * dxinv
        gx3 =  -a * dxinv
        gx4 =   a * dxinv
        gx5 =   b * dxinv
        gx6 =   c * dxinv

        if (optx.eq.0) then
          c = 1.0 / 90.0
        else if (optx.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = -(3.0*w**2 - 16.0*Sin(w/2.0)**2 + Sin(w)**2) / &
               (12.0*(-15.0*Sin(w/2.0)**2 + 6.0*Sin(w)**2 -  &
               Sin(3.0*w/2.0)**2))
        end if
        
        a = (4.0 + 45.0 * c) / 3.0
        b = -(1.0 + 72.0 * c) / 12.0

        dx1 =  c * dxsinv
        dx2 =  b * dxsinv
        dx3 =  a * dxsinv
        dx4 = -2.0 * ( a + b + c ) * dxsinv
        dx5 =  a * dxsinv
        dx6 =  b * dxsinv
        dx7 =  c * dxsinv

!.... seven point stencil in y

        if (opty.eq.0) then
          c = 1.0 / 60.0
        else if (opty.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gy1 =  -c * dyinv
        gy2 =  -b * dyinv
        gy3 =  -a * dyinv
        gy4 =   a * dyinv
        gy5 =   b * dyinv
        gy6 =   c * dyinv

        if (opty.eq.0) then
          c = 1.0 / 90.0
        else if (opty.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = -(3.0*w**2 - 16.0*Sin(w/2.0)**2 + Sin(w)**2) / &
               (12.0*(-15.0*Sin(w/2.0)**2 + 6.0*Sin(w)**2 -  &
               Sin(3.0*w/2.0)**2))
        end if
        
        a = (4.0 + 45.0 * c) / 3.0
        b = -(1.0 + 72.0 * c) / 12.0

        dy1 =  c * dysinv
        dy2 =  b * dysinv
        dy3 =  a * dysinv
        dy4 = -2.0 * ( a + b + c ) * dysinv
        dy5 =  a * dysinv
        dy6 =  b * dysinv
        dy7 =  c * dysinv

!=============================================================================!
!.... compute the second derivative in x
!=============================================================================!

        if (xper) then

          g11v(1,:)         = ( dx1 * v(1,nx-3,:)   + &
                                dx2 * v(1,nx-2,:)   + &
                                dx3 * v(1,nx-1,:)   + &
                                dx4 * v(1,1,:)      + &
                                dx5 * v(1,2,:)   + &
                                dx6 * v(1,3,:)   + &
                                dx7 * v(1,4,:)      ) 
  
          g11v(2,:)         = ( dx1 * v(1,nx-2,:)   + &
                                dx2 * v(1,nx-1,:)   + &
                                dx3 * v(1,1,:)      + &
                                dx4 * v(1,2,:)      + &
                                dx5 * v(1,3,:)      + &
                                dx6 * v(1,4,:)      + &
                                dx7 * v(1,5,:)      ) 
  
          g11v(3,:)         = ( dx1 * v(1,nx-1,:)   + &
                                dx2 * v(1,1,:)      + &
                                dx3 * v(1,2,:)      + &
                                dx4 * v(1,3,:)      + &
                                dx5 * v(1,4,:)      + &
                                dx6 * v(1,5,:)      + &
                                dx7 * v(1,6,:)      ) 
  
          g11v(nx-2,:)      = ( dx1 * v(1,nx-5,:)   + &
                                dx2 * v(1,nx-4,:)   + &
                                dx3 * v(1,nx-3,:)   + &
                                dx4 * v(1,nx-2,:)   + &
                                dx5 * v(1,nx-1,:)   + &
                                dx6 * v(1,1,:)      + &
                                dx7 * v(1,2,:)      ) 
  
          g11v(nx-1,:)      = ( dx1 * v(1,nx-4,:)   + &
                                dx2 * v(1,nx-3,:)   + &
                                dx3 * v(1,nx-2,:)   + &
                                dx4 * v(1,nx-1,:)   + &
                                dx5 * v(1,1,:)      + &
                                dx6 * v(1,2,:)      + &
                                dx7 * v(1,3,:)   ) 
  
          g11v(nx,:)        = g11v(1,:)
          
        else
        
          g11v(1,:)         = ( dd1 * v(1,1,:) + &
                                dd2 * v(1,2,:) + &
                                dd3 * v(1,3,:) + &
                                dd4 * v(1,4,:) + &
                                dd5 * v(1,5,:) ) * dxsinv
  
          g11v(2,:)         = ( db1 * v(1,1,:) + &
                                db2 * v(1,2,:) + &
                                db3 * v(1,3,:) + &
                                db4 * v(1,4,:) + &
                                db5 * v(1,5,:) ) * dxsinv
  
          g11v(3,:)         = ( da1 * v(1,1,:) + &
                                da2 * v(1,2,:) + &
                                da3 * v(1,3,:) + &
                                da4 * v(1,4,:) + &
                                da5 * v(1,5,:) ) * dxsinv
  
          g11v(nx-2,:)      = ( da1 * v(1,nx-4,:) + &
                                da2 * v(1,nx-3,:) + &
                                da3 * v(1,nx-2,:) + &
                                da4 * v(1,nx-1,:) + &
                                da5 * v(1,nx,:)   ) * dxsinv
  
          g11v(nx-1,:)      = ( db1 * v(1,nx,:)   + &
                                db2 * v(1,nx-1,:) + &
                                db3 * v(1,nx-2,:) + &
                                db4 * v(1,nx-3,:) + &
                                db5 * v(1,nx-4,:) ) * dxsinv
  
          g11v(nx,:)        = ( dd1 * v(1,nx,:)   + &
                                dd2 * v(1,nx-1,:) + &
                                dd3 * v(1,nx-2,:) + &
                                dd4 * v(1,nx-3,:) + &
                                dd5 * v(1,nx-4,:) ) * dxsinv

        end if
        
!.... implement symmetry conditions

        if (lsym) then
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
            
            g11v(1,idof)    = ( isign * dx1 * v(1,4,idof)       + &
                                isign * dx2 * v(1,3,idof)       + &
                                isign * dx3 * v(1,2,idof)       + &
                                dx4 * v(1,1,idof)               + &
                                dx5 * v(1,2,idof)               + &
                                dx6 * v(1,3,idof)               + &
                                dx7 * v(1,4,idof)               ) 
    
            g11v(2,idof)    = ( isign * dx1 * v(1,3,idof)       + &
                                isign * dx2 * v(1,2,idof)       + &
                                dx3 * v(1,1,idof)               + &
                                dx4 * v(1,2,idof)               + &
                                dx5 * v(1,3,idof)               + &
                                dx6 * v(1,4,idof)               + &
                                dx7 * v(1,5,idof)               ) 
    
            g11v(3,idof)    = ( isign * dx1 * v(1,2,idof)       + &
                                dx2 * v(1,1,idof)               + &
                                dx3 * v(1,2,idof)               + &
                                dx4 * v(1,3,idof)               + &
                                dx5 * v(1,4,idof)               + &
                                dx6 * v(1,5,idof)               + &
                                dx7 * v(1,6,idof)               ) 
          end do
        end if
        
!.... interior

        g11v(4:nx-3,:)   = ( dx1 * v(1,1:nx-6,:)   + &
                             dx2 * v(1,2:nx-5,:)   + &
                             dx3 * v(1,3:nx-4,:)   + &
                             dx4 * v(1,4:nx-3,:)   + &
                             dx5 * v(1,5:nx-2,:)   + &
                             dx6 * v(1,6:nx-1,:)   + &
                             dx7 * v(1,7:nx  ,:)   ) 

!=============================================================================!
!.... compute the second derivative in y
!=============================================================================!

          g22v(:,:)        =  ( dd1 * v(1,:,:) + &
                                dd2 * v(2,:,:) + &
                                dd3 * v(3,:,:) + &
                                dd4 * v(4,:,:) + &
                                dd5 * v(5,:,:) ) * dysinv
    
!=============================================================================!
!.... compute the cross derivative
!=============================================================================!

        if (xper) then

          g12v(1,:)         = ( gx1 * g2v(nx-3,:)       + &
                                gx2 * g2v(nx-2,:)       + &
                                gx3 * g2v(nx-1,:)       + &
                                gx4 * g2v(2,:)          + &
                                gx5 * g2v(3,:)          + &
                                gx6 * g2v(4,:)  ) 
  
          g12v(2,:)         = ( gx1 * g2v(nx-2,:)       + &
                                gx2 * g2v(nx-1,:)       + &
                                gx3 * g2v(1,:)          + &
                                gx4 * g2v(3,:)          + &
                                gx5 * g2v(4,:)          + &
                                gx6 * g2v(5,:)  ) 
  
          g12v(3,:)         = ( gx1 * g2v(nx-1,:)       + &
                                gx2 * g2v(1,:)          + &
                                gx3 * g2v(2,:)          + &
                                gx4 * g2v(4,:)          + &
                                gx5 * g2v(5,:)          + &
                                gx6 * g2v(6,:)  ) 
  
          g12v(nx-2,:)      = ( gx1 * g2v(nx-5,:)       + &
                                gx2 * g2v(nx-4,:)       + &
                                gx3 * g2v(nx-3,:)       + &
                                gx4 * g2v(nx-1,:)       + &
                                gx5 * g2v(1,:)          + &
                                gx6 * g2v(2,:)  ) 
  
          g12v(nx-1,:)      = ( gx1 * g2v(nx-4,:)       + &
                                gx2 * g2v(nx-3,:)       + &
                                gx3 * g2v(nx-2,:)       + &
                                gx4 * g2v(1,:)          + &
                                gx5 * g2v(2,:)          + &
                                gx6 * g2v(3,:)  ) 
  
          g12v(nx,:)     = g12v(1,:)
          
        else
        
          g12v(1,:)     = ( gc1 * g2v(1,:)  + &
                            gc2 * g2v(2,:)  + &
                            gc3 * g2v(3,:)  + &
                            gc4 * g2v(4,:)  + &
                            gc5 * g2v(5,:)  ) * dxinv
  
          g12v(2,:)     = ( gb1 * g2v(1,:)  + &
                            gb2 * g2v(2,:)  + &
                            gb3 * g2v(3,:)  + &
                            gb4 * g2v(4,:)  + &
                            gb5 * g2v(5,:)  ) * dxinv
  
          g12v(3,:)     = ( ga1 * g2v(1,:)  + &
                            ga2 * g2v(2,:)  + &
                            ga3 * g2v(4,:)  + &
                            ga4 * g2v(5,:)  ) * dxinv
          
          g12v(nx-2,:)  = ( ga1 * g2v(nx-4,:)  + &
                            ga2 * g2v(nx-3,:)  + &
                            ga3 * g2v(nx-1,:)  + &
                            ga4 * g2v(nx  ,:)  ) * dxinv
  
          g12v(nx-1,:) = -( gb1 * g2v(nx  ,:)  + &
                            gb2 * g2v(nx-1,:)  + &
                            gb3 * g2v(nx-2,:)  + &
                            gb4 * g2v(nx-3,:)  + &
                            gb5 * g2v(nx-4,:)  ) * dxinv
  
          g12v(nx,:)   = -( gc1 * g2v(nx  ,:)  + &
                            gc2 * g2v(nx-1,:)  + &
                            gc3 * g2v(nx-2,:)  + &
                            gc4 * g2v(nx-3,:)  + &
                            gc5 * g2v(nx-4,:)  ) * dxinv

        end if

!.... implement the symmetric conditions

        if (lsym) then
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            g12v(1,:)       = ( isign * gx1 * g2v(4,:)  + &
                                isign * gx2 * g2v(3,:)  + &
                                isign * gx3 * g2v(2,:)  + &
                                gx4 * g2v(2,:)          + &
                                gx5 * g2v(3,:)          + &
                                gx6 * g2v(4,:)          ) 
    
            g12v(2,:)       = ( isign * gx1 * g2v(3,:)  + &
                                isign * gx2 * g2v(2,:)  + &
                                gx3 * g2v(1,:)          + &
                                gx4 * g2v(3,:)          + &
                                gx5 * g2v(4,:)          + &
                                gx6 * g2v(5,:)          ) 
    
            g12v(3,:)       = ( isign * gx1 * g2v(2,:)  + &
                                gx2 * g2v(1,:)          + &
                                gx3 * g2v(2,:)          + &
                                gx4 * g2v(4,:)          + &
                                gx5 * g2v(5,:)          + &
                                gx6 * g2v(6,:)          ) 
          end do
        end if
        
!.... interior

        g12v(4:nx-3,:) = ( gx1 * g2v(1:nx-6,:)  + &
                           gx2 * g2v(2:nx-5,:)  + &
                           gx3 * g2v(3:nx-4,:)  + &
                           gx4 * g2v(5:nx-2,:)  + &
                           gx5 * g2v(6:nx-1,:)  + &
                           gx6 * g2v(7:nx  ,:)  ) 

        return
        end

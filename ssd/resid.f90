!=============================================================================!
        subroutine resid
!  
!  Computes the residual for the compressible Navier-Stokes solver. 
!  
!=============================================================================!
        use stuff
        use material
        implicit none
!=============================================================================!

!.... local variables
        
        integer ix, idof, jdof
        
        real vint(ny,ndof,nx), g1v(ny,ndof,nx),  g2v(ny,ndof,nx)
        real g11v(ny,ndof,nx), g12v(ny,ndof,nx), g22v(ny,ndof,nx)
        
        character*80 name
!=============================================================================!

!       write(*,*) 'Computing residual'

!.... compute v at the midpoint in time

!       if (iter .eq. 1) then
!         vint = v
!       else
!         vint = pt5 * ( v + vold )
!       endif

!.... compute first derivatives
        
!       call grad(ndof, nx, ny, vint, g1v, g2v, dxi, deta)

!.... compute second derivatives

!       call grad2(ndof, nx, ny, vint, g11v, g12v, g22v, dxi, deta)

!.... transform the gradients to physical space

!       do ix = 1, nx
!         do idof = 1, ndof
!           g22v(:,idof,ix) = g22v(:,idof,ix) * m1**2 + g2v(:,idof,ix) * m2
!           g2v (:,idof,ix) = g2v (:,idof,ix) * m1
!           g12v(:,idof,ix) = g12v(:,idof,ix) * m1
!         end do
!       end do

!.... initialize the residual

!       r = zero
        
!.... form the residual

!       if (iter .eq. 1) then
!         do ix = 1, nx
!           do idof = 1, ndof
!             do jdof = 1, ndof
!               r(:,idof,ix) = r(:,idof,ix) + delt * (                &
!                               A(:,idof,jdof,ix) *  g1v(:,jdof,ix) + &
!                               B(:,idof,jdof,ix) *  g2v(:,jdof,ix) + &
!                               D(:,idof,jdof,ix) * vint(:,jdof,ix) - &
!                             Vxx(:,idof,jdof,ix) * g11v(:,jdof,ix) - &
!                             Vxy(:,idof,jdof,ix) * g12v(:,jdof,ix) - &
!                             Vyy(:,idof,jdof,ix) * g22v(:,jdof,ix) )
!             end do
!           end do
!         end do
!       else
!         do ix = 1, nx
!           do idof = 1, ndof
!             do jdof = 1, ndof
!               r(:,idof,ix) = r(:,idof,ix) + delt * (                   &
!                               A(:,idof,jdof,ix) *   g1v(:,jdof,ix)   + &
!                               B(:,idof,jdof,ix) *   g2v(:,jdof,ix)   + &
!                               D(:,idof,jdof,ix) *  vint(:,jdof,ix)   - &
!                             Vxx(:,idof,jdof,ix) *  g11v(:,jdof,ix)   - &
!                             Vxy(:,idof,jdof,ix) *  g12v(:,jdof,ix)   - &
!                             Vyy(:,idof,jdof,ix) *  g22v(:,jdof,ix) ) + &
!                               G(:,idof,jdof,ix) * (   v(:,jdof,ix)   - &
!                                                    vold(:,jdof,ix) )
!             end do
!           end do
!         end do
!       end if
        
!.... calculate the residual statistics

!       call resstat(r)

!.... write out the disturbance field in Plot3d format

!       if (istep .eq. 1) then
!         write(*,*) 'Wrote residual to Plot3d file'
!         name = 'r.dat'
!         call wdata(name, r, nx, ny, nz, ndof, Ma, 0.0, Re, 0.0)
!       end if
        
        return
        end

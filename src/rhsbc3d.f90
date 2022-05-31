!=============================================================================!
        subroutine rhsBC3D(rl, vl, vml) 
!  
!  Satisfy the boundary conditions for the 3D residual
!  
!=============================================================================!
        use global
        use stencil
        use bump
        implicit none

        complex :: rl(ndof,nx,ny), vl(ndof,nx,ny)
        real    :: vml(ndof,nx,ny)
        
!.... local variables

        real :: tmp
        
        real, allocatable :: rhom(:), cm(:), tm(:), um(:)
        complex, allocatable :: c1(:), c2(:), c3(:), c4(:)
        complex, allocatable :: rho(:), u1(:), u2(:), u3(:), t(:), p(:)
        
        complex :: ub(nx)
        
        integer :: i, j

!.... forcing parameters

        real :: a, d, kk

!!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (yper) then
        
          rl(:,:,ny) = zero
        
        else                    ! yper

!=============================================================================!
!       W a l l
!=============================================================================!
        if (Navier) then
        
!.... no-slip boundary condition

!!$       if (wall.eq.3) then
!!$         rl(2,is:ie,1) = -(vl(2,is:ie,1) - u1w)
!!$         rl(3,is:ie,1) = -(vl(3,is:ie,1) - u2w)
!!$         rl(4,is:ie,1) = -(vl(4,is:ie,1) - u3w)
!!$       else
            rl(2,is:ie,1) = -(vl(2,is:ie,1) - zero)
            rl(3,is:ie,1) = -(vl(3,is:ie,1) - zero)
            rl(4,is:ie,1) = -(vl(4,is:ie,1) - zero)
!!$       end if

!.... density boundary condition

          if (wall.eq.1) then
            rl(1,is:ie,1) = gc1 * vl(1,is:ie,1) + gc2 * vl(1,is:ie,2) + &
                        gc3 * vl(1,is:ie,3) + gc4 * vl(1,is:ie,4) + &
                        gc5 * vl(1,is:ie,5)
          end if
          
!.... linear extrapolation of rho at the wall

          if (wall.eq.2) then
            rl(1,is:ie,1) = vl(1,is:ie,1) - two * vl(1,is:ie,2) + vl(1,is:ie,3)
          end if

!.... Third order extrapolation to the wall for density
!.... (in computational space)

          if (wall.eq.3) then
            rl(1,is:ie,1) = vl(1,is:ie,1) - (4.0*vl(1,is:ie,2) - 6.0*vl(1,is:ie,3) + &
                                     4.0*vl(1,is:ie,4) -     vl(1,is:ie,5) )
          end if

!.... isothermal wall

          if (wallt.eq.0) then
!           if (wall.eq.3) then
!             rl(ndof,is:ie,1) = -(vl(ndof,is:ie,1) - tw)
!           else
              rl(ndof,is:ie,1) = zero
!           end if
          end if

!.... adiabatic boundary condition

          if (wallt.eq.1) then
!           if (wall.eq.3) then
!             rl(ndof,is:ie,1) = gc1 * vl(ndof,is:ie,1) + gc2 * vl(ndof,is:ie,2) + &
!                            gc3 * vl(ndof,is:ie,3) + gc4 * vl(ndof,is:ie,4) + &
!                            gc5 * vl(ndof,is:ie,5) - twp
!           else
              rl(ndof,is:ie,1) = gc1 * vl(ndof,is:ie,1) + gc2 * vl(ndof,is:ie,2) + &
                             gc3 * vl(ndof,is:ie,3) + gc4 * vl(ndof,is:ie,4) + &
                             gc5 * vl(ndof,is:ie,5)
!           end if
          end if

        else     ! inviscid
        
!....     rl(2,is:ie,1) = wall tangent momentum
!....     rl(3,is:ie,1) = wall normal velocity
        
          rl(2,is:ie,1) =  ( bnb(is:ie,2) * rl(2,is:ie,1) - bnb(is:ie,1) * rl(3,is:ie,1) )        
          rl(3,is:ie,1) = -( bnb(is:ie,1) * vl(2,is:ie,1) + bnb(is:ie,2) * vl(3,is:ie,1) )

        end if   ! Navier
!=============================================================================!
!       T o p
!=============================================================================!

!.... freestream zero disturbance boundary conditions

        if (top.eq.0) then

          rl(:,:,ny) = zero

        else if (top.eq.1) then

!.... compute the characteristic amplitudes on the top boundary
          
          allocate( rhom(nx), tm(nx), cm(nx), um(nx), c3(nx) )

          rhom = vml(1,:,ny)
          tm   = vml(5,:,ny)
          cm   = sqrt( tm ) / Ma
          um   = vml(2,:,ny)

          d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

          do i = 1, nx
            kk = omega / (cm(i)+um(i))
            a  = omega**2 * d / (cm(i)+um(i))**3
            c3(i) = exp( -a * (x(ny,i) - x0) ) * exp( im * kk * x(ny,i) )
!           c3(i) = wamp(i) * exp( im * kk * x(ny,i) )
          end do

          rl(1,:,ny) = -(vl(1,:,ny) - pt5 * c3 / cm**2)
          rl(2,:,ny) = -(vl(2,:,ny) - c3 * pt5 / ( rhom * cm ))
          rl(3,:,ny) = -(vl(3,:,ny) - zero)
          rl(4,:,ny) = -(vl(4,:,ny) - zero)
          rl(5,:,ny) = -(vl(5,:,ny) - (gamma*Ma**2 * c3 * pt5 - &
                                       tm * pt5 * c3 / cm**2) / rhom)

          deallocate( rhom, tm, cm, um, c3 )

        end if
        
        end if                  ! yper

        if (xper) then
          rl(:,nx,:) = zero
        else                    ! xper
        
!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!
        if (left.eq.0) then             ! zero disturbance

          rl(:,1,:) = zero

        else if (left.eq.4) then        ! eigenfunction inflow

          !write(*,*) "rhsbc3d:  left.eq.4"
        
#if 0
          rl(:,1,:) = zero
#else
          allocate( rhom(ny), tm(ny), cm(ny) )
          allocate( c1(ny), c2(ny), c3(ny), c4(ny) )
          allocate( rho(ny), u1(ny), u2(ny), u3(ny), t(ny), p(ny) )

          rhom = vml(1,1,:)
          tm   = vml(5,1,:)
          cm   = one / Ma * sqrt( tm )

!.... compute incomming characteristics (eigenfunctions)

          rho = cmplx(rhor(1:ny), rhoi(1:ny))
          u1  = cmplx(  ur(1:ny),   ui(1:ny))
          u2  = cmplx(  vr(1:ny),   vi(1:ny))
          u3  = cmplx(  wr(1:ny),   wi(1:ny))
          t   = cmplx(  tr(1:ny),   ti(1:ny))

          rl(1,1,:) = -(vl(1,1,:) - rho)
          rl(2,1,:) = -(vl(2,1,:) -  u1)
          rl(3,1,:) = -(vl(3,1,:) -  u2)
          rl(4,1,:) = -(vl(4,1,:) -  u3)
          rl(5,1,:) = -(vl(5,1,:) -   t)

          deallocate( rhom, tm, cm )
          deallocate( c1, c2, c3, c4 )
          deallocate( rho, u1, u2, u3, t, p )
#endif
        
        else if (left.eq.5) then        ! acoustic wave inflow
          
          allocate( rhom(ny), tm(ny), cm(ny), um(ny), c3(ny) )

          rhom = vml(1,1,:)
          tm   = vml(5,1,:)
          cm   = sqrt( tm ) / Ma
          um   = vml(2,1,:)

          do j = 1, ny
            kk = omega / ( cm(j)+um(j) )
            c3(j) = exp( im * kk * x(j,1) )
          end do

          rl(1,1,:) = -(vl(1,1,:) - pt5 * c3 / cm**2)
          rl(2,1,:) = -(vl(2,1,:) - c3 * pt5 / ( rhom * cm ))
          rl(3,1,:) = -(vl(3,1,:) - zero)
          rl(4,1,:) = -(vl(4,1,:) - zero)
          rl(5,1,:) = -(vl(5,1,:) - (gamma*Ma**2 * c3 * pt5 - &
                                     tm * pt5 * c3 / cm**2) / rhom)

          deallocate( rhom, tm, cm, um, c3 )

        else if (left.eq.7) then           ! symmetry boundary

        end if

!=============================================================================!
!       R i g h t   B o u n d a r y
!=============================================================================!
        if (right.eq.0) then              ! zero disturbance

          rl(:,nx,:) = zero

        else if (right.eq.4) then
        
          rl(:,nx,:) = zero

        else if (right.eq.5) then

!.... compute the characteristic amplitudes on the right boundary
          
          allocate( rhom(ny), tm(ny), cm(ny), um(ny), c3(ny) )

          rhom = vml(1,nx,:)
          tm   = vml(5,nx,:)
          cm   = sqrt( tm ) / Ma
          um   = vml(2,nx,:)

          do j = 1, ny
            kk = omega / ( cm(j)+um(j) )
            c3(j) = exp( im * kk * x(j,nx) )
          end do

          rl(1,nx,:) = -(vl(1,nx,:) - pt5 * c3 / cm**2)
          rl(2,nx,:) = -(vl(2,nx,:) - c3 * pt5 / ( rhom * cm ))
          rl(3,nx,:) = -(vl(3,nx,:) - zero)
          rl(4,nx,:) = -(vl(4,nx,:) - zero)
          rl(5,nx,:) = -(vl(5,nx,:) - (gamma*Ma**2 * c3 * pt5 - &
                                     tm * pt5 * c3 / cm**2) / rhom)

          deallocate( rhom, tm, cm, um, c3 )

        else if (right.eq.7) then        ! symmetry boundary

        else if (right.eq.8) then        ! hold initial condition

          rl(:,nx,:) = zero

        else if (right .eq. 9) then      ! extrapolation 

          rl(:,nx,:) = -( vl(:,nx,:) - exp( two * log(vl(:,nx-1,:)) - &
                                                  log(vl(:,nx-2,:)) ) )

        end if

        end if                            ! xper
                        
        return
        end

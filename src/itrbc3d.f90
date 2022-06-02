!=============================================================================!
        subroutine itrBC3D(vl,vml) 
!  
!  Satisfy the boundary conditions for the 3D variables
!
!  Revised: 4-17-96
!=============================================================================!
        use global
        use stencil
        use bump
        implicit none

        complex :: vl(ndof,nx,ny)
        real :: vml(ndof,nx,ny)

!.... local variables

        real :: tmp
        
        real, allocatable :: rhom(:), cm(:), tm(:), um(:)
        complex, allocatable :: c1(:), c2(:), c3(:), c4(:)
        complex, allocatable :: rho(:), u1(:), u2(:), u3(:), t(:), p(:)
        
        complex :: ub(nx), vb(nx)
        
        integer :: i, j

!.... forcing parameters

        real :: a, d, kk

!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (yper) then
        
          vl(:,:,ny) = vl(:,:,1)
          
        else                    ! yper
        
!=============================================================================!
!       W a l l 
!=============================================================================!
        if (Navier) then

!.... no slip wall

!!$     if (wall.eq.3) then   ! bump
!!$       vl(2,:,1) = u1w
!!$       vl(3,:,1) = u2w
!!$       vl(4,:,1) = u3w
!!$     else
          vl(2,is:ie,1) = zero
          vl(3,is:ie,1) = zero
          vl(4,is:ie,1) = zero
!!$     end if

!.... set drho/dy = 0 which implies dp/deta = 0

        if (wall.eq.1) then
          vl(1,is:ie,1) = -( gc2 * vl(1,is:ie,2)  + &
                         gc3 * vl(1,is:ie,3)  + &
                         gc4 * vl(1,is:ie,4)  + &
                         gc5 * vl(1,is:ie,5)  ) / gc1
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          vl(1,is:ie,1) = two * vl(1,is:ie,2) - vl(1,is:ie,3)
        end if

!.... Third order extrapolation to the wall for density
!.... (in computational space)

        if (wall.eq.3) then
          vl(1,is:ie,1) = 4.0*vl(1,is:ie,2) - 6.0*vl(1,is:ie,3) + &
                      4.0*vl(1,is:ie,4) - vl(1,is:ie,5)
        end if

!.... isothermal wall

        if (wallt.eq.0) then
!         if (wall.eq.3) then   ! bump
!           vl(ndof,is:ie,1) = tw
!         else
            vl(ndof,is:ie,1) = zero
!         end if
        end if

!.... adiabatic temperature ( dT/deta = 0 )

        if (wallt.eq.1) then
!         if (wall.eq.3) then   ! bump
!           if (carp) then      ! Carpenter's stencil
!             vl(ndof,is:ie,1) = ( twp - ( gg2 * vl(ndof,is:ie,2)  + &
!                               gg3 * vl(ndof,is:ie,3)  + &
!                               gg4 * vl(ndof,is:ie,4)  + &
!                               gg5 * vl(ndof,is:ie,5)  + &
!                               gg6 * vl(ndof,is:ie,6)  ) ) / gg1
!           else
!             vl(ndof,is:ie,1) = ( twp - ( gc2 * vl(ndof,is:ie,2)  + &
!                               gc3 * vl(ndof,is:ie,3)  + &
!                               gc4 * vl(ndof,is:ie,4)  + &
!                               gc5 * vl(ndof,is:ie,5)  ) ) / gc1
!           end if        
!         else                  ! no bump
            if (carp) then      ! Carpenter's stencil
              vl(ndof,is:ie,1) = -( gg2 * vl(ndof,is:ie,2)  + &
                                gg3 * vl(ndof,is:ie,3)  + &
                                gg4 * vl(ndof,is:ie,4)  + &
                                gg5 * vl(ndof,is:ie,5)  + &
                                gg6 * vl(ndof,is:ie,6)  ) / gg1
            else
              vl(ndof,is:ie,1) = -( gc2 * vl(ndof,is:ie,2)  + &
                                gc3 * vl(ndof,is:ie,3)  + &
                                gc4 * vl(ndof,is:ie,4)  + &
                                gc5 * vl(ndof,is:ie,5)  ) / gc1
            end if
!         end if
        end if

        else    ! inviscid:  zero the wall-normal velocity

          ub(is:ie) = bnb(is:ie,2)*vl(2,is:ie,1) - bnb(is:ie,1)*vl(3,is:ie,1)
!         vb(is:ie) = bnb(is:ie,1)*vl(2,is:ie,1) + bnb(is:ie,2)*vl(3,is:ie,1)
          vb(is:ie) = zero
          
          vl(2,is:ie,1) =  bnb(is:ie,2) * ub(is:ie) + bnb(is:ie,1) * vb(is:ie)
          vl(3,is:ie,1) = -bnb(is:ie,1) * ub(is:ie) + bnb(is:ie,2) * vb(is:ie)

        end if  ! Navier
!=============================================================================!
!       T o p   B o u n d a r y
!=============================================================================!

!.... apply zero disturbance boundary condition

        if (top.eq.0) then
        
          vl(:,:,ny) = zero

!.... force an incomming acoustic wave

        else if (top.eq.1) then

!.... allocate room for the boundary amplitude

!         if (istep.eq.0) then
!           allocate( wamp(nx) )
!           open(66,file='amp.dat',form='formatted',status='unknown')
!           do i = 1, nx
!             read(66,*) tmp, wamp(i)
!           end do
!           close(66)
!         end if
        
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
            ! the first term is a viscous damping correction
            !c3(i) = exp(-a * (x(i,ny) - x0)) * exp(im * kk * x(i,ny))
            c3(i) = exp(im * kk * x(i,ny))
            !c3(i) = wamp(i) * exp( im * kk * x(i,ny) )
          end do

          vl(1,:,ny) = pt5 * c3 / cm**2
          vl(2,:,ny) = c3 * pt5 / ( rhom * cm )
          vl(3,:,ny) = zero
          vl(4,:,ny) = zero
          vl(5,:,ny) = (gamma*Ma**2*c3*pt5 - tm*pt5*c3/cm**2)/rhom

          deallocate( rhom, tm, cm, um, c3 )

        end if

        end if                  ! yper

        if (xper) then

          vl(:,nx,:) = vl(:,1,:)
          
        else                    ! xper

!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!

        if (left.eq.0) then             ! zero disturbance BC

          vl(1,1,:) = zero
          vl(2,1,:) = zero
          vl(3,1,:) = zero
          vl(4,1,:) = zero
          vl(5,1,:) = zero
        
        else if (left.eq.4) then        ! eigenfunction inflow disturbance

          !write(*,*) "itrbc3d:  left.eq.4"

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
#ifdef USE_TRANSIENT_EIGENFUNCTION 
          if (omega.eq.0) then
            rho = rho * exp(-im*lomega*time)
            u1  =  u1 * exp(-im*lomega*time)
            u2  =  u2 * exp(-im*lomega*time)
            u3  =  u3 * exp(-im*lomega*time)
            t   =   t * exp(-im*lomega*time)
          endif
#endif
#if 0

!.... SSC: turned this on for spatial TS wave case 5/30/22

          p   = ( rhom * t + tm * rho ) / (gamma * Ma**2)
                                  
          c1  = -cm**2 * rho + p                        ! entropy
          c2  =  rhom * cm * u2                         ! vorticity
          c3  =  rhom * cm * u1 + p                     ! right acoustic

#if 1

!.... compute outgoing characteristics
!.... SSC: try not setting left acoustic to zero

          rho = vl(1,1,:)
          u1  = vl(2,1,:)
          u2  = vl(3,1,:)
          u3  = vl(4,1,:)
          t   = vl(5,1,:)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

          c4  =  -rhom * cm * u1 + p                    ! left acoustic
#else
          c4  =  0                                      ! left acoustic
#endif

!.... update the boundary values

          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
          u2  = c2 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom
#endif
          vl(1,1,:) = rho
          vl(2,1,:) = u1
          vl(3,1,:) = u2
          vl(4,1,:) = u3
          vl(5,1,:) = t

          deallocate( rhom, tm, cm )
          deallocate( c1, c2, c3, c4 )
          deallocate( rho, u1, u2, u3, t, p )

        else if (left.eq.5) then        ! acoustic inflow disturbance
          
          allocate( rhom(ny), tm(ny), cm(ny), um(ny), c3(ny) )

          rhom = vml(1,1,:)
          tm   = vml(5,1,:)
          cm   = sqrt( tm ) / Ma
          um   = vml(2,1,:)

          do j = 1, ny
            kk = omega / ( cm(j)+um(j) )
            c3(j) = exp( im * kk * x(j,1) )
          end do

          vl(1,1,:) = pt5 * c3 / cm**2
          vl(2,1,:) = c3 * pt5 / ( rhom * cm )
          vl(3,1,:) = zero
          vl(4,1,:) = zero
          vl(5,1,:) = (gamma*Ma**2 * c3 * pt5 - tm * pt5 * c3 / cm**2) / rhom

          deallocate( rhom, tm, cm, um, c3 )

        else if (left.eq.7) then        ! symmetry boundary

        end if
!=============================================================================!
!       R i g h t   B o u n d a r y
!=============================================================================!

        if (right.eq.0) then            ! zero disturbance BC

          vl(1,nx,:) = zero
          vl(2,nx,:) = zero
          vl(3,nx,:) = zero
          vl(4,nx,:) = zero
          vl(5,nx,:) = zero

        else if (right.eq.4) then

!.... SSC: this is currently hardwired for the Ch.4 thesis spatial TS wave

!         ac = (2.2804739410500E-001,-6.5163146761218E-003)
!         ac = (-2.8831962908130E-001,-1.3854663671636E-002)

          allocate( rho(ny), u1(ny), u2(ny), u3(ny), t(ny), p(ny) )

!.... compute incomming characteristics (eigenfunctions)

          rho = cmplx(rhor(:), rhoi(:)) * exp(im * lalpha * x(nx,:))
          u1  = cmplx(  ur(:),   ui(:)) * exp(im * lalpha * x(nx,:))
          u2  = cmplx(  vr(:),   vi(:)) * exp(im * lalpha * x(nx,:))
          u3  = cmplx(  wr(:),   wi(:)) * exp(im * lalpha * x(nx,:))
          t   = cmplx(  tr(:),   ti(:)) * exp(im * lalpha * x(nx,:))

          vl(1,nx,:) = rho
          vl(2,nx,:) = u1
          vl(3,nx,:) = u2
          vl(4,nx,:) = u3
          vl(5,nx,:) = t

          deallocate( rho, u1, u2, u3, t, p )

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

          vl(1,nx,:) = pt5 * c3 / cm**2
          vl(2,nx,:) = c3 * pt5 / ( rhom * cm )
          vl(3,nx,:) = zero
          vl(4,nx,:) = zero
          vl(5,nx,:) = (gamma*Ma**2 * c3 * pt5 - tm * pt5 * c3 / cm**2) / rhom

          deallocate( rhom, tm, cm, um, c3 )

        else if (right .eq. 9) then      ! extrapolation 

          vl(:,nx,:) = exp( two * log(vl(:,nx-1,:)) - log(vl(:,nx-2,:)) )

        end if
        
        end if                  ! xper

        return
        end

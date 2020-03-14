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

        complex :: vl(ny,nx,ndof)
        real :: vml(ny,nx,ndof)

!.... local variables

        real :: tmp
        
        real, allocatable :: rhom(:), cm(:), tm(:), um(:)
        complex, allocatable :: c1(:), c2(:), c3(:), c4(:)
        complex, allocatable :: rho(:), u1(:), u2(:), u3(:), t(:), p(:)
        
        complex :: ub(nx), vb(nx)
        
        integer :: i, j

!.... forcing parameters

        real :: a, d, kk
        complex :: ac

!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (yper) then
        
          vl(ny,:,:) = vl(1,:,:)
          
        else                    ! yper
        
!=============================================================================!
!       W a l l 
!=============================================================================!
        if (Navier) then

!.... set drho/dy = 0 which implies dp/deta = 0

        if (wall.eq.1) then
          vl(1,:,1) = -( gc2 * vl(2,:,1)  + &
                         gc3 * vl(3,:,1)  + &
                         gc4 * vl(4,:,1)  + &
                         gc5 * vl(5,:,1)  ) / gc1
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          vl(1,:,1) = two * vl(2,:,1) - vl(3,:,1)
        end if

!.... no slip wall

        if (wall.eq.3) then   ! bump
          vl(1,:,2) = u1w
          vl(1,:,3) = u2w
          vl(1,:,4) = u3w
        else
          vl(1,:,2) = zero
          vl(1,:,3) = zero
          vl(1,:,4) = zero
        end if

!.... isothermal wall

        if (wallt.eq.0) then
          if (wall.eq.3) then   ! bump
            vl(1,:,ndof) = tw
          else
            vl(1,:,ndof) = zero
          end if
        end if

!.... adiabatic temperature ( dT/deta = 0 )

        if (wallt.eq.1) then
          if (wall.eq.3) then   ! bump
            if (impl.eq.0) then ! Carpenter's stencil
              vl(1,:,ndof) = ( twp - ( gg2 * vl(2,:,ndof)  + &
                                gg3 * vl(3,:,ndof)  + &
                                gg4 * vl(4,:,ndof)  + &
                                gg5 * vl(5,:,ndof)  + &
                                gg6 * vl(6,:,ndof)  ) ) / gg1
            else
              vl(1,:,ndof) = ( twp - ( gc2 * vl(2,:,ndof)  + &
                                gc3 * vl(3,:,ndof)  + &
                                gc4 * vl(4,:,ndof)  + &
                                gc5 * vl(5,:,ndof)  ) ) / gc1
            end if        
          else                  ! no bump
            if (impl.eq.0) then ! Carpenter's stencil
              vl(1,:,ndof) = -( gg2 * vl(2,:,ndof)  + &
                                gg3 * vl(3,:,ndof)  + &
                                gg4 * vl(4,:,ndof)  + &
                                gg5 * vl(5,:,ndof)  + &
                                gg6 * vl(6,:,ndof)  ) / gg1
            else
              vl(1,:,ndof) = -( gc2 * vl(2,:,ndof)  + &
                                gc3 * vl(3,:,ndof)  + &
                                gc4 * vl(4,:,ndof)  + &
                                gc5 * vl(5,:,ndof)  ) / gc1
            end if
          end if
        end if

        else    ! inviscid:  zero the wall-normal velocity

          ub(:) = bnb(:,2) * vl(1,:,2) - bnb(:,1) * vl(1,:,3)
!         vb(:) = bnb(:,1) * vl(1,:,2) + bnb(:,2) * vl(1,:,3)
          vb(:) = zero
          
          vl(1,:,2) =  bnb(:,2) * ub + bnb(:,1) * vb
          vl(1,:,3) = -bnb(:,1) * ub + bnb(:,2) * vb

        end if  ! Navier
!=============================================================================!
!       T o p   B o u n d a r y
!=============================================================================!

!.... apply zero disturbance boundary condition

        if (top.eq.0) then
        
          vl(ny,:,:) = zero

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

          rhom = vml(ny,:,1)
          tm   = vml(ny,:,5)
          cm   = sqrt( tm ) / Ma
          um   = vml(ny,:,2)

          d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

          do i = 1, nx
            kk = omega / (cm(i)+um(i))
            a  = omega**2 * d / (cm(i)+um(i))**3
            c3(i) = exp( -a * (x(ny,i) - x0) ) * exp( im * kk * x(ny,i) )
!           c3(i) = wamp(i) * exp( im * kk * x(ny,i) )
          end do

          vl(ny,:,1) = pt5 * c3 / cm**2
          vl(ny,:,2) = c3 * pt5 / ( rhom * cm )
          vl(ny,:,3) = zero
          vl(ny,:,4) = zero
          vl(ny,:,5) = (gamma*Ma**2 * c3 * pt5 - tm * pt5 * c3 / cm**2) / rhom

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

          vl(:,1,1) = zero
          vl(:,1,2) = zero
          vl(:,1,3) = zero
          vl(:,1,4) = zero
          vl(:,1,5) = zero
        
        else if (left.eq.4) then        ! eigenfunction inflow disturbance

          allocate( rhom(ny), tm(ny), cm(ny) )
          allocate( c1(ny), c2(ny), c3(ny), c4(ny) )
          allocate( rho(ny), u1(ny), u2(ny), u3(ny), t(ny), p(ny) )

          rhom = vml(:,1,1)
          tm   = vml(:,1,5)
          cm   = one / Ma * sqrt( tm )

!.... compute incomming characteristics (eigenfunctions)

          rho = cmplx(rhor(1:ny), rhoi(1:ny))
          u1  = cmplx(  ur(1:ny),   ui(1:ny))
          u2  = cmplx(  vr(1:ny),   vi(1:ny))
          u3  = cmplx(  wr(1:ny),   wi(1:ny))
          t   = cmplx(  tr(1:ny),   ti(1:ny))

!         p   = ( rhom * t + tm * rho ) / (gamma * Ma**2)
                                  
!         c1  = -cm**2 * rho + p                        ! entropy
!         c2  =  rhom * cm * u2                         ! vorticity
!         c3  =  rhom * cm * u1 + p                     ! right acoustic
!         c4  =  0                                      ! left acoustic

!.... update the boundary values

!         rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
!         u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
!         u2  = c2 / ( rhom * cm )
!         p   = ( c3 + c4 ) * pt5
!         t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          vl(:,1,1) = rho
          vl(:,1,2) = u1
          vl(:,1,3) = u2
          vl(:,1,4) = u3
          vl(:,1,5) = t

          deallocate( rhom, tm, cm )
          deallocate( c1, c2, c3, c4 )
          deallocate( rho, u1, u2, u3, t, p )

        else if (left.eq.5) then        ! acoustic inflow disturbance
          
          allocate( rhom(ny), tm(ny), cm(ny), um(ny), c3(ny) )

          rhom = vml(:,1,1)
          tm   = vml(:,1,5)
          cm   = sqrt( tm ) / Ma
          um   = vml(:,1,2)

          do j = 1, ny
            kk = omega / ( cm(j)+um(j) )
            c3(j) = exp( im * kk * x(j,1) )
          end do

          vl(:,1,1) = pt5 * c3 / cm**2
          vl(:,1,2) = c3 * pt5 / ( rhom * cm )
          vl(:,1,3) = zero
          vl(:,1,4) = zero
          vl(:,1,5) = (gamma*Ma**2 * c3 * pt5 - tm * pt5 * c3 / cm**2) / rhom

          deallocate( rhom, tm, cm, um, c3 )

        else if (left.eq.7) then        ! symmetry boundary

          vl(:,1,3) = zero

        end if
!=============================================================================!
!       R i g h t   B o u n d a r y
!=============================================================================!

        if (right.eq.0) then            ! zero disturbance BC

          vl(:,nx,1) = zero
          vl(:,nx,2) = zero
          vl(:,nx,3) = zero
          vl(:,nx,4) = zero
          vl(:,nx,5) = zero

        else if (right.eq.4) then

          ac = (2.2804739410500E-001,-6.5163146761218E-003)
!         ac = (-2.8831962908130E-001,-1.3854663671636E-002)

          allocate( rho(ny), u1(ny), u2(ny), u3(ny), t(ny), p(ny) )

!.... compute incomming characteristics (eigenfunctions)

          rho = cmplx(rhor(:), rhoi(:)) * exp(im * ac * x(:,nx))
          u1  = cmplx(  ur(:),   ui(:)) * exp(im * ac * x(:,nx))
          u2  = cmplx(  vr(:),   vi(:)) * exp(im * ac * x(:,nx))
          u3  = cmplx(  wr(:),   wi(:)) * exp(im * ac * x(:,nx))
          t   = cmplx(  tr(:),   ti(:)) * exp(im * ac * x(:,nx))

          vl(:,nx,1) = rho
          vl(:,nx,2) = u1
          vl(:,nx,3) = u2
          vl(:,nx,4) = u3
          vl(:,nx,5) = t

          deallocate( rho, u1, u2, u3, t, p )

        else if (right.eq.5) then

!.... compute the characteristic amplitudes on the right boundary
          
          allocate( rhom(ny), tm(ny), cm(ny), um(ny), c3(ny) )

          rhom = vml(:,nx,1)
          tm   = vml(:,nx,5)
          cm   = sqrt( tm ) / Ma
          um   = vml(:,nx,2)

          do j = 1, ny
            kk = omega / ( cm(j)+um(j) )
            c3(j) = exp( im * kk * x(j,nx) )
          end do

          vl(:,nx,1) = pt5 * c3 / cm**2
          vl(:,nx,2) = c3 * pt5 / ( rhom * cm )
          vl(:,nx,3) = zero
          vl(:,nx,4) = zero
          vl(:,nx,5) = (gamma*Ma**2 * c3 * pt5 - tm * pt5 * c3 / cm**2) / rhom

          deallocate( rhom, tm, cm, um, c3 )

        else if (right .eq. 9) then      ! extrapolation 

          vl(:,nx,:) = exp( two * log(vl(:,nx-1,:)) - log(vl(:,nx-2,:)) )

        end if
        
        end if                  ! xper

        return
        end

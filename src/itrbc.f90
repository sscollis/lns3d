!=============================================================================!
        subroutine itrBC(vl,vml) 
!  
!  Satisfy the boundary conditions 
!  
!  Partially revised for index change (nonlinear only)
!=============================================================================!
        use global
        use stencil
        use pot
        implicit none

        real :: vl(ndof,nx,ny), vml(ndof,nx,ny)
        real :: amp, tmp
        
        real, allocatable :: rho(:),  u1(:),  u2(:), u3(:), t(:), p(:)
        real, allocatable :: rhom(:), cm(:), tm(:), um(:), pnorm(:)
        real, allocatable :: c1(:),   c2(:), c3(:), c4(:)
        real, allocatable :: m1l(:),  m2l(:), dl(:), ul(:)

        real :: ub(nx), vb(nx)

        integer :: i, j, ij, jbl

!.... forcing parameters

        real :: a, d, kk, ra, rb

        real, external :: ramp
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (linear.eq.1) then

        if (yper) then
        
          vl(:,:,ny) = vl(:,:,1)
          
        else                    ! yper
        
          allocate( rhom(nx) ,tm(nx), cm(nx), um(nx), rho(nx), u1(nx), &
                    u2(nx), u3(nx), t(nx), p(nx), c1(nx), c2(nx),      &
                    c3(nx), c4(nx), m1l(nx), m2l(nx), dl(nx), ul(nx)   )
!=============================================================================!
!       W a l l 
!=============================================================================!

        if (Navier) then

!.... set drho/dy = 0 which implies dp/deta = 0

        if (wall.eq.1) then
          vl(1,:,1) = -( gc2 * vl(1,:,2)  + &
                         gc3 * vl(1,:,3)  + &
                         gc4 * vl(1,:,4)  + &
                         gc5 * vl(1,:,5)  ) / gc1
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          vl(1,:,1) = two * vl(1,:,2) - vl(1,:,3)
        end if

!.... no slip wall

        vl(2:ndof-1,:,1) = zero

!.... isothermal wall

        if (wallt.eq.0) vl(ndof,:,1) = zero

!.... adiabatic temperature ( dT/deta = 0 )

        if (wallt.eq.1) then

          if (carp) then        ! Carpenter's stencil
            vl(ndof,:,1) = -( gg2 * vl(ndof,:,2)  + &
                              gg3 * vl(ndof,:,3)  + &
                              gg4 * vl(ndof,:,4)  + &
                              gg5 * vl(ndof,:,5)  + &
                              gg6 * vl(ndof,:,6)  ) / gg1
          else
            vl(ndof,:,1) = -( gc2 * vl(ndof,:,2)  + &
                              gc3 * vl(ndof,:,3)  + &
                              gc4 * vl(ndof,:,4)  + &
                              gc5 * vl(ndof,:,5)  ) / gc1
          end if
        
        end if          ! wallt
        
        else            ! inviscid wall   
        
!         allocate( pnorm(nx) )
!         call lwallbc(vl, vml, rho, pnorm)
!         vl(1,:,1) = rho(:)
!         deallocate( pnorm )

!         do i = 1, nx
!           if ( abs(bnb(i,2)) .ge. abs(bnb(i,1)) ) then
!             vl(3,i,1) = -bnb(i,1) * vl(2,i,1) / bnb(i,2)
!           else
!             vl(2,i,1) = -bnb(i,2) * vl(3,i,1) / bnb(i,1)
!           end if
!         end do

          vl(3,:,1) = zero    ! assume a flat plate

        end if          ! Navier
!=============================================================================!
!       T o p   B o u n d a r y
!=============================================================================!

!.... apply zero disturbance boundary condition

        if (top.eq.0) then
        
          vl(:,:,ny) = zero
        
        else if (top.eq.1) then

!.... allocate room for the boundary amplitude

          if (useAmp) then
            if (istep.eq.0) then
              write(*,*) "Read ",amp_file
              allocate( wamp(nx) )
              open(66,file=amp_file,form='formatted',status='unknown')
              do i = 1, nx
                read(66,*) tmp, wamp(i)
              end do
              close(66)
            end if
          else
            if (istep.eq.0) then
              allocate( wamp(nx) )
              wamp = one
            end if
          endif
        
!.... compute the characteristic amplitudes on the top boundary
          
          rhom = vml(1,:,ny)
          um   = vml(2,:,ny)
          tm   = vml(5,:,ny)
          cm   = one / Ma * sqrt( tm )

!.... get the metrics along this boundary

          j = ny
          do i = 1, nx
            m1l(i)  = n1(i,j)
            m2l(i)  = n2(i,j)
          end do          
          
          dl = sqrt( m1l**2 + m2l**2)
          ul = vml(2,:,ny) * m1l + vml(3,:,ny) * m2l
          
          rho = vl(1,:,ny)
          u1  = vl(2,:,ny)
          u2  = vl(3,:,ny)
          u3  = vl(4,:,ny)
          t   = vl(5,:,ny)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

#if 0          

!.... subtract off the forcing wave

          if (.true.) then
          
            c1  = zero
            c2  = zero
          
!.... sin wave (accounting for viscosity)

            d = one / Re * pt5 / one * ( onept33 + gamma1 / Pr )
  
            do i = 1, nx
              kk = omega / (cm(i)+um(i))
              a  = omega**2 * d / (cm(i)+um(i))**3
              c3(i) = ramp((kk*(x0-x(i,ny)) + omega*time)/(two*pi)) * &
                      wamp(i) * cos( kk * x(i,ny) - omega * time ) * &
                      exp( -a * (x(i,ny) - x0) )
!             c3(i) = ramp((kk*(x0-x(i,ny)) + omega*time)/(two*pi)) * &
!                     wamp(i) * cos( kk * x(i,ny) - omega * time )
            end do

            c4  = zero
          
            rho = rho - ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
            u1  = u1 - ( c3 - c4 ) * pt5 / ( rhom * cm )
            u2  = u2 -  one / ( rhom * cm ) * c2
            p   = p - ( c3 + c4 ) * pt5
            t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          end if
          
!.... absorb boundary normal waves

          c1  = -cm**2 * rho + p
          c2  =  rhom * cm * ( -m2l * u1 + m1l * u2 )
          c3  =  rhom * cm * ( m1l/dl * u1 + m2l/dl * u2 ) + p
          c4  = -rhom * cm * ( m1l/dl * u1 + m2l/dl * u2 ) + p

          do i = 1, nx
            if (ul(i) .le. zero) then
              c1(i) = zero
              c2(i) = zero
            end if
            if (ul(i) + dl(i) * cm(i) .le. zero) c3(i) = zero
            if (ul(i) - dl(i) * cm(i) .le. zero) c4(i) = zero
          end do

          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = -m1l * m2l / ( dl**2 * rhom * cm) * c2 +        &
                ( c3 - c4 ) * m1l / dl * pt5 / ( rhom * cm )
          u2  =  m1l * m1l / ( dl**2 * rhom * cm ) * c2 +       &
                ( c3 - c4 ) * m2l / dl * pt5 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

!.... add back on the forcing wave

          if (.true.) then
          
            c1  = zero
            c2  = zero

!.... sin wave (accounting for viscosity)

            d = one / Re * pt5 / one * ( onept33 + gamma1 / Pr )
  
            do i = 1, nx
              kk = omega / (cm(i)+um(i))
              a  = omega**2 * d / (cm(i)+um(i))**3
              c3(i) = ramp((kk*(x0-x(i,ny)) + omega*time)/(two*pi)) * &
                      wamp(i) * cos( kk * x(i,ny) - omega * time ) * &
                      exp( -a * (x(i,ny) - x0) )
!             c3(i) = ramp((kk*(x0-x(i,ny)) + omega*time)/(two*pi)) * &
!                     wamp(i) * cos( kk * x(i,ny) - omega * time )
            end do
            
            c4  = zero
            
            rho = rho + ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
            u1  = u1 + ( c3 - c4 ) * pt5 / ( rhom * cm )
            u2  = u2 +  one / ( rhom * cm ) * c2
            p   = p + ( c3 + c4 ) * pt5
            t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          end if
#else
            c1  = zero                                    ! entropy
            c2  = zero                                    ! vorticity

!.... sin wave (accounting for viscosity)

            d = one / Re * pt5 / one * ( onept33 + gamma1 / Pr )

            do i = 1, nx
              kk = omega / (cm(i)+um(i))
              a  = omega**2 * d / (cm(i)+um(i))**3
              c3(i) = ramp((kk*(x0-x(i,ny)) + omega*time)/(two*pi)) * &
                      wamp(i) * cos( kk * x(i,ny) - omega * time ) * &
                      exp( -a * (x(i,ny) - x0) )
            end do

!           c4  = -rhom * cm * u1 + p                     ! left acoustic
            c4  = zero                                    ! left acoustic

            rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
            u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
            u2  = c2 / ( rhom * cm )
            u3  = zero
            p   = ( c3 + c4 ) * pt5
            t   = ( gamma * Ma**2 * p - tm * rho ) / rhom
#endif
          
!.... apply the boundary condition

          vl(1,:,ny) = rho
          vl(2,:,ny) = u1
          vl(3,:,ny) = u2
          vl(4,:,ny) = u3
          vl(5,:,ny) = t

        end if

        deallocate( rhom ,tm, cm, um, rho, u1, &
                    u2, u3, t, p, c1, c2,      &
                    c3, c4, m1l, m2l, dl, ul   )

        end if                  ! yper

        if (xper) then

          vl(:,nx,:) = vl(:,1,:)
          
        else                    ! xper

          allocate( rhom(ny) ,tm(ny), cm(ny), um(ny), rho(ny), u1(ny), &
                    u2(ny), u3(ny), t(ny), p(ny), c1(ny), c2(ny),      &
                    c3(ny), c4(ny), m1l(ny), m2l(ny), dl(ny), ul(ny)   )

!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!

        if (left.eq.0) then             ! zero disturbance BC

          vl(1,1,1:ny) = zero
          vl(2,1,1:ny) = zero
          vl(3,1,1:ny) = zero
          vl(4,1,1:ny) = zero
          vl(5,1,1:ny) = zero
        
        else if (left.eq.1) then        ! nonreflecting boundary

          rhom = vml(1,1,:)
          tm   = vml(5,1,:)
          cm   = sqrt( tm ) / Ma
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(2,1,:) * m1l + vml(3,1,:) * m2l
            
          rho = vl(1,1,:)
          u1  = vl(2,1,:)
          u2  = vl(3,1,:)
          u3  = vl(4,1,:)
          t   = vl(5,1,:)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

!.... set the boundary-layer edge

!         nbl = 40              ! set for the Stokes layer for omega=0.8
!         nbl = 1               ! set for the viscous layer on small grid

          vl(1,1,nbl+1:ny) = (gamma * Ma**2 * rhom(nbl+1:ny) * &
                              cm(nbl+1:ny) * u1(nbl+1:ny) - &
                              rhom(nbl+1:ny) * t(nbl+1:ny))/tm(nbl+1:ny)

!.... zeroth-order extrapolation in the viscous layer

!         vl(:,1,1:nbl) = vl(:,2,1:nbl)

        else if (left.eq.2) then           !.... symmetry boundary

!         vl(1,1,:) = -( gc2 * vl(1,2,:)  + &
!                        gc3 * vl(1,3,:)  + &
!                        gc4 * vl(1,4,:)  + &
!                        gc5 * vl(1,5,:)  ) / gc1

          vl(2,1,:) = -( gc2 * vl(2,2,:)  + &
                         gc3 * vl(2,3,:)  + &
                         gc4 * vl(2,4,:)  + &
                         gc5 * vl(2,5,:)  ) / gc1
                         
          vl(3,1,:) = zero

          vl(4,1,:) = -( gc2 * vl(4,2,:)  + &
                         gc3 * vl(4,3,:)  + &
                         gc4 * vl(4,4,:)  + &
                         gc5 * vl(4,5,:)  ) / gc1
                         
          vl(5,1,:) = -( gc2 * vl(5,2,:)  + &
                         gc3 * vl(5,3,:)  + &
                         gc4 * vl(5,4,:)  + &
                         gc5 * vl(5,5,:)  ) / gc1

        else if (left.eq.3) then           ! Giles nonreflecting boundary
        
          rhom = vml(1,1,:)
          tm   = vml(5,1,:)
          cm   = one / Ma * sqrt( tm )
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(2,1,:) * m1l + vml(3,1,:) * m2l
            
          rho = vl(1,1,:)
          u1  = vl(2,1,:)
          u2  = vl(3,1,:)
          u3  = vl(4,1,:)
          t   = vl(5,1,:)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )
          
          c1  = -cm**2 * rho + p
          c2  =  rhom * cm * ( -m2l * u1 + m1l * u2 )
          c3  =  rhom * cm * ( m1l/dl * u1 + m2l/dl * u2 ) + p
          c4  = -rhom * cm * ( m1l/dl * u1 + m2l/dl * u2 ) + p

          do j = 1, ny
            if (ul(j) .le. zero) then
              c1(j) = zero
              c2(j) = zero
            end if
            if (ul(j) + dl(j) * cm(j) .le. zero) c3(j) = zero
            if (ul(j) - dl(j) * cm(j) .le. zero) c4(j) = zero
          end do

          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = -m1l * m2l / ( dl**2 * rhom * cm) * c2 +        &
                ( c3 - c4 ) * m1l / dl * pt5 / ( rhom * cm )
          u2  =  m1l * m1l / ( dl**2 * rhom * cm ) * c2 +       &
                ( c3 - c4 ) * m2l / dl * pt5 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          do j = 1, ny

!.... set for the MSE scattering problem

!           ra = ramp( (y(j,1) - y(50,1)) / (y(70,1) - y(50,1)) )

!.... set for the MSE receptivity problem

!           ra = ramp( (y(j,1) - y(64,1)) / (y(84,1) - y(64,1)) )

            ra = one
            rb = one - ra

            vl(1,1,j) = rb * vl(1,1,j) + ra * rho(j)
            vl(2,1,j) = rb * vl(2,1,j) + ra * u1(j)
            vl(3,1,j) = rb * vl(3,1,j) + ra * u2(j)
            vl(4,1,j) = rb * vl(4,1,j) + ra * u3(j)
            vl(5,1,j) = rb * vl(5,1,j) + ra * t(j)

          end do

        else if (left.eq.4) then        ! eigenfunction inflow disturbance

          rhom = vml(1,1,:)
          tm   = vml(5,1,:)
          cm   = one / Ma * sqrt( tm )
          
          amp = ramp( (omega * time) / (two*pi) )

!.... compute incomming characteristics (eigenfunctions)

          rho = amp * (rhor(1:ny)*cos(omega*time) + &
                       rhoi(1:ny)*sin(omega*time))
          u1  = amp * (  ur(1:ny)*cos(omega*time) + &
                         ui(1:ny)*sin(omega*time))
          u2  = amp * (  vr(1:ny)*cos(omega*time) + &
                         vi(1:ny)*sin(omega*time))
          u3  = amp * (  wr(1:ny)*cos(omega*time) + &
                         wi(1:ny)*sin(omega*time))
          t   = amp * (  tr(1:ny)*cos(omega*time) + &
                         ti(1:ny)*sin(omega*time))

          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )
                                  
          c1  = -cm**2 * rho + p                        ! entropy
          c2  =  rhom * cm * u2                         ! vorticity
          c3  =  rhom * cm * u1 + p                     ! right acoustic

!.... compute outgoing characteristics

          rho = vl(1,1,:)
          u1  = vl(2,1,:)
          u2  = vl(3,1,:)
          u3  = vl(4,1,:)
          t   = vl(5,1,:)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

          c4  =  -rhom * cm * u1 + p                    ! left acoustic

!.... update the boundary values

          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
          u2  = c2 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          vl(1,1,:) = rho
          vl(2,1,:) = u1
          vl(3,1,:) = u2
          vl(4,1,:) = u3
          vl(5,1,:) = t

        else if (left.eq.5) then        ! acoustic inflow disturbance
        
          rhom = vml(1,1,:)
          um   = vml(2,1,:)
          tm   = vml(5,1,:)
          cm   = one / Ma * sqrt( tm )
          
          rho = vl(1,1,:)
          u1  = vl(2,1,:)
          u2  = vl(3,1,:)
          u3  = vl(4,1,:)
          t   = vl(5,1,:)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

          c1  = zero                                    ! entropy
          c2  = zero                                    ! vorticity

          d = one / Re * pt5 / one * ( onept33 + gamma1 / Pr )

          do j = 1, ny
            kk = omega / (cm(j)+um(j))
            a  = omega**2 * d / (cm(j)+um(j))**3
            c3(j) = ramp( (kk*(x0-x(1,j)) + omega * time)/(two*pi) ) * &
                    cos( kk * x(1,j) - omega * time )
!           c3(j) = cos( kk * x(1,j) - omega * time )
          end do

          c4  = -rhom * cm * u1 + p                     ! left acoustic
        
          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
          u2  = c2 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          vl(1,1,:) = rho
          vl(2,1,:) = u1
          vl(3,1,:) = u2
          vl(4,1,:) = u3
          vl(5,1,:) = t

        else if (left.eq.7) then           ! symmetry boundary

          vl(3,1,:) = zero

        end if
!=============================================================================!
!       R i g h t   B o u n d a r y
!=============================================================================!

        if (right.eq.0) then            ! zero disturbance BC

          vl(1,nx,1:ny) = zero
          vl(2,nx,1:ny) = zero
          vl(3,nx,1:ny) = zero
          vl(4,nx,1:ny) = zero
          vl(5,nx,1:ny) = zero

        else if (right.eq.1) then       ! nonreflecting boundary

          rhom = vml(1,nx,:)
          tm   = vml(5,nx,:)
          cm   = one / Ma * sqrt( tm )
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(2,nx,:) * m1l + vml(3,nx,:) * m2l
          
          rho = vl(1,nx,:)
          u1  = vl(2,nx,:)
          u2  = vl(3,nx,:)
          u3  = vl(4,nx,:)
          t   = vl(5,nx,:)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )
          
          vl(1,nx,nbl+1:ny) = (gamma * Ma**2 * rhom(nbl+1:ny) * &
                               cm(nbl+1:ny) * u1(nbl+1:ny) - &
                               rhom(nbl+1:ny) * t(nbl+1:ny)) / &
                               tm(nbl+1:ny)

!.... zeroth-order extrapolation in the viscous layer

!         vl(:,nx,1:nbl) = vl(:,nx-1,1:nbl)

        else if (right.eq.3) then       ! Giles' nonreflecting boundary

          rhom = vml(1,nx,:)
          tm   = vml(5,nx,:)
          cm   = one / Ma * sqrt( tm )
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(2,nx,:) * m1l + vml(3,nx,:) * m2l
          
          rho = vl(1,nx,:)
          u1  = vl(2,nx,:)
          u2  = vl(3,nx,:)
          u3  = vl(4,nx,:)
          t   = vl(5,nx,:)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )
          
          c1  = -cm**2 * rho + p
          c2  =  rhom * cm * ( -m2l * u1 + m1l * u2 )
          c3  =  rhom * cm * ( m1l/dl * u1 + m2l/dl * u2 ) + p
          c4  = -rhom * cm * ( m1l/dl * u1 + m2l/dl * u2 ) + p

          do j = 1, ny
            if (ul(j) .le. zero) then
              c1(j) = zero
              c2(j) = zero
            end if
            if (ul(j) + dl(j) * cm(j) .le. zero) c3(j) = zero
            if (ul(j) - dl(j) * cm(j) .le. zero) c4(j) = zero
          end do

          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = -m1l * m2l / ( dl**2 * rhom * cm) * c2 +        &
                ( c3 - c4 ) * m1l / dl * pt5 / ( rhom * cm )
          u2  =  m1l * m1l / ( dl**2 * rhom * cm ) * c2 +       &
                ( c3 - c4 ) * m2l / dl * pt5 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          do j = 1, ny
            vl(1,nx,j) = rho(j)
            vl(2,nx,j) = u1(j)
            vl(3,nx,j) = u2(j)
            vl(4,nx,j) = u3(j)
            vl(5,nx,j) = t(j)
          end do

        end if
        
        deallocate( rhom ,tm, cm, um, rho, u1, &
                    u2, u3, t, p, c1, c2,      &
                    c3, c4, m1l, m2l, dl, ul   )

        end if                  ! xper

        else                    ! linear = 0

!=============================================================================!
!       N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!

        if (yper) then
        
          vl(:,:,ny) = vl(:,:,1)
          
        else                    ! yper
        
!=============================================================================!
!       W a l l
!=============================================================================!

        if (Navier) then
        
!.... no slip wall

          vl(2:ndof-1,:,1)  = zero

!.... set drho/dy = 0 which implies dp/deta = 0

          if (wall.eq.1) then
            vl(1,:,1) = -( gc2 * vl(1,:,2)  + &
                           gc3 * vl(1,:,3)  + &
                           gc4 * vl(1,:,4)  + &
                           gc5 * vl(1,:,5)  ) / gc1
          end if

!.... Wall normal momentum equation to get rho at the wall 

          if (wall.eq.2) then
            allocate( p(nx), pnorm(nx) )
            call wallbc(vl, p, pnorm)
            vl(1,:,1) = p(:)                    ! p is really rho here
            deallocate( p, pnorm )
          end if

!.... Third order extrapolation to the wall for density
!.... (in computational space)

          if (wall.eq.3) then
            vl(1,:,1) = 4.0*vl(1,:,2) - 6.0*vl(1,:,3) + &
                        4.0*vl(1,:,4) - vl(1,:,5)
          end if

!.... isothermal wall:  set to adiabatic wall temperature where 
!.... the recovery factor is r = \sqrt(Pr).  This means that
!.... there will be non-zero heat transfer near the leading
!.... edge and immediately downstream (e.g. on the parabolic
!.... cylinder).
!.... NOTE:  This is valid for a laminar boundary layer on a 
!....        flat plate for 0.5 < Pr < 10.  Near a stagnation
!....        point (i.e. the leading-edge), r = 1.  For a 
!....        turbulent boundary layer, r = Pr^(1/3) for Pr
!....        near 1.0. See: www.thermopedia.com/content/291/
!....        and references therein.

          if (wallt.eq.0) then
            vl(ndof,:,1) = one + pt5 * gamma1 * Ma**2 * &
                           (one + tan(theta)**2) * sqrt(Pr)
          end if

!.... adiabatic temperature ( dT/deta = 0 )

          if (wallt.eq.1) then
            if (carp) then      ! Carpenter's stencil
              vl(ndof,:,1) = -( gg2 * vl(ndof,:,2)  + &
                                gg3 * vl(ndof,:,3)  + &
                                gg4 * vl(ndof,:,4)  + &
                                gg5 * vl(ndof,:,5)  + &
                                gg6 * vl(ndof,:,6)  ) / gg1
            else
              vl(ndof,:,1) = -( gc2 * vl(ndof,:,2)  + &
                                gc3 * vl(ndof,:,3)  + &
                                gc4 * vl(ndof,:,4)  + &
                                gc5 * vl(ndof,:,5)  ) / gc1
            end if
          end if
          
        else     ! inviscid:  zero the wall-normal velocity
        
!          if (wall.ne.4) then
            ub(:) = bnb(:,2) * vl(2,:,1) - bnb(:,1) * vl(3,:,1)
!           vb(:) = bnb(:,1) * vl(2,:,1) + bnb(:,2) * vl(3,:,1)
            vb(:) = zero

            vl(2,:,1) =  bnb(:,2) * ub + bnb(:,1) * vb
            vl(3,:,1) = -bnb(:,1) * ub + bnb(:,2) * vb
!          end if

          if (wall.eq.3) then
            do i=1, nx
              vl(1,i,1)=4.*vl(1,i,2)-6.*vl(1,i,3)+4.*vl(1,i,4)-vl(1,i,5)
              vl(5,i,1)=4.*vl(5,i,2)-6.*vl(5,i,3)+4.*vl(5,i,4)-vl(5,i,5)
            enddo
          end if

        end if   ! Navier
        
!=============================================================================!
!       T o p
!=============================================================================!
        if (Ma.gt.one) then

          vl(1,:,ny) = one
          vl(2,:,ny) = one
          vl(3,:,ny) = zero
          vl(4,:,ny) = zero
          vl(5,:,ny) = one

        else

          if (top.eq.0) then
            call Reimann( nx, vl(:,:,ny), vl(:,:,ny-1), vl(:,:,ny-2), &
                          bnt, x(:,ny), y(:,ny), rhobt, ubt, vbt, wbt, &
                          tbt, pbt, cbt )
          end if

          if (top.eq.3) then
            vl(1,:,ny) = one
            vl(2,:,ny) = zero
            vl(3,:,ny) = zero
            vl(4,:,ny) = zero
            vl(5,:,ny) = one
          end if          

        end if                  ! Ma
        
        end if                  ! yper

        if (xper) then

          vl(:,nx,:) = vl(:,1,:)
          
        else                    ! xper

!=============================================================================!
!       L e f t  B o u n d a r y
!=============================================================================!

!.... first-order Riemann

        if (left.eq.0) then
          call Reimann( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                        bnl, x(1,:), y(1,:), rhobl, ubl, vbl, wbl, &
                        tbl, pbl, cbl )
        end if

!.... Symmetry plane

        if (left.eq.2) then
          vl(1,1,:) = -( gc2 * vl(1,2,:)  + &
                         gc3 * vl(1,3,:)  + &
                         gc4 * vl(1,4,:)  + &
                         gc5 * vl(1,5,:)  ) / gc1

          vl(2,1,:) = -( gc2 * vl(2,2,:)  + &
                         gc3 * vl(2,3,:)  + &
                         gc4 * vl(2,4,:)  + &
                         gc5 * vl(2,5,:)  ) / gc1
                         
          vl(3,1,:) = zero

          vl(4,1,:) = -( gc2 * vl(4,2,:)  + &
                         gc3 * vl(4,3,:)  + &
                         gc4 * vl(4,4,:)  + &
                         gc5 * vl(4,5,:)  ) / gc1
                         
          vl(5,1,:) = -( gc2 * vl(5,2,:)  + &
                         gc3 * vl(5,3,:)  + &
                         gc4 * vl(5,4,:)  + &
                         gc5 * vl(5,5,:)  ) / gc1
        end if
        
        if (left.eq.7) then           ! symmetry boundary
          vl(3,1,:) = zero
        end if
        
        if (.false.) then
        
!.... zeroth order extrapolation in the viscous layers

        if (extrap.eq.0) then
          do j = 1, nbl
            vl(:,1,j) = vl(:,2,j)
          end do
        end if
        
!.... first-order extrapolation in the viscous layers

        if (extrap.eq.1) then
          do j = 1, nbl
            vl(:,1,j) = two * vl(:,2,j) - vl(:,3,j)
          end do
        end if
        
!.... second-order extrapolation in the viscous layers

        if (extrap.eq.2) then
          do j = 1, nbl
            vl(:,1,j) = three * vl(:,2,j) - three * vl(:,3,j) + vl(:,4,j)
          end do
        end if                  ! extrapolation type
        
        end if
!=============================================================================!
!       R i g h t  B o u n d a r y
!=============================================================================!

!.... First-order Riemann

        if (right.eq.0) then
          call Reimann( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                        bnr, x(nx,:), y(nx,:), rhobr, ubr, vbr, wbr, &
                        tbr, pbr, cbr )
        end if

!.... Symmetry boundary

        if (right.eq.7) then           ! symmetry boundary
          vl(3,nx,:) = zero
        end if
        
!.... find the boundary-layer edge (using delta_99 right now)
!.... It is wise to keep nbl fixed, since every change sets up a transient.
!       nbl = 70                ! for ny = 128
!       nbl = 35                ! for ny = 63
        if (update_nbl) then
          do jbl = 1, ny
            if (vl(2,nx-1,jbl) .gt. 0.99 ) exit
            !if ((vl(5,nx-1,jbl)-vl(5,nx-1,1))/(one-vl(5,nx-1,1)).gt.0.99) exit
          end do
          if (iter.eq.1.and.output_nbl) write(90,*) istep, lstep, jbl, nbl 
          nbl = jbl
        end if
        if (nbl.gt.ny) call error("itrbc$","nbl > ny$")

!.... zeroth order extrapolation in the viscous layers

        if (extrap.eq.0) then
          do j = 1, nbl
            vl(:,nx,j) = vl(:,nx-1,j)
          end do
        
!.... first-order extrapolation in the viscous layers

        else if (extrap.eq.1) then
          do j = 1, nbl
            vl(:,nx,j) = two * vl(:,nx-1,j) - vl(:,nx-2,j)
          end do

!.... second-order extrapolation in the viscous layers

        else if (extrap.eq.2) then
          do j = 1, nbl
            vl(:,nx,j) = three*vl(:,nx-1,j) - three*vl(:,nx-2,j) + vl(:,nx-3,j)
          end do
        end if                  ! extrapolation type

        end if                  ! xper

        end if                  ! linear
        
        return
        end

!=============================================================================!
        function ramp(t) 
!  
!  Ramp up the forcing amplitude from [0,1] over the t range [0,1]
!  
!=============================================================================!
        implicit none

        real :: ramp, t
!=============================================================================!
        if (t .le. 0.0) then
          ramp = 0.0
        else if (t .lt. 1.0) then
          ramp = 0.5 * ( 1.0 + tanh( 4.0 * ( 2.0 * t - 1.0 ) ) )
        else
          ramp = 1.0
        end if

        return
        end

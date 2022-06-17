!=============================================================================!
        subroutine itrBC(vl,vml) 
!  
!  Satisfy the boundary conditions 
!  
!=============================================================================!
        use global
        use stencil
        use pot
        implicit none

        real :: vl(ny,nx,ndof), vml(ny,nx,ndof)
        real :: amp, tmp
        
        real, allocatable :: rho(:),  u1(:),  u2(:), u3(:), t(:), p(:)
        real, allocatable :: rhom(:), cm(:), tm(:), um(:), pnorm(:)
        real, allocatable :: c1(:),   c2(:), c3(:), c4(:)
        real, allocatable :: m1l(:),  m2l(:), dl(:), ul(:)

        real :: ub(nx), vb(nx)

        integer :: i, j, ij

!.... forcing parameters

        real :: a, d, kk, ra, rb

        real, external :: ramp

!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!

        if (linear.eq.1) then

        if (yper) then
        
          vl(ny,:,:) = vl(1,:,:)
          
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

        vl(1,:,2:ndof-1) = zero

!.... isothermal wall

        if (wallt.eq.0) vl(1,:,ndof) = zero

!.... adiabatic temperature ( dT/deta = 0 )

        if (wallt.eq.1) then

          if (impl.eq.0) then   ! Carpenter's stencil
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
        
        end if          ! wallt
        
        else            ! inviscid wall   
        
!         allocate( pnorm(nx) )
!         call lwallbc(vl, vml, rho, pnorm)
!         vl(1,:,1) = rho(:)
!         deallocate( pnorm )

!         do i = 1, nx
!           if ( abs(bnb(i,2)) .ge. abs(bnb(i,1)) ) then
!             vl(1,i,3) = -bnb(i,1) * vl(1,i,2) / bnb(i,2)
!           else
!             vl(1,i,2) = -bnb(i,2) * vl(1,i,3) / bnb(i,1)
!           end if
!         end do

          vl(1,:,3) = zero    ! assume a flat plate

        end if          ! Navier

!=============================================================================!
!       T o p   B o u n d a r y
!=============================================================================!

!.... apply zero disturbance boundary condition

        if (top.eq.0) then
        
          vl(ny,:,:) = zero
        
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
          
          rhom = vml(ny,:,1)
          um   = vml(ny,:,2)
          tm   = vml(ny,:,5)
          cm   = one / Ma * sqrt( tm )

!.... get the metrics along this boundary

          j = ny
          do i = 1, nx
            ij = j + (i-1) * ny
            m1l(i)  = n1(ij)
            m2l(i)  = n2(ij)
          end do          
          
          dl = sqrt( m1l**2 + m2l**2)
          ul = vml(ny,:,2) * m1l + vml(ny,:,3) * m2l
          
          rho = vl(ny,:,1)
          u1  = vl(ny,:,2)
          u2  = vl(ny,:,3)
          u3  = vl(ny,:,4)
          t   = vl(ny,:,5)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )
          
!.... subtract off the forcing wave

          if (.true.) then
          
            c1  = zero
            c2  = zero
          
!.... sin wave (accounting for viscosity)

            d = one / Re * pt5 / one * ( onept33 + gamma1 / Pr )
  
            do i = 1, nx
              kk = omega / (cm(i)+um(i))
              a  = omega**2 * d / (cm(i)+um(i))**3
              c3(i) = ramp( (kk*(x0-x(ny,i)) + omega * time)/(two*pi) ) * &
                      cos( kk * x(ny,i) - omega * time ) * &
                      exp( -a * (x(ny,i) - x0) )
!             c3(i) = ramp( (kk*(x0-x(ny,i)) + omega * time)/(two*pi) ) * &
!                     wamp(i) * cos( kk * x(ny,i) - omega * time )
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
              c3(i) = ramp( (kk*(x0-x(ny,i)) + omega * time)/(two*pi) ) * &
                      cos( kk * x(ny,i) - omega * time ) * &
                      exp( -a * (x(ny,i) - x0) )
!             c3(i) = ramp( (kk*(x0-x(ny,i)) + omega * time)/(two*pi) ) * &
!                     wamp(i) * cos( kk * x(ny,i) - omega * time )
            end do
            
            c4  = zero
            
            rho = rho + ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
            u1  = u1 + ( c3 - c4 ) * pt5 / ( rhom * cm )
            u2  = u2 +  one / ( rhom * cm ) * c2
            p   = p + ( c3 + c4 ) * pt5
            t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          end if
          
!.... apply the boundary condition

          vl(ny,:,1) = rho
          vl(ny,:,2) = u1
          vl(ny,:,3) = u2
          vl(ny,:,4) = u3
          vl(ny,:,5) = t

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

          vl(1:ny,1,1) = zero
          vl(1:ny,1,2) = zero
          vl(1:ny,1,3) = zero
          vl(1:ny,1,4) = zero
          vl(1:ny,1,5) = zero
        
        else if (left.eq.1) then        ! nonreflecting boundary

          rhom = vml(:,1,1)
          tm   = vml(:,1,5)
          cm   = sqrt( tm ) / Ma
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(:,1,2) * m1l + vml(:,1,3) * m2l
            
          rho = vl(:,1,1)
          u1  = vl(:,1,2)
          u2  = vl(:,1,3)
          u3  = vl(:,1,4)
          t   = vl(:,1,5)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

!.... set the boundary-layer edge

!         nbl = 40              ! set for the Stokes layer for omega=0.8
          nbl = 1               ! set for the viscous layer on small grid

          vl(nbl+1:ny,1,1) = (gamma * Ma**2 * rhom(nbl+1:ny) * &
                              cm(nbl+1:ny) * u1(nbl+1:ny) - &
                              rhom(nbl+1:ny) * t(nbl+1:ny)) / tm(nbl+1:ny)

!.... zeroth-order extrapolation in the viscous layer

!         vl(1:nbl,1,:) = vl(1:nbl,2,:)

        else if (left.eq.2) then           !.... symmetry boundary

!         vl(:,1,1) = -( gc2 * vl(:,2,1)  + &
!                        gc3 * vl(:,3,1)  + &
!                        gc4 * vl(:,4,1)  + &
!                        gc5 * vl(:,5,1)  ) / gc1

          vl(:,1,2) = -( gc2 * vl(:,2,2)  + &
                         gc3 * vl(:,3,2)  + &
                         gc4 * vl(:,4,2)  + &
                         gc5 * vl(:,5,2)  ) / gc1
                         
          vl(:,1,3) = zero

          vl(:,1,4) = -( gc2 * vl(:,2,4)  + &
                         gc3 * vl(:,3,4)  + &
                         gc4 * vl(:,4,4)  + &
                         gc5 * vl(:,5,4)  ) / gc1
                         
          vl(:,1,5) = -( gc2 * vl(:,2,5)  + &
                         gc3 * vl(:,3,5)  + &
                         gc4 * vl(:,4,5)  + &
                         gc5 * vl(:,5,5)  ) / gc1

        else if (left.eq.3) then           ! Giles nonreflecting boundary
        
          rhom = vml(:,1,1)
          tm   = vml(:,1,5)
          cm   = one / Ma * sqrt( tm )
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(:,1,2) * m1l + vml(:,1,3) * m2l
            
          rho = vl(:,1,1)
          u1  = vl(:,1,2)
          u2  = vl(:,1,3)
          u3  = vl(:,1,4)
          t   = vl(:,1,5)
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

            vl(j,1,1) = rb * vl(j,1,1) + ra * rho(j)
            vl(j,1,2) = rb * vl(j,1,2) + ra * u1(j)
            vl(j,1,3) = rb * vl(j,1,3) + ra * u2(j)
            vl(j,1,4) = rb * vl(j,1,4) + ra * u3(j)
            vl(j,1,5) = rb * vl(j,1,5) + ra * t(j)

          end do

        else if (left.eq.4) then        ! eigenfunction inflow disturbance

          rhom = vml(:,1,1)
          tm   = vml(:,1,5)
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

          rho = vl(:,1,1)
          u1  = vl(:,1,2)
          u2  = vl(:,1,3)
          u3  = vl(:,1,4)
          t   = vl(:,1,5)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

          c4  =  -rhom * cm * u1 + p                    ! left acoustic

!.... update the boundary values

          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
          u2  = c2 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          vl(:,1,1) = rho
          vl(:,1,2) = u1
          vl(:,1,3) = u2
          vl(:,1,4) = u3
          vl(:,1,5) = t

        else if (left.eq.5) then        ! acoustic inflow disturbance
        
          rhom = vml(:,1,1)
          um   = vml(:,1,2)
          tm   = vml(:,1,5)
          cm   = one / Ma * sqrt( tm )
          
          rho = vl(:,1,1)
          u1  = vl(:,1,2)
          u2  = vl(:,1,3)
          u3  = vl(:,1,4)
          t   = vl(:,1,5)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )

          c1  = zero                                    ! entropy
          c2  = zero                                    ! vorticity

          do j = 1, ny
            kk = omega / (cm(j)+um(j))
            a  = omega**2 * d / (cm(j)+um(j))**3

            c3  = ramp( (kk*(x0-x(j,1)) + omega * time) / (two*pi) ) * &
                  cos( kk * x(j,1) - omega * time )

!           c3  = cos( kk * x(:,1) - omega * time )
          end do

          c4  = -rhom * cm * u1 + p                     ! left acoustic
        
          rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
          u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
          u2  = c2 / ( rhom * cm )
          p   = ( c3 + c4 ) * pt5
          t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

          vl(:,1,1) = rho
          vl(:,1,2) = u1
          vl(:,1,3) = u2
          vl(:,1,4) = u3
          vl(:,1,5) = t
        else if (left.eq.7) then           ! symmetry boundary
          vl(:,1,3) = zero
        end if

!=============================================================================!
!       R i g h t   B o u n d a r y
!=============================================================================!

        if (right.eq.0) then            ! zero disturbance BC

          vl(1:ny,nx,1) = zero
          vl(1:ny,nx,2) = zero
          vl(1:ny,nx,3) = zero
          vl(1:ny,nx,4) = zero
          vl(1:ny,nx,5) = zero

        else if (right.eq.1) then       ! nonreflecting boundary

          rhom = vml(:,nx,1)
          tm   = vml(:,nx,5)
          cm   = one / Ma * sqrt( tm )
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(:,nx,2) * m1l + vml(:,nx,3) * m2l
          
          rho = vl(:,nx,1)
          u1  = vl(:,nx,2)
          u2  = vl(:,nx,3)
          u3  = vl(:,nx,4)
          t   = vl(:,nx,5)
          p   = one/(gamma * Ma**2) * ( rhom * t + tm * rho )
          
          vl(nbl+1:ny,nx,1) = (gamma * Ma**2 * rhom(nbl+1:ny) * &
                               cm(nbl+1:ny) * u1(nbl+1:ny) - &
                               rhom(nbl+1:ny) * t(nbl+1:ny)) / tm(nbl+1:ny)

!.... zeroth-order extrapolation in the viscous layer

!         vl(1:nbl,nx,:) = vl(1:nbl,nx-1,:)

        else if (right.eq.3) then       ! Giles' nonreflecting boundary

          rhom = vml(:,nx,1)
          tm   = vml(:,nx,5)
          cm   = one / Ma * sqrt( tm )
          m1l  = one
          m2l  = zero
          dl   = sqrt( m1l**2 + m2l**2)
          ul   = vml(:,nx,2) * m1l + vml(:,nx,3) * m2l
          
          rho = vl(:,nx,1)
          u1  = vl(:,nx,2)
          u2  = vl(:,nx,3)
          u3  = vl(:,nx,4)
          t   = vl(:,nx,5)
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
            vl(j,nx,1) = rho(j)
            vl(j,nx,2) = u1(j)
            vl(j,nx,3) = u2(j)
            vl(j,nx,4) = u3(j)
            vl(j,nx,5) = t(j)
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

!.... no slip wall

          vl(1,:,2:ndof-1)  = zero

!.... isothermal wall:  set to adiabatic wall temperature where 
!....                   the recovery factor is sqrt(Pr)
!.... NOTE:  This is not correct for Pr .ne. 1

          if (wallt.eq.0) then
            vl(1,:,ndof) = one + pt5 * gamma1 * Ma**2 * &
                           (one + tan(alpha)**2) * sqrt(Pr)
          end if

!.... adiabatic temperature ( dT/deta = 0 )

          if (wallt.eq.1) then
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
          
!.... use the wall normal momentum equation to get rho at the wall 

          if (wall.eq.2) then
            allocate( p(nx), pnorm(nx) )
            call wallbc(vl, p, pnorm)
            vl(1,:,1) = p(:)                    ! p is really rho here
            deallocate( p, pnorm )
          end if

        else     ! inviscid:  zero the wall-normal velocity
        
          ub(:) = bnb(:,2) * vl(1,:,2) - bnb(:,1) * vl(1,:,3)
!         vb(:) = bnb(:,1) * vl(1,:,2) + bnb(:,2) * vl(1,:,3)
          vb(:) = zero

          vl(1,:,2) =  bnb(:,2) * ub + bnb(:,1) * vb
          vl(1,:,3) = -bnb(:,1) * ub + bnb(:,2) * vb

        end if   ! Navier
        
!=============================================================================!
!       T o p
!=============================================================================!

        if (Ma.gt.one) then

          vl(ny,:,1) = one
          vl(ny,:,2) = one
          vl(ny,:,3) = zero
          vl(ny,:,4) = zero
          vl(ny,:,5) = one

        else

          if (top.eq.0) then
            call Reimann( nx, vl(ny,:,:), vl(ny-1,:,:), vl(ny-2,:,:), &
                          bnt, x(ny,:), y(ny,:), rhobt, ubt, vbt, wbt, &
                          tbt, pbt, cbt )
          end if

          if (top.eq.3) then
            vl(ny,:,1) = one
            vl(ny,:,2) = zero
            vl(ny,:,3) = zero
            vl(ny,:,4) = zero
            vl(ny,:,5) = one
          end if          

        end if                  ! Ma
        
        end if                  ! yper

        if (xper) then

          vl(:,nx,:) = vl(:,1,:)
          
        else                    ! xper

!=============================================================================!
!       L e f t  B o u n d a r y
!=============================================================================!

        if (Ma.lt.one) then

!.... Symmetry plane

        if (left.eq.2) then
          vl(:,1,1) = -( gc2 * vl(:,2,1)  + &
                         gc3 * vl(:,3,1)  + &
                         gc4 * vl(:,4,1)  + &
                         gc5 * vl(:,5,1)  ) / gc1

          vl(:,1,2) = -( gc2 * vl(:,2,2)  + &
                         gc3 * vl(:,3,2)  + &
                         gc4 * vl(:,4,2)  + &
                         gc5 * vl(:,5,2)  ) / gc1
                         
          vl(:,1,3) = zero

          vl(:,1,4) = -( gc2 * vl(:,2,4)  + &
                         gc3 * vl(:,3,4)  + &
                         gc4 * vl(:,4,4)  + &
                         gc5 * vl(:,5,4)  ) / gc1
                         
          vl(:,1,5) = -( gc2 * vl(:,2,5)  + &
                         gc3 * vl(:,3,5)  + &
                         gc4 * vl(:,4,5)  + &
                         gc5 * vl(:,5,5)  ) / gc1
        end if
        
!.... first-order Riemann

        if (left.eq.0) then
          call Reimann( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                        bnl, x(:,1), y(:,1), rhobl, ubl, vbl, wbl, &
                        tbl, pbl, cbl )
        end if

        if (left.eq.7) then           ! symmetry boundary
          vl(:,1,3) = zero
        end if
        
        end if          ! Ma

        if (.false.) then
        
!.... zeroth order extrapolation in the viscous layers

        if (extrap.eq.0) then
          do j = 1, nbl
            vl(j,1,1) = vl(j,2,1)
            vl(j,1,2) = vl(j,2,2)
            vl(j,1,3) = vl(j,2,3)
            vl(j,1,4) = vl(j,2,4)
            vl(j,1,5) = vl(j,2,5)
          end do
        end if
        
!.... first-order extrapolation in the viscous layers

        if (extrap.eq.1) then
          do j = 1, nbl
            vl(j,1,:) = two * vl(j,2,:) - vl(j,3,:)
          end do
        end if
        
!.... second-order extrapolation in the viscous layers

        if (extrap.eq.2) then
          do j = 1, nbl
            vl(j,1,:) = three * vl(j,2,:) - three * vl(j,3,:) + &
                        vl(j,4,:)
          end do
        end if                  ! extrapolation type
        
        end if

!=============================================================================!
!       R i g h t  B o u n d a r y
!=============================================================================!

        if (Ma.lt.one) then

!.... First-order Riemann

        if (right.eq.0) then
          call Reimann( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                        bnr, x(:,nx), y(:,nx), rhobr, ubr, vbr, wbr, &
                        tbr, pbr, cbr )
        end if

!.... Symmetry boundary

        if (right.eq.7) then           ! symmetry boundary
          vl(:,nx,3) = zero
        end if
        
        end if          ! Ma

!.... find the boundary-layer edge (using delta_99 right now)
!.... It is wise to keep nbl fixed, since every change sets up a transient.

!       do nbl = 1, ny
!         if (vl(nbl,nx-1,2) .gt. 0.99 ) exit
!         if ((vl(nbl,nx-1,5)-vl(1,nx-1,5))/(one-vl(1,nx-1,5)) .gt. 0.99 ) exit
!       end do
!       nbl = 70                ! for ny = 128
        nbl = 35                ! for ny = 63
!       if (iter.eq.1) write(90,*) istep, nbl 

!.... zeroth order extrapolation in the viscous layers

        if (extrap.eq.0) then
          do j = 1, nbl
            vl(j,nx,:) = vl(j,nx-1,:)
          end do
        end if
        
!.... first-order extrapolation in the viscous layers

        if (extrap.eq.1) then
          do j = 1, nbl
            vl(j,nx,:) = two * vl(j,nx-1,:) - vl(j,nx-2,:)
          end do
        end if

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

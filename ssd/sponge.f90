!=============================================================================!
        subroutine sponge(spgl, spg2l)
!  
!  Computes the sponge term 
!  
!=============================================================================!
        use stuff
        use spg_mod
        implicit none
        
        real :: spgl(ny,nx), spg2l(ny,nx)
        integer :: i, j
        real :: rr

!       real, parameter    :: As = 0.2          ! Mahesh sponge
!       integer, parameter :: Ns = 3

!       real, parameter    :: xs = 8.0          ! for cyl
!       real, parameter    :: xt = 11.0         ! for cyl

!       real, parameter    :: xs = 0.88         ! for pcyl
!       real, parameter    :: xt = 1.00         ! for pcyl

!       real, parameter    :: xs = 10.99298818  ! for spatial TS wave
!       real, parameter    :: xt = 13.74123524  ! for spatial TS wave

        real, parameter    :: Ad = 20.0         ! Guo & Adams sponge
        real, parameter    :: bd = 1.3
        real, parameter    :: xd = 15.0
        integer, parameter :: nd = 4

!=============================================================================!
!                       C U R R E N T   S P O N G E
!=============================================================================!

!.... initialize the sponge

        spgl = zero
        if (ispg.ge.2) spg2l = zero

!.... general sponge

          open(9,file='sponge.dat')
          do i = 1, nx                 ! outflow sponge
            do j = 1, ny
              if ( xi(i) .ge. xs ) then
                spgl(j,i) = As * ((xi(i)-xs)/(xt-xs))**Ns
              end if
            end do
            write(9,"(3(1pe13.6,1x))") x(1,i), xi(i), spgl(1,i)
          end do
          close(9)

          if (ispg.ge.2) then          ! inflow sponge
            open(9,file='sponge2.dat')
            do j = 1, ny
              do i = 1, nx
                if ( eta(j) .ge. xs2 ) then
                  spg2l(j,i) = As2 * ((eta(j)-xs2)/(xt2-xs2))**Ns2
                end if
              end do
              write(9,"(3(1pe13.6,1x))") x(j,1), eta(j), spg2l(j,1)
            end do
            close(9)
          end if

        return

!=============================================================================!
!                    L I B R A R Y   O F   S P O N G E S
!=============================================================================!

!.... sponge for spatial TS wave

        do i = 1, nx
          do j = 1, ny
            if ( x(j,i) .ge. xs ) then
              spgl(j,i) = As * ((x(j,i)-xs)/(xt-xs))**Ns
            end if
          end do
        end do

!.... sponge for cylinder scattering (Guo & Adams)

        do i = 1, nx
          do j = 1, ny
            rr = sqrt( x(j,i)**2 + y(j,i)**2 )
!           rr = y(j,i)
            if ( rr .gt. xs2 .and. rr .lt. xt2 ) then
              spg2l(j,i) = As2 * 10.0**( -bd**nd * ((xd-rr)/(xd-xs2))**nd)
            else
              spg2l(j,i) = zero
            end if
          end do
        end do
        open(9,file='sponge.dat')
        i = 1
        do j = 1, ny
          write(9,"(3(1pe13.6,1x))") sqrt( x(j,i)**2 + y(j,i)**2 ), spg2l(j,i)
        end do
        close(9)

!.... sponge for wave problem (Guo & Adams)

        do i = 1, nx
          do j = 1, ny
            rr = x(j,i)
            if ( rr .gt. xs .and. rr .lt. xt ) then
              spgl(j,i) = As * 10.0**( -bd**nd * ((xd-rr)/(xd-xs2))**nd)
            else
              spgl(j,i) = zero
            end if
          end do
        end do
        open(8,file='sponge.dat')
        j = 1
        do i = 1, nx
          write(8,"(3(1pe13.6,1x))") x(j,i), spgl(j,i)
        end do
        close(8)

        close(9)

!.... sponge for cylinder scattering (Mahesh)

        do i = 1, nx
          do j = 1, ny
            rr = sqrt( x(j,i)**2 + y(j,i)**2 )
            if ( rr .ge. xs2 .and. rr .le. xt2 ) then
              spg2l(j,i) = As2 * ((rr-xs2)/(xt2-xs2))**Ns2
            else
              spg2l(j,i) = zero
            end if
          end do
        end do

!.... sponge for mean flow and parabolic cylinder

        if (.false.) then
          open(9,file='sponge.dat')
          do i = 1, nx                 ! outflow sponge
            do j = 1, ny
              if ( xi(i) .ge. xs ) then
                spgl(j,i) = As * ((xi(i)-xs)/(xt-xs))**Ns
              end if
            end do
            write(9,"(3(1pe13.6,1x))") x(1,i), xi(i), spgl(1,i)
          end do
          close(9)

          if (ispg.ge.2) then          ! inflow sponge
            open(9,file='sponge2.dat')
            do j = 1, ny
              do i = 1, nx
                if ( eta(j) .ge. xs2 ) then
                  spg2l(j,i) = As2 * ((eta(j)-xs2)/(xt2-xs2))**Ns2
                end if
              end do
              write(9,"(3(1pe13.6,1x))") x(j,1), eta(j), spg2l(j,1)
            end do
            close(9)
          end if
        end if

!.... sponge for the wave problem

        if (.false.) then

          do i = 1, nx
            do j = 1, ny
              if ( x(j,i) .ge. 20.0 .and. x(j,i) .le. 40.0 ) then
                spgl(j,i) = As * ((x(j,i)-20.0)/(40.0-20.0))**Ns
              end if
            end do
          end do
            
        end if

!.... sponge for the MSE problem  (This is what I last used)

        if (.false.) then

!.... inflow sponge

            do j = 1, ny
              if ( eta(j) .ge. xs ) then
                spgl(j,:) = As * ((eta(j)-xs)/(xt-xs))**Ns
              end if
            end do

!.... outflow sponge

            do i = 1, nx
              do j = 1, ny
                if ( x(j,i) .ge. 20.0 ) then
                  spg2l(j,i) = As * ((x(j,i)-20.0)/(40.0-20.0))**Ns
                end if
              end do
            end do

        end if

!.... sponge for the spatial instability problem  (This is what I last used)

        if (.false.) then

            do i = 1, nx
              do j = 1, ny
                if ( x(j,i) .ge. xs ) then
                  spgl(j,i) = As * ((x(j,i)-xs)/(xt-xs))**Ns
                end if
              end do
            end do

        end if

!.... sponge for the MSE problem (Guo & Adams)

        if (.false.) then

!.... inflow sponge

            do j = 1, ny
              if ( eta(j) .ge. xs .and. eta(j) .le. xt ) then
                spgl(j,:) = Ad * 10.0**(-(bd**nd)*((xt-eta(j))/(xt-xs))**nd)
              else
                spgl(j,:) = zero
              end if
!              write(9,"(2(1pe13.6,1x))") x(j,1), spgl(j,1)
            end do

!.... outflow sponge

!           do i = 1, nx
!             do j = 1, ny
!               if ( x(j,i) .ge. 35.0 .and. x(j,i) .le. 44.0 ) then
!                 spgl(j,i) = spgl(j,i) + Ad * 10.0**(-(bd**nd)* &
!                             ((44.0-x(j,i))/(44.0-35.0))**nd)
!               end if
!             end do
!           end do

        end if

!.... sponge for the cylinder scattering problem (Mahesh sponge)

        if (.false.) then

        do i = 1, nx
          do j = 1, ny
            rr = sqrt( x(j,i)**2 + y(j,i)**2 )

            if ( rr .gt. xs .and. rr .lt. xt ) then
              spgl(j,i) = As * ((rr-xs)/(xt-xs))**Ns
            else
              spgl(j,i) = zero
            end if

          end do
        end do

        end if
        
        if (.false.) then

!.... Adams (CTR)

        do i = 1, nx
          do j = 1, ny
            rr = sqrt( x(j,i)**2 + y(j,i)**2 )

            if ( rr .gt. xs .and. rr .lt. xt ) then
              spgl(j,i) = As * real(Ns+1) * real(Ns + 2)*       &
                          (rr - xs)**Ns * (xt - rr) /           &
                          (xt - xs)**(Ns+2)
            else
              spgl(j,i) = zero
            end if

          end do
        end do

        end if
!
!.... sponge for the spatial instabilities
!
!       do i = 1, nx
!         do j = 1, ny
!           if (x(j,i) .gt. xs .and. x(j,i) .lt. xt) then
!             spgl(j,i) = As*real(Ns+1)*real(Ns+2)*             &
!                         (x(j,i)-xs)**Ns*(xt-x(j,i))           &
!                         /(xt-xs)**(Ns+2)
!           else
!             spgl(j,i) = zero
!           end if
!
!         end do
!       end do

        return
        end

        module sponge

        real    :: As, xs, xt
        integer :: Ns

        real    :: As2, xs2, xt2
        integer :: Ns2

        integer :: ic_start=0

        real, allocatable :: vic(:,:,:)
        !$sgi distribute vic(*,*,block)

        complex, allocatable :: cvic(:,:,:)
        !$sgi distribute cvic(*,*,block)

        end module sponge

!=============================================================================!
        subroutine init_sponge
!
!  Computes the sponge term
!
!  Revised: 12-15-00   Switched to i,j indices [SSC]
!=============================================================================!
        use global
        use sponge
        implicit none

        integer :: i, j, ier
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

        character(80) :: code='init_sponge$'

!=============================================================================!
!                       C U R R E N T   S P O N G E
!=============================================================================!

        if (ispg.eq.0) return

        allocate( spg(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for spg$')

!.... initialize the sponge

        !$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            spg(i,j) = zero
          end do
        end do

#if 1

!.... Sponge on outflow (right) boundary

        !$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            if ( xi(i) .ge. xs ) then
              spg(i,j) = As * ((xi(i)-xs)/(xt-xs))**Ns
            end if
          end do
        end do
        j = 1
        open(9,file='sponge.dat')
        do i = 1, nx
          write(9,"(3(1pe13.6,1x))") x(i,j), xi(i), spg(i,j)
        end do
        close(9)

#else

!.... Sponge near top boundary

        !$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            if ( eta(j) .ge. xs ) then
              spg(i,j) = As * ((eta(j)-xs)/(xt-xs))**Ns
            end if
          end do
        end do
        i = 1
        open(9,file='sponge.dat')
        do j = 1, ny
          write(9,"(3(1pe13.6,1x))") x(i,j), eta(j), spg(i,j)
        end do
        close(9)

#endif

        return

!=============================================================================!
!                    L I B R A R Y   O F   S P O N G E S
!=============================================================================!

        if (ispg .ge. 2) then
          allocate( spg2(nx,ny), STAT=ier)
        end if

        !$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            spg2(i,j) = zero
          end do
        end do

!.... general sponge

        !$omp parallel do private(i,j)
        do j = 1, ny                 ! outflow sponge
          do i = 1, nx
            if ( xi(i) .ge. xs ) then
              spg(i,j) = As * ((xi(i)-xs)/(xt-xs))**Ns
            end if
          end do
        end do
        j = 1
        open(9,file='sponge.dat')
        do i = 1, nx
          write(9,"(3(1pe13.6,1x))") x(i,j), xi(i), spg(i,j)
        end do
        close(9)

        if (ispg.ge.2) then          ! inflow sponge
          !$omp parallel do private(i,j)
          do j = 1, ny
            do i = 1, nx
              if ( eta(j) .ge. xs2 ) then
                spg2(i,j) = As2 * ((eta(j)-xs2)/(xt2-xs2))**Ns2
              end if
            end do
          end do
          i = 1
          open(9,file='sponge2.dat')
          do j = 1, ny
            write(9,"(3(1pe13.6,1x))") x(i,j), eta(j), spg2(i,j)
          end do
          close(9)
        end if

!.... WARNING:  The indicies need switching [SSC]

!.... sponge for spatial TS wave

        do i = 1, nx
          do j = 1, ny
            if ( x(j,i) .ge. xs ) then
              spg(j,i) = As * ((x(j,i)-xs)/(xt-xs))**Ns
            end if
          end do
        end do

!.... sponge for cylinder scattering (Guo & Adams)

        do i = 1, nx
          do j = 1, ny
            rr = sqrt( x(j,i)**2 + y(j,i)**2 )
!           rr = y(j,i)
            if ( rr .gt. xs2 .and. rr .lt. xt2 ) then
              spg2(j,i) = As2 * 10.0**( -bd**nd * ((xd-rr)/(xd-xs2))**nd)
            else
              spg2(j,i) = zero
            end if
          end do
        end do
        open(9,file='sponge.dat')
        i = 1
        do j = 1, ny
          write(9,"(3(1pe13.6,1x))") sqrt( x(j,i)**2 + y(j,i)**2 ), spg2(j,i)
        end do
        close(9)

!.... sponge for wave problem (Guo & Adams)

        do i = 1, nx
          do j = 1, ny
            rr = x(j,i)
            if ( rr .gt. xs .and. rr .lt. xt ) then
              spg(j,i) = As * 10.0**( -bd**nd * ((xd-rr)/(xd-xs2))**nd)
            else
              spg(j,i) = zero
            end if
          end do
        end do
        open(8,file='sponge.dat')
        j = 1
        do i = 1, nx
          write(8,"(3(1pe13.6,1x))") x(j,i), spg(j,i)
        end do
        close(8)

        close(9)

!.... sponge for cylinder scattering (Mahesh)

        do i = 1, nx
          do j = 1, ny
            rr = sqrt( x(j,i)**2 + y(j,i)**2 )
            if ( rr .ge. xs2 .and. rr .le. xt2 ) then
              spg2(j,i) = As2 * ((rr-xs2)/(xt2-xs2))**Ns2
            else
              spg2(j,i) = zero
            end if
          end do
        end do

!.... sponge for mean flow and parabolic cylinder

        if (.false.) then
          open(9,file='sponge.dat')
          do i = 1, nx                 ! outflow sponge
            do j = 1, ny
              if ( xi(i) .ge. xs ) then
                spg(j,i) = As * ((xi(i)-xs)/(xt-xs))**Ns
              end if
            end do
            write(9,"(3(1pe13.6,1x))") x(1,i), xi(i), spg(1,i)
          end do
          close(9)

          if (ispg.ge.2) then          ! inflow sponge
            open(9,file='sponge2.dat')
            do j = 1, ny
              do i = 1, nx
                if ( eta(j) .ge. xs2 ) then
                  spg2(j,i) = As2 * ((eta(j)-xs2)/(xt2-xs2))**Ns2
                end if
              end do
              write(9,"(3(1pe13.6,1x))") x(j,1), eta(j), spg2(j,1)
            end do
            close(9)
          end if
        end if

!.... sponge for the wave problem

        if (.false.) then

          do i = 1, nx
            do j = 1, ny
              if ( x(j,i) .ge. 20.0 .and. x(j,i) .le. 40.0 ) then
                spg(j,i) = As * ((x(j,i)-20.0)/(40.0-20.0))**Ns
              end if
            end do
          end do

        end if

!.... sponge for the MSE problem  (This is what I last used)

        if (.false.) then

!.... inflow sponge

            do j = 1, ny
              if ( eta(j) .ge. xs ) then
                spg(j,:) = As * ((eta(j)-xs)/(xt-xs))**Ns
              end if
            end do

!.... outflow sponge

            do i = 1, nx
              do j = 1, ny
                if ( x(j,i) .ge. 20.0 ) then
                  spg2(j,i) = As * ((x(j,i)-20.0)/(40.0-20.0))**Ns
                end if
              end do
            end do

        end if

!.... sponge for the spatial instability problem  (This is what I last used)

        if (.false.) then

            do i = 1, nx
              do j = 1, ny
                if ( x(j,i) .ge. xs ) then
                  spg(j,i) = As * ((x(j,i)-xs)/(xt-xs))**Ns
                end if
              end do
            end do

        end if

!.... sponge for the MSE problem (Guo & Adams)

        if (.false.) then

!.... inflow sponge

            do j = 1, ny
              if ( eta(j) .ge. xs .and. eta(j) .le. xt ) then
                spg(j,:) = Ad * 10.0**(-(bd**nd)*((xt-eta(j))/(xt-xs))**nd)
              else
                spg(j,:) = zero
              end if
!              write(9,"(2(1pe13.6,1x))") x(j,1), spg(j,1)
            end do

!.... outflow sponge

!           do i = 1, nx
!             do j = 1, ny
!               if ( x(j,i) .ge. 35.0 .and. x(j,i) .le. 44.0 ) then
!                 spg(j,i) = spg(j,i) + Ad * 10.0**(-(bd**nd)* &
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
              spg(j,i) = As * ((rr-xs)/(xt-xs))**Ns
            else
              spg(j,i) = zero
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
              spg(j,i) = As * real(Ns+1) * real(Ns + 2)*        &
                          (rr - xs)**Ns * (xt - rr) /           &
                          (xt - xs)**(Ns+2)
            else
              spg(j,i) = zero
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
!             spg(j,i) = As*real(Ns+1)*real(Ns+2)*              &
!                         (x(j,i)-xs)**Ns*(xt-x(j,i))           &
!                         /(xt-xs)**(Ns+2)
!           else
!             spg(j,i) = zero
!           end if
!
!         end do
!       end do

        return
        end

!=============================================================================!
        subroutine spg_it( rl, vl, spgl)
!
!       This routine reads in the initial condition when first called, and
!       sponges to the initial condition.
!
!=============================================================================!
        use global
        use sponge
        implicit none

        real :: rl(ndof,nx,ny), vl(ndof,nx,ny), spgl(nx,ny)
        !$sgi distribute rl(*,*,block), vl(*,*,block), spgl(*,block)

        integer :: i, j

        real :: rtmp
        integer :: itmp
!=============================================================================!
        if (ic_start .eq. 0) then
          allocate( vic(ndof,nx,ny) )
          ic_start = 1
          open(10,file='output.R.0',form='unformatted',status='old')
          read(10) itmp, rtmp, itmp, itmp, itmp, itmp, &
                   rtmp, rtmp, rtmp, rtmp, rtmp
          read(10) vic
          close(10)
        end if

        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            rl(:,i,j) = rl(:,i,j) + spgl(i,j) * ( vl(:,i,j) - vic(:,i,j) )
          end do
        end do

        return
        end

!=============================================================================!
        subroutine cspg_it( rl, vl, spgl, spg2l )
!
!       This routine reads in the initial condition when first called, and
!       sponges to the initial condition.
!
!=============================================================================!
        use global
        use sponge
        implicit none

        complex :: rl(ndof,nx,ny), vl(ndof,nx,ny)
        real    :: spgl(nx,ny), spg2l(nx,ny)
        integer :: i, j

        complex, parameter :: ac=(2.2804739410500E-001,-6.5163146761218E-003)
!       complex, parameter :: ac=(-2.8831962908130E-001,-1.3854663671636E-002)

        real :: rtmp
        integer :: itmp
!=============================================================================!
        if (ic_start .eq. 0) then
          allocate( cvic(ndof,nx,ny) )
          ic_start = 1
!         open(10,file='output.R.0',form='unformatted',status='old')
!         read(10) itmp, rtmp, itmp, itmp, itmp, itmp, &
!                  rtmp, rtmp, rtmp, rtmp, rtmp
!         read(10) vic
!         close(10)
          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              cvic(1,i,j) = cmplx(rhor(j), rhoi(j)) * exp(im * ac * x(i,j))
              cvic(2,i,j) = cmplx(  ur(j),   ui(j)) * exp(im * ac * x(i,j))
              cvic(3,i,j) = cmplx(  vr(j),   vi(j)) * exp(im * ac * x(i,j))
              cvic(4,i,j) = cmplx(  wr(j),   wi(j)) * exp(im * ac * x(i,j))
              cvic(5,i,j) = cmplx(  tr(j),   ti(j)) * exp(im * ac * x(i,j))
            end do
          end do
        end if

        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            rl(:,i,j) = rl(:,i,j) + (spgl(i,j) + spg2l(i,j)) * &
                        ( vl(:,i,j) - cvic(:,i,j) )
          end do
        end do

        return
        end

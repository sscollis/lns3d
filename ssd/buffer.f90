!=============================================================================!
        subroutine buffer()
!=============================================================================!
        use stuff
        use buff_stuff
        implicit none
        
        integer :: i, j

        real    :: xs = 0.8
        real    :: xt = 1.0

        real    :: aad = 2.57
        integer :: nd  = 4

        real    :: xis, xit, etas, etat

        integer :: ier
!=============================================================================!
        allocate( buff(ny,nx), STAT=ier)
        if (ier .ne. 0) call error('buffer$','Insufficient Memory for buff$')

!.... initialize

        buff = zero

!.... set damping zone

        if (ibuff .eq. 0) then          !.... everywhere

          buff = one

        else if (ibuff .eq. 1) then     !.... top boundary

          xs = 0.8
          xt = 1.0

          do i = 1, nx
            do j = 1, ny
              if ( eta(j) .ge. xs ) then
                buff(j,i) = one - 10.0**( -(aad**nd) * ( (eta(j) - xs)/ &
                            (xt - xs) )**nd )
              end if
            end do
          end do

        else if (ibuff .eq. 2) then     !.... Wall boundary

          xs = 0.1
          xt = 0.0

          do i = 1, nx
            do j = 1, ny
              if ( eta(j) .le. xs ) then
                buff(j,i) = one - 10.0**( -(aad**nd) * ( (xs - eta(j))/ &
                            (xs - xt) )**nd )
              end if
            end do
          end do

        else if (ibuff .eq. 3) then     !.... outflow boundary

          xs = 0.8
          xt = 1.0

          do i = 1, nx
            do j = 1, ny
              if ( xi(i) .ge. xs ) then
                buff(j,i) = one - 10.0**( -(aad**nd) * ( (xi(i) - xs)/ &
                            (xt - xs) )**nd )
              end if
            end do
          end do

        else if (ibuff .eq. 4) then     !.... top and outflow boundaries

          xs = 0.8
          xt = 1.0

          do i = 1, nx
            do j = 1, ny
              if ( eta(j) .ge. xs .and. xi(i) .lt. xs ) then
                buff(j,i) = one - 10.0**( -(aad**nd) * ( (eta(j) - xs)/ &
                            (xt - xs) )**nd )
              else if ( eta(j) .lt. xs .and. xi(i) .ge. xs ) then
                buff(j,i) = one - 10.0**( -(aad**nd) * ( (xi(i) - xs)/ &
                            (xt - xs) )**nd )
              else if ( eta(j) .ge. xs .and. xi(i) .ge. xs ) then
                buff(j,i) = one - ( &
                            10.0**( -(aad**nd) * ( (eta(j) - xs)/ &
                            (xt - xs) )**nd ) * &
                            10.0**( -(aad**nd) * ( (xi(i) - xs)/ &
                            (xt - xs) )**nd ) )
              end if
            end do
          end do

        end if
        
!.... diagnostic

!       open(20,file='buffer.dat')
!       do i = 1, nx
!         do j = 1, ny
!           write(20,"(3(1pe13.6,1x))") x(j,i), eta(j), buff(j,i)
!         end do
!       end do
!       close(20)

        return
        end 

        module buffer

        real, allocatable :: buff(:,:)
        !$sgi distribute buff(*,block)

        integer :: ibuff = 0
        
        end module buffer

!=============================================================================!
        subroutine init_buffer
!=============================================================================!
        use global
        use buffer
        implicit none
        
        integer :: i, j

        real    :: xs = 0.8
        real    :: xt = 1.0

        real    :: aad = 2.57
        integer :: nd  = 4

        real    :: xis, xit, etas, etat

        integer :: ier
        character*80 :: code='init_buffer$'
!=============================================================================!
        if (eps_e .eq. zero) return

        allocate( buff(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for buff$')

!.... initialize

        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            buff(i,j) = zero
          end do
        end do

!.... set damping zone

        if (ibuff .eq. 0) then          !.... everywhere

!         write(*,*) 'ibuff == 0'

          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              buff(i,j) = one
            end do
          end do

        else if (ibuff .eq. 1) then     !.... top boundary

          xs = 0.8
          xt = 1.0

          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              if ( eta(j) .ge. xs ) then
                buff(i,j) = one - 10.0**( -(aad**nd) * ( (eta(j) - xs)/ &
                            (xt - xs) )**nd )
              end if
            end do
          end do

        else if (ibuff .eq. 2) then     !.... Wall boundary

          xs = 0.1
          xt = 0.0

          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              if ( eta(j) .le. xs ) then
                buff(i,j) = one - 10.0**( -(aad**nd) * ( (xs - eta(j))/ &
                            (xs - xt) )**nd )
              end if
            end do
          end do

        else if (ibuff .eq. 3) then     !.... outflow boundary

          xs = 0.8
          xt = 1.0

          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              if ( xi(i) .ge. xs ) then
                buff(i,j) = one - 10.0**( -(aad**nd) * ( (xi(i) - xs)/ &
                            (xt - xs) )**nd )
              end if
            end do
          end do

        else if (ibuff .eq. 4) then     !.... top and outflow boundaries

          xs = 0.8
          xt = 1.0

          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              if ( eta(j) .ge. xs .and. xi(i) .lt. xs ) then
                buff(i,j) = one - 10.0**( -(aad**nd) * ( (eta(j) - xs)/ &
                            (xt - xs) )**nd )
              else if ( eta(j) .lt. xs .and. xi(i) .ge. xs ) then
                buff(i,j) = one - 10.0**( -(aad**nd) * ( (xi(i) - xs)/ &
                            (xt - xs) )**nd )
              else if ( eta(j) .ge. xs .and. xi(i) .ge. xs ) then
                buff(i,j) = one - ( &
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
!           write(20,"(3(1pe13.6,1x))") x(i,j), eta(j), buff(i,j)
!         end do
!       end do
!       close(20)

        return
        end 

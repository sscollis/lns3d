!=============================================================================!
        subroutine init_damper
!  
!  Computes the damping term for parabolizing the equations near the 
!  outflow boundary.
!  
!  Not used right now:  6-20-96
!=============================================================================!
        use global
        implicit none
        
        real    :: xt
        real    :: xs           ! start of buffer
        real    :: xe           ! end of buffer
        real    :: ab = 1.25
        integer :: nb = 4
        
        integer :: i, j, ier
        character*80 :: code='init_damper$'
!=============================================================================!

        allocate( damp(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for damp$')

        xs = 0.8
        xe = 1.0

        do i = 1, nx
          do j = 1, ny
            if (xi(i).ge.xs) then
              xt = (xi(i)-xs)/(xe - xs)
              damp(i,j) = 10.0**( -(ab**nb) * xt**nb )
            else
              damp(i,j) = one
            end if
          end do
!         write(99,"(2(1pe13.6,1x))") x(i,1), damp(i,1)
        end do
        
        if (idamp.eq.0) then                    ! turn it off
          damp = one
        else if (idamp.eq.1) then               ! damp only the last node
          damp = one
          damp(nx,:) = zero
        else if (idamp.ne.2) then
          call error(code,'Illegal value for idamp$')
        end if

        return
        end

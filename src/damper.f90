!=============================================================================!
        subroutine init_damper
!  
!  Computes the damping term for parabolizing the equations near the 
!  outflow boundary.
!  
!  Not used right now:  6-20-96
!  Revised:             6-07-22
!=============================================================================!
        use global
        implicit none
        
        real    :: xt
        real    :: xs           ! start of buffer
        real    :: xe           ! end of buffer
        real    :: ab = 1.25
        integer :: nb = 4
        
        integer :: i, j, ier
        character(80) :: code='init_damper$'
!=============================================================================!

        allocate( damp(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for damp$')

        if (idamp.eq.0) then                    ! turn it off
          damp = one
        else if (idamp.eq.1) then               ! damp only the last node
          damp = one
          damp(nx,:) = zero
        else if (idamp.eq.2) then               ! damp the last 20% in xi
          xs = 0.8
          xe = 1.0
          do j = 1, ny
            do i = 1, nx
              if (xi(i).ge.xs) then
                xt = (xi(i)-xs)/(xe - xs)
                damp(i,j) = 10.0**( -(ab**nb) * xt**nb )
              else
                damp(i,j) = one
              end if
              if (j.eq.1) write(99,"(2(1pe13.6,1x))") x(i,j), damp(i,j)
            end do
          end do
        else 
          call error(code,'Illegal value for idamp$')
        end if

        return
        end subroutine init_damper

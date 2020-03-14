!=============================================================================!
        subroutine damper(dampl)
!  
!  Computes the damping term for parabolizing the equations near the 
!  outflow boundary.
!  
!  Not used right now:  6-20-96
!=============================================================================!
        use global
        implicit none
        
        real    :: dampl(ny,nx)
        
        real    :: xt
        real    :: xs           ! start of buffer
        real    :: xe           ! end of buffer
        real    :: ab = 1.25
        integer :: nb = 4
        
        integer :: i, j
        character*80 :: code='ImpDrv$'
!=============================================================================!

        xs = 0.8
        xe = 1.0

        do i = 1, nx
          do j = 1, ny
            if (xi(i).ge.xs) then
              xt = (xi(i)-xs)/(xe - xs)
              dampl(j,i) = 10.0**( -(ab**nb) * xt**nb )
            else
              dampl(j,i) = one
            end if
          end do
!         write(99,"(2(1pe13.6,1x))") x(1,i), dampl(1,i)
        end do
        
        if (idamp.eq.0) then                    ! turn it off
          dampl = one
        else if (idamp.eq.1) then               ! damp only the last node
          dampl = one
          dampl(:,nx) = zero
        else if (idamp.ne.2) then
          call error(code,'Illegal value for idamp$')
        end if

        return
        end

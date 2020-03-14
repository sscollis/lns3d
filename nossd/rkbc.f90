!=============================================================================!
        subroutine rkbc(rl)
!  
!       Satisfy boundary condition on explicit RHS
!
!=============================================================================!
        use global
        use local
        implicit none

        real :: rl(ny,nx,ndof)
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (linear.eq.1) then

        if (yper) then
          rl(ny,:,:) = zero
        else                    ! yper

          if (Navier) then
          
!.... density boundary condition

          if (wall.eq.1 .or. wall.eq.2) rl(1,:,1) = zero

!.... no-slip boundary condition

          rl(1,:,2:ndof-1) = zero

!.... temperature boundary condition

          if (wallt.eq.0 .or. wallt.eq.1) rl(1,:,ndof) = zero

!.... zero freestrean disturbance or Riemann boundary

          if (top.eq.0) rl(ny,:,:) = zero
        
          else            ! Inviscid

            rl(1,:,3) = zero         ! assume a flat plate

          end if          ! Navier

        end if                  ! yper
        
        if (xper) then
          rl(:,nx,:) = zero
        else                    ! xper
          if (left.eq.0) then
            rl(:,1,:)  = zero
          else if (left.eq.1) then
            rl(:,1,1) = zero
          else if (left.eq.2) then
            rl(:,1,2:5) = zero
          else if (left.eq.7) then
            rl(:,1,3) = zero
          end if

          if (right.eq.0) then
            rl(:,nx,:) = zero
          else if (right.eq.1) then
            rl(:,nx,1) = zero
          end if
        end if                  ! xper

        else                    ! linear = 0    
!=============================================================================!
!       N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (yper) then
          rl(ny,:,:) = zero
        else                    ! yper

!.... density boundary condition

          if (wall.eq.1 .or. wall.eq.2) rl(1,:,1) = zero

!.... no-slip boundary condition

          rl(1,:,2:ndof-1) = zero

!.... adiabatic boundary condition

          if (wallt.eq.1) rl(1,:,ndof) = zero

!.... Riemann boundary

          if (top.eq.0) rl(ny,:,:) = zero

        end if                  ! yper

        if (xper) then
          rl(:,nx,:) = zero
        else                    ! xper

!.... Left and Right Riemann boundary conditions

          if (left.eq.0) then
            rl(:,1,:) = zero
          else if (left.eq.2) then
            rl(:,1,:) = zero
          else if (left.eq.7) then
            rl(:,1,3) = zero
          end if

          if (right.eq.0) then
            rl(:,nx,:) = zero
          end if
          
        end if                  ! xper
        
        end if                  ! linear 
        
        return
        end

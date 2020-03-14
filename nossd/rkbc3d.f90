!=============================================================================!
        subroutine rkbc3d(rl)
!  
!       Satisfy boundary condition on explicit RHS
!
!=============================================================================!
        use global
        implicit none

        complex :: rl(ny,nx,ndof)
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (yper) then
          rl(ny,:,:) = zero
        else                    ! yper

!.... density boundary condition

          if (wall.eq.1 .or. wall.eq.2) rl(1,:,1) = zero

!.... no-slip boundary condition

          rl(1,:,2:ndof-1) = zero

!.... temperature boundary condition

          if (wallt.eq.0 .or. wallt.eq.1) rl(1,:,ndof) = zero

!.... zero freestrean disturbance

          if (top.eq.0 .or. top.eq.1) rl(ny,:,:) = zero
        
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
          else if (left.eq.5) then
            rl(:,1,2:5) = zero
          else if (left.eq.7) then
            rl(:,1,3) = zero
          end if

          if (right.eq.0) then
            rl(:,nx,:) = zero
          else if (right.eq.1) then
            rl(:,nx,1) = zero
          else if (right.eq.8) then
            rl(:,nx,:) = zero
          end if
        end if                  ! xper

        return
        end

!=============================================================================!
        subroutine checkbc(vl)
!=============================================================================!
        use stuff
        real :: vl(ny,nx,ndof)

        real, parameter :: gc1 = -2.083333333333333333333E+00
        real, parameter :: gc2 =  4.000000000000000000000E+00
        real, parameter :: gc3 = -3.000000000000000000000E+00
        real, parameter :: gc4 =  1.333333333333333333333E+00
        real, parameter :: gc5 = -2.500000000000000000000E-01
!=============================================================================!

!.... WARNING

        return

        if (.true.) then

!         vl(1,:,1) = zero              ! zero disturbance
        
          vl(1,:,1) = vl(2,:,1)         ! first-order difference

!         vl(1,:,1)      = -( gc2 * vl(2,:,1)  + &
!                             gc3 * vl(3,:,1)  + &
!                             gc4 * vl(4,:,1)  + &
!                             gc5 * vl(5,:,1)  ) / gc1
        end if

!.... no slip wall

          vl(1,:,2:ndof-1)  = zero

!.... adiabatic temperature ( dT/deta = 0 )

        if (.true.) then

!         vl(1,:,ndof) = zero           ! zero disturbance

          vl(1,:,ndof) = vl(2,:,ndof)   ! first-order difference

!         vl(1,:,ndof)   = -( gc2 * vl(2,:,ndof)  + &
!                             gc3 * vl(3,:,ndof)  + &
!                             gc4 * vl(4,:,ndof)  + &
!                             gc5 * vl(5,:,ndof)  ) / gc1
        end if
        
        return
        end

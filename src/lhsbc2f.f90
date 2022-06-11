!=============================================================================!
        subroutine lhsbc2f (mat, vl)
!  
!  Correct the LHS for boundary conditions in the \eta direction
!
!  This version supports the 4th order LHS
!
!  Revised: 10-17-95
!
!=============================================================================!
        use global
        use stencil
        use pot
        implicit none
        
        real :: mat(5,ndof,ndof,nx,ny), vl(ndof,nx,ny)

        real :: detainv, detasinv

        real :: matl(5,2)

!.... local variables for Riemann invariants

        real :: m1l(nx), m2l(nx), R1(nx), R2(nx), fact(nx), vint(nx,ndof)
        real :: rho(nx), p(nx), fact1(nx), fact2(nx), vn(nx), vt(nx), pnorm(nx)
        real :: ub(nx), vb(nx), tb(nx), rhob(nx), pb(nx), cb(nx)
        real :: xil(nx), etal(nx)
        real :: b, rbeta
        
        integer :: i, j
!=============================================================================!
        detainv  = one / deta
        detasinv = one / deta**2
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (linear.eq.1) then

        if (.not. yper) then

        if (Navier) then
        
!.... wall boundary condition

        if (wall.eq.1) then
          mat(:,1,:,:,1) = zero
          mat(3,1,1,:,1) = -gc1 
          mat(4,1,1,:,1) = -gc2
          mat(5,1,1,:,1) = -gc3
          mat(1,1,1,:,1) = -gc4
          mat(2,1,1,:,1) = -gc5
        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        mat(3,2,2,:,1)   = one
        mat(3,3,3,:,1)   = one
        mat(3,4,4,:,1)   = one

!.... isothermal

        if (wallt.eq.0) then
          mat(:,5,:,:,1) = zero
          mat(3,5,5,:,1) = one
        end if

!.... adiabatic

        if (wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          mat(3,5,5,:,1) = -gc1
          mat(4,5,5,:,1) = -gc2
          mat(5,5,5,:,1) = -gc3
          mat(1,5,5,:,1) = -gc4
          mat(2,5,5,:,1) = -gc5
        end if

        else                    ! inviscid wall
          call lwallbc(vl, vm, p, pnorm)
          mat(:,1,:,:,1) = zero
          mat(3,1,1,:,1) = -gc1 * detainv
          mat(4,1,1,:,1) = -gc2 * detainv
          mat(5,1,1,:,1) = -gc3 * detainv
          mat(1,1,1,:,1) = -gc4 * detainv
          mat(2,1,1,:,1) = -gc5 * detainv

          do i = 1, nx
            if ( abs(bnb(i,2)) .ge. abs(bnb(i,1)) .and. .false. ) then
              mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                               bnb(i,1) * mat(:,3,:,i,1)
              mat(:,3,:,i,1) = zero
              mat(3,3,2,i,1) = bnb(i,1)
              mat(3,3,3,i,1) = bnb(i,2)
            else
              mat(:,3,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                               bnb(i,1) * mat(:,3,:,i,1)
              mat(:,2,:,i,1) = zero
              mat(3,2,2,i,1) = bnb(i,1)
              mat(3,2,3,i,1) = bnb(i,2)
            end if
          end do
        end if                  ! Navier

!.... freestream zero disturbance boundary condition

        if (top.eq.0.or.top.eq.1) then
          mat(:,:,:,:,ny) = zero
          mat(3,1,1,:,ny) = one
          mat(3,2,2,:,ny) = one
          mat(3,3,3,:,ny) = one
          mat(3,4,4,:,ny) = one
          mat(3,5,5,:,ny) = one
        end if
        
        end if          ! yper

!=============================================================================!

        if (.not. xper) then

!.... Left boundary

          if (left.eq.0.or.left.eq.4) then  ! zero disturbance BC

            mat(:,:,:,1,:) = zero
            mat(3,1,1,1,:) = one
            mat(3,2,2,1,:) = one
            mat(3,3,3,1,:) = one
            mat(3,4,4,1,:) = one
            mat(3,5,5,1,:) = one

          else if (left.eq.1) then      ! nonreflecting BC

            mat(:,1,:,1,nbl+1:ny) = zero
            mat(3,1,1,1,nbl+1:ny) = one
            if (.false.) then
              mat(:,:,:,1,1:nbl) = zero
              mat(3,1,1,1,1:nbl) = one
              mat(3,2,2,1,1:nbl) = one
              mat(3,3,3,1,1:nbl) = one
              mat(3,4,4,1,1:nbl) = one
              mat(3,5,5,1,1:nbl) = one
            end if

          else if (left.eq.2) then      ! symmetry plane

            mat(:,2:5,:,1,:) = zero
!           mat(3,1,1,1,:) = one
            mat(3,2,2,1,:) = one
            mat(3,3,3,1,:) = one
            mat(3,4,4,1,:) = one
            mat(3,5,5,1,:) = one

          else if (left.eq.7) then      ! symmetry plane

!           mat(:,3,:,1,:) = zero
!           mat(3,3,3,1,:) = one

          end if

!.... Right boundary

          if (right.eq.0) then          ! zero disturbance

            mat(:,:,:,nx,:) = zero
            mat(3,1,1,nx,:) = one
            mat(3,2,2,nx,:) = one
            mat(3,3,3,nx,:) = one
            mat(3,4,4,nx,:) = one
            mat(3,5,5,nx,:) = one
        
          else if (right.eq.1) then     ! nonreflecting BC

            mat(:,1,:,nx,nbl+1:ny) = zero
            mat(3,1,1,nx,nbl+1:ny) = one
            if (.false.) then
              mat(:,:,:,nx,1:nbl) = zero
              mat(3,1,1,nx,1:nbl) = one
              mat(3,2,2,nx,1:nbl) = one
              mat(3,3,3,nx,1:nbl) = one
              mat(3,4,4,nx,1:nbl) = one
              mat(3,5,5,nx,1:nbl) = one
            end if

          end if

        end if          ! xper

!=============================================================================!
!       N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        else            ! nonlinear

        if (.not. yper) then

        if (Navier) then
        
!.... wall boundary condition

        if (wall.eq.1) then
          mat(:,1,:,:,1) = zero
          mat(3,1,1,:,1) = -gc1
          mat(4,1,1,:,1) = -gc2
          mat(5,1,1,:,1) = -gc3
          mat(1,1,1,:,1) = -gc4
          mat(2,1,1,:,1) = -gc5
        end if

!.... partial tangent for normal momentum equation

        if (wall.eq.2) then
          call wallbc(vl, p, pnorm)
          mat(:,1,:,:,1) = zero
          mat(3,1,1,:,1) = -gc1 * detainv
          mat(3,1,5,:,1) = -gamma * Ma**2 * Pnorm / vl(5,:,1)**2
          mat(4,1,1,:,1) = -gc2 * detainv
          mat(5,1,1,:,1) = -gc3 * detainv
          mat(1,1,1,:,1) = -gc4 * detainv
          mat(2,1,1,:,1) = -gc5 * detainv
        end if
        
!.... Third-order extrapolation of rho at the wall

        if (wall.eq.3) then
          mat(:,1,:,:,1) =  zero
          mat(3,1,1,:,1) = -one
          mat(4,1,1,:,1) =  4.0
          mat(5,1,1,:,1) = -6.0
          mat(1,1,1,:,1) =  4.0
          mat(2,1,1,:,1) = -one
        end if

!.... Lele-Poinsot

        if (wall.eq.4) then

          ! LHS must be modified

        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        mat(3,2,2,:,1) = one
        mat(3,3,3,:,1) = one
        mat(3,4,4,:,1) = one

!.... isothermal

        if (wallt.eq.0) then
          mat(:,5,:,:,1) = zero
          mat(3,5,5,:,1) = -one
        end if

!.... adiabatic

        if (wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          mat(3,5,5,:,1) = -gc1         ! fourth-order difference
          mat(4,5,5,:,1) = -gc2
          mat(5,5,5,:,1) = -gc3
          mat(1,5,5,:,1) = -gc4
          mat(2,5,5,:,1) = -gc5
        end if
        
        else      ! inviscid wall

          do i = 1, nx
            matl(1,1) = mat(3,1,2,i,1) * bnb(i,2) - mat(3,1,3,i,1) * bnb(i,1)
            matl(1,2) = mat(3,1,2,i,1) * bnb(i,1) + mat(3,1,3,i,1) * bnb(i,2)
            matl(2,1) = mat(3,2,2,i,1) * bnb(i,2) - mat(3,2,3,i,1) * bnb(i,1)
            matl(2,2) = mat(3,2,2,i,1) * bnb(i,1) + mat(3,2,3,i,1) * bnb(i,2)
            matl(3,1) = mat(3,3,2,i,1) * bnb(i,2) - mat(3,3,3,i,1) * bnb(i,1)
            matl(3,2) = mat(3,3,2,i,1) * bnb(i,1) + mat(3,3,3,i,1) * bnb(i,2)
            matl(4,1) = mat(3,4,2,i,1) * bnb(i,2) - mat(3,4,3,i,1) * bnb(i,1)
            matl(4,2) = mat(3,4,2,i,1) * bnb(i,1) + mat(3,4,3,i,1) * bnb(i,2)
            matl(5,1) = mat(3,5,2,i,1) * bnb(i,2) - mat(3,5,3,i,1) * bnb(i,1)
            matl(5,2) = mat(3,5,2,i,1) * bnb(i,1) + mat(3,5,3,i,1) * bnb(i,2)
            mat(3,:,2:3,i,1) = matl
            
            matl(1,1) = mat(2,1,2,i,2) * bnb(i,2) - mat(2,1,3,i,2) * bnb(i,1)
            matl(1,2) = mat(2,1,2,i,2) * bnb(i,1) + mat(2,1,3,i,2) * bnb(i,2)
            matl(2,1) = mat(2,2,2,i,2) * bnb(i,2) - mat(2,2,3,i,2) * bnb(i,1)
            matl(2,2) = mat(2,2,2,i,2) * bnb(i,1) + mat(2,2,3,i,2) * bnb(i,2)
            matl(3,1) = mat(2,3,2,i,2) * bnb(i,2) - mat(2,3,3,i,2) * bnb(i,1)
            matl(3,2) = mat(2,3,2,i,2) * bnb(i,1) + mat(2,3,3,i,2) * bnb(i,2)
            matl(4,1) = mat(2,4,2,i,2) * bnb(i,2) - mat(2,4,3,i,2) * bnb(i,1)
            matl(4,2) = mat(2,4,2,i,2) * bnb(i,1) + mat(2,4,3,i,2) * bnb(i,2)
            matl(5,1) = mat(2,5,2,i,2) * bnb(i,2) - mat(2,5,3,i,2) * bnb(i,1)
            matl(5,2) = mat(2,5,2,i,2) * bnb(i,1) + mat(2,5,3,i,2) * bnb(i,2)
            mat(2,:,2:3,i,2) = matl
            
            matl(1,1) = mat(1,1,2,i,3) * bnb(i,2) - mat(1,1,3,i,3) * bnb(i,1)
            matl(1,2) = mat(1,1,2,i,3) * bnb(i,1) + mat(1,1,3,i,3) * bnb(i,2)
            matl(2,1) = mat(1,2,2,i,3) * bnb(i,2) - mat(1,2,3,i,3) * bnb(i,1)
            matl(2,2) = mat(1,2,2,i,3) * bnb(i,1) + mat(1,2,3,i,3) * bnb(i,2)
            matl(3,1) = mat(1,3,2,i,3) * bnb(i,2) - mat(1,3,3,i,3) * bnb(i,1)
            matl(3,2) = mat(1,3,2,i,3) * bnb(i,1) + mat(1,3,3,i,3) * bnb(i,2)
            matl(4,1) = mat(1,4,2,i,3) * bnb(i,2) - mat(1,4,3,i,3) * bnb(i,1)
            matl(4,2) = mat(1,4,2,i,3) * bnb(i,1) + mat(1,4,3,i,3) * bnb(i,2)
            matl(5,1) = mat(1,5,2,i,3) * bnb(i,2) - mat(1,5,3,i,3) * bnb(i,1)
            matl(5,2) = mat(1,5,2,i,3) * bnb(i,1) + mat(1,5,3,i,3) * bnb(i,2)
            mat(1,:,2:3,i,3) = matl

            mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                             bnb(i,1) * mat(:,3,:,i,1)
            mat(:,3,:,i,1) = bnb(i,1) * mat(:,2,:,i,1) + &
                             bnb(i,2) * mat(:,3,:,i,1)

!.... constrain the normal velocity to be zero

            mat(:,3,:,i,1) = zero
            mat(3,3,3,i,1) = one

          end do

        end if                  ! Navier
!=============================================================================!
!       T o p
!=============================================================================!
        if (Ma.lt.one) then

        if (top.eq.0) then
          mat(:,:,:,:,ny) = zero
          call ReimannLHS( nx, vl(:,:,ny), vl(:,:,ny-1), vl(:,:,ny-2), &
                           bnt, x(:,ny), y(:,ny), mat(3,:,:,:,ny), &
                           mat(2,:,:,:,ny), mat(1,:,:,:,ny), &
                           rhobt, ubt, vbt, wbt, tbt, pbt, cbt )
        end if

        if (top.eq.3) then
          mat(:,:,:,:,ny) = zero
          mat(3,1,1,:,ny) = one
          mat(3,2,2,:,ny) = one
          mat(3,3,3,:,ny) = one
          mat(3,4,4,:,ny) = one
          mat(3,5,5,:,ny) = one
        end if

        if (.false.) then
        
!.... Riemann invariants far-field boundary condition (zero-order)
        
        j = ny
        do i = 1, nx
          m1l(i)  = n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
          m2l(i)  = n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
        end do    

!.... compute Vn and Vt
        
        vn(:) = m1l(:) * vl(2,:,ny) + m2l(:) * vl(3,:,ny)
        vt(:) = m2l(:) * vl(2,:,ny) - m1l(:) * vl(3,:,ny)
          
!.... perform extrapolation (zero-order)

        vint(:,:) = vl(:,:,ny-1)

!.... compute the boundary values from potential theory

        b = one / pi
        if (Ma.lt.one) then
          rbeta = sqrt(one - Ma**2)
        else
          rbeta = sqrt(Ma**2 - one)
        end if
        ub(:)   = one + (x(:,ny)-b)/pi/rbeta**2/((x(:,ny)-b)**2 + &
                  rbeta**2*y(:,ny)**2)
        vb(:)   = rbeta*y(:,ny)/pi/((x(:,ny)-b)**2+rbeta**2*y(:,ny)**2)
        tb(:)   = one - pt5*gamma1*Ma**2*( ub(:)**2 + vb(:)**2 - one )
        rhob(:) = tb(:)**(one/gamma1)
        pb(:)   = rhob(:) * tb(:) / ( gamma * Ma**2 )
        cb(:)   = sqrt(tb(:))/Ma

        R1(:) = m1l(:) * ub(:) + m2l(:) * vb(:) - two / Ma / gamma1
        R2(:) = m1l(:) * vint(:,2) + m2l(:) * vint(:,3) + &
                two * sqrt(vint(:,5)) / Ma / gamma1
                
        mat(:,:,:,:,ny) = zero          ! initialize
        
        do i = 1, nx

        if (vn(i) .le. zero) then
          rho(i) = rhob(i)
          p(i)   = pb(i)
         
          vn(i)  = pt5 * ( R1(i) + R2(i) )
          vt(i)  = m2l(i) * ub(i) - m1l(i) * vb(i)

          fact(i) = (rho(i)**gamma/(gamma*p(i))* &
                    (gamma1/four)**2)**(one/gamma1) * &
                    two/gamma1 * (R2(i) - R1(i))**((two-gamma1)/gamma1)
          mat(3,1,1,i,ny) = -one
          mat(2,1,2,i,ny) = (fact(i) * m1l(i))
          mat(2,1,3,i,ny) = (fact(i) * m2l(i))
          mat(2,1,5,i,ny) = (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
          
          fact(i) = pt5 * m1l(i) / (m1l(i)**2 + m2l(i)**2)
          mat(3,2,2,i,ny) = -one
          mat(2,2,2,i,ny) = (fact(i) * m1l(i))
          mat(2,2,3,i,ny) = (fact(i) * m2l(i))
          mat(2,2,5,i,ny) = (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
          fact(i) = pt5 * m2l(i) / (m1l(i)**2 + m2l(i)**2)
          mat(3,3,3,i,ny) = -one
          mat(2,3,2,i,ny) = (fact(i) * m1l(i))
          mat(2,3,3,i,ny) = (fact(i) * m2l(i))
          mat(2,3,5,i,ny) = (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
          mat(3,4,4,i,ny) = one
  
          fact(i) = (Ma * gamma1 / four)**2 * two * ( R2(i) - R1(i) )
          mat(3,5,5,i,ny) = -one
          mat(2,5,2,i,ny) = (fact(i) * m1l(i))
          mat(2,5,3,i,ny) = (fact(i) * m2l(i))
          mat(2,5,5,i,ny) = (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))

        else
        
          rho(i) = vint(i,1)
          p(i)   = vint(i,1) * vint(i,5) / (gamma * Ma**2)

          vn(i)  = pt5 * ( R1(i) + R2(i) )
          vt(i)  = m2l(i) * vint(i,2) - m1l(i) * vint(i,3)

          fact1(i) = (rho(i)**gamma/(gamma*p(i)))**(one/gamma1) * &
                     (gamma1/four)**(two/gamma1) * (two/gamma1) * &
                     (R2(i) - R1(i))**((two-gamma1)/gamma1)
          fact2(i) = (gamma1/four * (R2(i) - R1(i)))**(two/gamma1) * &
                     one/(gamma*gamma1) * &
                     (rho(i)**gamma/(gamma*p(i)))**((two-gamma)/gamma1)
          mat(3,1,1,i,ny) = -one
          mat(2,1,2,i,ny) = (fact1(i) * m1l(i) + &
                            gamma1 * rho(i)**gamma1 / p(i))
          mat(2,1,3,i,ny) = (fact1(i) * m2l(i))
          mat(2,1,5,i,ny) = (fact1(i) / ( gamma1 * Ma * &
                            sqrt(vint(i,5)) ) - &
                            fact2(i) * rho(i)**(gamma+1)/(gamma*Ma**2*p(i)**2))
          
          fact1(i) = pt5 * m1l(i) / (m1l(i)**2 + m2l(i)**2)
          fact2(i) =       m2l(i) / (m1l(i)**2 + m2l(i)**2)
          mat(3,2,2,i,ny) = -one
          mat(2,2,2,i,ny) = (fact1(i) * m1l(i) + fact2(i) * m2l(i))
          mat(2,2,3,i,ny) = (fact1(i) * m2l(i) - fact2(i) * m1l(i))
          mat(2,2,5,i,ny) = (fact1(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
          fact1(i) = pt5 * m2l(i) / (m1l(i)**2 + m2l(i)**2)
          fact2(i) =       m1l(i) / (m1l(i)**2 + m2l(i)**2)
          mat(3,3,3,i,ny) = -one
          mat(2,3,2,i,ny) = (fact1(i) * m1l(i) - fact2(i) * m2l(i))
          mat(2,3,3,i,ny) = (fact1(i) * m2l(i) + fact2(i) * m1l(i))
          mat(2,3,5,i,ny) = (fact1(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
          mat(3,4,4,i,ny) = -one
          mat(2,4,4,i,ny) = one
  
          fact(i) = (Ma * gamma1 / four)**2 * two * ( R2(i) - R1(i) )
          mat(3,5,5,i,ny) = -one
          mat(2,5,2,i,ny) = (fact(i) * m1l(i))
          mat(2,5,3,i,ny) = (fact(i) * m2l(i))
          mat(2,5,5,i,ny) = (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
        
        end if

        end do

        end if
        
        else            ! Mach > 0
        
!.... if supersonic use freestream hard boundary condition

          mat(:,:,:,:,ny) = zero
          mat(3,1,1,:,ny) = one
          mat(3,2,2,:,ny) = one
          mat(3,3,3,:,ny) = one
          mat(3,4,4,:,ny) = one
          mat(3,5,5,:,ny) = one
        
        end if          ! Mach

        end if          ! yper
        
!=============================================================================!
!       L e f t   a n d   R i g h t   B o u n d a r i e s
!=============================================================================!

        if (.not. xper) then

!.... constant pressure

!         mat(:,1,:,1,:) = zero
!         mat(3,1,1,1,:) = -one
!         mat(3,1,5,1,:) = -one / vl(5,1,:)**2

!         mat(:,1,:,nx,:) = zero
!         mat(3,1,1,nx,:) = -one
!         mat(3,1,5,nx,:) = -one / vl(5,nx,:)**2

        if (Ma.lt.one) then

!.... Riemann invariants or extrapolation on the x-boundaries

          if (left.eq.0) then
            mat(:,:,:,1,:) = zero
            mat(3,1,1,1,:) = one
            mat(3,2,2,1,:) = one
            mat(3,3,3,1,:) = one
            mat(3,4,4,1,:) = one
            mat(3,5,5,1,:) = one
          end if

!.... symmetry plane

          if (left.eq.2) then
            mat(:,1:5,:,1,:) = zero
            mat(3,1,1,1,:) = one
            mat(3,2,2,1,:) = one
            mat(3,3,3,1,:) = one
            mat(3,4,4,1,:) = one
            mat(3,5,5,1,:) = one
          end if

          if (left.eq.7) then           ! symmetry plane
            mat(:,3,:,1,:) = zero
            mat(3,3,3,1,:) = one
          end if
          
          if (right.eq.0 .or. right.eq.8) then
            mat(:,:,:,nx,:) = zero
            mat(3,1,1,nx,:) = one
            mat(3,2,2,nx,:) = one
            mat(3,3,3,nx,:) = one
            mat(3,4,4,nx,:) = one
            mat(3,5,5,nx,:) = one
          end if
        
        else
        
          mat(:,:,:,1,1:nbl) = zero
          mat(3,1,1,1,1:nbl) = one
          mat(3,2,2,1,1:nbl) = one
          mat(3,3,3,1,1:nbl) = one
          mat(3,4,4,1,1:nbl) = one
          mat(3,5,5,1,1:nbl) = one

          mat(:,:,:,nx,1:nbl) = zero
          mat(3,1,1,nx,1:nbl) = one
          mat(3,2,2,nx,1:nbl) = one
          mat(3,3,3,nx,1:nbl) = one
          mat(3,4,4,nx,1:nbl) = one
          mat(3,5,5,nx,1:nbl) = one

        end if

        end if          ! xper

        end if          ! linear

        return
        end

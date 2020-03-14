!=============================================================================!
        subroutine genB(B, v)
!
!  Form the B matrix which multiplies U,y
!
!=============================================================================!
        use global
        use local
        implicit none
        
        real :: B(ny*nx,ndof,ndof), v(ny*nx,ndof)
        real :: fact1(ny*nx), fact2(ny*nx), gmsinv
!=============================================================================!

!.... Continuity equation

        B(:,1,1) = u2
        B(:,1,2) = zero
        B(:,1,3) = rho
        B(:,1,4) = zero
        B(:,1,5) = zero

!=============================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re
        gmsinv = one / (gamma * Ma**2)

        B(:,2,1) = zero
        B(:,2,2) = u2 - fact1 * g2mu
        B(:,2,3) = -fact1 * g1lm
        B(:,2,4) = zero
        B(:,2,5) = -fact1 * dmu * two * S(:,1,2)

!=============================================================================!

!.... Momentum equation -- x_2

        B(:,3,1) = rhoinv * t * gmsinv
        B(:,3,2) = -fact1 * g1mu
        B(:,3,3) = u2 - fact1 * g2lm - fact1 * two * g2mu
        B(:,3,4) = -fact1 * g3mu
        B(:,3,5) = gmsinv - fact1 * dlm * divu - &
                   fact1 * dmu * two * S(:,2,2)

!=============================================================================!

!.... Momentum equation -- x_3 

        B(:,4,1) = zero
        B(:,4,2) = zero
        B(:,4,3) = -fact1 * g3lm
        B(:,4,4) = u2 - fact1 * g2mu
        B(:,4,5) = -fact1 * dmu * two * S(:,3,2)

!=============================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        B(:,5,1) =  zero 
        B(:,5,2) = -four * fact2 * mu * S(:,1,2)
        B(:,5,3) =  gamma1 * t - fact2 * two * lm * divu - &
                    four * fact2 * mu * S(:,2,2)
        B(:,5,4) = -four * fact2 * mu * S(:,3,2)
        B(:,5,5) =  u2 - fact1 * (g2con + dcon * gt(:,2))

!.... correct for the Lele & Poinsot boundary conditions

        if (linear.eq.0) call genB_n(B, v)
        if (linear.eq.1) call genB_l(B, v)

        return
        end

!=============================================================================!
        subroutine bcB(B)
!
!  WARNING: not updated with new equations
!
!  Apply the viscous BC to outflow
!
!=============================================================================!
        use global
        use local
        implicit none
        
        real    :: B(ny*nx,ndof,ndof)
        real    :: fact1, fact2, fact3, gmsinv, ginv

        integer :: i, j, ji
!=============================================================================!
        fact1  = one / Re
        fact2  = one / Re
        gmsinv = one / (gamma * Ma**2)
!=============================================================================!
        write(*,*) 'WARNING: not updated with new equations'

        do j = 1, ny
          i = 1
          ji = j + (i-1)*ny

          B(ji,2,3) = zero

          B(ji,3,2) = zero

          i = nx
          ji = j + (i-1)*ny

          B(ji,2,3) = zero

          B(ji,3,2) = zero
        end do

        return
        end

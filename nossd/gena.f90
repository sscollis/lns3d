!=============================================================================!
        subroutine genA(A, v)
!
!  Form the A matrix which multiplies U,x
!
!=============================================================================!
        use global
        use local
        implicit none
        
        real    :: A(ny*nx,ndof,ndof), v(ny*nx,ndof)
        real    :: fact1(ny*nx), fact2(ny*nx), gmsinv
!=============================================================================!

!.... Continuity equation

        A(:,1,1) = u1
        A(:,1,2) = rho
        A(:,1,3) = zero
        A(:,1,4) = zero
        A(:,1,5) = zero

!=============================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re
        gmsinv = one / (gamma * Ma**2)

        A(:,2,1) = rhoinv * t * gmsinv
        A(:,2,2) = u1 - fact1 * g1lm - fact1 * two * g1mu
        A(:,2,3) = -fact1 * g2mu
        A(:,2,4) = -fact1 * g3mu
        A(:,2,5) = gmsinv - fact1 * dlm * divu - &
                   fact1 * dmu * two * S(:,1,1)

!=============================================================================!

!.... Momentum equation -- x_2

        A(:,3,1) = zero
        A(:,3,2) = -fact1 * g2lm
        A(:,3,3) = u1 - fact1 * g1mu
        A(:,3,4) = zero
        A(:,3,5) = -fact1 * dmu * two * S(:,2,1)

!=============================================================================!

!.... Momentum equation -- x_3

        A(:,4,1) = zero
        A(:,4,2) = -fact1 * g3lm
        A(:,4,3) = zero
        A(:,4,4) = u1 - fact1 * g1mu
        A(:,4,5) = -fact1 * dmu * two * S(:,3,1)

!=============================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        A(:,5,1) =  zero
        A(:,5,2) =  gamma1 * t - fact2 * two * lm * divu - &
                    four * fact2 * mu * S(:,1,1)
        A(:,5,3) = -four * fact2 * mu * S(:,2,1)
        A(:,5,4) = -four * fact2 * mu * S(:,3,1)
        A(:,5,5) = u1 - fact1 * (g1con + dcon * gt(:,1))

!.... correct for the Lele & Poinsot boundary conditions

        if (linear.eq.0) call genA_n(A, v)
        if (linear.eq.1) call genA_l(A, v)

        return
        end

!=============================================================================!
        subroutine bcA(A)
!
!  WARNING: Not updated to new equations
!
!  Apply the viscous BC to outflow
!
!=============================================================================!
        use global
        use local
        implicit none
        
        real    :: A(ny*nx,ndof,ndof)
        real    :: fact1, fact2, fact3, gmsinv, ginv

        integer :: i, j, ji
!=============================================================================!
        fact1  = one / Re
        fact2  = one / Re
        gmsinv = one / (gamma * Ma**2)
        ginv   = one / gamma
!=============================================================================!
        write(*,*) 'WARNING: Not updated to new equations'

        do j = 1, ny
          i = 1
          ji = j + (i-1)*ny

          A(ji,2,2) = rho(ji) * u1(ji)
          A(ji,2,5) = rho(ji) * gmsinv

          A(ji,3,3) = rho(ji) * u1(ji)
          A(ji,3,5) = zero

          A(ji,4,4) = rho(ji) * u1(ji)
          A(ji,4,5) = zero

          A(ji,5,5) = rho(ji) * u1(ji) * ginv

          i = nx
          ji = j + (i-1)*ny

          A(ji,2,2) = rho(ji) * u1(ji)
          A(ji,2,5) = rho(ji) * gmsinv

          A(ji,3,3) = rho(ji) * u1(ji)
          A(ji,3,5) = zero

          A(ji,4,4) = rho(ji) * u1(ji)
          A(ji,4,5) = zero

          A(ji,5,5) = rho(ji) * u1(ji) * ginv

        end do

        return
        end

!=============================================================================!
        subroutine genD(D, v)
!
!  Form the D matrix which multiplies U
!
!=============================================================================!
        use stuff
        use local
        implicit none
        
        real :: D(ny*nx,ndof,ndof), v(ny*nx,ndof)
        real :: fact1(ny*nx), fact2(ny*nx), gmsinv
!=============================================================================!

!.... Continuity equation

        D(:,1,1) = divu
        D(:,1,2) = grho(:,1)
        D(:,1,3) = grho(:,2)
        D(:,1,4) = grho(:,3)
        D(:,1,5) = zero

!=============================================================================!

        fact1 = rhoinv / Re
        gmsinv = one / (gamma * Ma**2)

!.... Momentum equation -- x_1

        if (linear.eq.1) then
          D(:,2,1) = rhoinv * ( u1 * gu(:,1,1) + u2 * gu(:,1,2) +       &
                                u3 * gu(:,1,3) + gmsinv * gt(:,1) )
        else
          D(:,2,1) = rhoinv * ( gmsinv * gt(:,1) - rhoinv * gp(:,1)  +  &
                     fact1 * ( g1lm * divu + lm * g1divu ) +            &
                     fact1 * two * ( g1mu * S(:,1,1) + &
                                     g2mu * S(:,1,2) + &
                                     g3mu * S(:,1,3) + &
                                       mu * S1jj ) )
        end if  

        D(:,2,2) = gu(:,1,1)
        D(:,2,3) = gu(:,1,2)
        D(:,2,4) = gu(:,1,3)
        D(:,2,5) = rhoinv * gmsinv * grho(:,1) -                &
                   fact1 * ( g1dlm * divu + dlm * g1divu ) -    &
                   fact1 * two * ( g1dmu * S(:,1,1) + &
                                   g2dmu * S(:,1,2) + &
                                   g3dmu * S(:,1,3) + &
                                     dmu * S1jj )

!=============================================================================!

!.... Momentum equation -- x_2

        if (linear.eq.1) then
          D(:,3,1) = rhoinv * ( u1 * gu(:,2,1) + u2 * gu(:,2,2) + &
                                u3 * gu(:,2,3) + gmsinv * gt(:,2) )
        else
          D(:,3,1) = rhoinv * ( gmsinv * gt(:,2) - rhoinv * gp(:,2) + &
                     fact1 * ( g2lm * divu + lm * g2divu ) + &
                     fact1 * two * ( g1mu * S(:,2,1) + &
                                     g2mu * S(:,2,2) + &
                                     g3mu * S(:,2,3) + &
                                       mu * S2jj  ) )
        end if  

        D(:,3,2) = gu(:,2,1)
        D(:,3,3) = gu(:,2,2)
        D(:,3,4) = gu(:,2,3)
        D(:,3,5) = rhoinv * gmsinv * grho(:,2)  - &
                   fact1 * ( g2dlm * divu + dlm * g2divu ) - &
                   fact1 * two * ( g1dmu * S(:,2,1) + &
                                   g2dmu * S(:,2,2) + &
                                   g3dmu * S(:,2,3) + &
                                     dmu * S2jj )

!=============================================================================!

!.... Momentum equation -- x_3

        if (linear.eq.1) then
          D(:,4,1) = rhoinv * ( u1 * gu(:,3,1) + u2 * gu(:,3,2) + &
                                u3 * gu(:,3,3) + gmsinv * gt(:,3) )
        else
          D(:,4,1) = rhoinv * ( gmsinv * gt(:,3) - rhoinv * gp(:,3) + &
                     fact1 * ( g3lm * divu + lm * g3divu ) + &
                     fact1 * two * ( g1mu * S(:,3,1) + &
                                     g2mu * S(:,3,2) + &
                                     g3mu * S(:,3,3) + &
                                       mu * S3jj ) )
        end if  

        D(:,4,2) = gu(:,3,1)
        D(:,4,3) = gu(:,3,2)
        D(:,4,4) = gu(:,3,3)
        D(:,4,5) = rhoinv * gmsinv * grho(:,3)  - &
                   fact1 * ( g3dlm * divu + dlm * g3divu ) - &
                   fact1 * two * ( g1dmu * S(:,3,1) + &
                                   g2dmu * S(:,3,2) + &
                                   g3dmu * S(:,3,3) + &
                                     dmu * S3jj )

!=============================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        if (linear .eq. 1) then
          D(:,5,1) = rhoinv * ( u1 * gt(:,1) + u2 * gt(:,2) + &
                                u3 * gt(:,3) + gamma1 * t * divu )
        else
          D(:,5,1) = fact1 * rhoinv * ( g1con * gt(:,1) + &
                                        g2con * gt(:,2) + &
                                        g3con * gt(:,3) + &
                                          con * Lapt ) + &
                     fact2 * rhoinv * lm * divu**2 + &
                     two * fact2 * rhoinv * mu * ( &
                        S(:,1,1)**2 + S(:,1,2)**2 + S(:,1,3)**2 + &
                        S(:,2,1)**2 + S(:,2,2)**2 + S(:,2,3)**2 + &
                        S(:,3,1)**2 + S(:,3,2)**2 + S(:,3,3)**2 )

        end if
        D(:,5,2) = gt(:,1)
        D(:,5,3) = gt(:,2)
        D(:,5,4) = gt(:,3)

!.... to get the linear verion, it is assumed that \rho,t = 0
!.... Note that it really isn't necessary, it would be okay to always
!.... use the nonlinear version.

        if (linear .eq. 1) then
          D(:,5,5) = -gamma1 * rhoinv * ( u1 * grho(:,1) + &
                      u2 * grho(:,2) + u3 * grho(:,3) )
        else
          D(:,5,5) = gamma1 * divu
        end if
        
        D(:,5,5) = D(:,5,5) - &
                   fact1 * ( g1dcon * gt(:,1) + &
                             g2dcon * gt(:,2) + &
                             g3dcon * gt(:,3) + &
                               dcon * Lapt ) - &
                   fact2 * dlm * divu**2 - &
                   two * fact2 * dmu * ( &
                        S(:,1,1)**2 + S(:,1,2)**2 + S(:,1,3)**2 + &
                        S(:,2,1)**2 + S(:,2,2)**2 + S(:,2,3)**2 + &
                        S(:,3,1)**2 + S(:,3,2)**2 + S(:,3,3)**2 )

!.... correct for the Lele & Poinsot boundary conditions

        if (linear.eq.0) call genD_n(D, v)
        if (linear.eq.1) call genD_l(D, v)

        return
        end

!=============================================================================!
        subroutine bcD(D)
!
!  WARNING: not updated with new equations
!
!  Apply the viscous BC to the outflow
!
!=============================================================================!
        use stuff
        use local
        implicit none
        
        real    :: D(ny*nx,ndof,ndof)
        real    :: fact1, fact2, fact3, gmsinv, ginv, g1ms

        integer :: i, j, ji
!=============================================================================!
        write(*,*) 'WARNING: not updated with new equations'

        do j = 1, ny
          i = 1
          ji = j + (i-1)*ny

          fact1 = one / Re
          fact2 = two / Re
          gmsinv = one / (gamma * Ma**2)

          D(ji,2,5) = grho(ji,1) * gmsinv - &
            fact2 * ( g2dmu(ji) * S(ji,1,2) + g3dmu(ji) * S(ji,1,3) + &
            dmu(ji) * pt5 * ( g22v(ji,2) + g12v(ji,3) ) )

          D(ji,3,5) = grho(ji,2) * gmsinv - &
            fact1 * ( g2dlm(ji) * divu(ji) + dlm(ji) * g2divu(ji) ) - &
            fact2 * ( g2dmu(ji) * S(ji,2,2) + g3dmu(ji) * S(ji,2,3) + &
            dmu(ji) * pt5 * g22v(ji,3) )

          D(ji,4,5) = grho(ji,3) * gmsinv - &
            fact1 * ( g3dlm(ji) * divu(ji) + dlm(ji) * g3divu(ji) ) - &
            fact2 * ( g2dmu(ji) * S(ji,3,2) + g3dmu(ji) * S(ji,3,3) + &
            dmu(ji) * pt5 * g22v(ji,4) )

          fact1 = one / (Pr * Re)
          fact2 = gamma1 * Ma**2 / Re
          fact3 = two * gamma1 * Ma**2 / Re
          ginv  = one / gamma
          g1ms  = gamma1 * Ma**2

          D(ji,5,5) = -gamma1 * ginv * ( u1(ji) * grho(ji,1) + &
                u2(ji) * grho(ji,2) + u3(ji) * grho(ji,3) ) - &
                fact1 * (g2dcon(ji) * gt(ji,2) + g3dcon(ji) * gt(ji,3) + &
                         dcon(ji) * g22v(ji,5) ) - &
                fact2 * dlm(ji) * divu(ji)**2 - &
                fact3 * dmu(ji) * ( &
                        S(ji,1,1)**2 + S(ji,1,2)**2 + S(ji,1,3)**2 + &
                        S(ji,2,1)**2 + S(ji,2,2)**2 + S(ji,2,3)**2 + &
                        S(ji,3,1)**2 + S(ji,3,2)**2 + S(ji,3,3)**2 )

          i = nx
          ji = j + (i-1)*ny

          fact1 = one / Re
          fact2 = two / Re
          gmsinv = one / (gamma * Ma**2)

          D(ji,2,5) = grho(ji,1) * gmsinv - &
            fact2 * ( g2dmu(ji) * S(ji,1,2) + g3dmu(ji) * S(ji,1,3) + &
            dmu(ji) * pt5 * ( g22v(ji,2) + g12v(ji,3) ) )

          D(ji,3,5) = grho(ji,2) * gmsinv - &
            fact1 * ( g2dlm(ji) * divu(ji) + dlm(ji) * g2divu(ji) ) - &
            fact2 * ( g2dmu(ji) * S(ji,2,2) + g3dmu(ji) * S(ji,2,3) + &
            dmu(ji) * pt5 * g22v(ji,3) )

          D(ji,4,5) = grho(ji,3) * gmsinv - &
            fact1 * ( g3dlm(ji) * divu(ji) + dlm(ji) * g3divu(ji) ) -   &
            fact2 * ( g2dmu(ji) * S(ji,3,2) + g3dmu(ji) * S(ji,3,3) +   &
            dmu(ji) * pt5 * g22v(ji,4) )

          fact1 = one / (Pr * Re)
          fact2 = gamma1 * Ma**2 / Re
          fact3 = two * gamma1 * Ma**2 / Re
          ginv  = one / gamma
          g1ms  = gamma1 * Ma**2

          D(ji,5,5) = -gamma1 * ginv * ( u1(ji) * grho(ji,1) + &
                u2(ji) * grho(ji,2) + u3(ji) * grho(ji,3) ) - &
                fact1 * (g2dcon(ji) * gt(ji,2) + g3dcon(ji) * gt(ji,3) + &
                         dcon(ji) * g22v(ji,5) ) - &
                fact2 * dlm(ji) * divu(ji)**2 - &
                fact3 * dmu(ji) * ( &
                        S(ji,1,1)**2 + S(ji,1,2)**2 + S(ji,1,3)**2 + &
                        S(ji,2,1)**2 + S(ji,2,2)**2 + S(ji,2,3)**2 + &
                        S(ji,3,1)**2 + S(ji,3,2)**2 + S(ji,3,3)**2 )

        end do

        return
        end

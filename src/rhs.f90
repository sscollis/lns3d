!=============================================================================!
        subroutine rhs(rl, vl) 
!  
!       Generates the RHS for the compressible 2D-3C, N-S equations 
!
!       Written: 5-10-95
!       Revised: 6-19-95
!
!=============================================================================!
        use global
        use local
        implicit none

        integer :: ij
        real :: rl(ndof,nx*ny), vl(ndof,nx*ny)
        real :: fact1, fact2, fact3
!=============================================================================!

        fact1 = one / Re
        fact2 = gamma / ( Pr * Re )
        fact3 = gamma * gamma1 * Ma**2 / Re

!$doacross local(ij)
        do ij = 1, nx*ny

!.... continuity

        rl(1,ij) = rho(ij) * divu(ij) + grho(1,ij) * u1(ij) + &
                   grho(2,ij) * u2(ij)

!.... momentum

        rl(2,ij) = u1(ij) * gu(1,1,ij) + u2(ij) * gu(1,2,ij) +           &
                   rhoinv(ij) * gp(1,ij) -                               &
                   fact1 * ( g1lm(ij) * divu(ij) + lm(ij) * g1divu(ij) + &
                   two * (g1mu(ij) * S(1,1,ij) + g2mu(ij) * S(1,2,ij)  + &
                                   mu(ij) * S1jj(ij) ) ) * rhoinv(ij)

        rl(3,ij) = u1(ij) * gu(2,1,ij) + u2(ij) * gu(2,2,ij) +           &
                   rhoinv(ij) * gp(2,ij) -                               &
                   fact1 * ( g2lm(ij) * divu(ij) + lm(ij) * g2divu(ij) + &
                   two * (g1mu(ij) * S(2,1,ij) + g2mu(ij) * S(2,2,ij) +  &
                                   mu(ij) * S2jj(ij) ) ) * rhoinv(ij)

        rl(4,ij) = u1(ij) * gu(3,1,ij) + u2(ij) * gu(3,2,ij) -           &
                   fact1 * ( two * (g1mu(ij) * S(3,1,ij) +               &
                            g2mu(ij) * S(3,2,ij) + &
                                   mu(ij) * S3jj(ij) ) ) * rhoinv(ij)

!.... temperature

        rl(5,ij) = u1(ij) * gt(1,ij) + u2(ij) * gt(2,ij) +                &
                  gamma1 * t(ij) * divu(ij) - rhoinv(ij) * (              &
                  fact2 * ( g1con(ij) * gt(1,ij) + g2con(ij) * gt(2,ij) + &
                            con(ij) * Lapt(ij) ) +                        &
                  fact3 * ( lm(ij) * divu(ij)**2 + two * mu(ij) * (       &
                        S(1,1,ij)**2 + S(1,2,ij)**2 + S(1,3,ij)**2 +      &
                        S(2,1,ij)**2 + S(2,2,ij)**2 + S(2,3,ij)**2 +      &
                        S(3,1,ij)**2 + S(3,2,ij)**2 ) ) )

        end do

!.... correct the outflow boundaries using Lele & Poinsot BC's

!       call rhs_l(rl,vl)
        
!.... standard sponge

        if (ispg.gt.0) then
          if (ispg.eq.1) then
            call spg_it(rl,vl,spg)
          else if (ispg.eq.2) then
            !$doacross local(ij)
            do ij = 1, nx*ny
              rl(1,ij) = rl(1,ij) + (spg(ij) + spg2(ij)) * ( vl(1,ij) - one )
              rl(2,ij) = rl(2,ij) + (spg(ij) + spg2(ij)) * ( vl(2,ij) )
              rl(3,ij) = rl(3,ij) + (spg(ij) + spg2(ij)) * ( vl(3,ij) )
              rl(4,ij) = rl(4,ij) + (spg(ij) + spg2(ij)) * ( vl(4,ij) )
              rl(5,ij) = rl(5,ij) + (spg(ij) + spg2(ij)) * ( vl(5,ij) - one )
            end do
          else
            call error('rhs_p$','ispg > 2 is not supported$')
          end if
        end if

!.... explicit smoother

!       if (eps_e.ne.zero) call smoother( rl, vl )
        
        return
        end

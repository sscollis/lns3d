!=============================================================================!
        subroutine rk(vl, rl, vml, xl, yl, dtl)
!  
!  Computes the RK3 RHS for the linearized compressible Navier-Stokes solver. 
!  Updated to use the hat matrices
!
!  Revised: 10-25-95
!=============================================================================!
        use global
        use local2d
        implicit none
!=============================================================================!

        real :: vl(ndof,nx,ny), rl(ndof,nx,ny), vml(ndof,nx,ny)
        real :: xl(nx,ny), yl(nx,ny), dtl(nx,ny)

!.... local variables
        
        integer :: i, j, idof, jdof
        
!.... mean flow variables

        real :: cm, c3, um
        real :: a, d, kk

        character*80 :: name, code='rk$'

        real, external :: ramp
!=============================================================================!

!.... compute first derivatives
        
        call grad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... compute second derivatives

        call grad2(ndof, nx, ny, vl, g1v, g11v, g12v, g22v, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall temperature boundary condition

        call gradbc( g1v, g2v, g11v, g12v, g22v )

        !$omp parallel do private (i,cm,um,d,kk,a,c3)
        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!=============================================================================!
!.... form the RHS
!=============================================================================!

!.... U,xi term

        rl(1,i,j) =       - Ah(1,1,i,j) * g1v(1,i,j)    &
                          - Ah(1,2,i,j) * g1v(2,i,j)    &
                          - Ah(1,3,i,j) * g1v(3,i,j)    &
                          - Ah(1,4,i,j) * g1v(4,i,j)    &
                          - Ah(1,5,i,j) * g1v(5,i,j)    

        rl(2,i,j) =       - Ah(2,1,i,j) * g1v(1,i,j)    &
                          - Ah(2,2,i,j) * g1v(2,i,j)    &
                          - Ah(2,3,i,j) * g1v(3,i,j)    &
                          - Ah(2,4,i,j) * g1v(4,i,j)    &
                          - Ah(2,5,i,j) * g1v(5,i,j)

        rl(3,i,j) =       - Ah(3,1,i,j) * g1v(1,i,j)    &
                          - Ah(3,2,i,j) * g1v(2,i,j)    &
                          - Ah(3,3,i,j) * g1v(3,i,j)    &
                          - Ah(3,4,i,j) * g1v(4,i,j)    &
                          - Ah(3,5,i,j) * g1v(5,i,j)

        rl(4,i,j) =       - Ah(4,1,i,j) * g1v(1,i,j)    &
                          - Ah(4,2,i,j) * g1v(2,i,j)    &
                          - Ah(4,3,i,j) * g1v(3,i,j)    &
                          - Ah(4,4,i,j) * g1v(4,i,j)    &
                          - Ah(4,5,i,j) * g1v(5,i,j)

        rl(5,i,j) =       - Ah(5,1,i,j) * g1v(1,i,j)    &
                          - Ah(5,2,i,j) * g1v(2,i,j)    &
                          - Ah(5,3,i,j) * g1v(3,i,j)    &
                          - Ah(5,4,i,j) * g1v(4,i,j)    &
                          - Ah(5,5,i,j) * g1v(5,i,j)

!.... U,eta term

        rl(1,i,j) = rl(1,i,j) - Bh(1,1,i,j) * g2v(1,i,j)        &
                          - Bh(1,2,i,j) * g2v(2,i,j)    &
                          - Bh(1,3,i,j) * g2v(3,i,j)    &
                          - Bh(1,4,i,j) * g2v(4,i,j)    &
                          - Bh(1,5,i,j) * g2v(5,i,j)    

        rl(2,i,j) = rl(2,i,j) - Bh(2,1,i,j) * g2v(1,i,j)        &
                          - Bh(2,2,i,j) * g2v(2,i,j)    &
                          - Bh(2,3,i,j) * g2v(3,i,j)    &
                          - Bh(2,4,i,j) * g2v(4,i,j)    &
                          - Bh(2,5,i,j) * g2v(5,i,j)

        rl(3,i,j) = rl(3,i,j) - Bh(3,1,i,j) * g2v(1,i,j)        &
                          - Bh(3,2,i,j) * g2v(2,i,j)    &
                          - Bh(3,3,i,j) * g2v(3,i,j)    &
                          - Bh(3,4,i,j) * g2v(4,i,j)    &
                          - Bh(3,5,i,j) * g2v(5,i,j)

        rl(4,i,j) = rl(4,i,j) - Bh(4,1,i,j) * g2v(1,i,j)        &
                          - Bh(4,2,i,j) * g2v(2,i,j)    &
                          - Bh(4,3,i,j) * g2v(3,i,j)    &
                          - Bh(4,4,i,j) * g2v(4,i,j)    &
                          - Bh(4,5,i,j) * g2v(5,i,j)

        rl(5,i,j) = rl(5,i,j) - Bh(5,1,i,j) * g2v(1,i,j)        &
                          - Bh(5,2,i,j) * g2v(2,i,j)    &
                          - Bh(5,3,i,j) * g2v(3,i,j)    &
                          - Bh(5,4,i,j) * g2v(4,i,j)    &
                          - Bh(5,5,i,j) * g2v(5,i,j)

!.... U term

        rl(1,i,j) = rl(1,i,j) - Dh(1,1,i,j) * vl(1,i,j) &
                          - Dh(1,2,i,j) * vl(2,i,j)     &
                          - Dh(1,3,i,j) * vl(3,i,j)     &
                          - Dh(1,4,i,j) * vl(4,i,j)     &
                          - Dh(1,5,i,j) * vl(5,i,j)

        rl(2,i,j) = rl(2,i,j) - Dh(2,1,i,j) * vl(1,i,j) &
                          - Dh(2,2,i,j) * vl(2,i,j)     &
                          - Dh(2,3,i,j) * vl(3,i,j)     &
                          - Dh(2,4,i,j) * vl(4,i,j)     &
                          - Dh(2,5,i,j) * vl(5,i,j)

        rl(3,i,j) = rl(3,i,j) - Dh(3,1,i,j) * vl(1,i,j) &
                          - Dh(3,2,i,j) * vl(2,i,j)     &
                          - Dh(3,3,i,j) * vl(3,i,j)     &
                          - Dh(3,4,i,j) * vl(4,i,j)     &
                          - Dh(3,5,i,j) * vl(5,i,j)

        rl(4,i,j) = rl(4,i,j) - Dh(4,1,i,j) * vl(1,i,j) &
                          - Dh(4,2,i,j) * vl(2,i,j)     &
                          - Dh(4,3,i,j) * vl(3,i,j)     &
                          - Dh(4,4,i,j) * vl(4,i,j)     &
                          - Dh(4,5,i,j) * vl(5,i,j)

        rl(5,i,j) = rl(5,i,j) - Dh(5,1,i,j) * vl(1,i,j) &
                          - Dh(5,2,i,j) * vl(2,i,j)     &
                          - Dh(5,3,i,j) * vl(3,i,j)     &
                          - Dh(5,4,i,j) * vl(4,i,j)     &
                          - Dh(5,5,i,j) * vl(5,i,j)

!.... U,\xi\xi term

        rl(2,i,j) = rl(2,i,j) + Vh11(1,i,j) * g11v(2,i,j)       &
                              + Vh11(5,i,j) * g11v(3,i,j)
        rl(3,i,j) = rl(3,i,j) + Vh11(6,i,j) * g11v(2,i,j)       &
                              + Vh11(2,i,j) * g11v(3,i,j)
        rl(4,i,j) = rl(4,i,j) + Vh11(3,i,j) * g11v(4,i,j)
        rl(5,i,j) = rl(5,i,j) + Vh11(4,i,j) * g11v(5,i,j)

!.... U,\xi\eta term

        rl(2,i,j) = rl(2,i,j) + Vh12(1,i,j) * g12v(2,i,j)       &
                              + Vh12(5,i,j) * g12v(3,i,j)
        rl(3,i,j) = rl(3,i,j) + Vh12(6,i,j) * g12v(2,i,j)       &
                              + Vh12(2,i,j) * g12v(3,i,j)
        rl(4,i,j) = rl(4,i,j) + Vh12(3,i,j) * g12v(4,i,j)
        rl(5,i,j) = rl(5,i,j) + Vh12(4,i,j) * g12v(5,i,j)

!.... U,\eta\eta term

        rl(2,i,j) = rl(2,i,j) + Vh22(1,i,j) * g22v(2,i,j)       &
                              + Vh22(5,i,j) * g22v(3,i,j)
        rl(3,i,j) = rl(3,i,j) + Vh22(6,i,j) * g22v(2,i,j)       &
                              + Vh22(2,i,j) * g22v(3,i,j)
        rl(4,i,j) = rl(4,i,j) + Vh22(3,i,j) * g22v(4,i,j)
        rl(5,i,j) = rl(5,i,j) + Vh22(4,i,j) * g22v(5,i,j)

!.... Sponge term

        if (ispg .eq. 1) then

!.... standard sponge

          rl(1,i,j) = rl(1,i,j) - spg(i,j) * vl(1,i,j)
          rl(2,i,j) = rl(2,i,j) - spg(i,j) * vl(2,i,j)
          rl(3,i,j) = rl(3,i,j) - spg(i,j) * vl(3,i,j)
          rl(4,i,j) = rl(4,i,j) - spg(i,j) * vl(4,i,j)
          rl(5,i,j) = rl(5,i,j) - spg(i,j) * vl(5,i,j)
        
        else if (ispg .eq. 2) then

          rl(1,i,j) = rl(1,i,j) - (spg(i,j) + spg2(i,j)) * vl(1,i,j)
          rl(2,i,j) = rl(2,i,j) - (spg(i,j) + spg2(i,j)) * vl(2,i,j)
          rl(3,i,j) = rl(3,i,j) - (spg(i,j) + spg2(i,j)) * vl(3,i,j)
          rl(4,i,j) = rl(4,i,j) - (spg(i,j) + spg2(i,j)) * vl(4,i,j)
          rl(5,i,j) = rl(5,i,j) - (spg(i,j) + spg2(i,j)) * vl(5,i,j)

        else if (ispg .eq. 3) then

!.... first do the outflow sponge

          rl(1,i,j) = rl(1,i,j) - spg(i,j) * vl(1,i,j)
          rl(2,i,j) = rl(2,i,j) - spg(i,j) * vl(2,i,j)
          rl(3,i,j) = rl(3,i,j) - spg(i,j) * vl(3,i,j)
          rl(4,i,j) = rl(4,i,j) - spg(i,j) * vl(4,i,j)
          rl(5,i,j) = rl(5,i,j) - spg(i,j) * vl(5,i,j)  

!.... subtract off the forcing wave in the inflow sponge

          cm = sqrt( vml(5,i,j) ) / Ma
          um = vml(2,i,j)

          d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

          kk = omega / (cm+um)
          a  = omega**2 * d / (cm+um)**3
          c3 = ramp( (kk*(x0-xl(i,j)) + omega * time)/(two*pi) ) * &
               cos( kk * xl(i,j) - omega * time ) * &
               exp( -a * (xl(i,j) - x0) )
!         c3 = ramp( (kk*(x0-xl(i,j)) + omega * time)/(two*pi) ) * &
!              wamp(i) * cos( kk * xl(i,j) - omega * time )

          rl(1,i,j) = rl(1,i,j) - spg2(i,j) * ( vl(1,i,j) - &
                      pt5 * c3 / cm**2 )
          rl(2,i,j) = rl(2,i,j) - spg2(i,j) * ( vl(2,i,j) - &
                      c3 * pt5 / ( vml(1,i,j) * cm ) )
          rl(3,i,j) = rl(3,i,j) - spg2(i,j) * ( vl(3,i,j) )
          rl(4,i,j) = rl(4,i,j) - spg2(i,j) * ( vl(4,i,j) )
          rl(5,i,j) = rl(5,i,j) - spg2(i,j) * ( vl(5,i,j) - &
                      (gamma*Ma**2 * c3 * pt5 - &
                      vml(5,i,j) * pt5 * c3 / cm**2) / vml(1,i,j) )
        end if

!.... multiply by the time step

        do idof = 1, ndof
          rl(idof,i,j) = dtl(i,j) * rl(idof,i,j)
        end do
        
        end do loop_i
        end do loop_j

!.... set boundary conditions

        call rkbc(rl)

        return
        end

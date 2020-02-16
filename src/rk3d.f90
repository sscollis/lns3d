!=============================================================================!
        subroutine rk3d(vl, rl, vml, xl, yl, dtl)
!  
!  Computes the RK3 RHS for the linearized three-dimensional compressible 
!  Navier-Stokes solver.  Updated to use the hat matrices
!
!  NOTE: the sign is opposite from lrhs!
!
!  Revised: 4-22-96
!=============================================================================!
        use global
        use local2d
        use local3d
        implicit none
!=============================================================================!
        complex :: vl(ndof,nx,ny), rl(ndof,nx,ny)
        real    :: vml(ndof,nx,ny)
        real    :: xl(nx,ny), yl(nx,ny), dtl(nx,ny)

!.... local variables
        
        integer :: i, j, idof, jdof
        
!.... mean flow variables

        real :: cm, um
        real :: a, d, kk
        complex :: c3

        character*80 :: name, code='rk3d$'
!=============================================================================!

!.... compute first derivatives
        
        call cgrad( ndof, nx, ny, vl, c1v, c2v, dxi, deta, optx, opty, &
                    xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... compute second derivatives

        call cgrad2( ndof, nx, ny, vl, c1v, c11v, c12v, c22v, dxi, deta, &
                     optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall temperature boundary condition

        if (Navier) call cgradbc( c1v, c2v, c11v, c12v, c22v )

        !$omp parallel do private (i,cm,um,d,kk,a,c3)
        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!=============================================================================!
!.... form the RHS
!=============================================================================!

!.... U,xi term

        rl(1,i,j) =       - Ah(1,1,i,j) * c1v(1,i,j)    &
                          - Ah(1,2,i,j) * c1v(2,i,j)    &
                          - Ah(1,3,i,j) * c1v(3,i,j)    &
                          - Ah(1,4,i,j) * c1v(4,i,j)    &
                          - Ah(1,5,i,j) * c1v(5,i,j)    

        rl(2,i,j) =       - Ah(2,1,i,j) * c1v(1,i,j)    &
                          - Ah(2,2,i,j) * c1v(2,i,j)    &
                          - Ah(2,3,i,j) * c1v(3,i,j)    &
                          - Ah(2,4,i,j) * c1v(4,i,j)    &
                          - Ah(2,5,i,j) * c1v(5,i,j)

        rl(3,i,j) =       - Ah(3,1,i,j) * c1v(1,i,j)    &
                          - Ah(3,2,i,j) * c1v(2,i,j)    &
                          - Ah(3,3,i,j) * c1v(3,i,j)    &
                          - Ah(3,4,i,j) * c1v(4,i,j)    &
                          - Ah(3,5,i,j) * c1v(5,i,j)

        rl(4,i,j) =       - Ah(4,1,i,j) * c1v(1,i,j)    &
                          - Ah(4,2,i,j) * c1v(2,i,j)    &
                          - Ah(4,3,i,j) * c1v(3,i,j)    &
                          - Ah(4,4,i,j) * c1v(4,i,j)    &
                          - Ah(4,5,i,j) * c1v(5,i,j)

        rl(5,i,j) =       - Ah(5,1,i,j) * c1v(1,i,j)    &
                          - Ah(5,2,i,j) * c1v(2,i,j)    &
                          - Ah(5,3,i,j) * c1v(3,i,j)    &
                          - Ah(5,4,i,j) * c1v(4,i,j)    &
                          - Ah(5,5,i,j) * c1v(5,i,j)

!.... U,eta term

        rl(1,i,j) = rl(1,i,j) - Bh(1,1,i,j) * c2v(1,i,j)        &
                          - Bh(1,2,i,j) * c2v(2,i,j)    &
                          - Bh(1,3,i,j) * c2v(3,i,j)    &
                          - Bh(1,4,i,j) * c2v(4,i,j)    &
                          - Bh(1,5,i,j) * c2v(5,i,j)    

        rl(2,i,j) = rl(2,i,j) - Bh(2,1,i,j) * c2v(1,i,j)        &
                          - Bh(2,2,i,j) * c2v(2,i,j)    &
                          - Bh(2,3,i,j) * c2v(3,i,j)    &
                          - Bh(2,4,i,j) * c2v(4,i,j)    &
                          - Bh(2,5,i,j) * c2v(5,i,j)

        rl(3,i,j) = rl(3,i,j) - Bh(3,1,i,j) * c2v(1,i,j)        &
                          - Bh(3,2,i,j) * c2v(2,i,j)    &
                          - Bh(3,3,i,j) * c2v(3,i,j)    &
                          - Bh(3,4,i,j) * c2v(4,i,j)    &
                          - Bh(3,5,i,j) * c2v(5,i,j)

        rl(4,i,j) = rl(4,i,j) - Bh(4,1,i,j) * c2v(1,i,j)        &
                          - Bh(4,2,i,j) * c2v(2,i,j)    &
                          - Bh(4,3,i,j) * c2v(3,i,j)    &
                          - Bh(4,4,i,j) * c2v(4,i,j)    &
                          - Bh(4,5,i,j) * c2v(5,i,j)

        rl(5,i,j) = rl(5,i,j) - Bh(5,1,i,j) * c2v(1,i,j)        &
                          - Bh(5,2,i,j) * c2v(2,i,j)    &
                          - Bh(5,3,i,j) * c2v(3,i,j)    &
                          - Bh(5,4,i,j) * c2v(4,i,j)    &
                          - Bh(5,5,i,j) * c2v(5,i,j)

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

        rl(2,i,j) = rl(2,i,j) + Vh11(1,i,j) * c11v(2,i,j)       &
                              + Vh11(5,i,j) * c11v(3,i,j)
        rl(3,i,j) = rl(3,i,j) + Vh11(6,i,j) * c11v(2,i,j)       &
                              + Vh11(2,i,j) * c11v(3,i,j)
        rl(4,i,j) = rl(4,i,j) + Vh11(3,i,j) * c11v(4,i,j)
        rl(5,i,j) = rl(5,i,j) + Vh11(4,i,j) * c11v(5,i,j)

!.... U,\xi\eta term

        rl(2,i,j) = rl(2,i,j) + Vh12(1,i,j) * c12v(2,i,j)       &
                              + Vh12(5,i,j) * c12v(3,i,j)
        rl(3,i,j) = rl(3,i,j) + Vh12(6,i,j) * c12v(2,i,j)       &
                              + Vh12(2,i,j) * c12v(3,i,j)
        rl(4,i,j) = rl(4,i,j) + Vh12(3,i,j) * c12v(4,i,j)
        rl(5,i,j) = rl(5,i,j) + Vh12(4,i,j) * c12v(5,i,j)

!.... U,\eta\eta term

        rl(2,i,j) = rl(2,i,j) + Vh22(1,i,j) * c22v(2,i,j)       &
                              + Vh22(5,i,j) * c22v(3,i,j)
        rl(3,i,j) = rl(3,i,j) + Vh22(6,i,j) * c22v(2,i,j)       &
                              + Vh22(2,i,j) * c22v(3,i,j)
        rl(4,i,j) = rl(4,i,j) + Vh22(3,i,j) * c22v(4,i,j)
        rl(5,i,j) = rl(5,i,j) + Vh22(4,i,j) * c22v(5,i,j)

!.... 3D terms

!.... \hat{A}_i, \hat{B}_i terms

        rl(2,i,j) = rl(2,i,j) + im * ( ABhi(1,i,j) * c1v(4,i,j) + &
                                       ABhi(3,i,j) * c2v(4,i,j) )
        rl(3,i,j) = rl(3,i,j) + im * ( ABhi(2,i,j) * c1v(4,i,j) + &
                                       ABhi(4,i,j) * c2v(4,i,j) )
        rl(4,i,j) = rl(4,i,j) + im * ( ABhi(1,i,j) * c1v(2,i,j) + &
                                       ABhi(2,i,j) * c1v(3,i,j) + &
                    ABhi(3,i,j) * c2v(2,i,j) + ABhi(4,i,j) * c2v(3,i,j) )
        
!.... \hat{D}_i term

        rl(1,i,j) = rl(1,i,j) - im * ( Dhi(1,1,i,j) * vl(1,i,j) &
                                 + Dhi(1,2,i,j) * vl(2,i,j)     &
                                 + Dhi(1,3,i,j) * vl(3,i,j)     &
                                 + Dhi(1,4,i,j) * vl(4,i,j)     &
                                 + Dhi(1,5,i,j) * vl(5,i,j) )

        rl(2,i,j) = rl(2,i,j) - im * ( Dhi(2,1,i,j) * vl(1,i,j) &
                                 + Dhi(2,2,i,j) * vl(2,i,j)     &
                                 + Dhi(2,3,i,j) * vl(3,i,j)     &
                                 + Dhi(2,4,i,j) * vl(4,i,j)     &
                                 + Dhi(2,5,i,j) * vl(5,i,j) )

        rl(3,i,j) = rl(3,i,j) - im * ( Dhi(3,1,i,j) * vl(1,i,j) &
                                 + Dhi(3,2,i,j) * vl(2,i,j)     &
                                 + Dhi(3,3,i,j) * vl(3,i,j)     &
                                 + Dhi(3,4,i,j) * vl(4,i,j)     &
                                 + Dhi(3,5,i,j) * vl(5,i,j) )

        rl(4,i,j) = rl(4,i,j) - im * ( Dhi(4,1,i,j) * vl(1,i,j) &
                                 + Dhi(4,2,i,j) * vl(2,i,j)     &
                                 + Dhi(4,3,i,j) * vl(3,i,j)     &
                                 + Dhi(4,4,i,j) * vl(4,i,j)     &
                                 + Dhi(4,5,i,j) * vl(5,i,j) )

        rl(5,i,j) = rl(5,i,j) - im * ( Dhi(5,1,i,j) * vl(1,i,j) &
                                 + Dhi(5,2,i,j) * vl(2,i,j)     &
                                 + Dhi(5,3,i,j) * vl(3,i,j)     &
                                 + Dhi(5,4,i,j) * vl(4,i,j)     &
                                 + Dhi(5,5,i,j) * vl(5,i,j) )

!.... Sponge term

        if (ispg .eq. 1) then           !.... standard sponge

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

        else if (ispg .eq. 3) then      !.... double sponge

!.... now do the outflow sponge

          rl(1,i,j) = rl(1,i,j) - spg2(i,j) * vl(1,i,j)
          rl(2,i,j) = rl(2,i,j) - spg2(i,j) * vl(2,i,j)
          rl(3,i,j) = rl(3,i,j) - spg2(i,j) * vl(3,i,j)
          rl(4,i,j) = rl(4,i,j) - spg2(i,j) * vl(4,i,j)
          rl(5,i,j) = rl(5,i,j) - spg2(i,j) * vl(5,i,j) 

!.... the inflow sponge

          cm = sqrt( vml(5,i,j) ) / Ma
          um = vml(2,i,j)

          d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

          kk = omega / (cm+um)
          a  = omega**2 * d / (cm+um)**3
          c3 = exp( -a * (xl(i,j) - x0) ) * exp( im * kk * xl(i,j) )
!         c3 = wamp(i) * exp( im * kk * xl(i,j) )

          rl(1,i,j) = rl(1,i,j) - spg(i,j) * ( vl(1,i,j) - &
                      pt5 * c3/cm**2 )
          rl(2,i,j) = rl(2,i,j) - spg(i,j) * ( vl(2,i,j) - &
                      c3 * pt5 / ( vml(1,i,j) * cm ) )
          rl(3,i,j) = rl(3,i,j) - spg(i,j) * ( vl(3,i,j) )
          rl(4,i,j) = rl(4,i,j) - spg(i,j) * ( vl(4,i,j) )
          rl(5,i,j) = rl(5,i,j) - spg(i,j) * ( vl(5,i,j) - &
                      (gamma*Ma**2 * c3 * pt5 - &
                      vml(5,i,j) * pt5 * c3 / cm**2) / vml(1,i,j) )

        else if (ispg .eq. 4) then              

          call error(code,'ispg = 4 is not working$')
          call cspg_it( rl, vl, spg, spg2 )

        end if

        end do loop_i
        end do loop_j

!.... explicit smoother

        if (eps_e .ne. zero) call smoother3D( rl, vl )

!=============================================================================!

!.... multiply by the time step

        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              rl(idof,i,j) = dtl(i,j) * rl(idof,i,j)
            end do
          end do
        end do
        
        call rkbc3D(rl)

        return
        end

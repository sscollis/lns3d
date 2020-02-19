!=============================================================================!
        subroutine rk(vl, rl, vml, xl, yl, dtl)
!  
!  Computes the RK3 RHS for the linearized compressible Navier-Stokes solver. 
!  Updated to use the hat matrices
!
!  Revised: 10-25-95
!=============================================================================!
        use stuff
        use local
        implicit none
!=============================================================================!

        real :: vl(ny*nx,ndof), rl(ny*nx,ndof), vml(ny*nx,ndof)
        real :: xl(ny*nx), yl(ny*nx), dtl(ny*nx)

!.... local variables
        
        integer :: i, j, ij, idof, jdof
        
!.... mean flow variables

        real :: cm(ny*nx), c3(ny*nx), um(ny*nx)
        real :: a, d, kk

        character*80 name
        integer :: lrec, ier, istat

        real, external :: ramp

        logical, parameter :: carp = .false.
!=============================================================================!

!.... compute first derivatives
        
        call grad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... compute second derivatives

        call grad2(ndof, nx, ny, vl, g1v, g11v, g12v, g22v, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall temperature boundary condition

        call gradbc( g1v, g2v, g11v, g12v, g22v )

!=============================================================================!
!.... form the RHS
!=============================================================================!

!.... U,xi term

        lrec = ny * nx * ndof * ndof
        call readdr(MATRIX, Q1, lrec, 1, ier)           ! \hat{A}
        call waitdr(MATRIX, istat, ier)
        call readdr(MATRIX, Q2, lrec, 2, ier)           ! \hat{B}

        rl(:,1) =         - Q1(:,1,1) * g1v(:,1)        &
                          - Q1(:,1,2) * g1v(:,2)        &
                          - Q1(:,1,3) * g1v(:,3)        &
                          - Q1(:,1,4) * g1v(:,4)        &
                          - Q1(:,1,5) * g1v(:,5)        

        rl(:,2) =         - Q1(:,2,1) * g1v(:,1)        &
                          - Q1(:,2,2) * g1v(:,2)        &
                          - Q1(:,2,3) * g1v(:,3)        &
                          - Q1(:,2,4) * g1v(:,4)        &
                          - Q1(:,2,5) * g1v(:,5)

        rl(:,3) =         - Q1(:,3,1) * g1v(:,1)        &
                          - Q1(:,3,2) * g1v(:,2)        &
                          - Q1(:,3,3) * g1v(:,3)        &
                          - Q1(:,3,4) * g1v(:,4)        &
                          - Q1(:,3,5) * g1v(:,5)

        rl(:,4) =         - Q1(:,4,1) * g1v(:,1)        &
                          - Q1(:,4,2) * g1v(:,2)        &
                          - Q1(:,4,3) * g1v(:,3)        &
                          - Q1(:,4,4) * g1v(:,4)        &
                          - Q1(:,4,5) * g1v(:,5)

        rl(:,5) =         - Q1(:,5,1) * g1v(:,1)        &
                          - Q1(:,5,2) * g1v(:,2)        &
                          - Q1(:,5,3) * g1v(:,3)        &
                          - Q1(:,5,4) * g1v(:,4)        &
                          - Q1(:,5,5) * g1v(:,5)

!.... U,eta term

        call waitdr(MATRIX, istat, ier)
        call readdr(MATRIX, Q1, lrec, 3, ier)           ! \hat{D}

        rl(:,1) = rl(:,1) - Q2(:,1,1) * g2v(:,1)        &
                          - Q2(:,1,2) * g2v(:,2)        &
                          - Q2(:,1,3) * g2v(:,3)        &
                          - Q2(:,1,4) * g2v(:,4)        &
                          - Q2(:,1,5) * g2v(:,5)        

        rl(:,2) = rl(:,2) - Q2(:,2,1) * g2v(:,1)        &
                          - Q2(:,2,2) * g2v(:,2)        &
                          - Q2(:,2,3) * g2v(:,3)        &
                          - Q2(:,2,4) * g2v(:,4)        &
                          - Q2(:,2,5) * g2v(:,5)

        rl(:,3) = rl(:,3) - Q2(:,3,1) * g2v(:,1)        &
                          - Q2(:,3,2) * g2v(:,2)        &
                          - Q2(:,3,3) * g2v(:,3)        &
                          - Q2(:,3,4) * g2v(:,4)        &
                          - Q2(:,3,5) * g2v(:,5)

        rl(:,4) = rl(:,4) - Q2(:,4,1) * g2v(:,1)        &
                          - Q2(:,4,2) * g2v(:,2)        &
                          - Q2(:,4,3) * g2v(:,3)        &
                          - Q2(:,4,4) * g2v(:,4)        &
                          - Q2(:,4,5) * g2v(:,5)

        rl(:,5) = rl(:,5) - Q2(:,5,1) * g2v(:,1)        &
                          - Q2(:,5,2) * g2v(:,2)        &
                          - Q2(:,5,3) * g2v(:,3)        &
                          - Q2(:,5,4) * g2v(:,4)        &
                          - Q2(:,5,5) * g2v(:,5)

!.... U term

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh1, lrec, 4, ier)          ! \hat{V}_{\xi\xi}

        rl(:,1) = rl(:,1) - Q1(:,1,1) * vl(:,1) &
                          - Q1(:,1,2) * vl(:,2) &
                          - Q1(:,1,3) * vl(:,3) &
                          - Q1(:,1,4) * vl(:,4) &
                          - Q1(:,1,5) * vl(:,5)

        rl(:,2) = rl(:,2) - Q1(:,2,1) * vl(:,1) &
                          - Q1(:,2,2) * vl(:,2) &
                          - Q1(:,2,3) * vl(:,3) &
                          - Q1(:,2,4) * vl(:,4) &
                          - Q1(:,2,5) * vl(:,5)

        rl(:,3) = rl(:,3) - Q1(:,3,1) * vl(:,1) &
                          - Q1(:,3,2) * vl(:,2) &
                          - Q1(:,3,3) * vl(:,3) &
                          - Q1(:,3,4) * vl(:,4) &
                          - Q1(:,3,5) * vl(:,5)

        rl(:,4) = rl(:,4) - Q1(:,4,1) * vl(:,1) &
                          - Q1(:,4,2) * vl(:,2) &
                          - Q1(:,4,3) * vl(:,3) &
                          - Q1(:,4,4) * vl(:,4) &
                          - Q1(:,4,5) * vl(:,5)

        rl(:,5) = rl(:,5) - Q1(:,5,1) * vl(:,1) &
                          - Q1(:,5,2) * vl(:,2) &
                          - Q1(:,5,3) * vl(:,3) &
                          - Q1(:,5,4) * vl(:,4) &
                          - Q1(:,5,5) * vl(:,5)

!.... U,\xi\xi term

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh2, lrec, 5, ier)          ! \hat{V}_{\xi\eta}

        rl(:,2) = rl(:,2) + Vh1(:,1) * g11v(:,2)        &
                          + Vh1(:,5) * g11v(:,3)
        rl(:,3) = rl(:,3) + Vh1(:,6) * g11v(:,2)        &
                          + Vh1(:,2) * g11v(:,3)
        rl(:,4) = rl(:,4) + Vh1(:,3) * g11v(:,4)
        rl(:,5) = rl(:,5) + Vh1(:,4) * g11v(:,5)

!.... U,\xi\eta term

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh1, lrec, 6, ier)          ! \hat{V}_{\eta\eta}

        rl(:,2) = rl(:,2) + Vh2(:,1) * g12v(:,2)        &
                          + Vh2(:,5) * g12v(:,3)
        rl(:,3) = rl(:,3) + Vh2(:,6) * g12v(:,2)        &
                          + Vh2(:,2) * g12v(:,3)
        rl(:,4) = rl(:,4) + Vh2(:,3) * g12v(:,4)
        rl(:,5) = rl(:,5) + Vh2(:,4) * g12v(:,5)

!.... U,\eta\eta term

        call waitdr(MATRIX, istat, ier)

        rl(:,2) = rl(:,2) + Vh1(:,1) * g22v(:,2)        &
                          + Vh1(:,5) * g22v(:,3)
        rl(:,3) = rl(:,3) + Vh1(:,6) * g22v(:,2)        &
                          + Vh1(:,2) * g22v(:,3)
        rl(:,4) = rl(:,4) + Vh1(:,3) * g22v(:,4)
        rl(:,5) = rl(:,5) + Vh1(:,4) * g22v(:,5)

!.... Sponge term

        if (ispg .eq. 1) then

!.... standard sponge

          rl(:,1) = rl(:,1) - spg(:) * vl(:,1)
          rl(:,2) = rl(:,2) - spg(:) * vl(:,2)
          rl(:,3) = rl(:,3) - spg(:) * vl(:,3)
          rl(:,4) = rl(:,4) - spg(:) * vl(:,4)
          rl(:,5) = rl(:,5) - spg(:) * vl(:,5)
        
        else if (ispg .eq. 2) then

          rl(:,1) = rl(:,1) - (spg(:) + spg2(:)) * vl(:,1)
          rl(:,2) = rl(:,2) - (spg(:) + spg2(:)) * vl(:,2)
          rl(:,3) = rl(:,3) - (spg(:) + spg2(:)) * vl(:,3)
          rl(:,4) = rl(:,4) - (spg(:) + spg2(:)) * vl(:,4)
          rl(:,5) = rl(:,5) - (spg(:) + spg2(:)) * vl(:,5)

        else if (ispg .eq. 3) then

!.... first do the outflow sponge

          rl(:,1) = rl(:,1) - spg(:) * vl(:,1)
          rl(:,2) = rl(:,2) - spg(:) * vl(:,2)
          rl(:,3) = rl(:,3) - spg(:) * vl(:,3)
          rl(:,4) = rl(:,4) - spg(:) * vl(:,4)
          rl(:,5) = rl(:,5) - spg(:) * vl(:,5)  

!.... subtract off the forcing wave in the inflow sponge

          cm = sqrt( vml(:,5) ) / Ma
          um = vml(:,2)

          d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

          do i = 1, nx
            do j = 1, ny
              ij = j + (i-1)*ny
              kk = omega / (cm(ij)+um(ij))
              a  = omega**2 * d / (cm(ij)+um(ij))**3
              c3(ij) = ramp( (kk*(x0-xl(ij)) + omega * time)/(two*pi) ) * &
                       cos( kk * xl(ij) - omega * time ) * &
                       exp( -a * (xl(ij) - x0) )
!             c3(ij) = ramp( (kk*(x0-xl(ij)) + omega * time)/(two*pi) ) * &
!                      wamp(i) * cos( kk * xl(ij) - omega * time )
            end do
          end do

          do ij = 1, nx*ny
            rl(ij,1) = rl(ij,1) - spg2(ij) * ( vl(ij,1) - &
                                  pt5 * c3(ij) / cm(ij)**2 )
            rl(ij,2) = rl(ij,2) - spg2(ij) * ( vl(ij,2) - &
                                  c3(ij) * pt5 / ( vml(ij,1) * cm(ij) ) )
            rl(ij,3) = rl(ij,3) - spg2(ij) * ( vl(ij,3) )
            rl(ij,4) = rl(ij,4) - spg2(ij) * ( vl(ij,4) )
            rl(ij,5) = rl(ij,5) - spg2(ij) * ( vl(ij,5) - &
                                  (gamma*Ma**2 * c3(ij) * pt5 - &
                                  vml(ij,5) * pt5 * c3(ij) / cm(ij)**2) / &
                                  vml(ij,1) )
          end do

        end if

!.... multiply by the time step

        do idof = 1, ndof
          rl(:,idof) = dtl * rl(:,idof)
        end do
        
        call rkbc(rl)

!.... output RHS statistics

        if (mod(lstep,itout).eq.0 .and. iter.eq.1) then
          call resstat(rl, dtl)
        end if
        
        return
        end

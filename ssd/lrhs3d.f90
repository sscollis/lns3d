!=============================================================================!
        subroutine lrhs3D(vl, rl, vml, xl, yl)
!  
!  Computes the RK3 RHS for the linearized three-dimensional compressible 
!  Navier-Stokes solver using Implicit time advancement.
!
!  Revised: 4-22-96
!=============================================================================!
        use stuff
        use local
        use local3d
        implicit none
!=============================================================================!
        complex :: vl(ny*nx,ndof), rl(ny*nx,ndof)
        real    :: vml(ny*nx,ndof)
        real    :: xl(ny*nx), yl(ny*nx)

!.... local variables
        
        integer :: i, j, ij, idof, jdof
        
!.... mean flow variables

        real :: cm(ny*nx), um(ny*nx)
        real :: a, d, kk
        complex :: c3(ny*nx)

        character*80 name
        integer :: lrec, ier, istat
!=============================================================================!

!.... compute first derivatives
        
        call cgrad( ndof, nx, ny, vl, c1v, c2v, dxi, deta, optx, opty, &
                    xper, yper, lsym, rsym, bsym, tsym, .false. )
        
!.... compute second derivatives

        call cgrad2( ndof, nx, ny, vl, c1v, c11v, c12v, c22v, dxi, deta, &
                     optx, opty, xper, yper, lsym, rsym, bsym, tsym, .false. )

!.... enforce adiabatic wall temperature boundary condition

        if (Navier) call cgradbc( c1v, c2v, c11v, c12v, c22v )

!=============================================================================!
!.... form the RHS
!=============================================================================!

!.... U,xi term

        lrec = ny * nx * ndof * ndof
        call readdr(MATRIX, Q1, lrec, 1, ier)           ! \hat{A}
        call waitdr(MATRIX, istat, ier)
        call readdr(MATRIX, Q2, lrec, 2, ier)           ! \hat{B}

        rl(:,1) =         + Q1(:,1,1) * c1v(:,1)        &
                          + Q1(:,1,2) * c1v(:,2)        &
                          + Q1(:,1,3) * c1v(:,3)        &
                          + Q1(:,1,4) * c1v(:,4)        &
                          + Q1(:,1,5) * c1v(:,5)        

        rl(:,2) =         + Q1(:,2,1) * c1v(:,1)        &
                          + Q1(:,2,2) * c1v(:,2)        &
                          + Q1(:,2,3) * c1v(:,3)        &
                          + Q1(:,2,4) * c1v(:,4)        &
                          + Q1(:,2,5) * c1v(:,5)

        rl(:,3) =         + Q1(:,3,1) * c1v(:,1)        &
                          + Q1(:,3,2) * c1v(:,2)        &
                          + Q1(:,3,3) * c1v(:,3)        &
                          + Q1(:,3,4) * c1v(:,4)        &
                          + Q1(:,3,5) * c1v(:,5)

        rl(:,4) =         + Q1(:,4,1) * c1v(:,1)        &
                          + Q1(:,4,2) * c1v(:,2)        &
                          + Q1(:,4,3) * c1v(:,3)        &
                          + Q1(:,4,4) * c1v(:,4)        &
                          + Q1(:,4,5) * c1v(:,5)

        rl(:,5) =         + Q1(:,5,1) * c1v(:,1)        &
                          + Q1(:,5,2) * c1v(:,2)        &
                          + Q1(:,5,3) * c1v(:,3)        &
                          + Q1(:,5,4) * c1v(:,4)        &
                          + Q1(:,5,5) * c1v(:,5)

!.... U,eta term

        call waitdr(MATRIX, istat, ier)
        call readdr(MATRIX, Q1, lrec, 3, ier)           ! \hat{D}

        rl(:,1) = rl(:,1) + Q2(:,1,1) * c2v(:,1)        &
                          + Q2(:,1,2) * c2v(:,2)        &
                          + Q2(:,1,3) * c2v(:,3)        &
                          + Q2(:,1,4) * c2v(:,4)        &
                          + Q2(:,1,5) * c2v(:,5)        

        rl(:,2) = rl(:,2) + Q2(:,2,1) * c2v(:,1)        &
                          + Q2(:,2,2) * c2v(:,2)        &
                          + Q2(:,2,3) * c2v(:,3)        &
                          + Q2(:,2,4) * c2v(:,4)        &
                          + Q2(:,2,5) * c2v(:,5)

        rl(:,3) = rl(:,3) + Q2(:,3,1) * c2v(:,1)        &
                          + Q2(:,3,2) * c2v(:,2)        &
                          + Q2(:,3,3) * c2v(:,3)        &
                          + Q2(:,3,4) * c2v(:,4)        &
                          + Q2(:,3,5) * c2v(:,5)

        rl(:,4) = rl(:,4) + Q2(:,4,1) * c2v(:,1)        &
                          + Q2(:,4,2) * c2v(:,2)        &
                          + Q2(:,4,3) * c2v(:,3)        &
                          + Q2(:,4,4) * c2v(:,4)        &
                          + Q2(:,4,5) * c2v(:,5)

        rl(:,5) = rl(:,5) + Q2(:,5,1) * c2v(:,1)        &
                          + Q2(:,5,2) * c2v(:,2)        &
                          + Q2(:,5,3) * c2v(:,3)        &
                          + Q2(:,5,4) * c2v(:,4)        &
                          + Q2(:,5,5) * c2v(:,5)

!.... U term

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh1, lrec, 4, ier)          ! \hat{V}_{\xi\xi}

        rl(:,1) = rl(:,1) + Q1(:,1,1) * vl(:,1) &
                          + Q1(:,1,2) * vl(:,2) &
                          + Q1(:,1,3) * vl(:,3) &
                          + Q1(:,1,4) * vl(:,4) &
                          + Q1(:,1,5) * vl(:,5)

        rl(:,2) = rl(:,2) + Q1(:,2,1) * vl(:,1) &
                          + Q1(:,2,2) * vl(:,2) &
                          + Q1(:,2,3) * vl(:,3) &
                          + Q1(:,2,4) * vl(:,4) &
                          + Q1(:,2,5) * vl(:,5)

        rl(:,3) = rl(:,3) + Q1(:,3,1) * vl(:,1) &
                          + Q1(:,3,2) * vl(:,2) &
                          + Q1(:,3,3) * vl(:,3) &
                          + Q1(:,3,4) * vl(:,4) &
                          + Q1(:,3,5) * vl(:,5)

        rl(:,4) = rl(:,4) + Q1(:,4,1) * vl(:,1) &
                          + Q1(:,4,2) * vl(:,2) &
                          + Q1(:,4,3) * vl(:,3) &
                          + Q1(:,4,4) * vl(:,4) &
                          + Q1(:,4,5) * vl(:,5)

        rl(:,5) = rl(:,5) + Q1(:,5,1) * vl(:,1) &
                          + Q1(:,5,2) * vl(:,2) &
                          + Q1(:,5,3) * vl(:,3) &
                          + Q1(:,5,4) * vl(:,4) &
                          + Q1(:,5,5) * vl(:,5)

!.... U,\xi\xi term

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh2, lrec, 5, ier)          ! \hat{V}_{\xi\eta}

        rl(:,2) = rl(:,2) - Vh1(:,1) * c11v(:,2)        &
                          - Vh1(:,5) * c11v(:,3)
        rl(:,3) = rl(:,3) - Vh1(:,6) * c11v(:,2)        &
                          - Vh1(:,2) * c11v(:,3)
        rl(:,4) = rl(:,4) - Vh1(:,3) * c11v(:,4)
        rl(:,5) = rl(:,5) - Vh1(:,4) * c11v(:,5)

!.... U,\xi\eta term

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh1, lrec, 6, ier)          ! \hat{V}_{\eta\eta}

        rl(:,2) = rl(:,2) - Vh2(:,1) * c12v(:,2)        &
                          - Vh2(:,5) * c12v(:,3)
        rl(:,3) = rl(:,3) - Vh2(:,6) * c12v(:,2)        &
                          - Vh2(:,2) * c12v(:,3)
        rl(:,4) = rl(:,4) - Vh2(:,3) * c12v(:,4)
        rl(:,5) = rl(:,5) - Vh2(:,4) * c12v(:,5)

!.... U,\eta\eta term

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh2, lrec, 7, ier)          ! \hat{A}_i, \hat{B}_i

        rl(:,2) = rl(:,2) - Vh1(:,1) * c22v(:,2)        &
                          - Vh1(:,5) * c22v(:,3)
        rl(:,3) = rl(:,3) - Vh1(:,6) * c22v(:,2)        &
                          - Vh1(:,2) * c22v(:,3)
        rl(:,4) = rl(:,4) - Vh1(:,3) * c22v(:,4)
        rl(:,5) = rl(:,5) - Vh1(:,4) * c22v(:,5)

!.... 3D terms

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * ndof * ndof
        call readdr(MATRIX, Q1, lrec, 8, ier)           ! \hat{D}_i

!.... \hat{A}_i, \hat{B}_i terms

        rl(:,2) = rl(:,2) - im * ( Vh2(:,1) * c1v(:,4) + Vh2(:,3) * c2v(:,4) )
        rl(:,3) = rl(:,3) - im * ( Vh2(:,2) * c1v(:,4) + Vh2(:,4) * c2v(:,4) )
        rl(:,4) = rl(:,4) - im * ( Vh2(:,1) * c1v(:,2) + Vh2(:,2) * c1v(:,3) &
                                 + Vh2(:,3) * c2v(:,2) + Vh2(:,4) * c2v(:,3) )
        
!.... \hat{D}_i term

        call waitdr(MATRIX, istat, ier)

        rl(:,1) = rl(:,1) + im * ( Q1(:,1,1) * vl(:,1)  &
                                 + Q1(:,1,2) * vl(:,2)  &
                                 + Q1(:,1,3) * vl(:,3)  &
                                 + Q1(:,1,4) * vl(:,4)  &
                                 + Q1(:,1,5) * vl(:,5) )

        rl(:,2) = rl(:,2) + im * ( Q1(:,2,1) * vl(:,1)  &
                                 + Q1(:,2,2) * vl(:,2)  &
                                 + Q1(:,2,3) * vl(:,3)  &
                                 + Q1(:,2,4) * vl(:,4)  &
                                 + Q1(:,2,5) * vl(:,5) )

        rl(:,3) = rl(:,3) + im * ( Q1(:,3,1) * vl(:,1)  &
                                 + Q1(:,3,2) * vl(:,2)  &
                                 + Q1(:,3,3) * vl(:,3)  &
                                 + Q1(:,3,4) * vl(:,4)  &
                                 + Q1(:,3,5) * vl(:,5) )

        rl(:,4) = rl(:,4) + im * ( Q1(:,4,1) * vl(:,1)  &
                                 + Q1(:,4,2) * vl(:,2)  &
                                 + Q1(:,4,3) * vl(:,3)  &
                                 + Q1(:,4,4) * vl(:,4)  &
                                 + Q1(:,4,5) * vl(:,5) )

        rl(:,5) = rl(:,5) + im * ( Q1(:,5,1) * vl(:,1)  &
                                 + Q1(:,5,2) * vl(:,2)  &
                                 + Q1(:,5,3) * vl(:,3)  &
                                 + Q1(:,5,4) * vl(:,4)  &
                                 + Q1(:,5,5) * vl(:,5) )
!.... Sponge term

        if (ispg .eq. 1) then           !.... standard sponge

          rl(:,1) = rl(:,1) + spg(:) * vl(:,1)
          rl(:,2) = rl(:,2) + spg(:) * vl(:,2)
          rl(:,3) = rl(:,3) + spg(:) * vl(:,3)
          rl(:,4) = rl(:,4) + spg(:) * vl(:,4)
          rl(:,5) = rl(:,5) + spg(:) * vl(:,5)
        
        else if (ispg .eq. 2) then              
        
          rl(:,1) = rl(:,1) + (spg(:) + spg2(:)) * vl(:,1)
          rl(:,2) = rl(:,2) + (spg(:) + spg2(:)) * vl(:,2)
          rl(:,3) = rl(:,3) + (spg(:) + spg2(:)) * vl(:,3)
          rl(:,4) = rl(:,4) + (spg(:) + spg2(:)) * vl(:,4)
          rl(:,5) = rl(:,5) + (spg(:) + spg2(:)) * vl(:,5)

        else if (ispg .eq. 3) then

!.... the outflow sponge

          rl(:,1) = rl(:,1) + spg(:) * vl(:,1)
          rl(:,2) = rl(:,2) + spg(:) * vl(:,2)
          rl(:,3) = rl(:,3) + spg(:) * vl(:,3)
          rl(:,4) = rl(:,4) + spg(:) * vl(:,4)
          rl(:,5) = rl(:,5) + spg(:) * vl(:,5)  

!.... the inflow sponge

          cm = sqrt( vml(:,5) ) / Ma
          um = vml(:,2)

          d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

          do i = 1, nx
            do j = 1, ny
              ij = j + (i-1)*ny
              kk = omega / (cm(ij)+um(ij))
              a  = omega**2 * d / (cm(ij)+um(ij))**3
              c3(ij) = exp( -a * (xl(ij) - x0) ) * exp( im * kk * xl(ij) )
!             c3(ij) = wamp(i) * exp( im * kk * xl(ij) )
            end do
          end do

          do ij = 1, nx*ny
            rl(ij,1) = rl(ij,1) + spg2(ij) * ( vl(ij,1) - &
                       pt5 * c3(ij) / cm(ij)**2 )
            rl(ij,2) = rl(ij,2) + spg2(ij) * ( vl(ij,2) - &
                       c3(ij) * pt5 / ( vml(ij,1) * cm(ij) ) )
            rl(ij,3) = rl(ij,3) + spg2(ij) * ( vl(ij,3) )
            rl(ij,4) = rl(ij,4) + spg2(ij) * ( vl(ij,4) )
            rl(ij,5) = rl(ij,5) + spg2(ij) * ( vl(ij,5) - &
                       (gamma*Ma**2 * c3(ij) * pt5 - &
                       vml(ij,5) * pt5 * c3(ij) / cm(ij)**2) / vml(ij,1) )
          end do

        else if (ispg .eq. 4) then              

          call cspg_it( rl, vl, spg, spg2 )

        end if

!.... explicit smoother

        if (eps_e .ne. zero) call smoother3D( rl, vl )

        return
        end

!=============================================================================!
        module cic

        integer :: ic_start=0
        complex, allocatable :: vic(:,:,:)
        
        end module cic

!=============================================================================!
        subroutine cspg_it( rl, vl, spgl, spg2l )
!
!       This routine reads in the initial condition when first called, and
!       sponges to the initial condition.
!
!=============================================================================!
        use cic
        use stuff
        use local
        implicit none
        
        complex :: rl(ny,nx,ndof), vl(ny,nx,ndof)
        real    :: spgl(ny,nx), spg2l(ny,nx)
        integer :: i

        complex, parameter :: ac=(2.2804739410500E-001,-6.5163146761218E-003)
!       complex, parameter :: ac=(-2.8831962908130E-001,-1.3854663671636E-002)

        real :: rtmp
        integer :: itmp
!=============================================================================!
        if (ic_start .eq. 0) then
          allocate( vic(ny,nx,ndof) )
          ic_start = 1

!         open(10,file='output.R.0',form='unformatted',status='old')
!         read(10) itmp, rtmp, itmp, itmp, itmp, itmp, &
!                  rtmp, rtmp, rtmp, rtmp, rtmp
!         read(10) vic
!         close(10)
          
          do i = 1, nx
            vic(:,i,1) = cmplx(rhor(:), rhoi(:)) * exp(im * ac * x(:,i))
            vic(:,i,2) = cmplx(  ur(:),   ui(:)) * exp(im * ac * x(:,i))
            vic(:,i,3) = cmplx(  vr(:),   vi(:)) * exp(im * ac * x(:,i))
            vic(:,i,4) = cmplx(  wr(:),   wi(:)) * exp(im * ac * x(:,i))
            vic(:,i,5) = cmplx(  tr(:),   ti(:)) * exp(im * ac * x(:,i))
          end do
        end if

        do i = 1, nx
          rl(:,i,1) = rl(:,i,1) + (spgl(:,i) + spg2l(:,i)) * &
                      ( vl(:,i,1) - vic(:,i,1) )
          rl(:,i,2) = rl(:,i,2) + (spgl(:,i) + spg2l(:,i)) * &
                      ( vl(:,i,2) - vic(:,i,2) )
          rl(:,i,3) = rl(:,i,3) + (spgl(:,i) + spg2l(:,i)) * &
                      ( vl(:,i,3) - vic(:,i,3) )
          rl(:,i,4) = rl(:,i,4) + (spgl(:,i) + spg2l(:,i)) * &
                      ( vl(:,i,4) - vic(:,i,4) )
          rl(:,i,5) = rl(:,i,5) + (spgl(:,i) + spg2l(:,i)) * &
                      ( vl(:,i,5) - vic(:,i,5) )
        end do
        
        return
        end

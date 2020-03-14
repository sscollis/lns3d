!=============================================================================!
        subroutine genmtrx(vl, kzl, comp) 
!  
!       Generates the linearized coefficient matrices for the
!       compressible NS equations 
!
! Input:
!   vl:    The field for which the matrices are computed
!   kzl:   The local spanwise wave number
!   comp:  Make the complex matrices
!
! Output:
!   matrices:  Stored on the SSD
!
!       Revised: 5-10-95
!
!       1) added compact storage for G and Vij
!       2) improved G.Q product
!       3) now compute gradients of p, mu, con, lm using chain rule
!          instead of calling grad individually
!       4) Jump around C matrix if 2D
!
!       Revised: 5-23-95
!
!       1) removed mean flow designation to make for general
!       2) now do allocation external to the routine
!       3) merged with conmtrx
!
!       Revised: 9-9-97
!
!       1) removed SSD for use on non-Cray computers
!
!=============================================================================!
        use global
        use local
        use pot
        use material
        implicit none
        
        real    :: vl(ny*nx,ndof), kzl
        logical :: comp
!=============================================================================!
        integer ier, i, j, idof, ji
        character*80 name
        real :: fact
        logical :: carp

        integer :: lrec, istat

        real :: Q1(ny*nx,ndof,ndof), Q2(ny*nx,ndof,ndof), Vij(ny*nx,4)
        
        real :: gpr(ny*nx,2)
!=============================================================================!
        carp = .false.
        if (impl.eq.0) carp = .true.

!.... Compute first derivatives of field in the mapped space

        call grad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... Compute second derivatives of field
        
        call grad2(ndof, nx, ny, vl, g1v, g11v, g12v, g22v, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall BC on the field

        if (Navier) call gradbc( g1v, g2v, g11v, g12v, g22v )
        
!.... damp the second derivatives in the viscous terms and 
!.... compute the pressure gradient on the outflow plane

        if (idamp.eq.1) then
          do idof = 1, ndof
            g11v(:,idof) = damp(:) * g11v(:,idof)
          end do
          if (linear.eq.0) call calcp( vl, g2v, m1, m2, n1, n2, gpr )
        end if
        
!.... transform the gradients to physical space

        do idof = 1, ndof

          g1vl  = g1v(:,idof)*m1(:) + g2v(:,idof)*n1(:)
          g2vl  = g1v(:,idof)*m2(:) + g2v(:,idof)*n2(:)

          g11vl = g11v(:,idof)       * m1m1(:)  + &
                  two * g12v(:,idof) * m1n1(:)  + &
                  g22v(:,idof)       * n1n1(:)  + &
                  g1v(:,idof)        * m11(:)   + &
                  g2v(:,idof)        * n11(:)

          g12vl = g11v(:,idof)       * m1m2(:)  + &
                  g12v(:,idof)       * m1n2(:)  + &
                  g12v(:,idof)       * m2n1(:)  + &
                  g22v(:,idof)       * n1n2(:)  + &
                  g1v(:,idof)        * m12(:)   + &
                  g2v(:,idof)        * n12(:)

          g22vl = g11v(:,idof)       * m2m2(:)  + &
                  two * g12v(:,idof) * m2n2(:)  + &
                  g22v(:,idof)       * n2n2(:)  + &
                  g1v(:,idof)        * m22(:)   + &
                  g2v(:,idof)        * n22(:)

          g1v(:,idof)  = g1vl
          g2v(:,idof)  = g2vl
          g11v(:,idof) = g11vl
          g12v(:,idof) = g12vl
          g22v(:,idof) = g22vl

        end do
        
!=============================================================================!
!.... setup local variables
!=============================================================================!

        rho     = vl(:,1)
        u1      = vl(:,2)
        u2      = vl(:,3)
        u3      = vl(:,4)
        t       = vl(:,5)
        rhoinv  = one / rho

!.... compute pressure

        fact = one / (gamma * Ma**2)
        p    = fact * rho * t

!.... initialize gradients of mean field (uniform in span)

        grho(:,1) = g1v(:,1)
        grho(:,2) = g2v(:,1)
        grho(:,3) = zero

        gu(:,1,1) = g1v(:,2)
        gu(:,1,2) = g2v(:,2)
        gu(:,1,3) = zero
        
        gu(:,2,1) = g1v(:,3)
        gu(:,2,2) = g2v(:,3)
        gu(:,2,3) = zero

        gu(:,3,1) = g1v(:,4)
        gu(:,3,2) = g2v(:,4)
        gu(:,3,3) = zero

        gt(:,1)   = g1v(:,5)
        gt(:,2)   = g2v(:,5)
        gt(:,3)   = zero

!.... compute divergence
  
        divu = gu(:,1,1) + gu(:,2,2) + gu(:,3,3)
            
!.... compute the gradient of the divergence of um

        g1divu = g11v(:,2) + g12v(:,3)
        g2divu = g12v(:,2) + g22v(:,3)
        g3divu = zero

!.... compute strain rate tensor derivative for the viscous terms

        S1jj = pt5 * ( g11v(:,2) + g11v(:,2) + g22v(:,2) + &
                       g12v(:,3) )
                
        S2jj = pt5 * ( g11v(:,3) + g12v(:,2) + g22v(:,3) + &
                       g22v(:,3) )
        
        S3jj = pt5 * ( g11v(:,4) + g22v(:,4) )

!.... compute Laplacian of Temperature for the heat conduction term

        Lapt = g11v(:,5) + g22v(:,5)
        
!.... compute gradient of mean pressure using chain-rule

        fact = one / (gamma * Ma**2)
        gp(:,1) = fact * ( grho(:,1) * t + rho * gt(:,1) )
        gp(:,2) = fact * ( grho(:,2) * t + rho * gt(:,2) )
        gp(:,3) = fact * ( grho(:,3) * t + rho * gt(:,3) )

!.... correct the pressure near the outflow plane

        if (idamp.eq.1 .and. linear.eq.0) then
          gp(:,1) = damp * gp(:,1) + (one-damp) * gpr(:,1)
          gp(:,2) = damp * gp(:,2) + (one-damp) * gpr(:,2)
        end if
        
!.... compute strain rate tensor for the viscous terms

        do i = 1, nsd
          do j = 1, nsd
            S(:,i,j) = pt5 * ( gu(:,i,j) + gu(:,j,i) )
          end do
        end do

!.... compute the material properties

        call getmat(t(:),    mu(:),   lm(:),    con(:),         &
                             dmu(:),  d2mu(:),  dlm(:),         &
                             d2lm(:), dcon(:),  d2con(:)        )

!.... compute gradients of viscosity using chain-rule

        g1mu = dmu * gt(:,1)
        g2mu = dmu * gt(:,2)
        g3mu = dmu * gt(:,3)

        g1dmu = d2mu * gt(:,1)
        g2dmu = d2mu * gt(:,2)
        g3dmu = d2mu * gt(:,3)

!.... compute gradients of conductivity using chain-rule

        g1con = dcon * gt(:,1)
        g2con = dcon * gt(:,2)
        g3con = dcon * gt(:,3)

        g1dcon = d2con * gt(:,1)
        g2dcon = d2con * gt(:,2)
        g3dcon = d2con * gt(:,3)

!.... compute gradients of bulk viscosity using chain-rule

        g1lm = dlm * gt(:,1)
        g2lm = dlm * gt(:,2)
        g3lm = dlm * gt(:,3)

        g1dlm = d2lm * gt(:,1)
        g2dlm = d2lm * gt(:,2)
        g3dlm = d2lm * gt(:,3)
        
        if (linear.eq.1 .and. wall.eq.3) &
          call genbump( g1v, g2v, g11v, g12v, g22v )
!=============================================================================!
!                       N O N L I N E A R   R H S
!=============================================================================!

!.... Everything you need to compute the nonlinear residual for vl is in
!.... the local module right now

!=============================================================================!
!               L I N E A R I Z E D   M A T R I C E S
!=============================================================================!
!
!..... Matrices needed for linear Navier-Stokes solver and for implicit
!      time advancement
!
!=============================================================================!
!.... compute the A matrix
!=============================================================================!

        call genA(Q1, vl)

!=============================================================================!
!.... compute the B matrix
!=============================================================================!

        call genB(Q2, vl)

!=============================================================================!
!.... compute the Vij matrix in compact form
!=============================================================================!

        call genVij(Vij)

!=============================================================================!
!.... Compute \hat{A}
!=============================================================================!

!.... A, B, Vxx, Vyy, Vxy terms

        Ah(:,1,1) = Q1(:,1,1) * m1(:) + Q2(:,1,1) * m2(:)
        Ah(:,1,2) = Q1(:,1,2) * m1(:)
        Ah(:,1,3) =                     Q2(:,1,3) * m2(:)
        Ah(:,1,4) = zero
        Ah(:,1,5) = zero

        Ah(:,2,1) = Q1(:,2,1) * m1(:)
        Ah(:,2,2) = Q1(:,2,2) * m1(:) + Q2(:,2,2) * m2(:) - &
                    Vij(:,1) * m11(:) - Vij(:,2) * m22(:)
        Ah(:,2,3) = Q1(:,2,3) * m1(:) + Q2(:,2,3) * m2(:) - &
                    Vij(:,4) * m12(:)
        Ah(:,2,4) = Q1(:,2,4) * m1(:)
        Ah(:,2,5) = Q1(:,2,5) * m1(:) + Q2(:,2,5) * m2(:)

        Ah(:,3,1) =                     Q2(:,3,1) * m2(:)
        Ah(:,3,2) = Q1(:,3,2) * m1(:) + Q2(:,3,2) * m2(:) - &
                    Vij(:,4) * m12(:)
        Ah(:,3,3) = Q1(:,3,3) * m1(:) + Q2(:,3,3) * m2(:) - &
                    Vij(:,2) * m11(:) - Vij(:,1) * m22(:)
        Ah(:,3,4) =                     Q2(:,3,4) * m2(:)
        Ah(:,3,5) = Q1(:,3,5) * m1(:) + Q2(:,3,5) * m2(:)

        Ah(:,4,1) = zero
        Ah(:,4,2) = Q1(:,4,2) * m1(:)
        Ah(:,4,3) =                     Q2(:,4,3) * m2(:)
        Ah(:,4,4) = Q1(:,4,4) * m1(:) + Q2(:,4,4) * m2(:) - &
                    Vij(:,2) * m11(:) - Vij(:,2) * m22(:)
        Ah(:,4,5) = Q1(:,4,5) * m1(:) + Q2(:,4,5) * m2(:)

        Ah(:,5,1) = Q1(:,5,1) * m1(:) + Q2(:,5,1) * m2(:)
        Ah(:,5,2) = Q1(:,5,2) * m1(:) + Q2(:,5,2) * m2(:)
        Ah(:,5,3) = Q1(:,5,3) * m1(:) + Q2(:,5,3) * m2(:)
        Ah(:,5,4) = Q1(:,5,4) * m1(:) + Q2(:,5,4) * m2(:)
        Ah(:,5,5) = Q1(:,5,5) * m1(:) + Q2(:,5,5) * m2(:) - &
                    Vij(:,3) * m11(:) - Vij(:,3) * m22(:)

!.... correct the matrix for the parabolized outflow

        if (idamp.eq.1 .and. linear.eq.0) then
          Ah(:,2,1) = damp * Ah(:,2,1)
          Ah(:,2,5) = Ah(:,2,5) - (one-damp) * m1 / (gamma * Ma**2)

          Ah(:,3,1) = damp * Ah(:,3,1)
          Ah(:,3,5) = Ah(:,3,5) - (one-damp) * m2 / (gamma * Ma**2)
        end if

!=============================================================================!
!.... Compute \hat{B}
!=============================================================================!

!.... A, B, Vxx, Vyy, Vxy term

        Bh(:,1,1) = Q1(:,1,1) * n1(:) + Q2(:,1,1) * n2(:)
        Bh(:,1,2) = Q1(:,1,2) * n1(:)
        Bh(:,1,3) =                     Q2(:,1,3) * n2(:)
        Bh(:,1,4) = zero
        Bh(:,1,5) = zero

        Bh(:,2,1) = Q1(:,2,1) * n1(:)
        Bh(:,2,2) = Q1(:,2,2) * n1(:) + Q2(:,2,2) * n2(:) - &
                    Vij(:,1) * n11(:) - Vij(:,2) * n22(:)
        Bh(:,2,3) = Q1(:,2,3) * n1(:) + Q2(:,2,3) * n2(:) - &
                    Vij(:,4) * n12(:)
        Bh(:,2,4) = Q1(:,2,4) * n1(:)
        Bh(:,2,5) = Q1(:,2,5) * n1(:) + Q2(:,2,5) * n2(:)

        Bh(:,3,1) =                     Q2(:,3,1) * n2(:)
        Bh(:,3,2) = Q1(:,3,2) * n1(:) + Q2(:,3,2) * n2(:) - &
                    Vij(:,4) * n12(:)
        Bh(:,3,3) = Q1(:,3,3) * n1(:) + Q2(:,3,3) * n2(:) - &
                    Vij(:,2) * n11(:) - Vij(:,1) * n22(:)
        Bh(:,3,4) =                     Q2(:,3,4) * n2(:)
        Bh(:,3,5) = Q1(:,3,5) * n1(:) + Q2(:,3,5) * n2(:)

        Bh(:,4,1) = zero
        Bh(:,4,2) = Q1(:,4,2) * n1(:)
        Bh(:,4,3) =                     Q2(:,4,3) * n2(:)
        Bh(:,4,4) = Q1(:,4,4) * n1(:) + Q2(:,4,4) * n2(:) - &
                    Vij(:,2) * n11(:) - Vij(:,2) * n22(:)
        Bh(:,4,5) = Q1(:,4,5) * n1(:) + Q2(:,4,5) * n2(:)

        Bh(:,5,1) = Q1(:,5,1) * n1(:) + Q2(:,5,1) * n2(:)
        Bh(:,5,2) = Q1(:,5,2) * n1(:) + Q2(:,5,2) * n2(:)
        Bh(:,5,3) = Q1(:,5,3) * n1(:) + Q2(:,5,3) * n2(:)
        Bh(:,5,4) = Q1(:,5,4) * n1(:) + Q2(:,5,4) * n2(:)
        Bh(:,5,5) = Q1(:,5,5) * n1(:) + Q2(:,5,5) * n2(:) - &
                    Vij(:,3) * n11(:) - Vij(:,3) * n22(:)

!.... correct for the adiabatic condition at the wall:  dT/dn = 0

        if ( (.not. yper) .and. wallt.eq.2) then
          j = 1
          do i = 1, nx
            ji = j + (i-1)*ny
            Bh(ji,1,5) = zero
            Bh(ji,2,5) = zero
            Bh(ji,3,5) = zero
            Bh(ji,4,5) = zero
            Bh(ji,5,5) = zero
          end do
        end if

!=============================================================================!
!.... compute the D matrix and \hat{D}
!=============================================================================!
        call genD(Dh, vl)

!.... Vzz term

        if (kzl .ne. zero) then
          Dh(:,2,2) = Dh(:,2,2) + kzl**2 * Vij(:,2)
          Dh(:,3,3) = Dh(:,3,3) + kzl**2 * Vij(:,2)
          Dh(:,4,4) = Dh(:,4,4) + kzl**2 * Vij(:,1)
          Dh(:,5,5) = Dh(:,5,5) + kzl**2 * Vij(:,3)
        end if

!=============================================================================!
!.... Compute \hat{V}_{\xi\xi} in compact form
!
!       Vh1(1) = V(2,2)
!       Vh1(2) = V(3,3)
!       Vh1(3) = V(4,4)
!       Vh1(4) = V(5,5)
!       Vh1(5) = V(2,3)
!       Vh1(6) = V(3,2)
!
!       Note:  Vh1(5) ussually equals Vh1(6) except on a boundary where
!              viscous BC's are applied.
!
!=============================================================================!

!.... Vxx term and Vyy term

        Vh11(:,1) = Vij(:,1) * m1m1(:) + Vij(:,2) * m2m2(:)
        Vh11(:,2) = Vij(:,2) * m1m1(:) + Vij(:,1) * m2m2(:)
        Vh11(:,3) = Vij(:,2) * m1m1(:) + Vij(:,2) * m2m2(:)
        Vh11(:,4) = Vij(:,3) * m1m1(:) + Vij(:,3) * m2m2(:)

!.... Vxy term

        Vh11(:,5) = Vij(:,4) * m1m2(:)
        Vh11(:,6) = Vij(:,4) * m1m2(:)

!.... damp the second derivatives in the viscous terms

        if (idamp.eq.1) then
          do i = 1, 6
            Vh11(:,i) = damp(:) * Vh11(:,i)
          end do
        end if

!=============================================================================!
!.... Compute \hat{V}_{\xi\eta}
!=============================================================================!

!.... Vxx term and Vyy

        Vh12(:,1) = two * ( Vij(:,1) * m1n1(:) + Vij(:,2) * m2n2(:) )
        Vh12(:,2) = two * ( Vij(:,2) * m1n1(:) + Vij(:,1) * m2n2(:) )
        Vh12(:,3) = two * ( Vij(:,2) * m1n1(:) + Vij(:,2) * m2n2(:) )
        Vh12(:,4) = two * ( Vij(:,3) * m1n1(:) + Vij(:,3) * m2n2(:) )

!.... Vxy term

        Vh12(:,5) = Vij(:,4) * ( m1n2(:) + m2n1(:) )
        Vh12(:,6) = Vij(:,4) * ( m1n2(:) + m2n1(:) )

!=============================================================================!
!.... Compute \hat{V}_{\eta\eta}
!=============================================================================!

!.... Vxx term and Vyy term

        Vh22(:,1) = Vij(:,1) * n1n1(:) + Vij(:,2) * n2n2(:)
        Vh22(:,2) = Vij(:,2) * n1n1(:) + Vij(:,1) * n2n2(:)
        Vh22(:,3) = Vij(:,2) * n1n1(:) + Vij(:,2) * n2n2(:)
        Vh22(:,4) = Vij(:,3) * n1n1(:) + Vij(:,3) * n2n2(:)

!.... Vxy term

        Vh22(:,5) = Vij(:,4) * n1n2(:)
        Vh22(:,6) = Vij(:,4) * n1n2(:)

!=============================================================================!
!.... Compute 3D terms
!=============================================================================!

!       if (kzl .ne. zero .or. omega .ne. zero) then
        if (comp) then

!.... \hat{A}_i and \hat{B}_i

          ABhi(:,1) = kzl * Vij(:,4) * m1(:) 
          ABhi(:,2) = kzl * Vij(:,4) * m2(:)
          ABhi(:,3) = kzl * Vij(:,4) * n1(:)
          ABhi(:,4) = kzl * Vij(:,4) * n2(:)

!.... \hat{D}_i

          call genC(Dhi, vl)

          Dhi = kzl * Dhi

!.... put in unsteady term

          Dhi(:,1,1) = Dhi(:,1,1) - omega
          Dhi(:,2,2) = Dhi(:,2,2) - omega
          Dhi(:,3,3) = Dhi(:,3,3) - omega
          Dhi(:,4,4) = Dhi(:,4,4) - omega
          Dhi(:,5,5) = Dhi(:,5,5) - omega

        end if  

        return
        end

!=============================================================================!
        subroutine calcp( v, g2v, m1l, m2l, n1l, n2l, gpr )
!
!=============================================================================!
        use global
        use pot
        implicit none
        
        real :: v(ny,nx,ndof), g2v(ny,nx,ndof)
        real :: m1l(ny,nx), m2l(ny,nx), n1l(ny,nx), n2l(ny,nx)
        real :: gpr(ny,nx,2)
!=============================================================================!
        gpr(:,:,1) = m1l(:,:) * g1p(:,:) + n1l(:,:) / (gamma * Ma**2) * &
                     ( g2v(:,:,1) * v(:,:,5) + v(:,:,1) * g2v(:,:,5) )

        gpr(:,:,2) = m2l(:,:) * g1p(:,:) + n2l(:,:) / (gamma * Ma**2) * &
                     ( g2v(:,:,1) * v(:,:,5) + v(:,:,1) * g2v(:,:,5) )

        return
        end

!=============================================================================!
        subroutine nrk(vl, rl, dtl) 
!  
!  Computes the RK RHS for the nonlinear compressible Navier-Stokes solver. 
!
!  Written: 9-14-95
!
!  Revised: 6-11-96     Cleaned up extraneous calculations
!=============================================================================!
        use stuff
        use local
        use material
        implicit none
        
        real :: vl(ny*nx,ndof), rl(ny*nx,ndof), dtl(ny*nx)
!=============================================================================!
        integer :: ier, i, j, idof, ji
        real    :: fact

        real :: gpr(ny*nx,2)
!=============================================================================!

!.... Compute first derivatives of field in the mapped space

        call grad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, .true.)
        
!.... Compute second derivatives of field
        
        call grad2(ndof, nx, ny, vl, g1v, g11v, g12v, g22v, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, .true.)

!.... enforce adiabatic wall BC on the field

        call gradbc( g1v, g2v, g11v, g12v, g22v )
        
!.... damp the second derivatives in the viscous terms and 
!.... compute the pressure gradient on the outflow plane

        if (idamp.eq.1) then
          do idof = 1, ndof
            g11v(:,idof) = damp(:) * g11v(:,idof)
          end do
          call calcp( vl, g2v, m1, m2, n1, n2, gpr )
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

!.... initialize gradients of mean field

        grho(:,1) = g1v(:,1)
        grho(:,2) = g2v(:,1)
!       grho(:,3) = zero

        gu(:,1,1) = g1v(:,2)
        gu(:,1,2) = g2v(:,2)
!       gu(:,1,3) = zero
        
        gu(:,2,1) = g1v(:,3)
        gu(:,2,2) = g2v(:,3)
!       gu(:,2,3) = zero

        gu(:,3,1) = g1v(:,4)
        gu(:,3,2) = g2v(:,4)
!       gu(:,3,3) = zero

        gt(:,1)   = g1v(:,5)
        gt(:,2)   = g2v(:,5)
!       gt(:,3)   = zero

!.... compute the gradient of the divergence of um

        g1divu = g11v(:,2) + g12v(:,3)
        g2divu = g12v(:,2) + g22v(:,3)
!       g3divu = zero

!.... compute strain rate tensor derivative for the viscous terms

        S1jj = pt5 * ( g11v(:,2) + g11v(:,2) + g22v(:,2) + &
                       g12v(:,3) )
                
        S2jj = pt5 * ( g11v(:,3) + g12v(:,2) + g22v(:,3) + &
                       g22v(:,3) )
        
        S3jj = pt5 * ( g11v(:,4) + g22v(:,4) )

!.... compute Laplacian of Temperature for the heat conduction term

        Lapt = g11v(:,5) + g22v(:,5)
        
!.... compute divergence
  
        divu = gu(:,1,1) + gu(:,2,2) ! + gu(:,3,3)
            
!.... compute gradient of mean pressure using chain-rule

        fact = one / (gamma * Ma**2)
        gp(:,1) = fact * ( grho(:,1) * t + rho * gt(:,1) )
        gp(:,2) = fact * ( grho(:,2) * t + rho * gt(:,2) )
!       gp(:,3) = fact * ( grho(:,3) * t + rho * gt(:,3) )

!.... compute strain rate tensor for the viscous terms

        S(:,1,1) = gu(:,1,1)
        S(:,2,1) = pt5 * ( gu(:,2,1) + gu(:,1,2) )
        S(:,3,1) = pt5 * gu(:,3,1)

        S(:,1,2) = pt5 * ( gu(:,1,2) + gu(:,2,1) )
        S(:,2,2) = gu(:,2,2)
        S(:,3,2) = pt5 * gu(:,3,2)

        S(:,1,3) = pt5 * gu(:,3,1) 
        S(:,2,3) = pt5 * gu(:,3,2)

!.... correct the pressure near the outflow plane

        if (idamp.eq.1) then
          gp(:,1) = damp * gp(:,1) + (one-damp) * gpr(:,1)
          gp(:,2) = damp * gp(:,2) + (one-damp) * gpr(:,2)
        end if
        
!.... compute the material properties

        call getmat(t(:),                               &
                    mu(:),   lm(:),    con(:),          &
                    dmu(:),  d2mu(:),  dlm(:),          &
                    d2lm(:), dcon(:),  d2con(:)         )
        
!.... compute gradients of viscosity using chain-rule

        g1mu = dmu * gt(:,1)
        g2mu = dmu * gt(:,2)
!       g3mu = dmu * gt(:,3)

        g1dmu = d2mu * gt(:,1)
        g2dmu = d2mu * gt(:,2)
!       g3dmu = d2mu * gt(:,3)

!.... compute gradients of conductivity using chain-rule

        g1con = dcon * gt(:,1)
        g2con = dcon * gt(:,2)
!       g3con = dcon * gt(:,3)

        g1dcon = d2con * gt(:,1)
        g2dcon = d2con * gt(:,2)
!       g3dcon = d2con * gt(:,3)

!.... compute gradients of bulk viscosity using chain-rule

        g1lm = dlm * gt(:,1)
        g2lm = dlm * gt(:,2)
!       g3lm = dlm * gt(:,3)

        g1dlm = d2lm * gt(:,1)
        g2dlm = d2lm * gt(:,2)
!       g3dlm = d2lm * gt(:,3)

!=============================================================================!
!                       N O N L I N E A R   R H S
!=============================================================================!

!.... Everything you need to compute the nonlinear RHS for vl is in
!.... the local module right now

        call rhs(rl, vl)
        
!.... multiply by the local time step.  
!.... Negate because RK is written as y' = f(y)

        do idof = 1, ndof
          rl(:,idof) = -dtl(:) * rl(:,idof)
        end do
        
        call rkbc(rl)

!.... output RHS statistics

        if (mod(lstep,itout).eq.0 .and. iter.eq.4) then
          call resstat(rl, dtl)
        end if
        
        return
        end
        

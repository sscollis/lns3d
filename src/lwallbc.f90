!=============================================================================!
        subroutine lwallbc(vl, vml, rhow, gnrhow) 
!=============================================================================!
        use global
        use material
        use stencil
        implicit none

        real :: vl(ny,nx,ndof), vml(ny,nx,ndof), rhow(nx), gnrhow(nx)
        
!.... local wall variables

        integer :: i, j, ij, idof, jdof
        
        real :: g1vm(nx,ndof), g2vm(nx,ndof)
        real :: g1v(nx,ndof), g2v(nx,ndof)
        real :: g11v(nx,ndof), g12v(nx,ndof), g22v(nx,ndof)
        real :: g1vl, g2vl, g11vl, g12vl, g22vl

        real :: g1un(nx), g2un(nx), g1unm(nx), g2unm(nx), gnrhom(nx)

!.... stuff for the wall normal pressure gradient
        
        real :: bn1(nx), bn2(nx)
!=============================================================================!

!.... Compute first derivatives of field in the mapped space

        call wgrad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                   xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... Compute second derivatives of field
        
        call wgrad2(ndof, nx, ny, vl, g2v, g11v, g12v, g22v, dxi, deta, &
                    optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall BC on the field

        if (Navier .and. (wallt.eq.1 .or. wallt.eq.2)) g2v(:,ndof) = zero
        

#if 0

!.... transform the gradients to physical space

        j = 1
        do idof = 1, ndof
          do i = 1, nx
            ij = j + (i-1)*ny
            g1vl  = g1v(i,idof)*m1(ij) + g2v(i,idof)*n1(ij)
            g2vl  = g1v(i,idof)*m2(ij) + g2v(i,idof)*n2(ij)
  
            g11vl = g11v(i,idof)       * m1m1(ij)       + &
                    two * g12v(i,idof) * m1n1(ij)       + &
                    g22v(i,idof)       * n1n1(ij)       + &
                    g1v(i,idof)        * m11(ij)        + &
                    g2v(i,idof)        * n11(ij)
  
            g12vl = g11v(i,idof)       * m1m2(ij)       + &
                    g12v(i,idof)       * m1n2(ij)       + &
                    g12v(i,idof)       * m2n1(ij)       + &
                    g22v(i,idof)       * n1n2(ij)       + &
                    g1v(i,idof)        * m12(ij)        + &
                    g2v(i,idof)        * n12(ij)
  
            g22vl = g11v(i,idof)       * m2m2(ij)       + &
                    two * g12v(i,idof) * m2n2(ij)       + &
                    g22v(i,idof)       * n2n2(ij)       + &
                    g1v(i,idof)        * m22(ij)        + &
                    g2v(i,idof)        * n22(ij)
  
            g1v(i,idof)  = g1vl
            g2v(i,idof)  = g2vl
            g11v(i,idof) = g11vl
            g12v(i,idof) = g12vl
            g22v(i,idof) = g22vl
          end do
        end do

!.... now compute some derivatives of the mean field

!.... Compute first derivatives of field in the mapped space

        call wgrad(ndof, nx, ny, vml, g1vm, g2vm, dxi, deta, optx, opty, &
                   xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall BC on the field

        if (Navier .and. (wallt.eq.1 .or. wallt.eq.2)) g2vm(:,ndof) = zero
        
!.... transform the gradients to physical space

        j = 1
        do idof = 1, ndof
          do i = 1, nx
            ij = j + (i-1)*ny
            g1vl  = g1vm(i,idof)*m1(ij) + g2vm(i,idof)*n1(ij)
            g2vl  = g1vm(i,idof)*m2(ij) + g2vm(i,idof)*n2(ij)
            g1vm(i,idof)  = g1vl
            g2vm(i,idof)  = g2vl
          end do
        end do

!=============================================================================!

!.... compute the boundary normal unit-vector

        j = 1
        do i = 1, nx
          ij = j + (i-1) * ny
          bn1(i)  = n1(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
          bn2(i)  = n2(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
        end do    

!.... compute the gradients of the normal velocity

        g1un(:) = bn1(:) * g1v(:,2) + bn2(:) * g1v(:,3)
        g2un(:) = bn1(:) * g2v(:,2) + bn2(:) * g2v(:,3)

        g1unm(:) = bn1(:) * g1vm(:,2) + bn2(:) * g1vm(:,3)
        g2unm(:) = bn1(:) * g2vm(:,2) + bn2(:) * g2vm(:,3)
                
        gnrhom(:) = bn1(:) * g1vm(:,1) + bn2(:) * g2vm(:,1)
        
!.... wall normal density gradient from the normal momentum equation
!.... Note that gnrhow is normalized by 1/sqrt(n1^2 + n2^2)
!.... WARNING: This assumes that T,n = 0 which is valid for 
!.... inviscid and/or adiabatic flows (WRONG)

        do i = 1, nx
          ij = j + (i-1) * ny
          gnrhow(i) = -gamma * Ma**2 * ( &
                      vml(j,i,1) * ( vml(j,i,2) * g1un(i) + &
                      vml(j,i,3) * g2un(i) + vl(j,i,2) * g1unm(i) + &
                      vl(j,i,3) * g2unm(i) ) + &
                      vl(j,i,1) * ( vml(j,i,2) * g1unm(i) + & 
                      vml(j,i,3) * g2unm(i) ) - gnrhom(i) * vl(j,i,5) &
                      ) / ( vml(j,i,5) * sqrt( n1(ij)**2 + n2(ij)**2 ) )
        end do

!.... now compute the wall density.  

        j = 1
        do i = 1, nx
          ij = j + (i-1) * ny
          if (carp) then
            rhow(i) = ( deta * gnrhow(i) -              &
                        gg2 * vl(2,i,1) -               &
                        gg3 * vl(3,i,1) -               &
                        gg4 * vl(4,i,1) -               &
                        gg5 * vl(5,i,1) -               &
                        gg6 * vl(6,i,1) ) / gg1
          else
            rhow(i) = ( deta * gnrhow(i) -              &
                        gc2 * vl(2,i,1) -               &
                        gc3 * vl(3,i,1) -               &
                        gc4 * vl(4,i,1) -               &
                        gc5 * vl(5,i,1) ) / gc1
          end if
        end do

#else

      write(*,*) "Need to fix lwall.f90"

#endif   

        return
        end

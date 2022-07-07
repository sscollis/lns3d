!=============================================================================!
        program pre
!
!  This is a pre processor for the LNS code.  It takes a plot3d grid 
!  file as input and computes the metric file needed for LNS.
!
!=============================================================================!
        use const
        implicit none
        
        integer :: i, j, k, ier, ians
        integer :: nx, ny, nz
        
        real, allocatable :: xy(:,:,:),  g1xy(:,:,:), g2xy(:,:,:)
        real, allocatable :: g11xy(:,:,:), g12xy(:,:,:), g22xy(:,:,:)
        real, allocatable :: g1m(:,:,:), g2m(:,:,:)
        real, allocatable :: g1n(:,:,:), g2n(:,:,:)
        real, allocatable :: xi(:),      eta(:),      IdetJ(:,:)
        real, allocatable :: m(:,:,:),   m11(:,:),    m12(:,:),     m22(:,:)
        real, allocatable :: n(:,:,:),   n11(:,:),    n12(:,:),     n22(:,:)
        real, allocatable :: q(:,:,:)
        
        real, allocatable :: f(:,:)
        real, allocatable :: g1(:,:), g2(:,:), g1l(:,:), g2l(:,:)
        real, allocatable :: g11(:,:),  g12(:,:),  g22(:,:)
        real, allocatable :: g11l(:,:), g12l(:,:), g22l(:,:)
        
        real :: dxi, deta, tmp, kappa

!.... stuff for B-spline

        integer           :: kxord, kyord
        real, allocatable :: x(:,:), y(:,:), xknot(:), yknot(:)
        real, allocatable :: bsx(:,:), bsy(:,:)
        character*1 :: ans
        logical :: bspline

        character*80 :: name, filen
        
        logical :: xper = .false., yper = .false.
        
!=============================================================================!
        write(*,"('Hardwire periodicity (1,0) ==> ',$)") 
        read(*,*) ians
        if (ians.eq.1) then
          xper = .true.
        else
          xper = .false.
        end if

!.... read the physical grid file

 10     filen = 'grid.dat'
        write (*,"('Enter unformatted grid file name [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a20)") name
        if (name(1:1) .ne. ' ') filen = name

!.... allocate memory for the physical grid

        open (unit=10, file=filen, form='unformatted', status='unknown')
        read(10) nx, ny, nz
        write(*,"('Nx = ',i4,', Ny = ',i4)") nx, ny
        allocate (xy(ny,2,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$', &
            'Insufficient Memory for phys. grid$')

        read(10) (((xy(j,1,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((xy(j,2,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((      tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

        do i = 1, nx
          write(*,*) i, xy(1,1,i)
        end do
        
!.... allocate memory for the computational grid

        allocate (xi(nx), eta(ny), STAT=ier)
        if (ier .ne. 0) call error('pre$', &
           'Insufficient Memory for comp. grid$')

!.... make the computational grid

        dxi  = one / float(nx-1)
        deta = one / float(ny-1)
        
        do i = 1, nx
          xi(i) = (i-1) * dxi
        end do
        
        do j = 1, ny
          eta(j) = (j-1) * deta
        end do
        
!.... compute the B-spline interpolant

!       write(*,"(/,'Use (B)-spline or (f)inite difference ==> ',$)")
!       read(*,*) ans
!       if (ans.eq.'B'.or.ans.eq.'b') then
!         bspline=.true.
!       else
!         bspline=.false.
!       endif

        bspline = .false.       !.... hardwired finite difference

        if (bspline) then

        write(*,"('Enter korder ==> ',$)")
        read(*,*) kxord
        kyord = kxord

        allocate (xknot(nx+kxord), yknot(ny+kyord), &
                  bsx(ny,nx), bsy(ny,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for B-spline$')

!       call BSNAK( nx, xi,  kxord, xknot)
!       call BSNAK( ny, eta, kyord, yknot)
!       call BS2IN( ny, eta, nx, xi, xy(:,1,:), ny, kyord, kxord, &
!                    yknot, xknot, bsx)
!       call BS2IN( ny, eta, nx, xi, xy(:,2,:), ny, kyord, kxord, &
!                    yknot, xknot, bsy)

        end if

!.... allocate space for first derivatives

        allocate (g1xy(ny,2,nx), g2xy(ny,2,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for gxy$')

        if (bspline) then

!.... compute first derivative metrics using B-splines

!       call BS2GD( 0, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g1xy(:,1,:), ny )

!       call BS2GD( 0, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsy, g1xy(:,2,:), ny )

!       call BS2GD( 1, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g2xy(:,1,:), ny )

!       call BS2GD( 1, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsy, g2xy(:,2,:), ny )

        else

!.... compute first derivative metrics using finite differences
        
        call grad(2, nx, ny, xy, g1xy, g2xy, dxi, deta, &
                  -1, -1, xper, yper)

        end if

!.... form the determinate of the jacobian

        allocate (IdetJ(ny,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for detJ$')

        IdetJ(:,:) = one / ( g1xy(:,1,:) * g2xy(:,2,:) - &
                             g2xy(:,1,:) * g1xy(:,2,:) )
        
!.... allocate memory for the grid metrics

        allocate( m(ny,2,nx),  m11(ny,nx), m12(ny,nx), m22(ny,nx), &
                  n(ny,2,nx),  n11(ny,nx), n12(ny,nx), n22(ny,nx), &
                  STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for metrics$')

!.... form the first derivative metrics
        
        m(:,1,:) =  g2xy(:,2,:) * IdetJ
        n(:,1,:) = -g1xy(:,2,:) * IdetJ
        m(:,2,:) = -g2xy(:,1,:) * IdetJ
        n(:,2,:) =  g1xy(:,1,:) * IdetJ

        allocate (g1m(ny,2,nx), g2m(ny,2,nx), &
                  g1n(ny,2,nx), g2n(ny,2,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for gm and gn$')

        if (bspline) then

!.... compute second derivative metrics using B-splines

!       call BS2IN( ny, eta, nx, xi, m(:,1,:), ny, kyord, kxord, &
!                    yknot, xknot, bsx)

!       call BS2GD( 0, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g1m(:,1,:), ny )
!       call BS2GD( 1, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g2m(:,1,:), ny )

!       call BS2IN( ny, eta, nx, xi, m(:,2,:), ny, kyord, kxord, &
!                    yknot, xknot, bsx)

!       call BS2GD( 0, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g1m(:,2,:), ny )
!       call BS2GD( 1, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g2m(:,2,:), ny )

!       call BS2IN( ny, eta, nx, xi, n(:,1,:), ny, kyord, kxord, &
!                    yknot, xknot, bsx)

!       call BS2GD( 0, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g1n(:,1,:), ny )
!       call BS2GD( 1, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g2n(:,1,:), ny )

!       call BS2IN( ny, eta, nx, xi, n(:,2,:), ny, kyord, kxord, &
!                    yknot, xknot, bsx)

!       call BS2GD( 0, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g1n(:,2,:), ny )
!       call BS2GD( 1, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g2n(:,2,:), ny )

        else

!.... compute second derivatives of x using finite differences

        allocate (g11xy(ny,2,nx), g12xy(ny,2,nx), g22xy(ny,2,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for gxy$')

        call grad(2, nx, ny, g1xy, g11xy, g12xy, dxi, deta, &
                  -1, -1, xper, yper)
        call grad(2, nx, ny, g2xy, g12xy, g22xy, dxi, deta, &
                  -1, -1, xper, yper)

!.... compute second derivative metrics using finite difference

        call grad(2, nx, ny, m, g1m, g2m, dxi, deta, &
                  -1, -1, xper, yper)

        call grad(2, nx, ny, n, g1n, g2n, dxi, deta, &
                  -1, -1, xper, yper)
        
        end if  

        m11 = m(:,1,:) * g1m(:,1,:) + n(:,1,:) * g2m(:,1,:)
        n11 = m(:,1,:) * g1n(:,1,:) + n(:,1,:) * g2n(:,1,:)

        m22 = m(:,2,:) * g1m(:,2,:) + n(:,2,:) * g2m(:,2,:)
        n22 = m(:,2,:) * g1n(:,2,:) + n(:,2,:) * g2n(:,2,:)
        
        m12 = m(:,1,:) * g1m(:,2,:) + n(:,1,:) * g2m(:,2,:)     
        n12 = m(:,1,:) * g1n(:,2,:) + n(:,1,:) * g2n(:,2,:)
                
        j = 1
        do i = 1, nx
          write(55,21) xy(j,1,i), m12(j,i), &
                      m(j,2,i)*g1m(j,1,i)+n(j,2,i)*g2m(j,1,i), &
                      m12(j,i) - (m(j,2,i)*g1m(j,1,i)+n(j,2,i)*g2m(j,1,i)), &
                      m(j,1,i)*g2m(j,1,i)+m(j,2,i)*g2m(j,2,i), &
                      n(j,1,i)*g1n(j,1,i)+n(j,2,i)*g1n(j,2,i)

          write(56,21) xy(j,1,i), n(j,1,i), n(j,2,i), n11(j,i), n12(j,i), &
                       n22(j,i)
 21       format(7(1pe14.7,1x))
        end do

        deallocate( g1m, g2m, g1n, g2n )

!.... write out the metric file

        filen = 'metric.new'
        write (*,"('Enter grid metric file name [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a20)") name
        if (name(1:1) .ne. ' ') filen = name
        open (unit=10, file=filen, form='unformatted', status='unknown')
        write(10) m(:,1,:), m(:,2,:), n(:,1,:), n(:,2,:), &
                  m11, m12, m22, n11, n12, n22
        close(10)

!.... write out plot3d files

        allocate (q(ny,5,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for q$')

        q(:,1,:) = m(:,1,:)
        q(:,2,:) = n(:,1,:)
        q(:,3,:) = m(:,2,:)
        q(:,4,:) = n(:,2,:)
        q(:,5,:) = IdetJ
        
!.... diagnostic

        j = 1
        do i = 1, nx
          write(54,21) xy(j,1,i), &
                     (-m(j,1,i)*m(j,2,i)*m11(j,i) + m(j,1,i)**2*m12(j,i) - &
                     m(j,2,i)**2*m12(j,i) + m(j,1,i)*m(j,2,i)*m22(j,i) ) * &
                     sqrt(n(j,1,i)**2 + n(j,2,i)**2) /                  &
                     sqrt(m(j,1,i)**2 + m(j,2,i)**2) /                  &
                     (m(j,1,i)*n(j,2,i)-m(j,2,i)*n(j,1,i)),             &
                     (n(j,1,i)*n(j,2,i)*n11(j,i) - n(j,1,i)**2*n12(j,i) + &
                     n(j,2,i)**2*n12(j,i) - n(j,1,i)*n(j,2,i)*n22(j,i) ) * &
                     sqrt(m(j,1,i)**2 + m(j,2,i)**2) /                  &
                     sqrt(n(j,1,i)**2 + n(j,2,i)**2) /                  &
                     (m(j,1,i)*n(j,2,i)-m(j,2,i)*n(j,1,i))
        end do

        name = 'm1.dat'
        call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)
                
        q(:,1,:) = m(:,1,:) - 2.0 * xy(:,2,:) / ( acos(-1.0) * &
                              ( (xy(:,1,:)-0.1)**2 + xy(:,2,:)**2 ) )
        q(:,2,:) = n(:,1,:) - (xy(:,1,:) - 0.1) / sqrt( (xy(:,1,:)-0.1)**2 + &
                               xy(:,2,:)**2 )
        q(:,3,:) = m(:,2,:) - 2.0 * (0.1 - xy(:,1,:)) / ( acos(-1.0) * &
                              ( (xy(:,1,:)-0.1)**2 + xy(:,2,:)**2 ) )
        q(:,4,:) = n(:,2,:) - xy(:,2,:) / sqrt( (xy(:,1,:)-0.1)**2 + &
                              xy(:,2,:)**2 )
        q(:,5,:) = IdetJ
        
        name = 'error.dat'
        call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

        q(:,1,:) = m11
        q(:,2,:) = m12
        q(:,3,:) = m12
        q(:,4,:) = m22
        q(:,5,:) = IdetJ
        
        name = 'm11.dat'
        call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

        q(:,1,:) = n11
        q(:,2,:) = n12
        q(:,3,:) = n12
        q(:,4,:) = n22
        q(:,5,:) = IdetJ
        
        name = 'n11.dat'
        call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

        q(:,1,:) = g11xy(:,1,:)
        q(:,2,:) = g11xy(:,2,:)
        q(:,3,:) = g12xy(:,1,:)
        q(:,4,:) = g12xy(:,2,:)
        q(:,5,:) = g22xy(:,1,:)
        
        name = 'grad.dat'
        call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

        deallocate( g1xy, g2xy)

!.... write out the derivatives on the body

!       write(*,"('Enter j ==> ',$)") 
!       read(*,*) j
!       do i = 1, nx
!         write(20,20) xi(i), m(j,1,i), m(j,2,i)
!         write(21,20) xi(i), n(j,1,i), n(j,2,i)
!         write(22,20) xi(i), m11(j,i), n11(j,i), m22(j,i), m12(j,i), n12(j,i)
!         write(23,20) xi(i), IdetJ(j,i), one / Idetj(j,i)
!       end do

!==============================================================================
!
! Test the derivative operators
!
!==============================================================================

!.... put a sine wave on the mesh

        allocate (f(ny,nx), g1l(ny,nx), g2l(ny,nx), &
                  g1(ny,nx), g2(ny,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for function$')

        kappa = two * pi / 2.5
        f(:,:) = sin( kappa * xy(:,1,:))

        if (bspline) then

!.... compute the gradient using B-splines

!       call BS2IN( ny, eta, nx, xi, f, ny, kyord, kxord, &
!                    yknot, xknot, bsx)

!       call BS2GD( 0, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g1l(:,:), ny )

!       call BS2GD( 1, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g2l(:,:), ny )

        else

!.... compute the gradient with finite differences

        call grad(1, nx, ny, f, g1l, g2l, dxi, deta, &
                  -1, -1, xper, yper)

        endif

!.... transform to physical coordinates

        g1 = g1l(:,:) * m(:,1,:) + g2l(:,:) * n(:,1,:)
        g2 = g1l(:,:) * m(:,2,:) + g2l(:,:) * n(:,2,:)

        q(:,1,:) = f
        q(:,2,:) = g1 / kappa
        q(:,3,:) = g2
        q(:,4,:) = abs(g1 - kappa * cos( kappa * xy(:,1,:) ) )
        q(:,5,:) = abs(g2 - zero)
        
        name = 'f.dat'
        call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

        deallocate( g1, g2)
        
!... allocate memory for second derivatives

        allocate (g11(ny,nx),  g12(ny,nx),  g22(ny,nx), &
                  g11l(ny,nx), g12l(ny,nx), g22l(ny,nx), STAT=ier)
        if (ier .ne. 0) call error('pre$','Insufficient Memory for g2$')

        if (bspline) then

!.... compute second derivatives using B-splines

!       call BS2GD( 0, 2, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g11l, ny )

!       call BS2GD( 1, 1, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g12l, ny )

!       call BS2GD( 2, 0, ny, eta, nx, xi, kyord, kxord, yknot, xknot, &
!                    ny, nx, bsx, g22l, ny )

        else

!.... compute second derivatives using finite differencs

        call grad2(1, nx, ny, f, g11l, g12l, g22l, dxi, deta, &
                   -1, -1, xper, yper)

        end if

        g11  =  g11l(:,:)        * m(:,1,:) * m(:,1,:)  + &
                two * g12l(:,:)  * m(:,1,:) * n(:,1,:)  + &
                g22l(:,:)        * n(:,1,:) * n(:,1,:)  + &
                g1l(:,:)         * m11(:,:)             + &
                g2l(:,:)         * n11(:,:)

        g12  =  g11l(:,:)        * m(:,1,:) * m(:,2,:)  + &
                g12l(:,:)        * m(:,1,:) * n(:,2,:)  + &
                g12l(:,:)        * m(:,2,:) * n(:,1,:)  + &
                g22l(:,:)        * n(:,1,:) * n(:,2,:)  + &
                g1l(:,:)         * m12(:,:)             + &
                g2l(:,:)         * n12(:,:)

        g22  =  g11l(:,:)        * m(:,2,:) * m(:,2,:)  + &
                two * g12l(:,:)  * m(:,2,:) * n(:,2,:)  + &
                g22l(:,:)        * n(:,2,:) * n(:,2,:)  + &
                g1l(:,:)         * m22(:,:)             + &
                g2l(:,:)         * n22(:,:)

        q(:,1,:) = g11 / kappa**2
        q(:,2,:) = g12
        q(:,3,:) = g22
        q(:,4,:) = abs(g11/kappa**2 + sin( kappa * xy(:,1,:) ) )
        q(:,5,:) = zero
        
        name = 'g2.dat'
        call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

        stop
        
 20     format(8(1pe13.6,1x))

        end

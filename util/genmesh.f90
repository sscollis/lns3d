!-----------------------------------------------------------------------------
        program genmesh
!
!  Generate a simple mesh for temporal instability waves and acoustic waves
!
!  11-14-99:  Switched i,j indices [SSC]
!-----------------------------------------------------------------------------
        use const
        implicit none

        integer :: i, j, k

        real, allocatable :: x(:), y(:), eta(:), m1(:), m2(:)
        real, allocatable :: m11(:), m12(:), m22(:), n1(:), n2(:)
        real, allocatable :: n11(:), n12(:), n22(:)

        real    :: dx, dy, dxi, deta
        integer :: nx, ny

        character(80) name

!.... default values are for temporal instability calculation

        real    :: ymin = 0.0, ymax = 8.0
        real    :: xmin = 0.0, xmax = 2.03589244590181

!.... parameters for algebraic mapping

        real :: yi, aa, bb

!.... parameters for hyperbolic tangent mapping

        real :: ds1, ds2, dd, xi
!-----------------------------------------------------------------------------
        write(*,*) 'GenMesh:  simple 2d mesh generator.'
        write(*,*) 'NOTE:  Indices are I,J'
        write(*,"('Enter nx, ny ==> ',$)")
        read(*,*) nx, ny
        write(*,"('Enter Yi (Yi=0 uniform, Yi>0 algebraic, Yi<0 Tanh) ==> ',$)")
        read(*,*) yi
        write(*,"('Enter xmin, xmax ==> ',$)")
        read(*,*) xmin, xmax
        write(*,"('Enter ymin, ymax ==> ',$)")
        read(*,*) ymin, ymax

!.... allocate storage

        allocate ( x(nx), y(ny), eta(ny), m1(nx), m2(nx), &
                   m11(nx), m12(nx), m22(nx), n1(ny), n2(ny), &
                   n11(ny), n12(ny), n22(ny) )

        dxi  = one / float(nx-1)
        deta = one / float(ny-1)

!.... algebraic grid or uniform y-grid

        if (yi .gt. zero) then
          write(*,*) 'Algebraic grid in y'
          aa = ymax * yi / ( ymax - two * yi )
          bb = one + aa / ymax
          do j = 1, ny
            eta(j)  = (j-1) * deta
            y(j)    = aa * eta(j) / (bb - eta(j))
            n1(j)   = zero
            n2(j)   = (bb - eta(j))**2 / (aa * bb)
            n11(j)  = zero
            n12(j)  = zero
            n22(j)  = -two * (bb - eta(j))**3 / (aa * bb)**2
          end do
        else if (yi .lt. zero) then
          write(*,*) 'Tanh grid in y'
          ds1 = 0.0005d0        ! set for R=2400 MSE, r=50
          ds2 = 5.0d0
          dd  = 5.36966703089523d0
          do j = 1, ny
            eta(j) = (j-1) * deta
            y(j) = ymax*(0.5 + 1.0/Tanh(dd/2.0)*Tanh(dd*(-0.5 +         &
                    eta(j)))/2.0)/(Sqrt(ds2/ds1) +                      &
                    (1.0 - Sqrt(ds2/ds1))*                              &
                    (0.5 + 1.0/Tanh(dd/2.0)*Tanh(dd*(-0.5 +             &
                    eta(j)))/2.0))

            n1(j) = zero
            n2(j) = Cosh(dd*(-0.5 + eta(j)))*1.0/Cosh(dd/2.0)*          &
                    (Sqrt(ds2/ds1)*Sinh(dd*(1.0 - eta(j))) +            &
                    Sinh(dd*eta(j)))**2/(dd*Sqrt(ds2/ds1)*ymax*         &
                    (Sinh(dd*(1.0 - eta(j))) + Sinh(dd*eta(j))))

            n11(j) = zero
            n12(j) = zero
            n22(j) = ds1*Cosh(dd*(-0.5 + eta(j)))*                      &
                  (Cosh(dd*(1.0 - 3.0*eta(j))) -                        &
                  Sqrt(ds2/ds1)*Cosh(dd*(2.0 - 3.0*eta(j))) +           &
                  Cosh(dd*(1.0 - eta(j))) -                             &
                  2.0*Sqrt(ds2/ds1)*Cosh(dd*(1.0 - eta(j))) +           &
                  2.0*Cosh(dd*eta(j)) -                                 &
                  Sqrt(ds2/ds1)*Cosh(dd*eta(j)))*1.0/Sinh(dd/2.0)*      &
                  1.0/Cosh(dd/2.0)**2*1.0/Cosh(dd/2.0 - dd*eta(j))**2*  &
                  (Sqrt(ds2/ds1)*Sinh(dd*(1.0 - eta(j))) +              &
                  Sinh(dd*eta(j)))**3/                                  &
                  (4.0*dd*ds2*ymax**2*                                  &
                  (Sinh(dd*(1.0 - eta(j))) + Sinh(dd*eta(j))))
          end do
        else
          write(*,*) 'Uniform grid in y'
          dy = (ymax-ymin)/float(ny-1)
          do j = 1, ny
            y(j)   = ymin + (j-1)*dy
            n1(j)  = zero
            n2(j)  = one / (ymax-ymin)
            n11(j) = zero
            n12(j) = zero
            n22(j) = zero
          end do
        end if

!.... make the uniform x-grid

        dx = (xmax-xmin)/float(nx-1)
        do i = 1, nx
          x(i) = xmin + (i-1)*dx
          m1(i)  = one / (xmax-xmin)
          m2(i)  = zero
          m11(i) = zero
          m12(i) = zero
          m22(i) = zero
        end do

!.... write out a coord.dat file

        open(10,file='coord.dat')
        do j = 1, ny
          write(10,10) j, 0, 0.0, y(j)
 10       format(2(i5,1x),2(1pe20.13,1x))
        end do
        close(10)

!.... write out the grid file

        open (unit=10, file='grid.dat', form='unformatted')
        write(10) nx, ny, 1
        write(10) ((x(i), i = 1, nx), j = 1, ny), &
                  ((y(j), i = 1, nx), j = 1, ny), &
                  (( 0.0, i = 1, nx), j = 1, ny)
        close(10)

        !.... write out the body file

        open (unit=10, file='body.dat', form='formatted')
        write(10,20) (x(i), x(i), i = 1, nx)
 20     format(2(1pe20.13,1x))
        close(10)

!.... write out the metric file

        open (unit=10, file='metric.dat', form='unformatted')
        write(10) (( m1(i), i = 1, nx), j = 1, ny), &
                  (( m2(i), i = 1, nx), j = 1, ny), &
                  (( n1(j), i = 1, nx), j = 1, ny), &
                  (( n2(j), i = 1, nx), j = 1, ny), &
                  ((m11(i), i = 1, nx), j = 1, ny), &
                  ((m12(i), i = 1, nx), j = 1, ny), &
                  ((m22(i), i = 1, nx), j = 1, ny), &
                  ((n11(j), i = 1, nx), j = 1, ny), &
                  ((n12(j), i = 1, nx), j = 1, ny), &
                  ((n22(j), i = 1, nx), j = 1, ny)
        close(10)

        stop
        end

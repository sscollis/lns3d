!============================================================================!
        module bspline

        integer           :: ns, kord
        real, allocatable :: knot(:)
        real, allocatable :: bs(:,:)

        real, external :: bsder

        real :: rhoe, ue, ve, we, te, pe, phie
        real :: alpha, Ma, gamma

        end module bspline
!============================================================================!
        program npost 
!  
!  Post processor for nonlinear calculations
!  See the usage for a description of the command-line arguments
!  
!  Author:  Scott Collis
!
!  Revised: 7-25-96     Original version
!           11-14-99    Switched indices
!           11-22-99    Added periodicity and symmetry flags
!           02-16-00    Changed file name convention
!============================================================================!
        use bspline
        use const
        implicit none

!.... flow data

        real, allocatable :: v(:,:,:), q(:,:,:), p(:,:)
        !$sgi distribute v(*,*,block), q(*,*,block), p(*,block)
        
        real, allocatable :: g1v(:,:,:), g2v(:,:,:)
        real, allocatable :: g11v(:,:,:), g12v(:,:,:), g22v(:,:,:)
        !$sgi distribute g1v(*,*,block), g2v(*,*,block)
        !$sgi distribute g11v(*,*,block), g12v(*,*,block), g22v(*,*,block)

        real :: g1vl, g2vl, g11vl, g12vl, g22vl

!.... mesh

        real, allocatable :: x(:,:), y(:,:), xi(:), eta(:), s(:), z(:)
        real :: dxi, deta, dz

!.... metrics

        real, allocatable :: m1(:,:),  m2(:,:),  n1(:,:),  n2(:,:),     &
                             m11(:,:), m12(:,:), m22(:,:),              &
                             n11(:,:), n12(:,:), n22(:,:)
        !$sgi distribute m1(*,block), m2(*,block), 
        !$sgi distribute m11(*,block), m12(*,block), m22(*,block)
        !$sgi distribute n1(*,block), n2(*,block), 
        !$sgi distribute n11(*,block), n12(*,block), n22(*,block)

        real :: m1l, m2l, bn1, bn2

!.... parameters

        real    :: Re, Pr, cv, Rgas, time, vinf, gamma1
        integer :: lstep, ny, nx, nz, ndof

        integer :: i, j, k, idof, itmp

        character(80) file, temp, base, base1, base2, fname, fname1, fname2
        integer :: iloc, iend, imin, imax, inc, iver
#ifndef __GFORTRAN__    
        integer, external :: iargc
#endif
        integer :: nx2, ny2, nz2

        real :: tmp, arc, alp, umag
        
        real :: b, beta, eps, Lz
        
        integer :: optx=-1, opty=-1
        logical :: xper=.false., yper=.false.
        logical :: lsym=.false., rsym=.false., bsym=.false., tsym=.false.
        logical :: carp=.false.

!.... variables for thickness calculation

        real :: rho, u1, u2, u3, t, delta, theta, us, un, ws, edge, errest
        real :: H0, Hinf, thetaz
        
        real, allocatable :: ynn(:), uss(:), unn(:), wss(:), uel(:)

        real, external :: rtsafe
        real, external :: dfun, thfun, thfunz
        external :: fedge, fwmax, fwinf

        real    :: wmax, wmloc, w1, w2, wiloc, uinf, winf, epsi
        integer :: jmax, ierr

!.... flags

        logical :: profile=.false., thick=.false., debug=.false.
        logical :: swept=.false., cons=.false., lfilter=.false.
        logical :: plot3d=.true., lst=.false., inviscid=.false.
        logical :: readcons=.false.
        logical :: metric_ji=.false., v_ji=.false.
        integer :: type=0
        integer, parameter :: mfile = 5
        integer :: narg, iarg, nfile=0, ifile(mfile)
        character(80) :: arg
        
!.... potential flow stuff and boundary layer stuff

        real :: c1, Ma1, Cp, phip, Re1, ReL, L, L1, L2, L3, &
                betah, beta2, beta3, H, rk, p0, pinf, dusds, cur, fact

!$      integer, external :: omp_get_num_threads, omp_get_thread_num
!$      integer, external :: omp_get_num_procs

!============================================================================!

!.... parse the argument list

        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          if (arg(1:1) .ne. '-') then
            nfile = nfile + 1
            if (nfile .gt. mfile) then
              write(*,*) '>> Error in argument list, too many file names'
              call exit(1)
            end if
            ifile(nfile) = iarg
          else
            select case (arg(1:3))
            case ('-p ')                ! profiles
              profile = .true.
            case ('-l ')                ! lst
              lst = .true.
            case ('-t ')                ! thickness
              thick = .true.
            case ('-d ')                ! debug
              debug = .true.
            case ('-c ')                ! conservation variables
              cons = .true.
            case ('-s ')                ! swept
              swept = .true.
            case ('-f ')                ! filter
              lfilter = .true.
            case ('-n ')                ! output no plot3d files
              plot3d = .false.
            case ('-i ')                ! inviscid flag
              inviscid = .true.
!           case ('-D ')                ! compute a second Plot3d file
!             select case (arg(1:3))
              case ('-D1')              ! 2nd derivative in x
                type = 1
              case ('-D2')              ! 2nd derivative in y
                type = 2
              case ('-D3')              ! xy derivative
                type = 3
              case ('-D4')              ! c, Mx, My, M, entropy
                type = 4
              case ('-D5')              ! dilitation, vorticity, etc
                type = 5
!             case default
!               type = 1
!             end select
            case ('-xp')
              xper = .true.
            case ('-yp')
              yper = .true.
            case ('-rs')
              rsym = .true.
            case ('-ls')
              lsym = .true.
            case ('-ts')
              tsym = .true.
            case ('-bs')
              bsym = .true.
            case ('-ca')
              carp = .true.
            case ('-rc')
              readcons = .true.
            case ('-ms')
              metric_ji = .true.
            case ('-vs')
              v_ji = .true.
            case ('-ji')
              metric_ji = .true.
              v_ji = .true.
            case ('-h')
              write(*,"('-----------------------------------------------')")
              write(*,"('Usage:  npost [options] [file1] [file2] [file3]')")
              write(*,"('-----------------------------------------------')")
              write(*,"('   -h:  this help')")
              write(*,"('-----------------------------------------------')")
              write(*,"('   -p:  output velocity profiles')")
              write(*,"('   -l:  output LST restart file and quit')")
              write(*,"('   -t:  compute BL thickness')")
              write(*,"('   -d:  debug i,j')")
              write(*,"('   -c:  output conservation variables')")
              write(*,"('   -s:  Swept wing -- output 3D files')")
              write(*,"('   -f:  filter the field in y')")
              write(*,"('   -n:  no plot3d files')")
              write(*,"('   -i:  treat as inviscid')")
              write(*,"('  -xp:  x-periodicity')")
              write(*,"('  -yp:  y-periodicity')")
              write(*,"('  -rs:  right symmetry')")
              write(*,"('  -ls:  left symmetry')")
              write(*,"('  -ts:  top symmetry')")
              write(*,"('  -bs:  bottom symmetry')")
              write(*,"('  -ca:  Carpenters stencil')")
              write(*,"('  -rc:  read conservative variables')")
              write(*,"('  -ms:  read metrics assuming ji format')")
              write(*,"('  -vs:  read restart assuming ji format')")
              write(*,"('  -ji:  read all assuming ji format')")
              write(*,"('-----------------------------------------------')")
              write(*,"('  -D1:  output xx derivative')")
              write(*,"('  -D2:  output yy derivative')")
              write(*,"('  -D3:  output xy derivative')")
              write(*,"('  -D4:  output c, Mx, My, M, entropy')")
              write(*,"('  -D5:  vorticity; dilitation; p,x; p,y; E')")
              write(*,"('-----------------------------------------------')")
              call exit(0)
            case default
              write(*,"('Argument ',i2,' ignored.')") iarg
            end select
          end if
        end do

!$omp parallel
!$      if (omp_get_thread_num() == 0) then
!$        write(*,*) 'Running on ',omp_get_num_threads(),' processor(s)'
!$        write(*,*) 'There are ',omp_get_num_procs(),' processor(s) available'
!$      end if
!$omp end parallel

!.... read the restart file

        if ( nfile .gt. 0 ) then
          call getarg(ifile(1),temp)
          file = temp
        else
          file = 'output.res'
          write (*,"('Enter file name [',a,']? ',$)") file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        end if
 10     open(unit=10, file=file, form='unformatted', status='old', err=20)
        goto 30
 20     write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
        write (*,"('Enter file name [',a,']? ',$)") file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 10
 30     continue
        
        read(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
        gamma1 = gamma - one
        Rgas = gamma1 * cv
        if (nz .ne. 1) then
          write(*,"(' Error:  nz <> 1 ')")
          call exit(1)
        end if
        allocate( v(ndof,nx,ny) )
        if (v_ji) then
          read(10) (((v(idof,i,j), j=1,ny), i=1,nx), idof=1,ndof)
          close(10)
          open(unit=10,file='restart.ij',form='unformatted',status='unknown')
          write(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
          write(10) v
          close(10)
        else
          read(10) v
          close(10)
        end if
        write(*,"('Read flow field for ',a)") file(1:index(file,' '))
        write(*,"('  Time = ',1pe10.3,'  step = ',i6)") time, lstep

!.... Convert to Primitive if required

        if (readcons) then
          fact = one / (gamma * Ma**2)
          do j = 1, ny
            do i = 1, nx
              v(5,i,j) = (v(5,i,j) - pt5*( v(2,i,j)**2 + v(3,i,j)**2 + &
                          v(4,i,j)**2) / v(1,i,j) ) * (gamma-one) / &
                          (fact * v(1,i,j))
              v(2,i,j) = v(2,i,j) / v(1,i,j)
              v(3,i,j) = v(3,i,j) / v(1,i,j)
              v(4,i,j) = v(4,i,j) / v(1,i,j)
            end do
          end do

!.... write out a primative restart file

          open(10,file='primitive.dat',form='unformatted')
          write(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
          write(10) v
          close(10)
        end if

!.... if swept then input the sweep angle (wish this were automatic)

        if ( maxval(abs(v(4,:,:))) .gt. 1e-14 ) then
          write(*,"('Enter sweep angle (deg) ==> ',$)")
          read(*,*) alpha
          alpha = alpha * pi / 180.0
        end if
        
!.... read in the grid file

        allocate( x(nx,ny), y(nx,ny), xi(nx), eta(ny), s(nx) )
        open(unit=10,file='grid.dat',form='unformatted',status='old')
        read(10) nx2, ny2, nz2
        if (nx2.ne.nx .or. ny2.ne.ny .or. nz2.ne.nz) then
          write(*,*) 'Grid and data file dimensions do not match'
          call exit(1)
        end if
        read(10) (((x(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... read in the body discription file

        s = x(:,1)
        open(unit=10,file='body.dat',form='formatted',status='old',err=100)
        do i = 1, nx
           read(10,*) tmp, s(i)
        end do
        close(10)
        goto 101
 100    continue

!.... if no body file is found, assume that its a parabolic cylinder

        write(*,*) ' WARNING:  assumming a parabolic cylinder'

!       if (x(1,1) .lt. zero) then
!         x(1,1) = zero         ! quick fix for roundoff error at origin
!         y(1,1) = zero         ! better to actually fix the mesh generator
!       end if

        s(:) = Sqrt(x(:,1) + two*x(:,1)**2)/Sqrt(two) + &
               Log(one + four*x(:,1) + two**(onept5)*Sqrt(x(:,1) + &
               two*x(:,1)**2)) / four
 101    continue

!.... make the xi grid
        
        dxi = one / float(nx-1)
        do i = 1, nx
          xi(i) = real(i-1) * dxi
        end do

!.... make the eta grid

        deta = one / float(ny-1)
        do j = 1, ny
          eta(j) = float(j-1) * deta
        end do

!.... filter the field?

        if (lfilter) then
          write(*,*) 'Filtering in Y'
          call filter_y( nx, ny, v(1,:,:) )
          call filter_y( nx, ny, v(2,:,:) )
          call filter_y( nx, ny, v(3,:,:) )
          call filter_y( nx, ny, v(4,:,:) )
          call filter_y( nx, ny, v(5,:,:) )
          write(*,*) 'Filtering in X'
          call filter_x( nx, ny, v(1,:,:) )
          call filter_x( nx, ny, v(2,:,:) )
          call filter_x( nx, ny, v(3,:,:) )
          call filter_x( nx, ny, v(4,:,:) )
          call filter_x( nx, ny, v(5,:,:) )
          open(unit=10,file='filter.dat',form='unformatted',status='unknown')
          write(10) lstep, time, nx, ny, nz, ndof, &
                    Re, Ma, Pr, gamma, cv
          write(10) v
          close(10)
        end if
        
!.... allocate storage for metrics

        allocate (m1(nx,ny),  m2(nx,ny),  n1(nx,ny),  n2(nx,ny), &
                  m11(nx,ny), m12(nx,ny), m22(nx,ny),            &
                  n11(nx,ny), n12(nx,ny), n22(nx,ny) )

!.... read in the metric file

        if (metric_ji) then
          open (unit=10,file='metric.dat',form='unformatted', status='old')
          read(10) ((m1(i,j), j=1,ny),i=1,nx), &
                   ((m2(i,j), j=1,ny),i=1,nx), &
                   ((n1(i,j), j=1,ny),i=1,nx), &
                   ((n2(i,j), j=1,ny),i=1,nx), &
                   ((m11(i,j), j=1,ny),i=1,nx), &
                   ((m12(i,j), j=1,ny),i=1,nx), &
                   ((m22(i,j), j=1,ny),i=1,nx), &
                   ((n11(i,j), j=1,ny),i=1,nx), &
                   ((n12(i,j), j=1,ny),i=1,nx), &
                   ((n22(i,j), j=1,ny),i=1,nx)
          close(10)
          open (unit=10,file='metric.ij',form='unformatted', status='unknown')
          write(10) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
          close(10)
        else
          open (unit=10,file='metric.dat',form='unformatted', status='old')
          read(10) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
          close(10)
        end if
        
!.... Compute first derivatives of field in the mapped space

        allocate( g1v(ndof,nx,ny), g2v(ndof,nx,ny) )
        call grad(ndof, nx, ny, v, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... Compute second derivatives of field
                
        allocate( g11v(ndof,nx,ny), g12v(ndof,nx,ny), g22v(ndof,nx,ny) )
        call grad2(ndof, nx, ny, v, g1v, g11v, g12v, g22v, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... transform the gradients to physical space

        !$omp parallel do private(j, idof, g1vl, g2vl, g11vl, g12vl, g22vl)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              g1vl  = g1v(idof,i,j)*m1(i,j) + g2v(idof,i,j)*n1(i,j)
              g2vl  = g1v(idof,i,j)*m2(i,j) + g2v(idof,i,j)*n2(i,j)
              
              g11vl = g11v(idof,i,j)       * m1(i,j)*m1(i,j)    + &
                      two * g12v(idof,i,j) * m1(i,j)*n1(i,j)    + &
                      g22v(idof,i,j)       * n1(i,j)*n1(i,j)    + &
                      g1v(idof,i,j)        * m11(i,j)           + &
                      g2v(idof,i,j)        * n11(i,j)
  
              g12vl = g11v(idof,i,j)       * m1(i,j)*m2(i,j)    + &
                      g12v(idof,i,j)       * m1(i,j)*n2(i,j)    + &
                      g12v(idof,i,j)       * m2(i,j)*n1(i,j)    + &
                      g22v(idof,i,j)       * n1(i,j)*n2(i,j)    + &
                      g1v(idof,i,j)        * m12(i,j)           + &
                      g2v(idof,i,j)        * n12(i,j)
  
              g22vl = g11v(idof,i,j)       * m2(i,j)*m2(i,j)    + &
                      two * g12v(idof,i,j) * m2(i,j)*n2(i,j)    + &
                      g22v(idof,i,j)       * n2(i,j)*n2(i,j)    + &
                      g1v(idof,i,j)        * m22(i,j)           + &
                      g2v(idof,i,j)        * n22(i,j)
  
              g1v(idof,i,j)  = g1vl
              g2v(idof,i,j)  = g2vl
              g11v(idof,i,j) = g11vl
              g12v(idof,i,j) = g12vl
              g22v(idof,i,j) = g22vl
            end do
          end do
        end do
        
!.... form the output variables

        allocate( q(ndof,nx,ny), p(nx,ny) )
        !$omp parallel do
        do j = 1, ny
          p(:,j) = v(1,:,j) * v(5,:,j) / (gamma * Ma**2)
        end do

!.... make a restart and grid file in the body normal coordinates

        if (lst) then
          write(*,*) ' WARNING: only valid for body-fitted mesh'
!         write(*,*) ' WARNING: parallel flow assumption'
          allocate ( ynn(ny), uss(ny), unn(ny) )
          do i = 1, nx
            bn1 = n1(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
            bn2 = n2(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )       
            do j = 1, ny
              uss(j) = bn2 * v(2,i,j) - bn1 * v(3,i,j)
              unn(j) = bn1 * v(2,i,j) + bn2 * v(3,i,j)
              q(1,i,j) = v(1,i,j)
              q(2,i,j) = uss(j)
              q(3,i,j) = unn(j)
              q(4,i,j) = v(4,i,j)
              q(5,i,j) = v(5,i,j)
            end do
          end do

          arc = zero
          do j = 1, ny
            ynn(j) = arc
            if (j.ne.ny) arc = arc + sqrt( (x(1,j+1)-x(1,j))**2 + &
                                           (y(1,j+1)-y(1,j))**2 )
          end do
          
          open(10,file='lstq.dat',form='unformatted')
          write(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
          write(10) q
          close(10)
          
          open(unit=10,file='lstx.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10)(((  s(i), i = 1, nx), j = 1, ny), k = 1, nz), &
                   (((ynn(j), i = 1, nx), j = 1, ny), k = 1, nz), &
                   (((  zero, i = 1, nx), j = 1, ny), k = 1, nz)
          close(10)
          
          deallocate( ynn, uss, unn )
          call exit(0)
        end if

!.... use either primitive variables or 
!.... conservative variables in Plot3d normalization

        if (cons) then
          q(1,:,:) = v(1,:,:)                                   ! rho
          q(2,:,:) = v(1,:,:) * v(2,:,:) * Ma                   ! rho u
          q(3,:,:) = v(1,:,:) * v(3,:,:) * Ma                   ! rho v
          q(4,:,:) = v(1,:,:) * v(4,:,:) * Ma                   ! rho w
          q(5,:,:) = v(1,:,:) * ( one / gamma1 / gamma * v(5,:,:) + pt5 * &
                     Ma**2 * (v(2,:,:)**2+v(3,:,:)**2+v(4,:,:)**2) ) ! total E
        else
          q(1,:,:) = v(1,:,:)                                   ! rho
          q(2,:,:) = v(2,:,:)                                   ! u-velocity
          q(3,:,:) = v(3,:,:)                                   ! v-velocity
          q(4,:,:) = v(4,:,:)                                   ! w-velocity
          pinf = one/(gamma*Ma**2)
          p0   = (one + pt5*gamma1*Ma**2)**(gamma/gamma1)
!         q(4,:,:) = (p(:,:)/pinf - one)/(p0 - one)             ! Cp (comp)
!         q(4,:,:) = two*(p(:,:)-pinf)                          ! Cp
!         q(4,:,:) = v(5,:,:) + pt5 * Ma**2 * (gamma - one) * &
!                    ( v(2,:,:)**2 + v(3,:,:)**2 + v(4,:,:)**2 )! H_0
          q(5,:,:) = v(5,:,:)                                   ! temperature
        end if

!============================================================================!

!.... output debug files

        if (debug) then

!.... write out quantities along the ray from the nose

        i = 1
        write(*,"('Enter i ==> ',$)")
        read(*,*) i
        bn1 = n1(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
        bn2 = n2(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
        allocate( ynn(ny), uss(ny), unn(ny) )
        arc = zero
        do j = 1, ny
          ynn(j) = arc
          if (j.ne.ny) arc = arc + sqrt( (x(i,j+1)-x(i,j))**2 + &
            (y(i,j+1)-y(i,j))**2 )
          uss(j) = bn2 * v(2,i,j) - bn1 * v(3,i,j)
          unn(j) = bn1 * v(2,i,j) + bn2 * v(3,i,j)
        end do
        do j = 1, ny
          write(11,"(7(1pe16.9,1x))") ynn(j), v(1,i,j), uss(j), unn(j), &
                                      v(4,i,j), v(5,i,j), p(i,j)
        end do
        deallocate( ynn, uss, unn )

!.... write out quantities along the far-field boundary

        j = 1
        write(*,"('Enter j ==> ',$)")
        read(*,*) j

        do i = 1, nx
          m1l = n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
          m2l = n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
          un  = v(2,i,j) * m1l + v(3,i,j) * m2l 
          write(12,"(8(1pe16.9,1x))") s(i), v(1,i,j), v(2,i,j), &
                                      v(3,i,j), v(4,i,j), v(5,i,j), p(i,j)
        end do

!       j = 1
        do i = 1, nx
          write(17,"(8(1pe16.9,1x))") x(i,j), g2v(1,i,j), g2v(2,i,j), &
                                      g2v(3,i,j), g2v(4,i,j), g2v(5,i,j)
        end do


!.... write out quantities along the outflow boundnary

        i = nx
        do j = 1, ny
          write(13,"(7(1pe16.9,1x))") y(i,j), v(1,i,j), v(2,i,j), v(3,i,j), &
                                      v(4,i,j), v(5,i,j), p(i,j)
        end do

!.... write out quantities along the 1 node in from the outflow boundnary

        i = nx-1
        do j = 1, ny
          write(14,"(7(1pe16.9,1x))") y(i,j), v(1,i,j), v(2,i,j), v(3,i,j), &
                                      v(4,i,j), v(5,i,j), p(i,j)
        end do

!.... write out quantities along the 10 nodes in from the outflow boundnary

        i = nx-10
        do j = 1, ny
          write(15,"(7(1pe16.9,1x))") y(i,j), v(1,i,j), v(2,i,j), v(3,i,j), &
                                      v(4,i,j), v(5,i,j), p(i,j)
        end do

!.... write out quantities along the 20 nodes in from the outflow boundnary

        i = nx-20
        do j = 1, ny
          write(16,"(7(1pe16.9,1x))") y(i,j), v(1,i,j), v(2,i,j), v(3,i,j), &
                                      v(4,i,j), v(5,i,j), p(i,j)
        end do

        end if
        
!============================================================================!

!.... extract profiles at user supplied locations

        if (profile) then

        write(*,"('E R R O R:  profile extraction not updated yet...')")
        call exit(1)

        write(*,"('Enter the x-stations for profiles [min,max,inc] ==> ',$)")
        read(*,*) imin, imax, inc
        base  = 'profile'
        base1 = 'first'
        base2 = 'second'
        if (imin.ne.0) then
        kord = 5
        ns = ny
        allocate ( knot(ns+kord), bs(ns,ndof), ynn(ny), uss(ny), unn(ny) )
        open(11,file='loc.dat')
        do i = imin, imax, inc
          call makename(base,i,fname)
          call makename(base1,i,fname1)
          call makename(base2,i,fname2)
          iloc = index(fname,' ')
          write(*,110) s(i), x(i,1), fname(1:iloc)
 110      format('Saving profile at s, x_b = ',2(1pe20.13,1x),'in file: ',a)
          write(11,"(i4,1x,2(1pe13.6,1x))") i, s(i), x(i,1)

          bn1 = n1(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
          bn2 = n2(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
          
          arc = zero
          do j = 1, ny
            ynn(j) = arc
            if (j.ne.ny) arc = arc + sqrt( (x(i,j+1)-x(i,j))**2 + &
                                           (y(i,j+1)-y(i,j))**2 )
            uss(j) = bn2 * v(2,i,j) - bn1 * v(3,i,j)
            unn(j) = bn1 * v(2,i,j) + bn2 * v(3,i,j)
          end do
          
          call BSNAK( ny, ynn, kord, knot)
          call BSINT( ny, ynn, v(:,i,1), kord, knot, bs(:,1) )
          call BSINT( ny, ynn,      uss, kord, knot, bs(:,2) )
          call BSINT( ny, ynn,      unn, kord, knot, bs(:,3) )
          call BSINT( ny, ynn, v(:,i,4), kord, knot, bs(:,4) )
          call BSINT( ny, ynn, v(:,i,5), kord, knot, bs(:,5) )

!.... compute the boundary layer edge location (two ways)

!.... Maximum in U_c

!         do j = 1, ny-1
!           if ( uss(j+1) .lt. uss(j) ) goto 700        ! rough estimate
!         end do
! 700     edge = RTSAFE( fedge, ynn(j-1), ynn(j+2), 1.0e-8 )

!.... 0.999 point in W_c

          do j = 1, ny-1
            if ( v(4,i,j) .gt. 0.999*tan(alpha) ) goto 700  ! rough estimate
          end do
 700      edge = RTSAFE( fedge, ynn(j-1), ynn(j+2), 1.0e-8 )

          rhoe = BSDER( 0, edge, kord, knot, ny, bs(:,1) )
          ue   = BSDER( 0, edge, kord, knot, ny, bs(:,2) )
          ve   = BSDER( 0, edge, kord, knot, ny, bs(:,3) )
          we   = BSDER( 0, edge, kord, knot, ny, bs(:,4) )
          te   = BSDER( 0, edge, kord, knot, ny, bs(:,5) )
          pe   = rhoe * te / (gamma * Ma**2)
          write(*,"(8(1pe13.6,1x))") x(i,1), edge, rhoe, ue, ve, we, te
          alp = atan2( we, ue )
          umag = sqrt(ue**2+we**2)

          open(10,file=fname)
          open(12,file=fname1)
          open(13,file=fname2)
          do j = 1, ny
            write(10,"(9(1pe17.9e3,1x))") ynn(j), v(1,i,j), uss(j), unn(j), &
               v(4,i,j), v(5,i,j), p(i,j), &
               (cos(alp) * uss(j) + sin(alp) * v(4,i,j)), &
               (-sin(alp) * uss(j) + cos(alp) * v(4,i,j))
            write(12,"(9(1pe17.9e3,1x))") ynn(j),      &
               g1v(1,i,j)*bn1 + g2v(1,i,j)*bn2,      &
               bn2*(g1v(2,i,j)*bn1+g2v(2,i,j)*bn2) - &
               bn1*(g1v(3,i,j)*bn1+g2v(3,i,j)*bn2),  &
               bn1*(g1v(2,i,j)*bn1+g2v(2,i,j)*bn2) + &
               bn2*(g1v(3,i,j)*bn1+g2v(3,i,j)*bn2),  &
               g1v(4,i,j)*bn1 + g2v(4,i,j)*bn2,      &
               g1v(5,i,j)*bn1 + g2v(5,i,j)*bn2
            write(13,"(9(1pe17.9e3,1x))") ynn(j),       &
               (g11v(1,i,j)*bn1+g12v(1,i,j)*bn2)*bn1 + & 
               (g12v(1,i,j)*bn1+g22v(1,i,j)*bn2)*bn2,  &
               bn2*((g11v(2,i,j)*bn1+g12v(2,i,j)*bn2)*bn1  + &
                    (g12v(2,i,j)*bn1+g22v(2,i,j)*bn2)*bn2) - &
               bn1*((g11v(3,i,j)*bn1+g12v(3,i,j)*bn2)*bn1  + &
                    (g12v(3,i,j)*bn1+g22v(3,i,j)*bn2)*bn2),  &
               bn1*((g11v(2,i,j)*bn1+g12v(2,i,j)*bn2)*bn1  + &
                    (g12v(2,i,j)*bn1+g22v(2,i,j)*bn2)*bn2) + &
               bn2*((g11v(3,i,j)*bn1+g12v(3,i,j)*bn2)*bn1  + &
                    (g12v(3,i,j)*bn1+g22v(3,i,j)*bn2)*bn2),  &
               (g11v(4,i,j)*bn1+g12v(4,i,j)*bn2)*bn1 + & 
               (g12v(4,i,j)*bn1+g22v(4,i,j)*bn2)*bn2,  &
               (g11v(5,i,j)*bn1+g12v(5,i,j)*bn2)*bn1 + & 
               (g12v(5,i,j)*bn1+g22v(5,i,j)*bn2)*bn2
          end do
          close(10)
          close(12)
          close(13)
        end do
        close(11)
        deallocate ( knot, bs, ynn, uss, unn )
        end if
        end if

!============================================================================!

!.... write out a 2D plot3d file

        if (plot3d) then
        
        iloc = index(file,'.R.')
        if (iloc.eq.0) iloc = index(file,'.r.')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'q'//file(iloc+2:iend)
        else 
          iloc = index(file,'.res')
          if (iloc .ne. 0) then
            temp = file(1:iloc)//'dat'
          else
            temp = file(1:iend)//'.q'
          end if
        end if

        if ( nfile .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)")  &
                  file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( nfile .eq. 1 ) then
          file = temp
        else if ( nfile .eq. 2 ) then
          call getarg(ifile(2),temp)
          file = temp
        end if
 40     open(unit=10, file=file, form='unformatted', status='unknown', err=50)
        goto 60
 50     write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
        write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)") &
              file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 40
 60     continue
              
        write(10,err=1010) nx, ny, nz
        write(10,err=1010) Ma, Pr, Re, time
        write(10,err=1010) ((((q(idof,i,j), i = 1, nx), j = 1, ny),  &
                            k = 1, nz), idof = 1, ndof)
        close(10,err=1010)

        end if
        
!============================================================================!

!.... output a supplemental Plot3d file?

        if (type.ne.0) then

        if (type.eq.1) then
          q(1,:,:) = g11v(1,:,:)
          q(2,:,:) = g11v(2,:,:)
          q(3,:,:) = g11v(3,:,:)
          q(4,:,:) = g11v(4,:,:)
          q(5,:,:) = g11v(5,:,:)
        else if (type.eq.2) then
          q(1,:,:) = g22v(1,:,:)
          q(2,:,:) = g22v(2,:,:)
          q(3,:,:) = g22v(3,:,:)
          q(4,:,:) = g22v(4,:,:)
          q(5,:,:) = g22v(5,:,:)
        else if (type.eq.3) then
          q(1,:,:) = g12v(1,:,:)
          q(2,:,:) = g12v(2,:,:)
          q(3,:,:) = g12v(3,:,:)
          q(4,:,:) = g12v(4,:,:)
          q(5,:,:) = g12v(5,:,:)
        else if (type.eq.4) then
          q(1,:,:) = sqrt( v(5,:,:) ) / Ma                         ! c
          q(2,:,:) = v(2,:,:) / q(1,:,:)                           ! u-Mach
          q(3,:,:) = v(3,:,:) / q(1,:,:)                           ! v-Mach
          q(4,:,:) = sqrt(v(2,:,:)**2 + v(3,:,:)**2) / q(1,:,:)    ! Mach
          q(5,:,:) = log( v(1,:,:) * v(5,:,:) / &
                     (gamma * Ma**2) / v(1,:,:) ** gamma )         ! entropy
        else if (type.eq.5) then
!!$       q(1,:,:) = g2v(1,:,:)
!!$       q(2,:,:) = g1v(2,:,:)
!!$       q(3,:,:) = g1v(3,:,:)
!!$       q(4,:,:) = g1v(4,:,:)
!!$       q(5,:,:) = g1v(5,:,:)

          q(1,:,:) = g1v(3,:,:) - g2v(2,:,:)                       ! omega_z
          q(2,:,:) = g1v(2,:,:) + g2v(3,:,:)                       ! divergence
          q(3,:,:) = g1v(1,:,:) * v(5,:,:) + v(1,:,:) * g1v(5,:,:) ! p_,x
          q(4,:,:) = g2v(1,:,:) * v(5,:,:) + v(1,:,:) * g2v(5,:,:) ! p_,y
          q(5,:,:) = pt5 * (v(2,:,:)**2+v(3,:,:)**2+v(4,:,:)**2)   ! Kinetic E

!         q(5,:,:) = v(1,:,:) * ( one / gamma1 / gamma * v(5,:,:) + pt5 * &
!                    Ma**2 * (v(2,:,:)**2+v(3,:,:)**2+v(4,:,:)**2) ) ! total E

!         q(5,:,:) = v(1,:,:) * v(5,:,:) / (gamma * Ma**2) + &     ! Stagnation
!                    pt5 * v(1,:,:) * ( v(2,:,:)**2 + &            ! pressure
!                    v(3,:,:)**2 + v(4,:,:)**2 )                   ! (incomp)
        end if
        
!.... write out a supplemental plot3d file

        iloc = index(file,'.q.')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'d'//file(iloc+2:iend)
        else
          iloc = index(file,'.dat')
          if (iloc .ne. 0) then
            temp = file(1:iloc)//'d.dat'
          else
            temp = file(1:iend)//'.d'
          end if
        end if

        if ( nfile .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)")  &
                  file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( nfile .eq. 1 .or. nfile .eq. 2) then
          file = temp
        else if ( nfile .eq. 3 ) then
          call getarg(ifile(3),temp)
          file = temp
        end if
 400    open(unit=10, file=file, form='unformatted', status='unknown', err=500)
        goto 600
 500    write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
        write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)") &
              file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 400
 600    continue
              
        write(10,err=1010) nx, ny, nz
        write(10,err=1010) Ma, Pr, Re, time
        write(10,err=1010) ((((q(idof,i,j), i = 1, nx), j = 1, ny),  &
                            k = 1, nz), idof = 1, ndof)
        close(10,err=1010)

        end if

!============================================================================!
!     Wall quantities
!============================================================================!

!.... write out the wall quantities:
!....   1       x station
!....   2       arc length
!....   3       wall vorticity
!....   4       -Cp
!....   5       dp/ds
!....   6       du_s/dn
!....   7       dt/dn
!....   8       us
!....   9       un
!....   10      t

        open(unit=10, file='wall.dat', form='formatted', status='unknown')
        do i = 1, nx
          m1l = n1(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
          m2l = n2(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )

          pinf = one/(gamma*Ma**2)
          p0   = (one + pt5*gamma1*Ma**2)**(gamma/gamma1)
!         Cp   = (p(i,1)/pinf - one)/(p0 - one)               ! Cp (comp)
          Cp   = two*(p(i,1) - pinf)                          ! Cp

          write(10,11) x(i,1), s(i),                                    &
                       g1v(3,i,1) - g2v(2,i,1),                         &
                       Cp,                                              &
                       (m2l * (g1v(1,i,1) *   v(5,i,1)  +               &
                                 v(1,i,1) * g1v(5,i,1)) -               &
                        m1l * (g2v(1,i,1) *   v(5,i,1)  +               &
                                 v(1,i,1) * g2v(5,i,1))) /              &
                          (gamma * Ma**2),                              &
                       (g1v(2,i,1) * m1l + g2v(2,i,1) * m2l) * m2l -    &
                       (g1v(3,i,1) * m1l + g2v(3,i,1) * m2l) * m1l,     &
                       g1v(5,i,1) * m1l + g2v(5,i,1) * m2l,             &
                       v(2,i,1)*m2l-v(3,i,1)*m1l,                       &
                       v(2,i,1)*m1l+v(3,i,1)*m2l,                       &
                       v(5,i,1)
        end do
        close(10)

!=============================================================================
!
!.... write out potential flow data
!
!=============================================================================
        if (Re.eq.zero .or. inviscid) then
        
        allocate ( uss(nx), unn(nx) )

        Re = 0

!.... This is useful when trying to estimate thicknesses at AT-line

!       write(*,"('Enter Effective Re ==> ',$)")
!       read(*,*) Re

        if (Re .ne. 0) then
        
        kord = 5
        ns = nx
        allocate ( knot(ns+kord), bs(ns,3) )
        
        call BSNAK( ns, s, kord, knot)
        call BSINT( ns, s, x(:,1), kord, knot, bs(:,1) )
        call BSINT( ns, s, y(:,1), kord, knot, bs(:,2) )

!.... compute the length-scale at the nose

        bn1 = n1(2,1) / sqrt( n1(2,1)**2 + n2(2,1)**2 )
        bn2 = n2(2,1) / sqrt( n1(2,1)**2 + n2(2,1)**2 )
        uss(1) = bn2 * v(2,2,1) - bn1 * v(3,2,1)
        L2 = sqrt(s(2) / uss(1) / Re)
        
        bn1 = n1(3,1) / sqrt( n1(3,1)**2 + n2(3,1)**2 )
        bn2 = n2(3,1) / sqrt( n1(3,1)**2 + n2(3,1)**2 )
        uss(1) = bn2 * v(2,3,1) - bn1 * v(3,3,1)
        L3 = sqrt(s(3) / uss(1) / Re)
        
        L1 = L2 + (L3 - L2)/(s(3)-s(2)) * (s(1)-s(2))

!.... compute the velocities
        
        do i = 1, nx
          bn1 = n1(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
          bn2 = n2(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
          
          uss(i) = bn2 * v(2,i,1) - bn1 * v(3,i,1)
          unn(i) = bn1 * v(2,i,1) + bn2 * v(3,i,1)
        end do
        
        call BSINT( ns, s, uss, kord, knot, bs(:,3) )

        open(unit=10, file='pot.dat', status='unknown')
        open(unit=14, file='ibl3d.dat',status='unknown')
        open(unit=15, file='beta.dat',status='unknown')
        do i = 1, nx
          bn1 = n1(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
          bn2 = n2(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )

          cur = sqrt( BSDER( 2, s(i), kord, knot, ns, bs(:,1) )**2 + &
                      BSDER( 2, s(i), kord, knot, ns, bs(:,2) )**2 )

          c1   = sqrt( v(5,i,1) ) / Ma
          ue   = v(2,i,1)
          ve   = v(3,i,1)
          we   = v(4,i,1)
          te   = v(5,i,1)
          U1   = sqrt( uss(i)**2 + v(4,i,1)**2 )
          Ma1  = U1 / c1 
          pinf = one/(gamma*Ma**2)
          p0   = (one + pt5*gamma1*Ma**2)**(gamma/gamma1)
!         Cp   = (p(i,1)/pinf - one)/(p0 - one)               ! Cp (comp)
          Cp   = two*( p(i,1) - pinf )
          if (alpha.ne.zero) then
            phip = (atan2( v(4,i,1), uss(i) ) - alpha) * 180.0 / pi
          else
            phip = zero
          end if
          Re1  = Re * U1 * s(i)
          if ( i .ne. 1 ) then
            L  = sqrt( s(i) / uss(i) / Re )
          else
            L  = L1
          end if
          ReL = Re * U1 * L

!.... compute the Hartree pressure gradient parameter
!....   Beta_H = s / u (\partial u / \partial s)

          if (i .eq. 1) then
            beta2 = s(2) * (uss(3)-uss(1)) / (s(3)-s(1)) / uss(2)
            beta3 = s(3) * (uss(4)-uss(2)) / (s(4)-s(2)) / uss(3)
            betah = beta2 + (beta3-beta2)/(s(3)-s(2)) * (s(1)-s(2))
            dusds = (uss(i+1)-uss(i)) / (s(i+1)-s(i))
          else if (i .eq. nx) then
            betah = s(i) * (uss(i)-uss(i-1)) / (s(i)-s(i-1)) / uss(i)
            dusds = (uss(i)-uss(i-1)) / (s(i)-s(i-1))
          else
            betah = s(i) * (uss(i+1)-uss(i-1)) / (s(i+1)-s(i-1)) / uss(i)
            dusds = (uss(i+1)-uss(i-1)) / (s(i+1)-s(i-1))
          end if

!.... more accurate Bspline method 

          if (i.eq.1) then
            betah = (s(i)+1.0e-6)*BSDER(1,s(i)+1.0e-6,kord,knot,ns,bs(:,3)) / &
                    BSDER(0,s(i)+1.0e-6,kord,knot,ns,bs(:,3))
          else
            betah = s(i) * BSDER(1,s(i),kord,knot,ns,bs(:,3)) / uss(i)
          end if

!.... pot.dat

          write(10,12) s(i), uss(i), unn(i), v(4,i,1), c1, &
                       Ma1, Cp, betah, phip, Re1, &
                       ReL, L

!.... ibl3d.dat for BL3D solver (Cartesian coordinates)

          write(14,12) x(i,1), zero, y(i,1), ue, we, ve, te, Cp

!.... beta.dat for NSBL5 solver

          write(15,12) s(i), uss(i), dusds, betah, cur, zero, x(i,1), &
                       y(i,1), bn1, bn2

        end do
        close(10); close(14); close(15)
        deallocate ( knot, bs )

        end if

        deallocate ( uss, unn )

        stop
        end if          ! Re = 0

!=============================================================================
!
!.... compute the boundary layer thickness parameters
!
!=============================================================================
        if (thick) then
        
        write(*,"('E R R O R:  thickness computation not updated yet...')")
        call exit(1)
        
!.... Setup the B-spline interpolation

        write(*,"(/,'Boundary layer thickness calculation . . .')")
        write(*,"(/,'W A R N I N G:  only valid on a body-fitted mesh!')")

        kord = 5
!       write(*,"('Enter korder (0 = terminate) ==> ',$)")
!       read(*,*) kord
!       if (kord.eq.0) call exit(0)

        ns = ny
        allocate ( knot(ny+kord), bs(ny,ndof), ynn(ny), &
                   uss(ny), unn(ny), wss(ny), uel(nx) )

!.... compute the boundary layer parameters

        open(10,file='delta.dat',status='unknown')
        open(12,file='delta2.dat',status='unknown')
        open(11,file='edge.dat',status='unknown')
        open(13,file='stat.dat',status='unknown')
        open(14,file='bl3d.dat',status='unknown')

        do i = 1, nx
          delta = zero
          theta = zero

          bn1 = n1(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )
          bn2 = n2(i,1) / sqrt( n1(i,1)**2 + n2(i,1)**2 )

          arc = zero
          do j = 1, ny
            ynn(j) = arc
            if (j.ne.ny) arc = arc + sqrt( (x(i,j+i)-x(i,j))**2 + &
                                           (y(i,j+i)-y(i,j))**2 )
            uss(j) = bn2 * v(2,i,j) - bn1 * v(3,i,j)
            unn(j) = bn1 * v(2,i,j) + bn2 * v(3,i,j)
          end do
          
          call BSNAK( ny, ynn, kord, knot)
          call BSINT( ny, ynn, v(:,i,1), kord, knot, bs(:,1) )
          call BSINT( ny, ynn,      uss, kord, knot, bs(:,2) )
          call BSINT( ny, ynn,      unn, kord, knot, bs(:,3) )
          call BSINT( ny, ynn, v(:,i,4), kord, knot, bs(:,4) )
          call BSINT( ny, ynn, v(:,i,5), kord, knot, bs(:,5) )

!.... compute the boundary layer edge location (three ways)

!.... Maximum in U_c

!         do j = 1, ny-1
!           if ( uss(j+1) .lt. uss(j) ) goto 70         ! rough estimate
!         end do
! 70      edge = RTSAFE( fedge, ynn(j-1), ynn(j+2), 1.0e-8 )

!.... 0.999 point in W_c

          do j = 1, ny-1
            if ( v(4,i,j) .gt. 0.999*tan(alpha) ) goto 70  ! rough estimate
          end do
  70      edge = RTSAFE( fedge, ynn(j-1), ynn(j+2), 1.0e-8 )

!.... 0.999 point in T

!         do j = 1, ny-1
!           if ( v(5,i,j) .lt. 0.999 ) goto 70  ! rough estimate
!         end do
! 70      edge = RTSAFE( fedge, ynn(j-1), ynn(j+2), 1.0e-8 )

          rhoe = BSDER( 0, edge, kord, knot, ny, bs(:,1) )
          ue   = BSDER( 0, edge, kord, knot, ny, bs(:,2) )
          ve   = BSDER( 0, edge, kord, knot, ny, bs(:,3) )
          we   = BSDER( 0, edge, kord, knot, ny, bs(:,4) )
          te   = BSDER( 0, edge, kord, knot, ny, bs(:,5) )
          uel(i) = ue
          pe   = rhoe * te / (gamma * Ma**2)
          pinf = one/(gamma*Ma**2)
          p0   = (one + pt5*gamma1*Ma**2)**(gamma/gamma1)
!         Cp   = (pe/pinf - one)/(p0 - one)               ! Cp (comp)
          Cp   = two*( pe - pinf )
          alp  = atan2( we, ue )
          phie = (alp - alpha) * 180.0 / pi

          write(*,"(9(1pe13.6,1x))")  s(i),edge,rhoe,ue,ve,we,te,pe,phie
          write(11,"(9(1pe13.6,1x))") s(i),edge,rhoe,ue,ve,we,te,pe,phie

!....     output a file for bl3d (Cartesian coordinates)

          write(14,12) x(i,1), zero, y(i,1), bn2*ue+bn1*ve, we, &
                       -bn1*ue+bn2*ve, te, Cp

          if (abs(s(i)).le.1.0e-8) goto 555

!.... switch to streamline coordinates, normalize by U1, and interpolate
!.... should U1 have un in it?

          U1 = sqrt(ue**2 + we**2)
          ue = one
          we = zero
!         write(*,"(9(1pe13.6,1x))") s(i), U1, ue, we, phie, cos(alp) * U1
          do j = 1, ny
            us =  cos(alp) * uss(j) + sin(alp) * v(4,i,j)
            ws = -sin(alp) * uss(j) + cos(alp) * v(4,i,j)
            uss(j) = us / U1
            wss(j) = ws / U1
          end do
          call BSINT( ny, ynn, uss, kord, knot, bs(:,2) )
          call BSINT( ny, ynn, wss, kord, knot, bs(:,4) )

!.... Integrate to get the boundary layer integral thicknesses

          call QDAG(  dfun, zero, edge, 1.0e-8, 1.0e-8, 2,  delta, errest)
          call QDAG( thfun, zero, edge, 1.0e-8, 1.0e-8, 2,  theta, errest)
!         call QDAG(thfunz, zero, edge, 1.0e-8, 1.0e-8, 2, thetaz, errest)

!.... Boundary layer statistics

          H   = theta / delta
          L   = sqrt( s(i) / (cos(alp)* U1) / Re )
          ReL = Re * U1 * L
          c1  = sqrt( Te ) / Ma
          Ma1 = U1 / c1 
          pinf = one/(gamma*Ma**2)
          p0   = (one + pt5*gamma1*Ma**2)**(gamma/gamma1)
!         Cp   = (pe/pinf - one)/(p0 - one)               ! Cp (comp)
          Cp   = two*( pe - pinf )

!.... determine the location of maximum crossflow

          do j = 1, ny-1
            if ( abs(wss(j+1)) .lt. abs(wss(j)) ) goto 75
          end do
 75       if (j.eq.ny) then
            write(*,*) 'Could not find maximum crossflow at i = ', i
            call exit(1)
          end if
          jmax = j
          wmloc = RTSAFE( fwmax, ynn(j-1), ynn(j+2), 1.0e-8 )
          wmax  = BSDER( 0, wmloc, kord, knot, ny, bs(:,4) )

!.... find the inflection point in the crossflow profile

          do j = jmax, ny-1
            w1 = BSDER( 2,   ynn(j), kord, knot, ny, bs(:,4) )
            w2 = BSDER( 2, ynn(j+1), kord, knot, ny, bs(:,4) )
            if ( (w1.gt.zero .and. w2.lt.zero) .or. &
                 (w2.gt.zero .and. w1.lt.zero) ) goto 80
          end do
 80       if (j.eq.ny) then
            write(*,*) 'Could not find crossflow inflection point at i = ', i
            call exit(1)
          end if
          wiloc = RTSAFE( fwinf, ynn(j-1), ynn(j+2), 1.0e-8 )
          uinf  = BSDER( 0, wiloc, kord, knot, ny, bs(:,2) )
          winf  = BSDER( 0, wiloc, kord, knot, ny, bs(:,4) )
          epsi  = atan2( winf, uinf ) * 180.0 / pi
          
!.... output the results

          rk = (two * x(i,1) + one)**(-onept5) * delta
          
!.... delta.dat

          write(10,"(20(1pe13.6,1x))") s(i), delta, theta, H, wmax, &
                                       ReL, L, Ma1, c1, Cp, &
                                       delta/L, theta/L, wmloc, wiloc, uinf, &
                                       winf, epsi, edge*wmax*Re, rk

!.... delta2.dat

          write(12,"(20(1pe13.6,1x))") s(i), delta, theta, wmax, U1, &
                                       Te, phie, Ma1, edge*wmax*Re

!.... stat.dat

          write(13,"(20(1pe13.6,1x))") s(i), delta, alp
          
!.... the last item is Poll's crossflow Reynolds number

!.... output a profile

          if(.false.) then
            do j = 1, ny
              rho = BSDER( 0, ynn(j), kord, knot, ny, bs(:,1) )
              us  = BSDER( 0, ynn(j), kord, knot, ny, bs(:,2) )
              un  = BSDER( 0, ynn(j), kord, knot, ny, bs(:,3) )
              u3  = BSDER( 0, ynn(j), kord, knot, ny, bs(:,4) )
              t   = BSDER( 0, ynn(j), kord, knot, ny, bs(:,5) )
              write(*,"(8(1pe13.6,1x))") ynn(j), rho, us, un, u3, t
            end do
          end if
555       continue        
        end do
        close(10); close(11); close(12); close(13); close(14)

!.... compute the Hartree pressure gradient parameter
!....
!....   m = (s / u) (\partial u / \partial s)
!....   Beta_H = 2 m / (m + 1)

        open(10,file='betah.dat')
        do i = 1, nx
          if (i .eq. 1) then
            beta2 = s(2) * (uel(3)-uel(1)) / (s(3)-s(1)) / uel(2)
            beta3 = s(3) * (uel(4)-uel(2)) / (s(4)-s(2)) / uel(3)
            betah = beta2 + (beta3-beta2)/(s(3)-s(2)) * (s(1)-s(2))
          else if (i .eq. nx) then
            betah = s(i) * (uel(i)-uel(i-1)) / (s(i)-s(i-1)) / uel(i)
          else
            betah = s(i) * (uel(i+1)-uel(i-1)) / (s(i+1)-s(i-1)) / uel(i)
          end if
          betah = two * betah / (betah + one)
          write(10,"(20(1pe13.6,1x))") s(i), betah
        end do
        close(10)

        deallocate (knot, bs, ynn, uss, unn, wss, uel)

        end if          ! thick

!=============================================================================
!
!.... write out a 3D grid and data file
!
!=============================================================================
        if (swept) then

!.... always use conservation variables for swept cases

          q(1,:,:) = v(1,:,:)
          q(2,:,:) = v(1,:,:) * v(2,:,:) * Ma
          q(3,:,:) = v(1,:,:) * v(3,:,:) * Ma
          q(4,:,:) = v(1,:,:) * v(4,:,:) * Ma
          q(5,:,:) = v(1,:,:) * ( one / gamma1 / gamma * v(5,:,:) + pt5 * &
                     Ma**2 * (v(2,:,:)**2+v(3,:,:)**2+v(4,:,:)**2) )  

          write(*,"('Enter nz, Lz ==> ',$)")
          read(*,*) nz, Lz

          vinf = sqrt( one + tan(alpha)**2 )

          write(10,err=1010) nx, ny, nz
          write(10,err=1010) Ma * vinf, alpha*180.0/pi, Re * vinf, time
          write(10,err=1010) ((((q(idof,i,j), i = 1, nx), j = 1, ny),  &
                             k = 1, nz), idof = 1, ndof)
          close(10,err=1010)

          dz = Lz/float(nz-1)
          allocate( z(nz) )
          do k = 1, nz
            z(k) = float(k-1)*dz
          end do
          open(unit=10,file='grid3d.dat',form='unformatted',status='unknown')
          write(10) nx, ny, nz
          write(10) (((x(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                    (((y(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                    (((  z(k), i = 1, nx), j = 1, ny), k = 1, nz)
          close(10)
        end if

        call exit(0)    

 11     format(20(1pe21.14,1x))
 12     format(20(1pe13.6,1x))

 1010   write(*,"('>> Error writing PLOT3D file')")
        call exit(1)

        end
!=============================================================================!
        function dfun(x)
!
!       Function to compute the dispacement thickness
!
!=============================================================================!
        use bspline

        real :: x, dfun, rho, u
!=============================================================================!
        rho = BSDER( 0, x, kord, knot, ns, bs(:,1) )
        u   = BSDER( 0, x, kord, knot, ns, bs(:,2) )

        dfun = 1.0 - rho * u / (rhoe * ue)
        
        return
        end 
!=============================================================================!
        function thfun(x)
!
!       Function to compute the momentum thickness
!
!=============================================================================!
        use bspline

        real :: x, thfun, rho, u
!=============================================================================!
        rho = BSDER( 0, x, kord, knot, ns, bs(:,1) )
        u   = BSDER( 0, x, kord, knot, ns, bs(:,2) )

        thfun = rho * u / (rhoe * ue) * ( 1.0 - u / ue )
        
        return
        end 
!=============================================================================!
        function thfunz(x)
!
!       Function to compute the momentum thickness in z
!
!=============================================================================!
        use bspline

        real :: x, thfun, rho, w
!=============================================================================!
        rho = BSDER( 0, x, kord, knot, ns, bs(:,1) )
        w   = BSDER( 0, x, kord, knot, ns, bs(:,4) )

        thfunz = rho * w / (rhoe * we) * ( 1.0 - w / we )
        
        return
        end 
!=============================================================================!
        subroutine fedge( x, g, d )
!
!       Function to find the boundary layer edge
!
!=============================================================================!
        use bspline

        real :: x, f, g, d
!=============================================================================!

!.... maximum in U_c

!       f = BSDER( 0, x, kord, knot, ns, bs(:,2) )
!       g = BSDER( 1, x, kord, knot, ns, bs(:,2) )
!       d = BSDER( 2, x, kord, knot, ns, bs(:,2) )

!.... Delta 0.999 in W_c

        g = BSDER( 0, x, kord, knot, ns, bs(:,4) ) - 0.999*tan(alpha)
        d = BSDER( 1, x, kord, knot, ns, bs(:,4) )

!.... Delta 0.999 in T
        
!       g = BSDER( 0, x, kord, knot, ns, bs(:,5) ) - 0.999
!       d = BSDER( 1, x, kord, knot, ns, bs(:,5) )

        return
        end 
!=============================================================================!
        subroutine fwmax( x, g, d )
!
!       Function to find the maximum crossflow velocity
!
!=============================================================================!
        use bspline

        real :: x, f, g, d
!=============================================================================!
        f = BSDER( 0, x, kord, knot, ns, bs(:,4) )
        g = BSDER( 1, x, kord, knot, ns, bs(:,4) )
        d = BSDER( 2, x, kord, knot, ns, bs(:,4) )

        return
        end 
!=============================================================================!
        subroutine fwinf( x, g, d )
!
!       Function to find the inflection in the crossflow velocity
!
!=============================================================================!
        use bspline

        real :: x, f, g, d
!=============================================================================!
        f = BSDER( 1, x, kord, knot, ns, bs(:,4) )
        g = BSDER( 2, x, kord, knot, ns, bs(:,4) )
        d = BSDER( 3, x, kord, knot, ns, bs(:,4) )

        return
        end 
!=============================================================================!
    subroutine makename(base,iver,fname)
!
!.... put a version number on the filename
!
!=============================================================================!
      character(80) base, fname

      length = index(base,' ')
      fname = base
      if (iver .lt. 10) then
        write(fname(length:80),"('.',i1)") iver
      else if (iver .lt. 100) then
        write(fname(length:80),"('.',i2)") iver
      else
        write(fname(length:80),"('.',i3)") iver
      end if

      return
      end

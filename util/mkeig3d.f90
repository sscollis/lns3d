!=============================================================================!
      program convert
!
!     Projects an eigenfunction file to a new grid
!
!     Revised: 5-9-96
!=============================================================================!
      use const
      implicit none

      real, allocatable    :: yc(:), ql(:)
      complex, allocatable :: q(:,:), q2(:,:), v(:,:,:)
      real, allocatable    :: x(:,:), y(:,:), s(:)

      complex :: u1, u2

!.... metrics

      real, allocatable :: m1(:,:),  m2(:,:),  n1(:,:),  n2(:,:),       &
                           m11(:,:), m12(:,:), m22(:,:),                &
                           n11(:,:), n12(:,:), n22(:,:)

      real, allocatable :: bn1(:), bn2(:)
        
      complex :: alpha, amp
      real    :: tmp, Lmap, eta, th, alphar, alphai, ampr, ampi, x0, r
      real    :: Ma, Re, Pr, gamma=1.4, cv=716.5
      integer :: i, j, k, m, idof
      integer :: nx, ny, nz, n, ndof=5
      character(80) :: efun

      logical :: switch_ij = .true.

!=============================================================================!
      write(*,"('Interpolates an eigenfunction from a Chebyshev grid')")
      write(*,"('  to an arbitrary grid in file grid.dat')")

!.... determine the number of points in the eigenfunction

      write(*,"('Enter Eigenfunction file name ==> ',$)")
      read(*,*) efun
      open (unit=10,file=efun,status='old')
      i = 1
 10   continue
         read(10,*,end=20) tmp
         i = i + 1
         goto 10
 20   continue
      close (10)
      n = i - 1

!.... read the eigenfunction for real

      allocate( yc(n), q(n,ndof), ql(2*ndof) )
      open (unit=10,file=efun,status='old')
      do i = 1, n
        read(10,*) yc(i), (ql(j),j=1,2*ndof)
        q(i,1) = ql(1) + im * ql(2)
        q(i,2) = ql(3) + im * ql(4)
        q(i,3) = ql(5) + im * ql(6)
        q(i,4) = ql(7) + im * ql(8)
        q(i,5) = ql(9) + im * ql(10)
      end do

!.... compute the complex Chebyshev transform

      do j = 1, ndof
         call CHEBYSHEV( q(1,j), n-1, 1 )
      end do

!.... read the grid file

      open(unit=10,file='grid.dat',form='unformatted',status='old')
      read(10) nx, ny, nz
      allocate( x(ny,nx), y(ny,nx), s(nx) )
      read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
               (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
               (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
      close(10)

!.... read in the body discription file

        s = x(1,:)
        open(unit=10,file='body.dat',form='formatted',status='old',err=100)
        do i = 1, nx
           read(10,*) tmp, s(i)
        end do
        close(10)
        goto 101
 100    continue

!.... if no body file is found, assume that its a parabolic cylinder

        write(*,*) 'WARNING:  assumming a parabolic cylinder'
        s(:) = Sqrt(x(1,:) + two*x(1,:)**2)/Sqrt(two) + &
               Log(one + four*x(1,:) + two**(onept5)*Sqrt(x(1,:) + &
               two*x(1,:)**2)) / four
 101    continue

!.... allocate storage for metrics

        allocate (m1(ny,nx),  m2(ny,nx),  n1(ny,nx),  n2(ny,nx), &
                  m11(ny,nx), m12(ny,nx), m22(ny,nx),            &
                  n11(ny,nx), n12(ny,nx), n22(ny,nx) )

!.... read in the metric file

        open (unit=10,file='metric.dat',form='unformatted', status='old')
        if (switch_ij) then
          read(10) (( m1(j,i), i=1,nx),j=1,ny), &
                   (( m2(j,i), i=1,nx),j=1,ny), &
                   (( n1(j,i), i=1,nx),j=1,ny), &
                   (( n2(j,i), i=1,nx),j=1,ny), &
                   ((m11(j,i), i=1,nx),j=1,ny), &
                   ((m12(j,i), i=1,nx),j=1,ny), &
                   ((m22(j,i), i=1,nx),j=1,ny), &
                   ((n11(j,i), i=1,nx),j=1,ny), &
                   ((n12(j,i), i=1,nx),j=1,ny), &
                   ((n22(j,i), i=1,nx),j=1,ny)
        else
          read(10) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
        endif
        close(10)
        
!.... compute the wall-normals

        allocate( bn1(nx), bn2(nx) )
        do i = 1, nx
          bn1(i) = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          bn2(i) = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )    
        end do

!.... write out the Chebyshev interpolant onto the new mesh

      write (*,"('Enter Yi ==> ',$)")
      read (*,*) Lmap                   ! called Lmap herein

      i = 1
      allocate( q2(ny,ndof) )
      q2 = zero
      do j = 1, ny
        r = sqrt( (x(j,i)-x(1,i))**2 + (y(j,i)-y(1,i))**2 )
        do idof = 1, ndof
          eta = (r - Lmap) / (r + Lmap)
!         eta = (y(j,1) - Lmap) / (y(j,1) + Lmap)
          th  = acos(eta)
          do m = 0, n-1
            q2(j,idof) = q2(j,idof) + q(m+1,idof) * COS(float(m)*th)
          end do
        end do
        u1 =  bn2(i) * q2(j,2) + bn1(i) * q2(j,3)
        u2 = -bn1(i) * q2(j,2) + bn2(i) * q2(j,3)
        q2(j,2) = u1
        q2(j,3) = u2
      end do

!.... make the inflow profile

      open(unit=33,file='inflow.dat',status='unknown')
      i = 1
      r = 0
      do j = 1, ny
        write(33,50) r, real(q2(j,1)), aimag(q2(j,1)), &
                        real(q2(j,2)), aimag(q2(j,2)), &
                        real(q2(j,3)), aimag(q2(j,3)), &
                        real(q2(j,4)), aimag(q2(j,4)), &
                        real(q2(j,5)), aimag(q2(j,5))
50      format(11(1pe20.13,1x))
        if (j.ne.ny) r = r + sqrt((x(j+1,i)-x(j,i))**2 + (y(j+1,i)-y(j,i))**2)
      end do
      
!.... write out on a mesh

      write(*,"('Enter Alpha_r, Alpha_i ==> ',$)") 
      read(*,*) alphar, alphai
      alpha = cmplx(alphar,alphai)

      write(*,"('Enter amp_r, amp_i ==> ',$)")
      read(*,*) ampr,ampi
      amp = cmplx(ampr,ampi)
      
      write(*,"('Enter s0 ==> ',$)")
      read(*,*) x0
      
      allocate( v(ny,nx,ndof) )

!      do idof = 1, ndof
!       do i = 1, nx
!         do j = 1, ny
!           v(j,i,idof) = amp * q2(j,idof) * exp( im * alpha * (x(j,i)-x0) )
!         end do
!       end do
!      end do

      do i = 1, nx
        q2 = zero
        r  = zero
        do j = 1, ny
          do idof = 1, ndof
            eta = (r - Lmap) / (r + Lmap)
            th  = acos(eta)
            do m = 0, n-1
              q2(j,idof) = q2(j,idof) + q(m+1,idof) * COS(float(m)*th)
            end do
          end do
          u1 =  bn2(i) * q2(j,2) + bn1(i) * q2(j,3)
          u2 = -bn1(i) * q2(j,2) + bn2(i) * q2(j,3)
          q2(j,2) = u1
          q2(j,3) = u2
          do idof = 1, ndof
            v(j,i,idof) = amp * q2(j,idof) * exp( im * alpha * (s(i) - x0) )
          end do
          if (j.ne.ny) r = r + sqrt((x(j+1,i)-x(j,i))**2+(y(j+1,i)-y(j,i))**2)
        end do
      end do

!.... write the restart file

      write(*,"('Enter Re, Ma, Pr ==> ',$)")
      read(*,*) Re, Ma, Pr
      open(unit=10,file='output.R.0',form='unformatted',status='unknown')
      write(10) 0, zero, nx, ny, nz, ndof, &
                Re, Ma, Pr, gamma, cv
      if (switch_ij) then
        write(10) (((v(j,i,idof), idof=1,ndof), i=1,nx), j=1,ny)
      else
        write(10) v
      endif
      close(10)

      stop
      end
            
!=============================================================================!
      subroutine CHEBYSHEV (Y, N, ISIGN)
!=============================================================================!
!
!     Take the Chebyshev transform and normalize
!
!=============================================================================!
      complex y(0:n)

      if (isign .ne. 1) then
        if (isign .ne. -1) then
          write (*,*) 'ERROR:  Invalid Isign in CHEBYSHEV'
          stop
        end if
      end if

      call COSFT3 (y,n,isign)

      return
      end

!=============================================================================!
      SUBROUTINE COSFT3 (YY,N,ISIGN)
!=============================================================================!
!
!     Do a brute force transform.  ISIGN = 1 does a forward transform.
!     ISIGN = -1 a backward transform.
!
!=============================================================================!
      integer N, ISIGN

      complex YY(0:N), TT(0:N)

      PI = DACOS(-1.0D0)

      DO I = 0, N
        TT(I) = YY(I)
      END DO

      if (isign .eq. 1) then
        DO I = 0, N
          YY(I) = TT(0)/2.D0 + TT(N)*COS(FLOAT(I)*PI)/2.D0
          DO M = 1, N-1
            YY(I) = YY(I) + TT(M)*COS(FLOAT(M)*FLOAT(I)*PI/FLOAT(N))
          END DO
          YY(I) = YY(I)*2./FLOAT(N)
        end do
        YY(0) = YY(0)/2.D0
        YY(N) = YY(N)/2.D0
      else 
        DO I = 0, N
          YY(I) = 0.0D0
          DO M = 0, N
            YY(I) = YY(I) + TT(M)*COS(FLOAT(M)*FLOAT(I)*PI/FLOAT(N))
          END DO
        END DO
      end if
      
      RETURN
      END

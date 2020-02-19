program ctest

! Test of the implicit solvers.
! Both have achieved solid performance on the C90
!
! cpenta1bc:  496 MFLOPS with 128 vector length
!
! cpenta2bc:  456 MFLOPS with 128 vector length
!
! WARNING:  Performance can drop to 55 MFLOPS in cpenta2bc if there are
!           bank conflicts.

      implicit none

      integer :: nx = 256, ny = 127, ndof = 5
      complex, allocatable :: mat(:,:,:,:,:), rhs(:,:,:)
      complex, allocatable :: mult(:), fact(:)

      integer :: i, j, idof, jdof, iter, niter = 100, ilhs

      write(*,"('Enter (nx, ny, niter) ==> ',$)")
      read(*,*) nx, ny, niter

      write(*,"('LHS (1 or 2) ==> ',$)")
      read(*,*) ilhs

      allocate( mat(ny,nx,ndof,ndof,5), rhs(ny,nx,ndof) )
      allocate( mult(max(ny,nx)), fact(max(ny,nx)) )

      do i = 1, nx
        do j = 1, ny
          mat(j,i,1,1,3) = 1.0
          mat(j,i,2,2,3) = 1.0
          mat(j,i,3,3,3) = 1.0
          mat(j,i,4,4,3) = 1.0
          mat(j,i,5,5,3) = 1.0

          rhs(j,i,1) = 1.0
          rhs(j,i,2) = 1.0
          rhs(j,i,3) = 1.0
          rhs(j,i,4) = 1.0
          rhs(j,i,5) = 1.0
        end do
      end do

      if (ilhs.eq.1) then
        call cpenta1bc( ny, nx, ndof, mat, rhs, mult, fact, 0)
        do iter = 1, niter
          call cpenta1bc( ny, nx, ndof, mat, rhs, mult, fact, 1)
        end do
      else
        call cpenta2bc( ny, nx, ndof, mat, rhs, mult, fact, 0)
        do iter = 1, niter
          call cpenta2bc( ny, nx, ndof, mat, rhs, mult, fact, 1)
        end do
      end if

      stop
end program ctest

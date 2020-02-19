program test

! Test of the implicit solvers.
! Both have achieved solid performance on the C90
!
! penta1bc:  419 MFLOPS (256,255)
!
! penta2bc:  360 MFLOPS (256,255)
!
! WARNING:  Performance can drop to 55 MFLOPS in penta2bc if there are
!           bank conflicts.

      implicit none

      integer :: nx = 256, ny = 127, ndof = 5
      real, allocatable :: mat(:,:,:,:,:), rhs(:,:,:)
      real, allocatable :: mult(:), fact(:)

      integer :: i, j, k, idof, jdof, iter, niter = 100, ilhs

      write(*,"('Enter (nx, ny, niter) ==> ',$)")
      read(*,*) nx, ny, niter

      write(*,"('LHS (1 or 2) ==> ',$)")
      read(*,*) ilhs

      allocate( mat(ny,nx,ndof,ndof,5), rhs(ny,nx,ndof) )
      allocate( mult(max(ny,nx)), fact(max(ny,nx)) )

      mat(:,:,1,1,3) = 1.0
      mat(:,:,2,2,3) = 1.0
      mat(:,:,3,3,3) = 1.0
      mat(:,:,4,4,3) = 1.0
      mat(:,:,5,5,3) = 1.0

      rhs(:,:,1) = 1.0
      rhs(:,:,2) = 1.0
      rhs(:,:,3) = 1.0
      rhs(:,:,4) = 1.0
      rhs(:,:,5) = 1.0
      
      if (ilhs.eq.1) then
        do iter = 1, niter
          call penta1bc( ny, nx, ndof, mat, rhs, mult, fact)
        end do
      else
        do iter = 1, niter
          call penta2bc( ny, nx, ndof, mat, rhs, mult, fact)
        end do
      end if

      stop
end program test

!=============================================================================!
        subroutine eigdrv3d()
!  
!  Driver for eigenvalue problem (LAPack)
!
!=============================================================================!
        use global
        use local2d
        use local3d
        implicit none
!=============================================================================!
        complex :: v(ndof,nx,ny)
        complex :: matx(5,ndof,ndof,nx,ny)
        complex :: mat(ndof,nx,ndof,nx)
        real    :: rmat(ndof,nx,ndof,nx)
        real    :: dtl(nx,ny)
        
!.... LAPack variables

        integer :: lwork, info
        real, allocatable :: rwork(:)
        complex, allocatable :: work(:)
        complex :: eval(ndof*nx), evec(ndof,nx,ndof*nx)
        real :: reval(ndof*nx), ieval(ndof*nx), revec(ndof,nx,ndof*nx)
        
        integer :: i, j, idof, jdof
        character*80 :: name, code='EigDrv$'
!=============================================================================!

        alfa = one
        dtl = one

        v = zero
        mat = zero
        matx = zero

        call lhs1f3D( matx, Ah, Dh, Dhi, Vh11, ABhi, dtl, v )

!.... Remove identity matrix

        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              matx(3,idof,idof,i,j) = matx(3,idof,idof,i,j) - one
            end do
          end do
        end do

!.... Convert to LAPack general matrix storage

        j = 1

        i = 1
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,nx-2) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,nx-1) = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  )  = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,i+1)  = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,i+2)  = matx(5,idof,jdof,i,j)
          end do
        end do
        
        i = 2
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,nx-1) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,i-1)  = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  )  = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,i+1)  = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,i+2)  = matx(5,idof,jdof,i,j)
          end do
        end do

        do i = 3, nx-2
          do idof = 1, ndof
            do jdof = 1, ndof
              mat(idof,i,jdof,i-2) = matx(1,idof,jdof,i,j)
              mat(idof,i,jdof,i-1) = matx(2,idof,jdof,i,j)
              mat(idof,i,jdof,i  ) = matx(3,idof,jdof,i,j)
              mat(idof,i,jdof,i+1) = matx(4,idof,jdof,i,j)
              mat(idof,i,jdof,i+2) = matx(5,idof,jdof,i,j)
            end do
          end do
        end do

        i = nx-2
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,i-2) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,i-1) = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  ) = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,i+1) = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,1)   = matx(5,idof,jdof,i,j)
          end do
        end do

        i = nx-1
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,i-2) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,i-1) = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  ) = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,1)   = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,2)   = matx(5,idof,jdof,i,j)
          end do
        end do
 
        i = nx
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,i  ) = one
          end do
        end do

        if (.false.) then

        write(*,*) 'Starting eigenvalue solver'
        lwork = max(1,3*(ndof*nx))
        allocate( work(lwork), rwork(2*(ndof*nx)) )
        call ZGEEV ( 'N', 'V', ndof*(nx-1), mat, ndof*nx, eval, &
                     evec, ndof*nx, evec, ndof*nx,          &
                     work, lwork, rwork, info )
        write(*,*) 'ZGEEV info = ',info
        deallocate( work, rwork )
        write(*,*) 'Finished eigenvalue solver'

        else

        rmat = real(mat)
        write(*,*) 'Starting eigenvalue solver'
        lwork = max(1,4*(ndof*nx))
        allocate( rwork(lwork) )
        call DGEEV ( 'N', 'V', ndof*(nx-1), rmat, ndof*nx, reval, ieval, &
                     revec, ndof*nx, revec, ndof*nx,          &
                     rwork, lwork, info )
        write(*,*) 'DGEEV info = ',info
        deallocate( rwork )
        write(*,*) 'Finished eigenvalue solver'

        do j = 1, ndof*(nx-1)
          eval(j) = cmplx(reval(j),ieval(j))
        end do

        do j = 2, ndof*(nx-1)
          if (eval(j) == conjg(eval(j-1))) then
            evec(:,:,j-1) = revec(:,:,j-1) + (zero,one)*revec(:,:,j)
            evec(:,:,j)   = revec(:,:,j-1) - (zero,one)*revec(:,:,j)
          end if
        end do

        end if

!.... Correct eigenvectors for periodic boundary conditions

        do j = 1, ndof*nx
          evec(:,nx,j) = evec(:,1,j)
        end do
        do idof = 1, ndof
          i = ndof*(nx-1) + idof
          eval(i) = zero
        end do

!.... Output the eigenvalues and eigenvectors

        open(11,file='eig.dat', form='unformatted', status='unknown')
        write(11) ndof*nx, 0.0, nx, 1, nz, ndof, &
                  Re, Ma, Pr, gamma, cv
        write(11) eval
        write(11) evec
        close(11)

        end subroutine eigdrv3d

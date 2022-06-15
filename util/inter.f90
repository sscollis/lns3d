!=============================================================================!
        program inter
!
!  This program interpolates from one computational mesh to another
!  using B-Spline basis function with the order set by the user.
!  Since this routine works in computational space, you cannot
!  change the mapping function between meshes.
!
!  Revised: 6/15/22
!
!=============================================================================!
        use const
        implicit none

!.... first grid and data

        integer :: nx1, ny1, nz1
        real, allocatable :: xy1(:,:,:), v1(:,:,:), xi1(:), eta1(:)
        real :: dxi1, deta1
        
!.... second grid and data

        integer :: nx2, ny2, nz2
        real, allocatable :: xy2(:,:,:), v2(:,:,:), xi2(:), eta2(:)
        real :: dxi2, deta2

!.... stuff for B-spline

        integer           :: kxord, kyord
        real, allocatable :: x(:,:), y(:,:), xknot(:), yknot(:)
        real, allocatable :: bsv(:,:,:)

!.... stuff for file I/O

        character(80) :: name, filen

!.... local vars

        real :: Re, Ma, Pr, gamma, cv
        integer :: i, j, k, idof, ndof=5, lstep
        real :: tmp, time
        logical, parameter :: useIJ=.true.
!=============================================================================!

!.... read the first grid file

        filen = 'grid.dat'
        write (*,"('Enter grid to interpolate from [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

!.... allocate memory for the physical grid

        open (unit=10, file=filen, form='unformatted', status='old')
        read(10) nx1, ny1, nz1
        write(*,"('Nx = ',i4,', Ny = ',i4)") nx1, ny1
        allocate (xy1(ny1,nx1,2))
        read(10) (((xy1(j,i,1), i = 1, nx1), j = 1, ny1), k = 1, nz1), &
                 (((xy1(j,i,1), i = 1, nx1), j = 1, ny1), k = 1, nz1), &
                 (((       tmp, i = 1, nx1), j = 1, ny1), k = 1, nz1)
        close(10)

!.... allocate memory for the computational grid

        allocate (xi1(nx1), eta1(ny1))

!.... make the computational grid

        dxi1  = two / float(nx1-1)
        deta1 = one / float(ny1-1)
        
        do i = 1, nx1
          xi1(i) = -one + (i-1) * dxi1
        end do
        
        do j = 1, ny1
          eta1(j) = (j-1) * deta1
        end do

!.... read in the first data file

        filen = 'output.dat'
 201    write (*,"('Enter data filename to interpolate from [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

        allocate(v1(ny1,nx1,ndof))
        open(unit=10, file=filen, form='unformatted', status='old')
        read(10) lstep, time, nx1, ny1, nz1, ndof, &
                 Re, Ma, Pr, gamma, cv
        if (useIJ) then
          read(10,err=201) (((v1(j,i,idof), idof=1,ndof), i=1,nx1), j=1,ny1)
        else
          read(10,err=201) v1
        endif
        close(10)

!.... read the second grid file

        filen = 'grid.dat'
        write (*,"('Enter grid to interpolate too [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

!.... allocate memory for the physical grid

        open (unit=10, file=filen, form='unformatted', status='old')
        read(10) nx2, ny2, nz2
        write(*,"('Nx = ',i4,', Ny = ',i4)") nx2, ny2
        allocate (xy2(ny2,nx2,2))
        read(10) (((xy2(j,i,1), i = 1, nx2), j = 1, ny2), k = 1, nz2), &
                 (((xy2(j,i,1), i = 1, nx2), j = 1, ny2), k = 1, nz2), &
                 (((       tmp, i = 1, nx2), j = 1, ny2), k = 1, nz2)
        close(10)

!.... allocate memory for the computational grid

        allocate (xi2(nx2), eta2(ny2))

!.... make the computational grid

        dxi2  = two / float(nx2-1)
        deta2 = one / float(ny2-1)
        
        do i = 1, nx2
          xi2(i) = -one + (i-1) * dxi2
        end do
        
        do j = 1, ny2
          eta2(j) = (j-1) * deta2
        end do

        allocate(v2(ny2,nx2,ndof))

!.... form the B-Spline Interpolant

        write(*,"('Enter korder ==> ',$)")
        read(*,*) kxord
        kyord = kxord

        allocate (xknot(nx1+kxord), yknot(ny1+kyord), bsv(ny1,nx1,ndof))

        call BSNAK( nx1, xi1,  kxord, xknot)
        call BSNAK( ny1, eta1, kyord, yknot)
        do idof = 1, ndof
          call BS2IN( ny1, eta1, nx1, xi1, v1(:,:,idof), ny1, kyord, kxord, &
                      yknot, xknot, bsv(:,:,idof))
        end do
        
!.... Interpolate onto the new mesh

        do idof = 1, ndof
          call BS2GD( 0, 0, ny2, eta2, nx2, xi2, kyord, kxord, yknot, xknot, &
                      ny1, nx1, bsv(:,:,idof), v2(:,:,idof), ny2 )
        end do

!.... write out the new data file

        filen = 'output.dat'
        write (*,"('Enter interpolated data filename [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

        open(unit=10, file=filen, form='unformatted', status='new')
        write(10) 0, 0, nx2, ny2, nz2, ndof, &
                  Re, Ma, Pr, gamma, cv
        if (useIJ) then
          write(10) (((v2(j,i,idof), idof=1,ndof), i=1,nx2), j=1,ny2)
        else
          write(10) v2
        endif
        close(10)

        end program inter

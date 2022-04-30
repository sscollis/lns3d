!=============================================================================!
        program conv
!  
!  Convert Fortran unformatted double PLOT3D file to a binary float 
!  PLOT3D file on the SGI
!
!  Default conversion on the qfile (or use -q) -- To translate a grid file 
!  use -g option
!
!  Author:   S. Scott Collis
!
!  Date:     1-4-97
!  
!=============================================================================!
        use const
        implicit none

!.... flow data

        real, allocatable :: q(:,:,:,:)
                
!.... mesh

        real, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)

!.... parameters

        integer :: lstep, nx, ny, nz, ndof=5
        real  :: Ma, Re, Pr, time

!.... local variables

        integer :: i, j, k, idof

        character*80 file, temp
        integer :: iloc, iend
#ifndef __GFORTRAN__
        integer, external :: iargc
#endif
        real :: tmp
        
        logical :: grid=.false., qfile=.true., twod=.false.
        
        integer, parameter :: mfile = 2
        integer :: narg, iarg, nfile=0, ifile(mfile)
        character*80 :: arg
!=============================================================================! 
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
            select case (arg(1:2))
            case ('-g')
              grid = .true.
              qfile = .false.
            case ('-2')
              twod = .true.
            case ('-h')
              write(*,"('-----------------------------------------------')")
              write(*,"('Usage:  conv [options] [file1] [file2]')")
              write(*,"('-----------------------------------------------')")
              write(*,"(' Convert an unformated double precision PLOT3D')")
              write(*,"(' file to a binary single precision PLOT3D file')")
              write(*,"('-----------------------------------------------')")
              write(*,"('   -h:  this help')")
              write(*,"('   -g:  convert a grid file')")
              write(*,"('   -2:  make a 2-D Plot3d file')")
!             write(*,"('   -q:  convert q data file')")
              write(*,"('-----------------------------------------------')")
              call exit(0)
            case default
              write(*,"('Argument ',i2,' ignored.')") iarg
            end select
          end if
        end do

!.... read in the unformatted PLOT3D grid file

        if (grid) then
        
        if ( nfile .gt. 0 ) then
          call getarg(ifile(1),temp)
          file = temp
        else
          file = 'grid.dat'
          write (*,"('Enter unformatted grid file name [',a,']? ',$)") &
            file(1:index(file,' ')-1)
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

        open(unit=10,file=file,form='unformatted',status='old')
        read(10) nx, ny, nz
!       write(*,*) nx,ny,nz
        allocate( x(nx,ny,nz), y(nx,ny,nz), z(nx,ny,nz) )
        read(10) (((x(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((z(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... write out the binary PLOT3D grid file

        iloc = index(file,'.dat')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'bin'
        else
          temp = file(1:iend)//'.bin'
        end if
        if ( nfile .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter binary PLOT3D file name [',a,']? ',$)")   &
                  file(1:index(file,' ')-1)
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( nfile .eq. 1 ) then
          file = temp
        else if ( nfile .eq. 2 ) then
          call getarg(ifile(2),temp)
          file = temp
        end if
        iloc = index(file,' ')
        file = file(1:(iloc-1))//char(0)

        if (twod) then
          call wgrid2d(x, y, nx, ny, nx, file)
        else
          call wgrid  (x, y, z, nx, ny, nz, nx, ny, file)
        end if

        end if
        
!.... read the unformatted Plot3d q-file

        if (qfile) then
        
        if ( nfile .gt. 0 ) then
          call getarg(ifile(1),temp)
          file = temp
        else
          file = 'output.q.1'
          write (*,"('Enter unformatted qfile name [',a,']? ',$)") &
            file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        end if
 50     open(unit=10, file=file, form='unformatted', status='old', err=60)
        goto 70
 60     write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
        write (*,"('Enter file name [',a,']? ',$)") file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 50
 70     continue
        
        read(10) nx, ny, nz
        read(10) Ma, Pr, Re, time
        allocate ( q(nx,ny,nz,ndof) )
        read(10) ((((q(i,j,k,idof), i = 1, nx), j = 1, ny),  &
                                    k = 1, nz), idof = 1, ndof)
        close(10)

!.... write out the binary Plot3d q-file

        iloc = index(file,'.dat')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'bin'
        else
          temp = file(1:iend)//'.bin'
        end if

        iloc = index(file,'.q.')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'b'//file(iloc+2:iend)
        else
          iloc = index(file,'.dat')
          iend = index(file,' ')-1
          if (iloc .ne. 0) then
            temp = file(1:iloc)//'bin'
          else
            temp = file(1:iend)//'.bin'
          end if
        end if
        if ( nfile .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter binary PLOT3D qfile name [',a,']? ',$)")  &
                  file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( nfile .eq. 1 ) then
          file = temp
        else if ( nfile .eq. 2 ) then
          call getarg(ifile(2),temp)
          file = temp
        end if
        iloc = index(file,' ')
        file = file(1:(iloc-1))//char(0)
        
        if (twod) then
          call wdata2d( Ma, Pr, Re, time, q(1,1,1,1), q(1,1,1,2), &
                        q(1,1,1,3), q(1,1,1,5), nx, ny, nx, file )
        else
          call wdata  ( Ma, Pr, Re, time, q(1,1,1,1), q(1,1,1,2), &
                        q(1,1,1,3), q(1,1,1,4), q(1,1,1,5), &
                        nx, ny, nz, nx, ny, file )
        end if

        end if 
        
        call exit(0)    

 1000   write(*,"('>> Error in file I/O')")
        call exit(1)

        end program conv

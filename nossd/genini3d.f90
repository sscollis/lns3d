!=============================================================================!
	subroutine genini3d(v, vold)
!  
!  Read the initial condition for the 3D code
!  
!=============================================================================!
	use global
	use material
	implicit none

	complex :: v(ny*nx*ndof), vold(ny*nx*ndof)
	integer :: nxi, nyi, nzi, ndofi
	real    :: Rei, Mai, Pri, gammai, cvi
	character*80 :: name, code='GenIni3D$'
	
	integer, external :: igetver
!=============================================================================!

!.... read the restart file

	iver = igetver( base, 'R'//char(0) )
	call addver(base, 'R'//char(0), iver, filen, lfile)
        write (*,10) filen(1:index(filen,' ')-1)
  10    format(' Restarting from file:  ',a)
  
	open(unit=10, file=filen, form='unformatted', status='old', &
	     err=1000)
	
	read(10,err=1100,end=1100) lstep, time, nxi, nyi, nzi, ndofi, &
	                           Rei, Mai, Pri, gammai, cvi

!.... check to make sure that the grid and restart files are compatible

	if (nxi .ne. nx .or. nyi .ne. ny .or. nzi .ne. nz .or.  &
	    ndofi .ne. ndof) then
	  call error(code,'Grid and data files are incompatible$')
	end if

!.... check flow parameters against the input file.  Issue a warning on
!.... on a mismatch, but continue execution using the input values.

	if (Rei .ne. Re) &
	  write(*,100) 'WARNING: Re does not match input file', Rei, Re
	if (Mai .ne. Ma) &
	  write(*,100) 'WARNING: Ma does not match input file', Mai, Ma
	if (Pri .ne. Pr) &
	  write(*,100) 'WARNING: Pr does not match input file', Pri, Pr
	if (gammai .ne. gamma) &
	  write(*,100) 'WARNING: Gamma does not match input file', gammai,gamma
	if (cvi .ne. cv) &
	  write(*,*) 'WARNING: Cv does not match input file', cvi, cv
	
	read(10,err=1200,end=1200) v

!.... read the previous field if available

	if (lstep .gt. 0 .and. impl.eq.2) read(10,err=1300,end=1300) vold
	
	close(10)
	
	return

100	format(a,1x,'(',1pe11.3E3,',',1pe11.3E3,')')
1000    call error(code,'Could not open restart file$')
1100    call error(code,'Error reading info from restart file$')
1200    call error(code,'Error reading V from restart file$')
1300    call error(code,'Error reading Vold from restart file$')

	end

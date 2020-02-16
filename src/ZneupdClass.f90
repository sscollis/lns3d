
MODULE ZneupdClass 
  use ZnaupdClass
  TYPE :: ZneupdIO
    LOGICAL :: rvec                               ! Compute Ritz?
    CHARACTER :: howmny*1                         ! Specifies the form of the basis 
    LOGICAL, POINTER, DIMENSION(:) :: select      ! dimension ncv  
    COMPLEX, POINTER, DIMENSION(:) :: D           ! NEV+1 
    COMPLEX, POINTER, DIMENSION(:,:) :: Z         ! N by NEV array
    INTEGER :: LDZ                                ! leading dimension of Z
    COMPLEX :: sigma         
    COMPLEX, POINTER, DIMENSION(:) :: workev      ! 2*NCV 
  END TYPE ZneupdIO 

CONTAINS

!... Setup for Standard Eigenvalue Problem

  SUBROUTINE StandardZneupdSetup(myZneupd, myZnaupd)
    IMPLICIT NONE

       TYPE(ZnaupdIO), INTENT(IN) :: myZnaupd
       TYPE(ZneupdIO), INTENT(INOUT) :: myZneupd 

       write(*,*) 'Beginning ZneupdClass Setup...'
       write(*,*)
       write(*,*)
       write(*,*) '**************************************************'
       write(*,*) 'Welcome and Enjoy Using and Re-using ZneupdClass.'
       write(*,*) '**************************************************'
       write(*,*)

       myZneupd%rvec = .true.        
       myZneupd%howmny = 'A' 
       myZneupd%ldz = myZnaupd%ldv 

       ALLOCATE( myZneupd%select(myZnaupd%ncv), &
                 myZneupd%D(myZnaupd%nev+1),    &
                 myZneupd%Z(myZnaupd%n,myZnaupd%nev), &
                 myZneupd%workev(2*myZnaupd%ncv) )
 
  END SUBROUTINE StandardZneupdSetup
     
  SUBROUTINE CallZneupd(myZneupd, myZnaupd)
   IMPLICIT NONE
       TYPE(ZnaupdIO), INTENT(INOUT) :: myZnaupd 
       TYPE(ZneupdIO), INTENT(INOUT) :: myZneupd
         CALL zneupd ( myZneupd%rvec, myZneupd%howmny, myZneupd%select, myZneupd%D, &
             myZnaupd%v, myZnaupd%ldv, myZneupd%sigma, myZneupd%workev,             &
             myZnaupd%bmat, myZnaupd%n, myZnaupd%which, myZnaupd%nev,               &
             myZnaupd%tol, myZnaupd%resid, myZnaupd%ncv, myZnaupd%v, myZnaupd%ldv,  &
             myZnaupd%iparam, myZnaupd%ipntr, myZnaupd%workd, myZnaupd%workl,       &  
             myZnaupd%lworkl, myZnaupd%rwork, myZnaupd%info)
  END SUBROUTINE CallZneupd

END MODULE  ZneupdClass

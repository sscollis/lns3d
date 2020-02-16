!-----------------------------------------------------------------------------
        subroutine error(name,msg)
!  
!       Generic error handler 
!  
!-----------------------------------------------------------------------------
        implicit none
        
        integer loc
        character*80 name, msg

        loc = index(name,'$')-1
        write(*,"(/,'Error in:  ',a)") name(1:loc)
        loc = index(msg,'$')-1
        write(*,"('  ',a)") msg(1:loc)
        
        call exit(1)
        
        stop
        end

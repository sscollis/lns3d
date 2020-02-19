!-----------------------------------------------------------------------------
        subroutine getmat(t, rmu, rlm, con, dmu, d2mu, dlm, d2lm, dcon, d2con)
                           
!.... Calculate the Navier-Stokes material properties
!.... Note: everything is nondimensional

        use stuff
        implicit none
        real t(:), rmu(:), rlm(:), con(:), dmu(:), d2mu(:)
        real dlm(:), d2lm(:), dcon(:), d2con(:)
        
        if (mattyp .eq. 0) then         ! constant visc
          rmu  = one
          dmu  = zero
          d2mu = zero
          dlm  = zero
          d2lm = zero
        else                            ! Sutherland's law
          rmu  = t**onept5*(one+S0)/(t+S0)
          dmu  = sqrt(t)*(one+S0)*(3.0*S0+t)*pt5/(S0+t)**2
          d2mu = (one+S0)*(3.0*S0**2-6.0*S0*t-t**2)*pt25/sqrt(t)/(S0+t)**3
        end if

        con   = rmu
        dcon  = dmu
        d2con = d2mu
        
        rlm  = -pt66 * rmu
        dlm  = -pt66 * dmu
        d2lm = -pt66 * d2mu
        
        return
        end

!-----------------------------------------------------------------------------
        subroutine sgetmat(t, rmu, rlm, con, dmu, d2mu, dlm, d2lm, dcon, d2con)
                           
!.... Calculate the Navier-Stokes material properties

        use stuff
        implicit none
        real t, rmu, rlm, con, dmu, d2mu
        real dlm, d2lm, dcon, d2con
        
        if (mattyp .eq. 0) then         ! constant visc
          rmu  = one
          dmu  = zero
          d2mu = zero
          dlm  = zero
          d2lm = zero
        else                            ! Sutherland's law
          rmu  = t**onept5*(one+S0)/(t+S0)
          dmu  = sqrt(t)*(one+S0)*(3.0*S0+t)*pt5/(S0+t)**2
          d2mu = (one+S0)*(3.0*S0**2-6.0*S0*t-t**2)*pt25/sqrt(t)/(S0+t)**3
        end if

        con = rmu
        dcon = dmu
        d2con = d2mu
        
        rlm = -pt66 * rmu
        dlm = -pt66 * dmu
        d2lm = -pt66 * d2mu
        
        return
        end


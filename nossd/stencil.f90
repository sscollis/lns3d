!=============================================================================!
        module stencil
!
!  Finite Difference coefficient of various orders
!
!=============================================================================!
!.... First derivatives
!=============================================================================!

!.... fourth order central difference ( 1 2 x 4 5 )

        real, parameter :: ga1 =  8.333333333333333333333E-02
        real, parameter :: ga2 = -6.666666666666666666667E-01
        real, parameter :: ga3 =  6.666666666666666666667E-01
        real, parameter :: ga4 = -8.333333333333333333333E-02

!.... fourth order one-sided ( x 2 3 4 5 ) 

        real, parameter :: gc1 = -2.083333333333333333333E+00
        real, parameter :: gc2 =  4.000000000000000000000E+00
        real, parameter :: gc3 = -3.000000000000000000000E+00
        real, parameter :: gc4 =  1.333333333333333333333E+00
        real, parameter :: gc5 = -2.500000000000000000000E-01

!.... fourth order biased difference ( 1 x 2 3 4 5 )

        real, parameter :: gb1 = -2.500000000000000000000E-01
        real, parameter :: gb2 = -8.333333333333333333333E-01
        real, parameter :: gb3 =  1.500000000000000000000E+00
        real, parameter :: gb4 = -5.000000000000000000000E-01
        real, parameter :: gb5 =  8.333333333333333333333E-02

!.... sixth order one-sided ( x 2 3 4 5 6 7 )

        real, parameter :: ge1 = -2.450000000000000000000E+00
        real, parameter :: ge2 =  6.000000000000000000000E+00
        real, parameter :: ge3 = -7.500000000000000000000E+00
        real, parameter :: ge4 =  6.666666666666666666667E+00
        real, parameter :: ge5 = -3.750000000000000000000E+00
        real, parameter :: ge6 =  1.200000000000000000000E+00
        real, parameter :: ge7 = -1.666666666666666666667E-01

!.... sixth order one-pt biased ( 1 x 3 4 5 6 7 ) 

        real, parameter :: gf1 = -1.666666666666666666667E-01
        real, parameter :: gf2 = -1.283333333333333333333E+00
        real, parameter :: gf3 =  2.500000000000000000000E+00 
        real, parameter :: gf4 = -1.666666666666666666667E+00
        real, parameter :: gf5 =  8.333333333333333333333E-01
        real, parameter :: gf6 = -2.500000000000000000000E-01
        real, parameter :: gf7 =  3.333333333333333333333E-02

!==============================================================================
!.... The following stencils are Carpenter's stable and third
!.... order accurate boundary treatment for the explicit fourth
!.... order interior scheme.
!==============================================================================

!.... third order one-sided [Carpenter] ( x 2 3 4 5 6) 
  
        real, parameter :: gg1 = -1.8760320556207377229
        real, parameter :: gg2 =  3.185383225577892867
        real, parameter :: gg3 = -1.8145456794375275725
        real, parameter :: gg4 =  0.5916582410526027442
        real, parameter :: gg5 = -0.10105206800050562464
        real, parameter :: gg6 =  0.014588336428275308769

!.... third order biased [Carpenter] ( 1 x 3 4 5 ) 
  
        real, parameter :: gh1 = -0.38425423267792540204
        real, parameter :: gh2 = -0.29063894776734868107
        real, parameter :: gh3 =  0.6717647845153154114
        real, parameter :: gh4 =  0.07108165983739987271
        real, parameter :: gh5 = -0.07363071876172424507
        real, parameter :: gh6 =  0.00567745485428304409

!.... third order biased [Carpenter] ( 1 2 x 4 5 ) 
  
        real, parameter :: gi1 =  0.18288527868682620658
        real, parameter :: gi2 = -1.0800147541745551643
        real, parameter :: gi3 =  0.6578728964966252582
        real, parameter :: gi4 =  0.17761704868919314564
        real, parameter :: gi5 =  0.0767798363958275586
        real, parameter :: gi6 = -0.015140306093917004671

!.... third order biased [Carpenter] ( 1 2 3 x 5 ) 
  
        real, parameter :: gj1 = -0.03418371033652918578
        real, parameter :: gj2 =  0.22482902574010312173
        real, parameter :: gj3 = -0.8908123329284539625
        real, parameter :: gj4 =  0.16529994771003501478
        real, parameter :: gj5 =  0.6134395520875252998
        real, parameter :: gj6 = -0.07857248227268028806

!=============================================================================!
!.... Second derivatives
!=============================================================================!

!.... fourth order central difference

        real, parameter :: da1 = -8.333333333333333333333E-02
        real, parameter :: da2 =  1.333333333333333333333E+00
        real, parameter :: da3 = -2.500000000000000000000E+00
        real, parameter :: da4 =  1.333333333333333333333E+00
        real, parameter :: da5 = -8.333333333333333333333E-02
        
!.... third order biased difference

        real, parameter :: db1 =  9.166666666666666666667E-01
        real, parameter :: db2 = -1.666666666666666666667E+00
        real, parameter :: db3 =  5.000000000000000000000E-01
        real, parameter :: db4 =  3.333333333333333333333E-01
        real, parameter :: db5 = -8.333333333333333333333E-02
        
!.... fourth order one-sided (assumes f'=0) [not smooth]

        real, parameter :: dc1 = -5.763888888888888888889E+00
        real, parameter :: dc2 =  8.000000000000000000000E+00
        real, parameter :: dc3 = -3.000000000000000000000E+00
        real, parameter :: dc4 =  8.888888888888888888889E-01
        real, parameter :: dc5 = -1.250000000000000000000E-01
        
!.... third order one-sided

        real, parameter :: dd1 =  2.916666666666666666667E+00
        real, parameter :: dd2 = -8.666666666666666666667E+00
        real, parameter :: dd3 =  9.500000000000000000000E+00
        real, parameter :: dd4 = -4.666666666666666666667E+00
        real, parameter :: dd5 =  9.166666666666666666667E-01
        
!=============================================================================!
!.... Fourth derivatives
!=============================================================================!

!.... second-order, fourth derivative

        real, parameter :: fb1 =  1.000000000000000E+00
        real, parameter :: fb2 = -4.000000000000000E+00
        real, parameter :: fb3 =  6.000000000000000E+00
        real, parameter :: fb4 = -4.000000000000000E+00
        real, parameter :: fb5 =  1.000000000000000E+00

!.... fourth-order, fourth derivative

        real, parameter :: fa1 = -1.666666666666667E-01
        real, parameter :: fa2 =  2.000000000000000E+00
        real, parameter :: fa3 = -6.500000000000000E+00
        real, parameter :: fa4 =  9.333333333333333E+00
        real, parameter :: fa5 = -6.500000000000000E+00
        real, parameter :: fa6 =  2.000000000000000E+00
        real, parameter :: fa7 = -1.666666666666667E-01
        
!.... sixth-order, fourth derivative

        real, parameter :: fc1 =  2.916666666666667E-02
        real, parameter :: fc2 = -4.000000000000000E-01
        real, parameter :: fc3 =  2.816666666666667E+00
        real, parameter :: fc4 = -8.133333333333333E+00
        real, parameter :: fc5 =  1.137500000000000E+01
        real, parameter :: fc6 = -8.133333333333333E+00
        real, parameter :: fc7 =  2.816666666666667E+00
        real, parameter :: fc8 = -4.000000000000000E-01
        real, parameter :: fc9 =  2.916666666666667E-02

        end module stencil


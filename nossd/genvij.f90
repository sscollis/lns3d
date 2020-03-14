!=============================================================================!
	subroutine genVij(Vij)
!
!  Form the Vij terms in compact storage.
!
!  Compact Storage Scheme:
!
!	Vxx = D[ 0, Vij(1), Vij(2), Vij(2), Vij(3)]
!	Vyy = D[ 0, Vij(2), Vij(1), Vij(2), Vij(3)]
!	Vzz = D[ 0, Vij(2), Vij(2), Vij(1), Vij(3)]
!	Vxy(2,3) = Vij(4)
!       Vxy(3,2) = Vij(4)
!	Vxz(2,4) = Vij(4)
!	Vxz(4,2) = Vij(4)
!	Vyz(3,4) = Vij(4)
!	Vyz(4,3) = Vij(4)
!
!=============================================================================!
	use global
	use local
	implicit none
	
	real :: Vij(ny*nx,4)

	real :: fact1, fact2
!=============================================================================!

	fact1 = one / Re
	fact2 = gamma / (Pr * Re)

	Vij(:,1) = fact1 * rhoinv * ( lm + two * mu )
	Vij(:,2) = fact1 * rhoinv * mu
	Vij(:,3) = fact2 * rhoinv * con
	Vij(:,4) = fact1 * rhoinv * ( lm + mu )

	return
	end

module spring_mod
	implicit none
	private
	real*8, parameter :: epsil = .1
	real*8, parameter :: TWOSQRT3 = 3.46410161513775d0
	real*8, parameter :: ERROR_TOL = 0.01
	integer, parameter :: MAX_ITER = 200
	integer, parameter :: PRINT_INTERVAL = 500
	real*8, dimension(:), allocatable :: X_aux, Y_aux
	real*8, dimension(:), allocatable :: R, D, b
	real*8, dimension(:), allocatable :: S, E
	public :: main
contains
	subroutine main(X, Y, inpoel, fixed)
		use PointNeighbor!, only: psup1, psup2
		real*8, dimension(:), intent(inout) :: X, Y
		integer, dimension(:,:), intent(in) :: inpoel
		logical, dimension(:), intent(in) :: fixed
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		integer, save :: npoin, nelem, counter
		logical, save :: isFirstCall = .true.
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		counter = counter + 1
		if(isFirstCall) then
			npoin = size(X,1)
			nelem = size(inpoel,2)
			counter = PRINT_INTERVAL
			allocate(b(npoin))
			allocate(X_aux(npoin), Y_aux(npoin))
			call buildSE()
			isFirstCall = .false.
		end if
		if(.not.allocated(psup1) .or. .not.allocated(psup2)) call getPsup(inpoel, nelem, npoin)
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		call buildRD()
		! b = 0
		call SpMV(R, psup1, psup2, X, b)
		b = b + D*x
		call jacobi(R, D, psup1, psup2, X, b)
		call SpMV(R, psup1, psup2, Y, b)
		b = b + D*y
		call jacobi(R, D, psup1, psup2, Y, b)

		! X_aux = X
		! call jacobi(S, E, psup1, psup2, X, X_aux)
		! Y_aux = Y
		! call jacobi(S, E, psup1, psup2, Y, Y_aux)
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(counter == PRINT_INTERVAL) then
			print*, 'Min Distorion Metric:', getMeshQuality()
			counter = 0
		end if
	contains
		subroutine buildRD()
			logical, save :: isFirstCall = .true.
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			integer :: ipoin, jpoin, ipsup
			real*8 :: dist, c
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(isFirstCall) then
				allocate(R(size(psup1)), D(npoin))
				isFirstCall = .false.
			end if
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ipoin, ipsup, jpoin, dist, c)
			do ipoin = 1, npoin
			D(ipoin) = 0
			do ipsup = psup2(ipoin) + 1, psup2(ipoin + 1)
			jpoin = psup1(ipsup)
			dist = norm2((/X(ipoin) - X(jpoin), Y(ipoin) - Y(jpoin)/))
			c = 1d0/(dist**2)
			R(ipsup) = -c
			D(ipoin) = D(ipoin) + c
			end do
			end do
			!$OMP END PARALLEL DO
		end subroutine buildRD

		subroutine buildSE()
			logical, save :: isFirstCall = .true.
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			integer :: ipoin, ipsup, n
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(isFirstCall) then
				allocate(S(size(psup1)), E(npoin))
				isFirstCall = .false.
			end if
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ipoin, ipsup)
			do ipoin = 1, npoin
			n = psup2(ipoin + 1) - psup2(ipoin)
			E(ipoin) = 1 + epsil*n
			do ipsup = psup2(ipoin) + 1, psup2(ipoin + 1)
			S(ipsup) = -epsil
			end do
			end do
			!$OMP END PARALLEL DO
		end subroutine buildSE

		subroutine jacobi(R, D, idx, ptr, x, b)
			real*8, dimension(:) :: R, D, x, b
			integer, dimension(:) :: idx, ptr
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			real*8, dimension(:), allocatable, save :: x_fix, y
			integer :: iter
			logical, save :: isFirstCall = .true.
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(isFirstCall) then
				allocate(x_fix(size(x)))
				allocate(y(size(x)))
				isFirstCall = .false.
			end if
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			where(fixed) x_fix = x
			do iter = 1, MAX_ITER
			call SpMV(R, idx, ptr, x, y)
			x = (b - y)/D
			where(fixed) x = x_fix
			end do
		end subroutine jacobi

		subroutine SpMV(A, idx, ptr, x, y)
			real*8, dimension(:) :: A, x, y
			integer, dimension(:) :: idx, ptr
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			real*8 :: dot
			integer :: i, j, n
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			n = size(x)
			!$OMP PARALLEL DO PRIVATE(i, j, dot)
			do i = 1, n
			dot = 0.d0
			do j = ptr(i) + 1, ptr(i + 1)
			dot = dot + A(j)*x(idx(j))
			end do
			y(i) = dot
			end do
			!$OMP END PARALLEL DO 
		end subroutine SpMV

		! real*8 function error(last)
		! 	logical, intent(in), optional :: last
		! 	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		! 	logical, save :: isFirstCall = .true.
		! 	real*8, dimension(:), allocatable, save :: distance_old
		! 	real*8, save :: error_last
		! 	integer :: ipoin
		! 	real*8 :: error0, distance
		! 	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		! 	if(isFirstCall) then 
		! 		allocate(distance_old(npoin))
		! 		distance_old(:) = tiny(0d0)
		! 		isFirstCall = .false.
		! 		if(present(last)) stop 'laplacianSmoothing_mod: No se dispone de error_last'
		! 	end if
		! 	if(present(last)) then
		! 		error = error_last
		! 		return 
		! 	end if
		! 	error = 0
		! 	do ipoin = 1, npoin
		! 	distance = norm2((/X(ipoin) - coord_new(1,ipoin),Y(ipoin) - coord_new(2,ipoin)/))
		! 	if(distance < tiny(0d0)) distance = tiny(0d0)
		! 	error0 = dabs(1 - distance/distance_old(ipoin))
		! 	if(error0 > error) error = error0
		! 	distance_old(ipoin) = distance
		! 	end do
		! 	error_last = error
		! end function error

		! subroutine displayStats()
		! 	print*, 'Estimacion del error:', error(last=.true.)
		! 	print*, '# iteraciones:', iter
		! end subroutine displayStats

		real*8 function getMeshQuality()
			!%%%%%%%%%%%%%%%%%%%% 
			integer :: ielem
			real*8 :: X_loc(3), Y_loc(3)
			real*8 :: mu1
			!%%%%%%%%%%%%%%%%%%%% 
			getMeshQuality = 1
			do ielem = 1, nelem
			X_loc(:) = X(inpoel(:,ielem))
			Y_loc(:) = Y(inpoel(:,ielem))
			mu1 = mu(X_loc, Y_loc)
			if(mu1 < getMeshQuality) getMeshQuality = mu1
			end do
		end function getMeshQuality

		pure real*8 function mu(X_loc, Y_loc)
			real*8, intent(in) :: X_loc(:), Y_loc(:)
			!%%%%%%%%%%%%%%%%%%%%
			real*8 :: l1, l2, l3, l, area
			!%%%%%%%%%%%%%%%%%%%%
			area = X_loc(2)*Y_loc(3)+X_loc(3)*Y_loc(1)+X_loc(1)*Y_loc(2) &
				-(X_loc(2)*Y_loc(1)+X_loc(3)*Y_loc(2)+X_loc(1)*Y_loc(3))
			l1 = (X_loc(3)-X_loc(2))**2 + (Y_loc(3)-Y_loc(2))**2
			l2 = (X_loc(1)-X_loc(3))**2 + (Y_loc(1)-Y_loc(3))**2
			l3 = (X_loc(2)-X_loc(1))**2 + (Y_loc(2)-Y_loc(1))**2
			l = l1 + l2 + l3
			mu = TWOSQRT3*area/l
		end function mu
	end subroutine main
end module spring_mod

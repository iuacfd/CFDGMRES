module laplacianSmoothing_mod
	implicit none
	private
	real*8, parameter :: e = .2d0
	real*8, parameter :: dt = 0.2d0
	real*8, parameter :: ERROR_TOL = 0.01
	integer, parameter :: MAX_ITER = 100000
	integer, parameter :: PRINT_INTERVAL = 10
	real*8, dimension(:,:), allocatable :: coord_new
	public :: main
contains
	subroutine main(X, Y, inpoel, fixed)
		use PointNeighbor!, only: psup1, psup2
		real*8, dimension(:), intent(inout) :: X, Y
		integer, dimension(:,:), intent(in) :: inpoel
		logical, dimension(:), intent(in) :: fixed
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		integer, save :: npoin, nelem, counter
		real*8 :: X_sum, Y_sum
		real*8 :: x_j, y_j, l
		integer :: iter, ipoin, ipsup, n
		logical, save :: isFirstCall = .true.
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		counter = counter + 1
		if(isFirstCall) then
			npoin = size(X,1)
			nelem = size(inpoel,2)
			allocate(coord_new(2,npoin))
			counter = PRINT_INTERVAL
			isFirstCall = .false.
		end if
		if(.not.allocated(psup1) .or. .not.allocated(psup2)) call getPsup(inpoel, nelem, npoin)
		do iter = 1, MAX_ITER
		!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(X, Y, coord_new, psup1, psup2, fixed, inpoel, npoin)
		do ipoin = 1, npoin
		if(.not.fixed(ipoin)) then
			X_sum = 0
			Y_sum = 0
			n = psup2(ipoin + 1) - psup2(ipoin)
			do ipsup = psup2(ipoin) + 1, psup2(ipoin + 1)
			x_j = X(psup1(ipsup))
			y_j = Y(psup1(ipsup))
			l = norm2((/x_j - X(ipoin), y_j - Y(ipoin)/))
			X_sum = X_sum + (x_j - X(ipoin))/(l)
			Y_sum = Y_sum + (y_j - Y(ipoin))/(l)
			end do
			coord_new(1,ipoin) = X(ipoin) + X_sum*dt/n
			coord_new(2,ipoin) = Y(ipoin) + Y_sum*dt/n
		else
			coord_new(1,ipoin) = X(ipoin)
			coord_new(2,ipoin) = Y(ipoin)
		end if
		end do
		!$OMP END PARALLEL DO
		if(error() < ERROR_TOL) then
			! UPDATE COORDINATES
			X = coord_new(1,:)
			Y = coord_new(2,:)
			exit
		end if
		! UPDATE COORDINATES
		X = coord_new(1,:)
		Y = coord_new(2,:)
		end do
		if(counter == PRINT_INTERVAL) then
			call displayStats()
			counter = 0
		end if
	contains
		real*8 function error(last)
			logical, intent(in), optional :: last
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			logical, save :: isFirstCall = .true.
			real*8, dimension(:), allocatable, save :: distance_old
			real*8, save :: error_last
			integer :: ipoin
			real*8 :: error0, distance
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(isFirstCall) then 
				allocate(distance_old(npoin))
				distance_old(:) = tiny(0d0)
				isFirstCall = .false.
				if(present(last)) stop 'laplacianSmoothing_mod: No se dispone de error_last'
			end if
			if(present(last)) then
				error = error_last
				return 
			end if
			error = 0
			do ipoin = 1, npoin
			distance = norm2((/X(ipoin) - coord_new(1,ipoin),Y(ipoin) - coord_new(2,ipoin)/))
			if(distance < tiny(0d0)) distance = tiny(0d0)
			error0 = dabs(1 - distance/distance_old(ipoin))
			if(error0 > error) error = error0
			distance_old(ipoin) = distance
			end do
			error_last = error
		end function error

		subroutine displayStats()
			print*, 'Estimacion del error:', error(last=.true.)
			print*, '# iteraciones:', iter
		end subroutine displayStats
	end subroutine
end module

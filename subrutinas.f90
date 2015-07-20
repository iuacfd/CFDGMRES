#define timer(func, store) call system_clock(start_t, rate); call func; call system_clock(end_t); store  = store + real(end_t - start_t)/real(rate);
module Mnormales
	integer, private :: m
	integer, dimension(:), allocatable, private :: n_ipoin
	real(8), dimension(:), allocatable, private :: n_x, n_y, num_l_x, num_l_y, den_l
contains
	subroutine normales
		!use mgeometria
		use meshdata, only: wall, nwall, x, y, npoin
		implicit none
		real(8) l_edge_x, l_edge_y, l_edge, l_x, l_y, l_nrm
		integer ielem, ipoin

		if(.not.allocated(n_ipoin)) allocate(n_ipoin(npoin))
		if(.not.allocated(n_x)) allocate(n_x(npoin))
		if(.not.allocated(n_y)) allocate(n_y(npoin))
		if(.not.allocated(num_l_x)) allocate(num_l_x(npoin))
		if(.not.allocated(num_l_y)) allocate(num_l_y(npoin))
		if(.not.allocated(den_l)) allocate(den_l(npoin))

		num_l_x = 0.d0 
		num_l_y = 0.d0 
		den_l = 0.d0 
		m = 0 

		!$omp parallel 
		!$omp do private(l_edge_x, l_edge_y, l_edge, ielem)
		do ielem = 1, nwall
		l_edge_x = y(wall(2, ielem)) - y(wall(1, ielem))
		l_edge_y = -(x(wall(2, ielem)) - x(wall(1, ielem)))
		l_edge = dsqrt(l_edge_x*l_edge_x + l_edge_y*l_edge_y)

		!$omp atomic
		num_l_x(wall(1, ielem)) = num_l_x(wall(1, ielem)) + l_edge_x
		!$omp atomic
		num_l_y(wall(1, ielem)) = num_l_y(wall(1, ielem)) + l_edge_y
		!$omp atomic
		den_l(wall(1, ielem)) = den_l(wall(1, ielem)) + l_edge
		!$omp atomic
		num_l_x(wall(2, ielem)) = num_l_x(wall(2, ielem)) + l_edge_x
		!$omp atomic
		num_l_y(wall(2, ielem)) = num_l_y(wall(2, ielem)) + l_edge_y
		!$omp atomic
		den_l(wall(2, ielem)) = den_l(wall(2, ielem)) + l_edge
		end do
		!$omp end do
		!$omp end parallel

		!paralelizacion de esta parte requiere omp critical en el nested if,
		!en la implementacion demostro que no rinde
		do ipoin = 1, npoin
		if (den_l(ipoin) > 1.d-6) then
			l_x = num_l_x(ipoin)/den_l(ipoin)
			l_y = num_l_y(ipoin)/den_l(ipoin)
			l_nrm = dsqrt(l_x*l_x + l_y*l_y)
			if (l_nrm > 0.2d0) then
				m = m + 1
				n_ipoin(m) = ipoin
				n_x(m) = l_x/l_nrm
				n_y(m) = l_y/l_nrm
			end if
		end if
		end do
	end subroutine normales

 subroutine normalvel(velocidadx,velocidady)
   !use mallocar
   !use mnormales
   use meshdata, only : npoin
   use mvelocidades
   implicit none
   integer i, ipoin 
   real(8) vx, vy, wx, wy, p
   real(8) velocidadx(npoin),velocidady(npoin)
   !$omp parallel do private(vx, vy, wx, wy, p, i, ipoin)
   do i = 1, m
      ipoin = n_ipoin(i)
      !original       vx = vel_x(ipoin); vy = vel_y(ipoin)
      vx=velocidadx(ipoin);vy=velocidady(ipoin) 
      wx = w_x(ipoin); wy = w_y(ipoin)
      p = -n_y(i)*(vx - wx) + n_x(i)*(vy - wy)
      velocidadx(ipoin) = -n_y(i)*p + wx
      velocidady(ipoin) = n_x(i)*p + wy
   end do
   !$omp end parallel do
 end subroutine normalvel
end module Mnormales

subroutine deriv(hmin) !<----- PARA QUE SE USA HMIN?
  !------------------------------------------
  ! Calcula dNx, dNy y area de cada elemento
  !------------------------------------------
  !use MGEOMETRIA
  use MeshData, only: X, Y, inpoel, area, HH, HHX, HHY, dNx, dNy, nelem
  implicit none
  integer ielem
  real(8) hmin
  real(8) X_loc(3), Y_loc(3), a

  !$OMP PARALLEL &
  !$OMP PRIVATE(X_loc, Y_loc, ielem, a)
  !$OMP DO
  do ielem = 1,nelem
     X_loc(:) = X(inpoel(:,ielem))
     Y_loc(:) = Y(inpoel(:,ielem))

     area(ielem) = dabs((X_loc(2)*Y_loc(3) + X_loc(3)*Y_loc(1) + X_loc(1)*Y_loc(2)&
          - (X_loc(2)*Y_loc(1) + X_loc(3)*Y_loc(2) + X_loc(1)*Y_loc(3)))*0.5d0)
     a = area(ielem)

     dNx(1, ielem) = (Y_loc(2) - Y_loc(3))/(2.d0*a)
     dNx(2, ielem) = (Y_loc(3) - Y_loc(1))/(2.d0*a)
     dNx(3, ielem) = (Y_loc(1) - Y_loc(2))/(2.d0*a)
     dNy(1, ielem) = (X_loc(3) - X_loc(2))/(2.d0*a)
     dNy(2, ielem) = (X_loc(1) - X_loc(3))/(2.d0*a)
     dNy(3, ielem) = (X_loc(2) - X_loc(1))/(2.d0*a)

     HH(ielem) = dsqrt(area(ielem))
     ! HHX(ielem) = abs(min( X_loc(3) - X_loc(2), X_loc(1) - X_loc(3), X_loc(2) - X_loc(1) ))
     ! HHY(ielem) = abs(min( Y_loc(3) - Y_loc(2), Y_loc(1) - Y_loc(3), Y_loc(2) - Y_loc(1) ))
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! hmin = minval(hh) 
  ! if(hmin > 1.d10) hmin = 1.d10
end subroutine deriv

subroutine MASAS()
	!----------------------------------------
	! Ensambla la matriz de masas condensada
	!----------------------------------------
	!use MGEOMETRIA
	use MeshData, only: M, inpoel, area, nelem, npoin, dNx, dNy
	implicit none
	integer ipoin, ielem

	!$OMP PARALLEL DO PRIVATE(ipoin)
	do ipoin = 1, npoin
	M(ipoin) = 0.d0
	end do
	!$OMP END PARALLEL DO

	!$OMP PARALLEL DO PRIVATE(ielem)
	do ielem = 1, nelem
	!$OMP ATOMIC
	M(inpoel(1, ielem)) = M(inpoel(1, ielem)) + area(ielem)/3.d0
	!$OMP ATOMIC
	M(inpoel(2, ielem)) = M(inpoel(2, ielem)) + area(ielem)/3.d0
	!$OMP ATOMIC
	M(inpoel(3, ielem)) = M(inpoel(3, ielem)) + area(ielem)/3.d0
	end do
	!$OMP END PARALLEL DO
end subroutine MASAS

subroutine deltat(dtmin, dt)
  !-------------------------------
  !		Calcula delta t
  !-------------------------------
  !use MALLOCAR
  !use DATOS_ENTRADA
  use MeshData, only: inpoel, nelem, hh
  use InputData, only: FSAFE, fr, gama, t_inf
  use MVELOCIDADES
  use MVARIABLES
  implicit real(8) (A-H,O-Z)
  real(8) DT(nelem), dtmin
  integer ielem, i, ipoin
  real(8) T_iel, vumax, vvmax, vu, vv, vel, dtelem, cota, vc

  DTMIN = 1.D20

  do ielem = 1, nelem
     VUMAX = 0.d0
     VVMAX = 0.d0

     T_iel = (T(inpoel(1, ielem)) + T(inpoel(2, ielem)) + T(inpoel(3, ielem)))/3.d0

     !GM = GAMA-1.d0 !(GAMM(inpoel(1, ielem))+GAMM(inpoel(2, ielem))+GAMM(inpoel(3, ielem)))/3.d0
     VC = DSQRT(GAMA*FR*T_iel)

     !CCCC ----> CALCULO DE LA MAXIMA VELOCIDAD ELEMENTAL
     do i = 1, 3
	ipoin = inpoel(i, ielem)
	VU = dabs(VEL_X(ipoin) - W_X(ipoin))
	VV = dabs(VEL_Y(ipoin) - W_Y(ipoin))

	if (VU.GT.VUMAX) VUMAX = VU
	if (VV.GT.VVMAX) VVMAX = VV
     end do

     VEL=(VUMAX**2.D0+VVMAX**2.D0)**.5D0         

     DTELEM=.5D0*HH(IELEM)/(VC+VEL)
     DT(IELEM)=DTELEM*FSAFE


     IF (DTELEM.LT.DTMIN) DTMIN=DTELEM*fsafe
  END DO

  COTA = 10.d0*DTMIN !EL VALOR 10 ESTA PUESTO A OJO #MODIFICAR SI ES NECESARIO#
  do ielem = 1, nelem

     if(DT(ielem) > COTA) DT(ielem) = COTA
  end do

  return
end subroutine deltat

subroutine CUARTO_ORDEN(U, U_n, FR, gamm)
	!---------------------------------------------
	!		Calculo de la proyeccion U_n
	!---------------------------------------------
	!use MALLOCAR
	!use MGEOMETRIA
	use MeshData
	implicit real(8) (A-H,O-Z)
	real(8) U(4,npoin), U_n(4,npoin), GAMM(npoin)

	real(8) A1(3,4,3), A2(3,4,3), Adv(4,3), U_loc(4,3), Ux(4), Uy(4), Nx(3), Ny(3), U_n_tmp(4,3)
	real(8) sp(3,3) !Funciones de forma
	real(8) c(3), vx(3), vy(3), V_sq(3), TEMP(3), e(3)
	INTEGER ipoin(3)

	U_n = 0.d0
	sp(:,1) = (/ .5d0, .5d0, 0.d0 /)
	sp(:,2) = (/ 0.d0, .5d0, .5d0 /)
	sp(:,3) = (/ .5d0, 0.d0, .5d0 /)

	!$OMP PARALLEL
	!$OMP DO PRIVATE(ielem, ipoin, gama, temp, fmu, Nx, Ny, Ux, Uy, ar,&
	!$OMP U_loc, vx, vy, v_sq, e, c, A1, A2, Adv, U_n_tmp, i)
	do ielem = 1,nelem
	U_n_tmp = 0.d0
	ipoin = inpoel(:, ielem)

	gama = (GAMM(ipoin(1))+GAMM(ipoin(2))+GAMM(ipoin(3)))/3.d0

	Nx = dNx(:, ielem)
	Ny = dNy(:, ielem)

	Ux(:) = U(:,ipoin(1))*Nx(1) + U(:,ipoin(2))*Nx(2) + U(:,ipoin(3))*Nx(3)
	Uy(:) = U(:,ipoin(1))*Ny(1) + U(:,ipoin(2))*Ny(2) + U(:,ipoin(3))*Ny(3)

	AR = area(ielem)/3.d0

	!CCCC  ----> INTEGRO LAS VARIABLES EN LOS PUNTOS DE GAUSS
	U_loc(:,1) = sp(1,1)*U(:,ipoin(1)) + sp(2,1)*U(:,ipoin(2)) + sp(3,1)*U(:,ipoin(3))
	U_loc(:,2) = sp(1,2)*U(:,ipoin(1)) + sp(2,2)*U(:,ipoin(2)) + sp(3,2)*U(:,ipoin(3))
	U_loc(:,3) = sp(1,3)*U(:,ipoin(1)) + sp(2,3)*U(:,ipoin(2)) + sp(3,3)*U(:,ipoin(3))

	!CCCC  ----> DEFINO VARIABLES PRIMITIVAS
	vx(:) = U_loc(2,:)/U_loc(1,:)
	vy(:) = U_loc(3,:)/U_loc(1,:)
	e(:) = U_loc(4,:)/U_loc(1,:)
	V_sq(:) = vx*vx+vy*vy
	temp(:) = (gama-1.d0)/FR*(e-.5d0*V_sq) !FR = CTE. UNIVERSAL DE LOS GASES
	c(:) = dsqrt(gama*fr*temp(:))

	! = = = GAUSS QUAD = ==
	do i = 1,3
	!A1(0,:,i) = (/ 0.d0, 1.d0, 0.d0, 0.d0 /)
	A1(1,:,i) = (/ (gama-1.d0)/2.d0*V_sq(i)-vx(i)*vx(i),&
		(3.d0-gama)*vx(i),&
		-(gama-1.d0)*vy(i),&
		(gama-1.d0) /)
	A1(2,:,i) = (/ -vx(i)*vy(i),&
		vy(i),&
		vx(i),&
		0.d0 /)
	A1(3,:,i) = (/ ((gama-1.d0)*V_sq(i)-gama*e(i))*vx(i),&
		gama*e(i)-(gama-1.d0)/2.d0*V_sq(i)-(gama-1.d0)*vx(i)*vx(i),&
		-(gama-1.d0)*vx(i)*vy(i),&
		gama*vx(i) /)
	!A2(0,:,i) = (/ 0.d0, 0.d0, 1.d0, 0.d0 /)
	A2(1,:,i) = (/ -vx(i)*vy(i),&
		vy(i),&
		vx(i),&
		0.d0 /)
	A2(2,:,i) = (/ (gama-1.d0)/2.d0*V_sq(i)-vy(i)*vy(i),&
		-(gama-1.d0)*vx(i),&
		(3.d0-gama)*vy(i),&
		(gama-1.d0) /)
	A2(3,:,i) = (/ ((gama-1.d0)*V_sq(i)-gama*e(i))*vy(i),&
		-(gama-1.d0)*vx(i)*vy(i),&
		gama*e(i)-(gama-1.d0)/2.d0*V_sq(i)-(gama-1.d0)*vy(i)*vy(i),&
		gama*vy(i) /)

	Adv(1,i) = Ux(2) + Uy(3)
	Adv(2:4,i) = A1(:,1,i)*Ux(1) + A1(:,2,i)*Ux(2) + A1(:,3,i)*Ux(3) + A1(:,4,i)*Ux(4)&
		+A2(:,1,i)*Uy(1) + A2(:,2,i)*Uy(2) + A2(:,3,i)*Uy(3) + A2(:,4,i)*Uy(4)
	end do
	! = = = end GAUSS QUAD = ==

	U_n_tmp(:,1) = Adv(:,1)*sp(1,1)*AR + Adv(:,2)*sp(1,2)*AR + Adv(:,3)*sp(1,3)*AR
	U_n_tmp(:,2) = Adv(:,1)*sp(2,1)*AR + Adv(:,2)*sp(2,2)*AR + Adv(:,3)*sp(2,3)*AR
	U_n_tmp(:,3) = Adv(:,1)*sp(3,1)*AR + Adv(:,2)*sp(3,2)*AR + Adv(:,3)*sp(3,3)*AR

	do i = 1,3
	!$OMP ATOMIC
	U_n(1,ipoin(i)) = U_n(1,ipoin(i)) + U_n_tmp(1,i)
	!$OMP ATOMIC
	U_n(2,ipoin(i)) = U_n(2,ipoin(i)) + U_n_tmp(2,i)
	!$OMP ATOMIC
	U_n(3,ipoin(i)) = U_n(3,ipoin(i)) + U_n_tmp(3,i)
	!$OMP ATOMIC
	U_n(4,ipoin(i)) = U_n(4,ipoin(i)) + U_n_tmp(4,i)
	end do

	end do
	!$OMP END DO

	!$OMP DO PRIVATE(i)
	do i= 1, npoin
	U_n(:, i) = -U_n(:, i)/M(i)
	end do

	!$OMP END DO
	!$OMP END PARALLEL
end subroutine CUARTO_ORDEN


SUBROUTINE ESTAB(U,T,GAMA,FR,fMU &
     ,DTMIN,RHOINF,TINF,UINF,VINF,GAMM)

  use MeshData
  use MVELOCIDADES
  use MESTABILIZACION
  use InputData, only: cte
  implicit real(8) (A-H,O-Z) 


  REAL(8) U(4,npoin),T(npoin),GAMM(npoin)


  !CCCC  ----> CTES DE CALCULO PARA SHOCK-CAPTURING      
  CC=DSQRT(GAMA*FR*TINF)
  VEL2=DSQRT(UINF*UINF+VINF*VINF)  

  DO IELEM=1,NELEM

     N1=inpoel(1,IELEM)
     N2=inpoel(2,IELEM)
     N3=inpoel(3,IELEM)
     GM=(GAMM(N1)+GAMM(N2)+GAMM(N3))/3.D0
     TAU=0.D0
     H_RGNE=0.D0
     H_RGN=0.D0
     H_JGN=0.D0

     !CCCC  ----> VARIABLES ELEMENTALES
     RHO_ELEM=(U(1,N1)+U(1,N2)+U(1,N3))/3.D0
     VX=(VEL_X(N1)+VEL_X(N2)+VEL_X(N3))/3.D0
     VY=(VEL_Y(N1)+VEL_Y(N2)+VEL_Y(N3))/3.D0

     !PARTES NUEVAS
     WX=(W_X(N1)+W_X(N2)+W_X(N3))/3.D0
     WY=(W_Y(N1)+W_Y(N2)+W_Y(N3))/3.D0
     VX=VX-WX ; VY=VY-WY
     VEL2=DSQRT(VX*VX+VY*VY)

     !CCCC  ----> DERIVADA DE RHO           
     DRX= U(1,N1)*DNX(1,IELEM)+U(1,N2)*DNX(2,IELEM)+U(1,N3)*DNX(3,IELEM)
     DRY= U(1,N1)*DNY(1,IELEM)+U(1,N2)*DNY(2,IELEM)+U(1,N3)*DNY(3,IELEM)
     DR2=DSQRT(DRX*DRX+DRY*DRY)+1.D-20
     !CCCC  ----> DERIVADA DE TEMPERATURA
     DTX= T(N1)*DNX(1,IELEM)+T(N2)*DNX(2,IELEM)+T(N3)*DNX(3,IELEM)
     DTY= T(N1)*DNY(1,IELEM)+T(N2)*DNY(2,IELEM)+T(N3)*DNY(3,IELEM)
     DT2=DSQRT(DTX*DTX+DTY*DTY)+1.D-20
     !CCCC  ----> DERIVADA DE LA VELOCIDAD 
     DUX= VEL2*DNX(1,IELEM)+VEL2*DNX(2,IELEM)+VEL2*DNX(3,IELEM)
     DUY= VEL2*DNY(1,IELEM)+VEL2*DNY(2,IELEM)+VEL2*DNY(3,IELEM)

     DU2=DSQRT(DUX*DUX+DUY*DUY)+1.D-20

     !CCCC  ----> VECTOR UNIDAD THETA
     RTX=DTX/DT2
     RTY=DTY/DT2
     !CCCC  ----> VECTOR UNIDAD J
     RJX=DRX/DR2
     RJY=DRY/DR2
     !CCCC  ----> VECTOR UNIDAD VELOCIDAD 
     RUX=DUX/DU2
     RUY=DUY/DU2  

     TEMP=(T(N1)+T(N2)+T(N3))/3.D0
     C=DSQRT(GM*FR*TEMP)

     !FMU= RMU*162.6/(TEMP-110.55)*(TEMP/273.15)**.75D0 !SUTHERLAND
     smu = 110.d0
     !fmu = 0.000d0!*(temp/tinf)**1.5d0*(tinf+smu)/(temp+smu)
     DO I = 1,3

	TERM_1 = DABS(VX*DNX(I,IELEM) + VY*DNY(I,IELEM))
	TERM_2 = DABS(RJX*DNX(I,IELEM) + RJY*DNY(I,IELEM))

	H_RGN1 = DABS(RTX*DNX(I,IELEM) + RTY*DNY(I,IELEM)) !CALCULO PARA ECU. ENERGIA
	H_RGN2 = DABS(RUX*DNX(I,IELEM) + RUY*DNY(I,IELEM)) !CALCULO PARA ECU. MOMENTO

	TAU = TAU + TERM_1 + TERM_2*C

	H_RGNE = H_RGNE + H_RGN1
	H_RGN = H_RGN + H_RGN2
	H_JGN = H_JGN + TERM_2
     END DO

     TAU = 1.D0/TAU

     if(H_JGN > 0.01d0) then
        H_JGN = 2.D0/H_JGN         
        TR1 = DR2*H_JGN/RHO_ELEM
        ZZZ = H_JGN/(2.D0*C)
        SHOC(IELEM) = (TR1 + TR1**2)*.5D0*C**2*ZZZ
     else
        SHOC(ielem) = 0.d0
     end if

     RESUMEN = 1.D0/TAU**2.D0 +(2.D0/DTMIN)**2.D0
     RRR = RESUMEN**(-.5D0)  
     T_SUGN1(IELEM) = RRR
     T_SUGN2(IELEM) = RRR
     T_SUGN3(IELEM) = RRR

     IF(FMU.NE.0.D0)THEN
        if(H_RGN > 0.01d0) then
           H_RGN = 2.D0/H_RGN
           TAU_SUNG3 = H_RGN**2.D0/(4.D0*FMU/RHOINF)
           T_SUGN2(IELEM) = (RESUMEN + 1.D0/TAU_SUNG3**2.D0)**(-.5D0)
        else
           T_SUGN2(IELEM) = RESUMEN**(-.5D0)
        end if

        if(H_RGNE > 0.01d0) then
           H_RGNE = 2.D0/H_RGNE
           TAU_SUNG3_E = H_RGNE**2.D0/(4.D0*FMU/RHOINF)
           T_SUGN3(IELEM) = (RESUMEN + 1.D0/TAU_SUNG3_E**2.D0)**(-.5D0)
        else
           T_SUGN3(IELEM) = RESUMEN**(-.5D0)
        end if
     END IF

  END DO

  RETURN
END SUBROUTINE ESTAB


! subroutine calcRHS(U, U_n, rhs, P, RMU, dtl, gamm)
! 	!RMU, FK, P, T_inf, FCV al pedos en Euler
! 	!use MALLOCAR
! 	!use MGEOMETRIA
! 	use Mvelocidades, only: W_x, W_y
! 	use MeshData
! 	use MESTABILIZACION
! 	use InputData, only: cte, T_inf, FK, FCv, FR, printFlag
! 	implicit real(8) (A-H,O-Z)
! 	real(8) U(4,npoin),P(npoin),U_n(4,npoin),T(npoin), dtl(nelem), gamm(npoin)
! 	real(8) rhs(4,npoin)
! 	real(8) A1(3,4,3), A2(3,4,3), Adv(4,3), AA1(4), AA2(4), Adv_phi(4,3), U_loc(4,3), Ux(4), Uy(4), Nx(3), Ny(3)
! 	real(8) tau(4), choq(3), phi_loc(4,3), rhs_tmp(4,3)
! 	real(8) sp(3,3) !Funciones de forma
! 	real(8) c(3), vx(3), vy(3), V_sq(3), temp(3), e(3)
! 	real*8 W_xp(3), W_yp(3)
! 	INTEGER ipoin(3)

! 	sp(:,1) = (/ .5d0, .5d0, 0.d0 /)
! 	sp(:,2) = (/ 0.d0, .5d0, .5d0 /)
! 	sp(:,3) = (/ .5d0, 0.d0, .5d0 /)

! 	!$omp parallel &
! 	!$omp private(ielem,ipoin,gama,temp,FMU,Nx,Ny,Ux,Uy,AR,HLONG,HLONGX,HLONGY,tau,ALFA_MU,i,&
! 	!$omp U_loc,phi_loc,vx,vy,V_sq,e,C,A1,A2,Adv,AA1,AA2,Adv_phi,choq,rhs_tmp)

! 	!$omp do
! 	do ielem = 1,nelem
! 	rhs_tmp = 0.d0
! 	ipoin = inpoel(:, ielem)

! 	gama = (GAMM(ipoin(1))+GAMM(ipoin(2))+GAMM(ipoin(3)))/3.d0
! 	!temp = (T(ipoin(1))+T(ipoin(2))+T(ipoin(3)))/3.d0
! 	!FMU= 1.716D-5*162.6/(temp-110.55)*(temp/273.15)**.75D0     !SUTHERLAND

! 	Nx = dNx(:, ielem)
! 	Ny = dNy(:, ielem)

! 	Ux(:) = U(:,ipoin(1))*Nx(1) + U(:,ipoin(2))*Nx(2) + U(:,ipoin(3))*Nx(3)
! 	Uy(:) = U(:,ipoin(1))*Ny(1) + U(:,ipoin(2))*Ny(2) + U(:,ipoin(3))*Ny(3)

! 	AR = area(ielem)*dtl(ielem)/3.d0
! 	!CCCC  ----> LONG. CARACTERISTICA
! 	HLONGX = HHX(ielem)
! 	HLONGY = HHY(ielem)
! 	HLONG = DSQRT(area(ielem))
! 	!CCCC  ----> ESTAB. TEZDUyAR
! 	tau(1) = T_SUGN1(ielem)
! 	tau(2) = T_SUGN2(ielem)
! 	tau(3) = T_SUGN2(ielem)
! 	tau(4) = T_SUGN3(ielem)
! 	ALFA_MU = SHOC(ielem)
! 	choq(1) = ALFA_MU*AR*CTE
! 	choq(2) = ALFA_MU*AR*CTE
! 	choq(3) = ALFA_MU*AR*CTE

! 	!CCCC  ----> INTEGRO LAS VARIABLES EN LOS PUNTOS DE GAUSS
! 	U_loc(:,1) = sp(1,1)*U(:,ipoin(1)) + sp(2,1)*U(:,ipoin(2)) + sp(3,1)*U(:,ipoin(3))
! 	U_loc(:,2) = sp(1,2)*U(:,ipoin(1)) + sp(2,2)*U(:,ipoin(2)) + sp(3,2)*U(:,ipoin(3))
! 	U_loc(:,3) = sp(1,3)*U(:,ipoin(1)) + sp(2,3)*U(:,ipoin(2)) + sp(3,3)*U(:,ipoin(3))

! 	phi_loc(:,1) = sp(1,1)*U_n(:,ipoin(1)) + sp(2,1)*U_n(:,ipoin(2)) + sp(3,1)*U_n(:,ipoin(3))
! 	phi_loc(:,2) = sp(1,2)*U_n(:,ipoin(1)) + sp(2,2)*U_n(:,ipoin(2)) + sp(3,2)*U_n(:,ipoin(3))
! 	phi_loc(:,3) = sp(1,3)*U_n(:,ipoin(1)) + sp(2,3)*U_n(:,ipoin(2)) + sp(3,3)*U_n(:,ipoin(3))

! 	W_xp(:) = 0.0
! 	W_yp(:) = 0.0
! 	! W_xp(1) = sp(1,1)*W_x(ipoin(1)) + sp(2,1)*W_x(ipoin(2)) + sp(3,1)*W_x(ipoin(3))
! 	! W_xp(2) = sp(1,2)*W_x(ipoin(1)) + sp(2,2)*W_x(ipoin(2)) + sp(3,2)*W_x(ipoin(3))
! 	! W_xp(3) = sp(1,3)*W_x(ipoin(1)) + sp(2,3)*W_x(ipoin(2)) + sp(3,3)*W_x(ipoin(3))
! 	! W_yp(1) = sp(1,1)*W_y(ipoin(1)) + sp(2,1)*W_y(ipoin(2)) + sp(3,1)*W_y(ipoin(3))
! 	! W_yp(2) = sp(1,2)*W_y(ipoin(1)) + sp(2,2)*W_y(ipoin(2)) + sp(3,2)*W_y(ipoin(3))
! 	! W_yp(3) = sp(1,3)*W_y(ipoin(1)) + sp(2,3)*W_y(ipoin(2)) + sp(3,3)*W_y(ipoin(3))

! 	!CCCC  ----> DEFINO VARIABLES PRIMITIVAS
! 	vx(:) = U_loc(2,:)/U_loc(1,:)
! 	vy(:) = U_loc(3,:)/U_loc(1,:)
! 	e(:) = U_loc(4,:)/U_loc(1,:)
! 	V_sq(:) = vx*vx+vy*vy
! 	temp(:) = (gama-1.d0)/FR*(e-.5d0*V_sq) !FR = CTE. UNIVERSAL DE LOS GASES
! 	c(:) = DSQRT(gama*FR*temp(:))

! 	! = = = GAUSS QUAD = ==
! 	do i = 1,3 
! 	!A1(0,:,i) = (/ 0.d0, 1.d0, 0.d0, 0.d0 /)
! 	A1(1,:,i) = (/ (gama-1.d0)/2.d0*V_sq(i)-vx(i)*vx(i),&
! 		(3.d0-gama)*vx(i) - W_xp(i),&
! 		-(gama-1.d0)*vy(i),&
! 		(gama-1.d0) /)
! 	A1(2,:,i) = (/ -vx(i)*vy(i),&
! 		vy(i),&
! 		vx(i) - W_xp(i),&
! 		0.d0 /)
! 	A1(3,:,i) = (/ ((gama-1.d0)*V_sq(i)-gama*e(i))*vx(i),&
! 		gama*e(i)-(gama-1.d0)/2.d0*V_sq(i)-(gama-1.d0)*vx(i)*vx(i),&
! 		-(gama-1.d0)*vx(i)*vy(i),&
! 		gama*vx(i) - W_xp(i)/)
! 	!A2(0,:,i) = (/ 0.d0, 0.d0, 1.d0, 0.d0 /)
! 	A2(1,:,i) = (/ -vx(i)*vy(i),&
! 		vy(i) - W_yp(i),&
! 		vx(i),&
! 		0.d0 /)
! 	A2(2,:,i) = (/ (gama-1.d0)/2.d0*V_sq(i)-vy(i)*vy(i),&
! 		-(gama-1.d0)*vx(i),&
! 		(3.d0-gama)*vy(i) - W_yp(i),&
! 		(gama-1.d0) /)
! 	A2(3,:,i) = (/ ((gama-1.d0)*V_sq(i)-gama*e(i))*vy(i),&
! 		-(gama-1.d0)*vx(i)*vy(i),&
! 		gama*e(i)-(gama-1.d0)/2.d0*V_sq(i)-(gama-1.d0)*vy(i)*vy(i),&
! 		gama*vy(i) - W_yp(i) /)

! 	Adv(1,i) = -W_xp(i)*Ux(1) + Ux(2) + -W_yp(i)*Uy(1) + Uy(3)
! 	Adv(2:4,i) = A1(:,1,i)*Ux(1) + A1(:,2,i)*Ux(2) + A1(:,3,i)*Ux(3) + A1(:,4,i)*Ux(4)&
! 		+A2(:,1,i)*Uy(1) + A2(:,2,i)*Uy(2) + A2(:,3,i)*Uy(3) + A2(:,4,i)*Uy(4)
! 	end do 
! 	! = = = end GAUSS QUAD = ==

! 	Adv_phi = Adv + phi_loc

! 	AA1(1) = (-W_xp(1)*Adv_phi(1,1) + Adv_phi(2,1) -W_xp(2)*Adv_phi(1,2) + Adv_phi(2,2) -W_xp(3)*Adv_phi(1,3) + Adv_phi(2,3))*AR
! 	AA1(2:4) = (A1(:,1,1)*Adv_phi(1,1) + A1(:,2,1)*Adv_phi(2,1) + A1(:,3,1)*Adv_phi(3,1) + A1(:,4,1)*Adv_phi(4,1) + &
! 		A1(:,1,2)*Adv_phi(1,2) + A1(:,2,2)*Adv_phi(2,2) + A1(:,3,2)*Adv_phi(3,2) + A1(:,4,2)*Adv_phi(4,2) + &
! 		A1(:,1,3)*Adv_phi(1,3) + A1(:,2,3)*Adv_phi(2,3) + A1(:,3,3)*Adv_phi(3,3) + A1(:,4,3)*Adv_phi(4,3))*AR
! 	AA2(1) = (-W_yp(1)*Adv_phi(1,1) + Adv_phi(3,1) -W_yp(2)*Adv_phi(1,2) + Adv_phi(3,2) -W_yp(3)*Adv_phi(1,3) + Adv_phi(3,3))*AR
! 	AA2(2:4) = (A2(:,1,1)*Adv_phi(1,1) + A2(:,2,1)*Adv_phi(2,1) + A2(:,3,1)*Adv_phi(3,1) + A2(:,4,1)*Adv_phi(4,1) + &
! 		A2(:,1,2)*Adv_phi(1,2) + A2(:,2,2)*Adv_phi(2,2) + A2(:,3,2)*Adv_phi(3,2) + A2(:,4,2)*Adv_phi(4,2) + &
! 		A2(:,1,3)*Adv_phi(1,3) + A2(:,2,3)*Adv_phi(2,3) + A2(:,3,3)*Adv_phi(3,3) + A2(:,4,3)*Adv_phi(4,3))*AR

! 	rhs_tmp(:,1) = (Nx(1)*AA1 + Ny(1)*AA2)*tau + &
! 		(Adv(:,1)*sp(1,1) + Adv(:,2)*sp(1,2) + Adv(:,3)*sp(1,3))*AR + 3*(Nx(1)*Ux + Ny(1)*Uy)*choq(1)
! 	rhs_tmp(:,2) = (Nx(2)*AA1 + Ny(2)*AA2)*tau + &
! 		(Adv(:,1)*sp(2,1) + Adv(:,2)*sp(2,2) + Adv(:,3)*sp(2,3))*AR + 3*(Nx(2)*Ux + Ny(2)*Uy)*choq(2)
! 	rhs_tmp(:,3) = (Nx(3)*AA1 + Ny(3)*AA2)*tau + &
! 		(Adv(:,1)*sp(3,1) + Adv(:,2)*sp(3,2) + Adv(:,3)*sp(3,3))*AR + 3*(Nx(3)*Ux + Ny(3)*Uy)*choq(3)

! 	do i = 1,3
! 	!$omp atomic
! 	rhs(1,ipoin(i)) = rhs(1,ipoin(i)) + rhs_tmp(1,i)
! 	!$omp atomic
! 	rhs(2,ipoin(i)) = rhs(2,ipoin(i)) + rhs_tmp(2,i)
! 	!$omp atomic
! 	rhs(3,ipoin(i)) = rhs(3,ipoin(i)) + rhs_tmp(3,i)
! 	!$omp atomic
! 	rhs(4,ipoin(i)) = rhs(4,ipoin(i)) + rhs_tmp(4,i)
! 	end do

! 	end do
! 	!$omp end do
! 	!$omp end parallel
! 	return
! end subroutine calcRHS

subroutine fixvel(velocidadx,velocidady)
  !use mallocar
  !use mvariabfix
  use mvelocidades
  use meshdata
  implicit none
  real(8) velocidadx(npoin),velocidady(npoin)
  integer i, j

  !$OMP PARALLEL DO PRIVATE(i, j)
  do i = 1, nfixv
     j = ifixv_node(i)
     velocidadx(j) = rfixv_valuex(i)+w_x(j)
     velocidady(j) = rfixv_valuey(i)+w_y(j)
  end do
  !$OMP END PARALLEL DO
end subroutine FIXVEL

subroutine FIX(velocidadx,velocidady)
  !use MALLOCAR
  !use MVARIABFIX
  ! use MVELOCIDADES
  use Inputdata, only: gama,fr
  use MVARIABLES
  use MeshData
  use VARIMPLICIT, only : sol 
  !use PointNeighbor, only: esup7
  implicit real(8) (A-H,O-Z)
  !real(8) GAMM(npoin)
  real(8) velocidadx(npoin),velocidady(npoin)
  !$OMP PARALLEL
  !$OMP DO PRIVATE(i)
  do i = 1, nfixrho
    sol((ifixrho_node(i)-1)*4+1) = rfixrho_value(i)
  end do
  !$OMP END DO

  !$OMP DO PRIVATE(i, j, GM)
  do i = 1, NFIXT
     j = IFIXT_NODE(i)
     GM = GAMa - 1.d0
     T(j) = RFIXT_VALUE(i)
     sol((j)*4) = (T(j)*FR/GM + .5d0*(VELocidadX(j)**2 + VELocidadY(j)**2))*sol((j-1)*4+1)
  end do

  !$OMP END DO
  !$OMP END PARALLEL
end subroutine FIX



subroutine setcondition

  use varimplicit, only: sol
  use MeshData, only: npoin 
  use Mnormales
  real(8) velocidadx(npoin),velocidady(npoin)


  !SEPARA RHO*U y RHO*V
  do i=1,npoin
     velocidadx(i)=sol((i-1)*4+2)/sol((i-1)*4+1)
     velocidady(i)=sol((i-1)*4+3)/sol((i-1)*4+1)
  end do

  !SETEAR CONDICIONES DENTRO DEL PGMRES
  !CCCC---------------------------------------CCCC
  !CCCC  ----> CONDICIONES DE CONTORNO <----  CCCC
  !CCCC---------------------------------------CCCC
  !CCCC----> VELOCIDADES IMPUESTAS
  !CCCC---------------------------
  call FIXVEL(velocidadx,velocidady)

  !CCCC----> CORRECCION DE LAS VELOCIDADES NORMALES
  !CCCC--------------------------------------------
  call NORMALVEL(velocidadx,velocidady) 

  !CCCC----> VALORES IMPUESTOS
  !CCCC-----------------------
  call FIX(velocidadx,velocidady)


  !JUNTAR RHO*U y RHO*V
  do i=1,npoin
     sol((i-1)*4+2)=velocidadx(i)*sol((i-1)*4+1)
     sol((i-1)*4+3)=velocidady(i)*sol((i-1)*4+1)
  end do

end subroutine setcondition


subroutine RK(DTMIN, NRK, BANDERA, GAMM, dtl)
  !use DATOS_REFINAMIENTO
  !use DATOS_ENTRADA
  !use MALLOCAR
  !use MVARIABFIX
  !use MGEOMETRIA
  !use MFUERZAS
  !use MMOVIMIENTO
  use calcRHS_mod
  use InputData
  use MNORMALES
  use MVELOCIDADES
  use MVARIABGEN
  use MeshData
  use MVARIABLES
  use MESTABILIZACION
  use TIMERS
  use PointNeighbor, only: esup1,esup2,esup3,esup4
  use varimplicit
  use implicit
  implicit real(8)(A-H,O-Z)
  integer BANDERA,dataval,niteracion
  real(8) GAMM(npoin),converg1,epsilon,converg
  real(8) dtl(nelem),ukold(4,npoin),convergloc
  !save uold
  !if (.not. allocated(uold)) allocate(uold(4,npoin))
  ! dataval=npoin*4
  ! DATA ((uold(I,J), I=1,4), J=1,npoin) /dataval*0./ 
  ! do IRK = 1,NRK
  !RK_FACT = 1.d0/(NRK + 1 - IRK)

  if (.not.allocated(alu)) allocate (alu(size(esup3,1)*5000))
  if (.not.allocated(jlu)) allocate (jlu(size(esup3,1)*5000))
  if (.not.allocated(ju)) allocate (ju(npoin*4))
  if (.not.allocated(vv)) allocate (vv(npoin*4,im+1))


!print*,uold

  converg1=10
  epsilon=1d-10

  niteracion=0

  DO while (converg1.gt.epsilon.and.niteracion.lt.30)
     ukold=U1
     vv=0.d0 
     jlu=0
     ju=0
     alu=0.d0
     converg=0.d0
     convergloc=0.d0
     niteracion=niteracion+1
     !Setear tolerancia para picard o max numero de iteraciones


     !CCCC  ----> SOLO CALCULO UNA VEZ EL TERMINO DE ESTABILIZACION
     !if(IRK.EQ.1)THEN
     !DEBUGGG
     timer(cuarto_orden(U1,UN,FR,gamm), cuarto_t)

     timer(estab(U1,T,GAMA,FR,fMU,DTMIN,RHO_inf,T_inf,U_inf,V_inf,GAMM), estab_t)



     timer(ncalcRHS( U1, UN, dNx, dNy, area, shoc, dtl, t_sugn1, t_sugn2, t_sugn3, nelem, npoin,ukold), calcrhs_t)
     ! timer(calcRHS(U, UN, RHS, P, RMU, dtl, gamm), calcrhs_t)
     ! call rhs_diffusive&
     ! (rhs, U, dnx, dny, area, dtl, T, gamm, fr, fmu, fk, fcv, T_inf, inpoel, npoin, nelem)

     !CCCC  ----> CALCULO DE LOS TERMINOS FUENTES
     !DEBUGGG
     !timer(FUENTE(dtl), fuente_t)


!!$ LLAMADO A SUBRUTINA DE GMRES

     call intimpli

!!!!!PASAR DE SOL A U1
     call setcondition
!!$       do i=1,28
!!$          print*,i
!!$          do j=1,4
!!$             print*,sol((i-1)*4+j)
!!$          end do
!!$       end do
!!$       pause
     do i=1,npoin
        do k=1,4
           u1(k,i)=sol((i-1)*4+k)
        end do
     end do
     !$OMP END PARALLEL DO

     !CCCC----------------------------------------------
     !CCCC----> PASA A LA VARIABLE PRIMARIA PARA APLICAR
     !CCCC----> LAS CONDICIONES DE CONTORNO								
     !CCCC----------------------------------------------

     if(NGAS.EQ.1) GO TO 112
     !$OMP PARALLEL DO PRIVATE(IPOIN, VEL2)
     do ipoin = 1, npoin
        RHO(ipoin) = U1(1,ipoin)
        VEL_X(ipoin) = U1(2,ipoin)/RHO(ipoin)
        VEL_Y(ipoin) = U1(3,ipoin)/RHO(ipoin)
        E(ipoin) = U1(4,ipoin)/RHO(ipoin)
        VEL2 = (VEL_X(ipoin)**2 + VEL_Y(ipoin)**2)
        P(ipoin) = RHO(ipoin)*(GAMM(ipoin)-1.d0)*(E(ipoin)-.5d0*VEL2)
        T(ipoin) = P(ipoin)/(RHO(ipoin)*FR)
        RMACH(ipoin) = DSQRT(VEL2/(T(ipoin)*GAMM(ipoin)*FR))
     end do
     !$OMP END PARALLEL DO

112  CONTINUE

     !CCCC----> CASO PARA AIRE EN EQUILIBRIO
     !CCCC----------------------------------------------
     if(NGAS.NE.0)THEN
        do ipoin = 1, npoin
           RHO(ipoin) = U1(1,ipoin)
           E(ipoin) = U1(4,ipoin)/RHO(ipoin)
           VEL_X(ipoin) = U1(2,ipoin)/RHO(ipoin)
           VEL_Y(ipoin) = U1(3,ipoin)/RHO(ipoin)
           VEL2 = VEL_X(ipoin)**2.d0 + VEL_Y(ipoin)**2.d0
           P(ipoin) = RHO(ipoin)*(GAMM(ipoin)-1.d0)*(E(ipoin)-.5d0*VEL2)
           T(ipoin) = P(ipoin)/(RHO(ipoin)*FR)
           if(IRK.EQ.NRK)THEN
              call TGAS(E(ipoin)-.5d0*VEL2,RHO(ipoin),PGAS,AGAS,TGASi,GAMI)
              GAMM(ipoin) = GAMI
              P(ipoin) = PGAS
              T(ipoin) = TGASi
              RMACH(ipoin) = DSQRT(VEL2)/AGAS
           end if
        end do
     end if



     !SETEAR CONVERGENCIA PARA PICARD
     do ipoin = 1, npoin
        convergloc=convergloc+(u1(4,ipoin)-ukold(4,ipoin))**2.d0
        converg=converg+(u1(4,ipoin))**2.d0
     end do

     converg1=dsqrt(convergloc/converg)
     !print*,converg1
     ukold=u1
  end do

  uold=u1
  !sprint*,uold
end subroutine RK

subroutine ADAMSB(DTMIN, NESTAB, GAMM, dtl)
	!use DATOS_REFINAMIENTO
	!use DATOS_ENTRADA
	!use MGEOMETRIA
	use calcRHS_mod
	use InputData
	use MVELOCIDADES
	use MVARIABGEN
	use MeshData
	use MVARIABLES
	use MESTABILIZACION
	use TIMERS
	use Mnormales
	implicit real(8) (A-H,O-Z)
	integer NESTAB
	real(8) GAMM(npoin)
	real(8) dtl(nelem)

	!CCCC  ----> SOLO CALCULO UNA VEZ EL TERMINO DE ESTABILIZACION
	if (NESTAB.EQ.4) NESTAB = 1
	if (NESTAB.EQ.2) THEN
		timer(cuarto_orden(U1,UN,FR,gamm), cuarto_t)
		timer(estab(U1,T,GAMA,FR,RMU,DTMIN,RHO_inf,T_inf,U_inf,V_inf, GAMM), estab_t)
	end if
	NESTAB = NESTAB + 1

	!$OMP PARALLEL DO PRIVATE(ipoin)
	do ipoin = 1, npoin
	RHS(:, ipoin) = 0.d0
	end do
	!$OMP END PARALLEL DO

!	timer( calcRHS(rhs, U, UN, dNx, dNy, area, shoc, dtl, t_sugn1, t_sugn2, t_sugn3, inpoel, nelem, npoin), calcrhs_t)
	! timer(calcRHS(U, UN, RHS, P, RMU, dtl, gamm), calcrhs_t)
	! call rhs_diffusive&
	! (rhs, U, dnx, dny, area, dtl, T, gamm, fr, fmu, fk, fcv, T_inf, inpoel, npoin, nelem)
	!CCCC  ----> CALCULO DE LOS TERMINOS FUENTES
	timer(FUENTE(dtl), fuente_t)
	!	call FUENTE(dtl)

	!CCCC ----> INTEGRADOR TEMPORAL
	!$OMP PARALLEL DO PRIVATE(ipoin, RL)
	do ipoin = 1, npoin
	RL = 24.d0*M(ipoin)
	U1(:, ipoin) = U(:, ipoin) - (55.d0*RHS(:, ipoin) - 59.d0*RHS1(:, ipoin) + 37.d0*RHS2(:, ipoin) - 9.d0*RHS3(:, ipoin))/RL
	end do
	!$OMP END PARALLEL DO

	!$OMP PARALLEL DO PRIVATE(ipoin)
	do ipoin = 1, npoin
	RHS3(:, ipoin) = RHS2(:, ipoin)
	RHS2(:, ipoin) = RHS1(:, ipoin)
	RHS1(:, ipoin) = RHS(:, ipoin)
	end do
	!$OMP END PARALLEL DO

	!CCCC----------------------------------------------
	!CCCC----> PASA A LA VARIABLE PRIMARIA PARA APLICAR
	!CCCC----> LAS CONDICIONES DE CONTORNO
	!CCCC----------------------------------------------
	if(NGAS.EQ.1) GO TO 111
	!$OMP PARALLEL DO PRIVATE(ipoin, vel2)
	do ipoin = 1, npoin
	RHO(ipoin) = U1(1,ipoin)
	VEL_X(ipoin) = U1(2,ipoin)/RHO(ipoin)
	VEL_Y(ipoin) = U1(3,ipoin)/RHO(ipoin)
	E(ipoin) = U1(4,ipoin)/RHO(ipoin)
	VEL2 = (VEL_X(ipoin)**2.d0+VEL_Y(ipoin)**2.d0)
	P(ipoin) = RHO(ipoin)*(GAMM(ipoin)-1.d0)*(E(ipoin)-.5d0*VEL2)
	T(ipoin) = P(ipoin)/(RHO(ipoin)*FR)
	RMACH(ipoin) = DSQRT(VEL2/(T(ipoin)*GAMM(ipoin)*FR))
	end do
	!$OMP END PARALLEL DO
	111 CONTINUE

	!CCCC----> CASO PARA AIRE EN EQUILIBRIO
	!CCCC----------------------------------------------
	if(NGAS.NE.0)THEN
		do ipoin = 1, npoin
		RHO(ipoin) = U1(1,ipoin)
		E(ipoin) = U1(4,ipoin)/RHO(ipoin)
		VEL_X(ipoin) = U1(2,ipoin)/RHO(ipoin)
		VEL_Y(ipoin) = U1(3,ipoin)/RHO(ipoin)
		VEL2 = VEL_X(ipoin)**2.d0+VEL_Y(ipoin)**2.d0
		P(ipoin) = RHO(ipoin)*(GAMM(ipoin)-1.d0)*(E(ipoin)-.5d0*VEL2)
		T(ipoin) = P(ipoin)/(RHO(ipoin)*FR)

		call TGAS(E(ipoin)-.5d0*VEL2,RHO(ipoin),PGAS,AGAS,TGASi,GAMI)
		GAMM(ipoin) = GAMI
		P(ipoin) = PGAS
		T(ipoin) = TGASi
		RMACH(ipoin) = DSQRT(VEL2)/AGAS

		end do
	end if

	!CCCC---------------------------------------CCCC
	!CCCC  ----> CONDICIONES DE CONTORNO <----  CCCC
	!CCCC---------------------------------------CCCC

	!CCCC----> VELOCIDADES IMPUESTAS
	!CCCC---------------------------
!	call FIXVEL

	!CCCC----> CORRECCION DE LAS VELOCIDADES NORMALES
	!CCCC--------------------------------------------
!	call NORMALVEL

	!CCCC----> VALORES IMPUESTOS
	!CCCC-----------------------
!	call FIX(FR,GAMM)

	!$OMP PARALLEL DO PRIVATE(ipoin)
	do ipoin = 1, npoin
	U1(1,ipoin) = RHO(ipoin)
	U1(2,ipoin) = VEL_X(ipoin)*RHO(ipoin)
	U1(3,ipoin) = VEL_Y(ipoin)*RHO(ipoin)
	U1(4,ipoin) = E(ipoin)*RHO(ipoin)
	end do
	!OMP END PARALLEL
end subroutine ADAMSB

subroutine FUENTE(dtl)
	!use MALLOCAR
	!use MGEOMETRIA
	!use DATOS_ENTRADA
	use MeshData
	use MVELOCIDADES
	use MVARIABGEN
	use InputData

	implicit real(8) (A-H, O-Z)

	real(8) rhs_tmp(4, 3), Ux(4), Uy(4), Nx(3), Ny(3), dtl(nelem)
	real(8) sp(3, 3)
	real(8) wx(3), wy(3)
	integer ipoin(3)

	sp(:, 1) = (/ .5d0, .5d0, 0.d0 /)
	sp(:, 2) = (/ 0.d0, .5d0, .5d0 /)
	sp(:, 3) = (/ .5d0, 0.d0, .5d0 /)

	!$OMP PARALLEL DO &
	!$OMP PRIVATE(ielem, Ux, Uy, AR, wx, wy, rhs_tmp, ipoin, i)

	do ielem = 1, nelem
	rhs_tmp = 0.d0
	ipoin = inpoel(:,ielem)

	Ux = U(:,ipoin(1))*dNx(1,ielem) + U(:,ipoin(2))*dNx(2,ielem) + U(:,ipoin(3))*dNx(3,ielem)
	Uy = U(:,ipoin(1))*dNy(1,ielem) + U(:,ipoin(2))*dNy(2,ielem) + U(:,ipoin(3))*dNy(3,ielem)

	AR = AREA(ielem)*dtl(ielem)/3.d0

	wx(1) = sum(sp(:,1)*w_x(ipoin(:)))
	wx(2) = sum(sp(:,2)*w_x(ipoin(:)))
	wx(3) = sum(sp(:,3)*w_x(ipoin(:)))
	wy(1) = sum(sp(:,1)*w_y(ipoin(:)))
	wy(2) = sum(sp(:,2)*w_y(ipoin(:)))
	wy(3) = sum(sp(:,3)*w_y(ipoin(:)))

	rhs_tmp(:,1) = - AR*(sp(1,1)*(Ux*wx(1) + Uy*wy(1)) + sp(1,2)*(Ux*wx(2) + Uy*wy(2)) + sp(1,3)*(Ux*wx(3) + Uy*wy(3)))
	rhs_tmp(:,2) = - AR*(sp(2,1)*(Ux*wx(1) + Uy*wy(1)) + sp(2,2)*(Ux*wx(2) + Uy*wy(2)) + sp(2,3)*(Ux*wx(3) + Uy*wy(3)))
	rhs_tmp(:,3) = - AR*(sp(3,1)*(Ux*wx(1) + Uy*wy(1)) + sp(3,2)*(Ux*wx(2) + Uy*wy(2)) + sp(3,3)*(Ux*wx(3) + Uy*wy(3)))

	do i = 1, 3
	!$OMP ATOMIC
	rhs(1,ipoin(i)) = rhs(1,ipoin(i)) + rhs_tmp(1,i)
	!$OMP ATOMIC
	rhs(2,ipoin(i)) = rhs(2,ipoin(i)) + rhs_tmp(2,i)
	!$OMP ATOMIC
	rhs(3,ipoin(i)) = rhs(3,ipoin(i)) + rhs_tmp(3,i)
	!$OMP ATOMIC
	rhs(4,ipoin(i)) = rhs(4,ipoin(i)) + rhs_tmp(4,i)
	end do

	end do
	!$OMP END PARALLEL DO
end subroutine FUENTE

real*8 function cmpMtx(m, n, a, b)
	implicit none
	integer m, n, i, j
	real*8 a(m,n), b(m,n)
	real*8, allocatable, dimension(:,:) :: c
	real*8 res

	allocate(c(m,n))
	c = a - b
	res = 0.d0
	do i = 1, m
	do j = 1, n
	res = res + c(i,j)**2
	end do
	end do
	cmpMtx = res
	deallocate(c)
end function cmpMtx

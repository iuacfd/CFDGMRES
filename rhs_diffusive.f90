subroutine rhs_diffusive&
(rhs, U, dnx, dny, area, dtl, T, gamm, fr, fmur, fkr, fcv, Tinf, inpoel, npoin, nelem)
	implicit none
	integer nelem, npoin, inpoel(3,nelem)
	real*8 rhs(4,npoin), U(4,npoin), dnx(3,nelem), dny(3,nelem)
	real*8 gamm(npoin), T(npoin), area(nelem), dtl(nelem)
	real*8 fr, fmur, fkr, fcv, Tinf,sk
	real*8 sp(3,3), Nx(3), Ny(3), Ux(4), Uy(4), U_p(4)
	real*8 RK11_21, RK11_22, RK11_31, RK11_33, RK11_41, RK11_42, RK11_43, RK11_44
	real*8 RK12_21, RK12_23, RK12_31, RK12_32, RK12_41, RK12_42, RK12_43 
	real*8 RK21_21, RK21_23, RK21_31, RK21_32, RK21_41, RK21_42, RK21_43
	real*8 RK22_21, RK22_22, RK22_31, RK22_33, RK22_41, RK22_42, RK22_43, RK22_44 
	real*8 gama, gm, ar, T1, FMU, FK, VX, VY, ET, RMOD2, TEMP
	real*8 rmu, rfk, c,smu
	integer p, ipoi1, ipoi2, ipoi3
	integer ielem
	real*8 fmu1, fmu23, fmu43, ru1
	real*8 vv2, vv3, vv4, vv6, vv7, vv8

	sp(:,1) = (/ .5d0, .5d0, 0.d0 /)
	sp(:,2) = (/ 0.d0, .5d0, .5d0 /)
    sp(:,3) = (/ .5d0, 0.d0, .5d0 /)

	do ielem = 1, nelem
		Nx = dNx(:,ielem)
		Ny = dNy(:,ielem)
		ipoi1 = inpoel(1,ielem)
		ipoi2 = inpoel(2,ielem)
		ipoi3 = inpoel(3,ielem)
        gama = (gamm(ipoi1) + gamm(ipoi2) + gamm(ipoi3))/3.d0
		gm = gama - 1.d0
        AR = area(ielem)*dtl(ielem)/3.d0
		Ux(:) = U(:,ipoi1)*Nx(1) + U(:,ipoi2)*Nx(2) + U(:,ipoi3)*Nx(3)
		Uy(:) = U(:,ipoi1)*Ny(1) + U(:,ipoi2)*Ny(2) + U(:,ipoi3)*Ny(3)
		T1 = (T(ipoi1) + T(ipoi2) + T(ipoi3))/3.D0
        smu=110.d0
	   fmu=0.017d0*(t1/tinf)**1.5d0*(tinf+smu)/(t1+smu)
		!FMU = RMU!*FMUR

		!FK = FKR!*RFK(T1,TINF)
	    sk=194.d0
	    FK=0.001d0*(t1/tinf)**1.5d0*(tinf+sk)/(t1+sk)
		do p = 1, 3
			U_p = sp(1,p)*U(:,ipoi1) + sp(2,p)*U(:,ipoi2) + sp(3,p)*U(:,ipoi3)

			!CCCC  ----> DEFINO VARIABLES PRIMITIVAS
			VX = U_p(2)/U_p(1)
			VY = U_p(3)/U_p(1)
			ET = U_p(4)/U_p(1)
			RMOD2 = VX*VX+VY*VY
			TEMP = GM/FR*(ET - .5D0*RMOD2) !FR=CTE. UNIVERSAL DE LOS GASES
			C = DSQRT(GAMA*FR*TEMP)
			!CCCC  ----> VARIABLES COMPACTAS
			FMU1 = FMU-(FK/FCV)   ! FK=CONDUCTIVIDAD TERMICA; FCV=CALOR ESPECIFICO A VOL. CTE
			FMU43 = 4.d0/3.d0*FMU ! FMU=VISCOSIDAD
			FMU23 = 2.d0/3.d0*FMU
			RU1 = AR/U_p(1)          ! RU1= FACTOR COMUN DEL JACOBIANO DE LAS MATRICES K

			!CCCC  ----> DEFINICION DE LAS MATRICES K11, K12, K21 Y K22

			!CCCC  ----> K11
			RK11_21 = -FMU43*VX
			RK11_22 = FMU43
			RK11_31 = -FMU*VY
			RK11_33 = FMU 
			RK11_41 = FK/FCV*(.5D0*RMOD2 - FCV*TEMP) - FMU*RMOD2 - FMU/3.D0*VX*VX 
			RK11_42 = (FMU/3.D0 + FMU1)*VX
			RK11_43 = FMU1*VY 
			RK11_44 = FK/FCV
			!  !CCCC  ----> K12
			RK12_21 = FMU23*VY 
			RK12_23 = -FMU23
			RK12_31 = -FMU*VX 
			RK12_32 = FMU
			RK12_41 = -FMU/3.D0*VX*VY 
			RK12_42 = FMU*VY 
			RK12_43 = -FMU23*VX
			!  !CCCC  ----> K21
			RK21_21 = -FMU*VY 
			RK21_23 = FMU
			RK21_31 = FMU23*VX 
			RK21_32 = -FMU23
			RK21_41 = -FMU/3.D0*VX*VY 
			RK21_42 = -FMU23*VY 
			RK21_43 = FMU*VX
			!  !CCCC  ----> K22
			RK22_21 = -FMU*VX 
			RK22_22 = FMU
			RK22_31 = -FMU43*VY 
			RK22_33 = FMU43
			RK22_41 = FK/FCV*(.5D0*RMOD2 - FCV*TEMP) - FMU*RMOD2 - FMU/3.D0*VY*VY 
			RK22_42 = FMU1*VX
			RK22_43 = (FMU/3.D0 + FMU1)*VY 
			RK22_44 = FK/FCV

			!CCCC  ----> TERMINOS VISCOSOS
			VV2 = RK11_21*Ux(1) + RK11_22*Ux(2) + RK12_21*Uy(1) + RK12_23*Uy(3)
			VV3 = RK11_31*Ux(1) + RK11_33*Ux(3) + RK12_31*Uy(1) + RK12_32*Uy(2)
			VV4 = RK11_41*Ux(1) + RK11_42*Ux(2) + RK11_43*Ux(3) + RK11_44*Ux(4) + RK12_41*Uy(1) + RK12_42*Uy(2) + RK12_43*Uy(3)

			VV6 = RK21_21*Ux(1) + RK21_23*Ux(3) + RK22_21*Uy(1) + RK22_22*Uy(2)
			VV7 = RK21_31*Ux(1) + RK21_32*Ux(2) + RK22_31*Uy(1) + RK22_33*Uy(3)
			VV8 = RK21_41*Ux(1) + RK21_42*Ux(2) + RK21_43*Ux(3) + RK22_41*Uy(1) + RK22_42*Uy(2) + RK22_43*Uy(3) + RK22_44*Uy(4)
			!  
			!  !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N1
			!  
			RHS(2,ipoi1) = RHS(2,ipoi1) + RU1*(Nx(1)*VV2 + Ny(1)*VV6)
			RHS(3,ipoi1) = RHS(3,ipoi1) + RU1*(Nx(1)*VV3 + Ny(1)*VV7) 
			RHS(4,ipoi1) = RHS(4,ipoi1) + RU1*(Nx(1)*VV4 + Ny(1)*VV8)
			!  !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N2
			!  
			RHS(2,ipoi2) = RHS(2,ipoi2) + RU1*(Nx(2)*VV2 + Ny(2)*VV6)
			RHS(3,ipoi2) = RHS(3,ipoi2) + RU1*(Nx(2)*VV3 + Ny(2)*VV7)
			RHS(4,ipoi2) = RHS(4,ipoi2) + RU1*(Nx(2)*VV4 + Ny(2)*VV8)
			!  !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N3
			!  
			RHS(2,ipoi3) = RHS(2,ipoi3) + RU1*(Nx(3)*VV2 + Ny(3)*VV6)
			RHS(3,ipoi3) = RHS(3,ipoi3) + RU1*(Nx(3)*VV3 + Ny(3)*VV7)
			RHS(4,ipoi3) = RHS(4,ipoi3) + RU1*(Nx(3)*VV4 + Ny(3)*VV8)    
		end do
	end do
	
end subroutine rhs_diffusive



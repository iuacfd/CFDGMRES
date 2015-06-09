 module calcRHS_mod
	implicit none
contains
  subroutine calcRHS(rhs, U, theta, dNx, dNy, area, shoc, dtl, t_sugn1, t_sugn2, t_sugn3, inpoel, nelem, npoin)
    use InputData, only: Cv => FCV, lambda_ref => FK, mu_ref => FMU, gamma0 => gama, T_inf, cte
    use Mvariables, only: T
    real*8, parameter, dimension(3,3) :: N = &
         reshape((/&
         0.d0, .5d0, .5d0, &
         .5d0, 0.d0, .5d0, &
         .5d0, .5d0, 0.d0 &
         /), (/3,3/))


    NGAUSS=3    !PTOS DE GAUSS DONDE VOY A INTEGRAR

    DO IELEM=1,NELEM

       N1=N(1,IELEM) ; N2=N(2,IELEM) ; N3=N(3,IELEM)

       GAMA=(GAMM(N1)+GAMM(N2)+GAMM(N3))/3.D0
       GM=GAMA-1.D0
       TEMP=(T(N1)+T(N2)+T(N3))/3.D0
       ! FMU= 1.716d-5*162.6/(TEMP-110.55)*(TEMP/273.15)**.75D0     !SUTHERLAND

       ! RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
       ! RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)

       RNx(1:3) = dNx(1:3,ielem)
       RNy(1:3) = dNy(1:3,ielem)

       !NO SE USAN LAS DERIVADAS DE LAS VARIABLES SINO DE LAS FUNCIONES DE FORMA

       AR=DTL(IELEM)*AREA(IELEM)/3.D0

       !CCCC  ----> ESTAB. TEZDUYAR
       TAU1=T_SUGN1(IELEM)
       TAU2=T_SUGN2(IELEM)
       TAU3=T_SUGN3(IELEM)
       ALFA_MU=SHOC(IELEM)

       DO J=1,NGAUSS

          !CCCC  ----> INTEGRO LAS VARIABLES EN LOS PUNTOS DE GAUSS
          U_k(1:4) = &
               N(1,k)*U(1:4,ipoi1) +&
               N(2,k)*U(1:4,ipoi2) +&
               N(3,k)*U(1:4,ipoi3)
          theta_k(1:4) = &
               N(1,k)*theta(1:4,ipoi1) +&
               N(2,k)*theta(1:4,ipoi2) +&
               N(3,k)*theta(1:4,ipoi3)
          !CCCC  ----> DEFINO VARIABLES PRIMITIVAS
          rho = U_k(1)
          vx = U_k(2)/rho
          vy = U_k(3)/rho
          et = U_k(4)/rho
          rmod2 = v1**2 + v2**2

          TEMP=GM/FR*(ET-.5D0*RMOD2) !FR=CTE. UNIVERSAL DE LOS GASES

          C=DSQRT(GAMA*FR*TEMP)

          !CCCC  ----> VARIABLES COMPACTAS
          FMU1=FMU-(FK/FCV)   ! FK=CONDUCTIVIDAD TERMICA; FCV=CALOR ESPECIFICO A VOL. CTE
          FMU43=4.D0/3.D0*FMU ! FMU=VISCOSIDAD
          FMU23=2.D0/3.D0*FMU
          RU1=1.d0/U1           ! RU1= FACTOR COMUN DEL JACOBIANO DE LAS MATRICES K

          !CCCC  ----> DEFINICION DE LAS MATRICES A1 Y A2

          !CCCC  ----> A1
          A1_21= GM/2.D0*RMOD2-VX*VX ; A1_22=(3.D0-GAMA)*VX ; A1_23=-GM*VY ; A1_24=GM
          A1_31= -VX*VY ; A1_32=VY ; A1_33=VX
          A1_41=(GM*RMOD2-GAMA*ET)*VX ; A1_42=GAMA*ET-GM/2.D0*RMOD2-GM*VX*VX ; A1_43=-GM*VX*VY ; A1_44=GAMA*VX
          !CCCC  ----> A2
          A2_21=-VX*VY ; A2_22=VY ; A2_23=VX
          A2_31=GM/2.D0*RMOD2-VY*VY ; A2_32=-GM*VX ; A2_33=(3.D0-GAMA)*VY ; A2_34=GM
          A2_41=(GM*RMOD2-GAMA*ET)*VY ; A2_42=-GM*VX*VY ; A2_43=GAMA*ET-GM/2.D0*RMOD2-GM*VY*VY ; A2_44=GAMA*VY 

          !CCCC----------------------------------------
          !CCCC  ----> NO CALCULO LOS TERMINOS VISCOSOS CUANDO MU=0.0
          IF(FMU.GT.0.D0)THEN
             !CCCC  ----> TERMINOS VISCOSOS
             !Matriz de derivadas de funciones de forma de 8*12
             do l=1,3
                do m=1,4

                   nn(m,m+4*(l-1))=rnx(l)
                   nn(m+4,m+4*(l-1))=rny(l)

                end do
             end do
             !CCCC  ----> DEFINICION DE LAS MATRICES K11, K12, K21 Y K22
             !CCCC  ----> K11
             RKMAT(2,1)=-FMU43*VX ; RKMAT(2,2)=FMU43
             RKMAT(3,1)=-FMU*VY ; RKMAT(3,3)=FMU 
             RKMAT(4,1)=FK/FCV*(.5D0*RMOD2-FCV*TEMP)-FMU*RMOD2-FMU/3.D0*VX*VX ; RKMAT(4,2)=(FMU/3.D0+FMU1)*VX
             RKMAT(4,3)=FMU1*VY ; RKMAT(4,4)=FK/FCV
             !CCCC  ----> K12
             RKMAT(2,5)=FMU23*VY ; RKMAT(2,7)=-FMU23
             RKMAT(3,5)=-FMU*VX ; RKMAT(3,6)=FMU
             RKMAT(4,5)=-FMU/3.D0*VX*VY ; RKMAT(4,6)=FMU*VY ; RKMAT(4,7)=-FMU23*VX
             !CCCC  ----> K21
             RKMAT(6,1)=-FMU*VY ; RKMAT(6,3)=FMU
             RKMAT(7,1)=FMU23*VX ; RKMAT(7,2)=-FMU23
             RKMAT(8,1)=-FMU/3.D0*VX*VY ; RKMAT(8,2)=-FMU23*VY ; RKMAT(8,3)=FMU*VX
             !CCCC  ----> K22
             RKMAT(6,5)=-FMU*VX ; RKMAT(6,6)=FMU
             RKMAT(7,5)=-FMU43*VY ; RKMAT(7,7)=FMU43
             RKMAT(8,5)=FK/FCV*(.5D0*RMOD2-FCV*TEMP)-FMU*RMOD2-FMU/3.D0*VY*VY ; RKMAT(8,6)=FMU1*VX
             RKMAT(8,7)=(FMU/3.D0+FMU1)*VY ; RKMAT(8,8)=FK/FCV
             !dgemm= ver https://software.intel.com/es-es/node/520775#AE8380B9-CAC8-4C57-9AF3-2EAAC6ACFC1B
             !Primero calculo rkmat*dn y eso da kdn y luego hago dnt*kdn
             call dgemm('n','n',8,12,8,1.d0,RKMAT,8,nn,8,1.d0,kdn,8)
             call dgemm('t','n',12,12,8,1.d0,nn,8,kdn,12,1.d0,dntkdn,12)

          END IF

          !CCCC-----------------------------------------------------------------------------------------------------------------
          !CCCC----> TERMINOS DE ESTABILIZACION Y CAPTURA DE CHOQUE                        
          ARR1=AR*TAU1
          ARR2=AR*TAU2
          ARR3=AR*TAU3

          !CCCC  ----> Captura de Choque

          CHOQ=alfa_mu*CTE!*ar 
          !CHOQ2=alfa_mu*AR*CTE 
          !CHOQ3=alfa_mu*AR*CTE 
          !CCCC----------------------------------------------------------------------------------------------------
          !CCCC---->************************************<----!CCCC 
          !CCCC---->  CALCULO DE MUa Y SUS COMPONENTES  <----!CCCC
          !CCCC---->************************************<----!CCCC


          !CCCC----------------------------------------------------------------------------------------------------
          !CCCC  ----> MULTIPLICO POR PARTES PARA SIMPLIFICAR EL ASUNTO


          !CCCC  ----> 'A'*dN
          AdN=0.d0
          do l=1,3

             EAUXA12(l)= rnx(l) 
             EAUXA13(l)= rny(1)

             EAUXA21(l)= A1_21*rnx(l) + A2_21*rny(l)  
             EAUXA22(l)= A1_22*rnx(l) + A2_22*rny(l)
             EAUXA23(l)= A1_23*rnx(l) + A2_23*rny(l)
             EAUXA24(l)= A1_24*rnx(l)

             EAUXA31(l)= A1_31*rnx(l) + A2_31*rny(l)   
             EAUXA32(l)= A1_32*rnx(l) + A2_32*rny(l)
             EAUXA33(l)= A1_33*rnx(l) + A2_33*rny(l)
             EAUXA34(l)=                A2_34*rny(l)

             EAUXA41(l)= A1_41*rnx(l) + A2_41*rny(l)   
             EAUXA42(l)= A1_42*rnx(l) + A2_42*rny(l)
             EAUXA43(l)= A1_43*rnx(l) + A2_43*rny(l)
             EAUXA44(l)= A1_44*rnx(l) + A2_44*rny(l)

             AdN(1,2+4*(l-1))=EAUXA12(l)
             AdN(1,3+4*(l-1))=EAUXA13(l)
             AdN(2,1+4*(l-1))=EAUXA21(l)
             AdN(2,2+4*(l-1))=EAUXA22(l)
             AdN(2,3+4*(l-1))=EAUXA23(l)
             AdN(2,4+4*(l-1))=EAUXA24(l)
             AdN(3,1+4*(l-1))=EAUXA31(l)
             AdN(3,2+4*(l-1))=EAUXA32(l)
             AdN(3,3+4*(l-1))=EAUXA33(l)
             AdN(3,4+4*(l-1))=EAUXA34(l)
             AdN(4,1+4*(l-1))=EAUXA41(l)
             AdN(4,2+4*(l-1))=EAUXA42(l)
             AdN(4,3+4*(l-1))=EAUXA43(l)
             AdN(4,4+4*(l-1))=EAUXA44(l)
          end do

          !CCCC  ----> dNt*'A't*'A'*dN = ('A'*dN)t*'A'*dN
          call mkl_set_num_threads(1)
          !Solestab=solucion para el termino de la estabilizacion (12*12)
          !dgemm= ver https://software.intel.com/es-es/node/520775#AE8380B9-CAC8-4C57-9AF3-2EAAC6ACFC1B
          call dgemm('t','n',12,12,4,1.d0,AdN,12,AdN,12,0.d0,Solestab,12 )

          !PARA MULTIPLICAR POR TAU A LAS FILAS
          do m=1,12
             do i=1,3

                Solestab(1+4*(i-1),m)=tau1*Solestab(1+4*(i-1),m)
                Solestab(2+4*(i-1):3+4*(i-1),m)=tau2*Solestab(2+4*(i-1):3+4*(i-1),m)
                Solestab(4+4*(i-1),m)=tau3*Solestab(4+4*(i-1),m)

             end do
          end do

          !CCCC  ----> Nt*'A'*dN
          do l=1,4
             do m=1,3

                nt(l+(m-1)*4,l)=n(m,j)

             end do
          end do

          call dgemm('n','n',12,12,4,1.d0,nt,12,AdN,12,0.d0,SolNtAdN,12 )

          !CCCC  ----> dNt*dN
          do l=1,4
             do m=1,3

                soldNtdN(l+(m-1)*4,l)   = rnx(1)*rnx(m) + rny(1)*rny(m)
                soldNtdN(l+(m-1)*4,4+l) = rnx(2)*rnx(m) + rny(2)*rny(m)
                soldNtdN(l+(m-1)*4,8+l) = rnx(3)*rnx(m) + rny(3)*rny(m)

             end do
          end do

          !CCCC  ----> Nt*N
          do l=1,4
             do m=1,3

                solNtN(l+(m-1)*4,l)   = N(1,j)*N(m,j) 
                solNtN(l+(m-1)*4,4+l) = N(2,j)*N(m,j) 
                solNtN(l+(m-1)*4,8+l) = N(3,j)*N(m,j) 

             end do
          end do

          !SUMA DE TODOS LOS TERMINOS
          LHS=LHS+solNtN*AREA(IELEM)/3.D0+AR*(SolNtAdN+soldNtdN*CHOQ+dntkdn*ru1+Solestab)

          !CALCULAR RHS QUE ES NT*N*U VIEJO+DELTAT*tau*dNt*Ai*TITA (VECTOR)

          !!TRANSPONE AdN y lo multiplica por el vector Theta_k (theta integrados en el elemento)
          call dgemv('t', 4, 12, 1.d0, AdN, 4,theta_k , 1, 1.d0, dNAtheta, 1)
          !CALCULA NtN*Uviejo (da como resultado un vector)
          call dgemv('n', 12, 12, 1.d0, solNtN, 4,Uold , 1, 1.d0, NtNUold, 1)

          do i=1,3
             RHS(1+4*(i-1))= RHS(1+4*(i-1))+NtNUold(1+4*(i-1))*AREA(IELEM)/3.D0+AR*dNAtheta(1+4*(i-1))*tau1
             RHS(2+4*(i-1))= RHS(2+4*(i-1))+NtNUold(1+4*(2-1))*AREA(IELEM)/3.D0+AR*dNAtheta(2+4*(i-1))*tau2
             RHS(3+4*(i-1))= RHS(3+4*(i-1))+NtNUold(1+4*(3-1))*AREA(IELEM)/3.D0+AR*dNAtheta(3+4*(i-1))*tau2
             RHS(4+4*(i-1))= RHS(4+4*(i-1))+NtNUold(1+4*(4-1))*AREA(IELEM)/3.D0+AR*dNAtheta(4+4*(i-1))*tau3
          end do
       END DO

       call lrhsvector(lhs,rhs,ipoi1,ipoi2,ipoi3)

    END DO

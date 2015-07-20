module calcRHS_mod
  implicit none
contains
  subroutine ncalcRHS(U, theta, dNx, dNy, area, shoc, dtl, t_sugn1, t_sugn2, t_sugn3, nelem, npoin,uold)
    use Meshdata, only: inpoel
    use InputData, only: FCV, FK, mu_ref => FMU, gamm => gama, T_inf, cte,fr
    use Mvariables, only: T
    use PointNeighbor, only: esup6,esup7,esup3,esup4
    real*8, parameter, dimension(3,3) :: N = &
         reshape((/&
         0.d0, .5d0, .5d0, &
         .5d0, 0.d0, .5d0, &
         .5d0, .5d0, 0.d0 &
         /), (/3,3/))
    real*8, parameter, dimension(12,12) :: solntn  = &
         reshape((/ 0.5d0, 0.d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, &
         0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, &
         0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, &
         0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.25d0, 0.0d0, &
         0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, &
         0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, &
         0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, &
         0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, &
         0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, &
         0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, &
         0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.25d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0 /), shape(solntn))
    real*8, dimension(3) ::RNx, RNy
    integer(4) ipoi1,ipoi2,ipoi3,nelem,npoin,ngauss,ielem,   j1
    real(8) gama,gm,temp,dnx(3,nelem),dny(3,nelem),uold(4,npoin),uold1(12),u(4,npoin),theta(4,npoin)
    real*8, intent(in), dimension(nelem) :: area, dtl, shoc, t_sugn1, t_sugn2, t_sugn3
    real(8) FMU1, FMU43, FMU23, RU1,rho,vx,vy,et,rmod2,tau1,tau2,tau3,alfa_mu,nn(8,12),fmu,sum
    real(8) A1(4,4),A2(4,4),adnaux(4,12),dntat(12,4), dntataux(12,4),A(4,8)
    real*8, dimension(4) :: U_k, theta_k
    real(8) RKMAT(8,8),kdN(8,12),dNtkdN(12,12),ar,arr1,arr2,arr3,AdN(4,12),solestab(12,12),nt(12,4)
    real(8) SolNtAdN(12,12),soldNtdN(12,12),LHS(12,12),dNAtheta(12),NtNUold(12),rhs(12)
    integer(4) i,j,l,m
    real(8) choq,result1
    integer(4) auxesup4,cantcol,auxcol,p
    real(8) adns(72,72)

    esup6=0.d0
    esup7=0.d0
    open(101, FILE='Matriz_completa', STATUS='UNKNOWN')
    open(102, FILE='RHS', STATUS='UNKNOWN')




    NGAUSS=3    !PTOS DE GAUSS DONDE VOY A INTEGRAR

    DO IELEM=1,NELEM
       LHS=0.d0
       RHS=0.d0
       uold1=0.d0
       nn=0.d0

       ipoi1 = inpoel(1,ielem)
       ipoi2 = inpoel(2,ielem)
       ipoi3 = inpoel(3,ielem)

       do i =1,3       
          do j=1,4
             uold1(j+(i-1)*4)=uold(j,inpoel(i,ielem))
          end do
       end do

       GAMA=gamm!(GAMM(ipoi1)+GAMM(ipoi2)+GAMM(ipoi3))/3.D0
       GM=GAMA-1.D0
       TEMP=(T(ipoi1)+T(ipoi2)+T(ipoi3))/3.D0
       FMU= mu_ref!1.716d-5*162.6/(TEMP-110.55)*(TEMP/273.15)**.75D0     !SUTHERLAND

       ! RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
       ! RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)

       RNx(1:3) = dNx(1:3,ielem)
       RNy(1:3) = dNy(1:3,ielem)

       !NO SE USAN LAS DERIVADAS DE LAS VARIABLES SINO DE LAS FUNCIONES DE FORMA

       AR=DTL(IELEM)*AREA(IELEM)/3.D0
       ! AR=AREA(IELEM)/3.D0*.0002
       !CCCC  ----> ESTAB. TEZDUYAR
       TAU1=T_SUGN1(IELEM)*100.d0
       TAU2=T_SUGN2(IELEM)*100.d0
       TAU3=T_SUGN3(IELEM)*100.d0

       ALFA_MU=SHOC(IELEM)
       !Matriz de derivadas de funciones de forma de 8*12
       i=0
       do l=1,3
          do m=1,4
             nn(m,m+i)=rnx(l)
             nn(m+4,m+i)=rny(l)
          end do
          i=i+4
       end do

       DO J=1,NGAUSS
          Nt=0.d0
          rkmat=0.d0
          KdN=0.d0
          dNtKdN=0.d0
          AdN=0.d0
          solestab=0.d0
          solNtAdN=0.d0
          soldNtdN=0.d0
          dNAtheta=0.d0
          NtNUold=0.d0
          !CCCC  ----> INTEGRO LAS VARIABLES EN LOS PUNTOS DE GAUSS
          U_k(1:4) = &
               N(1,j)*U(1:4,ipoi1) +&
               N(2,j)*U(1:4,ipoi2) +&
               N(3,j)*U(1:4,ipoi3)
          theta_k(1:4) = &
               N(1,j)*theta(1:4,ipoi1) +&
               N(2,j)*theta(1:4,ipoi2) +&
               N(3,j)*theta(1:4,ipoi3)

          !CCCC  ----> DEFINO VARIABLES PRIMITIVAS
          rho = U_k(1)
          vx = U_k(2)/rho
          vy = U_k(3)/rho
          et = U_k(4)/rho
          rmod2 = vx**2 + vy**2

          TEMP=GM/FR*(ET-.5D0*RMOD2) !FR=CTE. UNIVERSAL DE LOS GASES

          !         C=DSQRT(GAMA*FR*TEMP)

          !CCCC  ----> VARIABLES COMPACTAS
          FMU1=FMU-(FK/FCV)   ! FK=CONDUCTIVIDAD TERMICA; FCV=CALOR ESPECIFICO A VOL. CTE
          FMU43=4.D0/3.D0*FMU ! FMU=VISCOSIDAD
          FMU23=2.D0/3.D0*FMU
          RU1=1.d0/U_k(1)          ! RU1= FACTOR COMUN DEL JACOBIANO DE LAS MATRICES K

          !CCCC  ----> DEFINICION DE LAS MATRICES A1 Y A2
          A1=0.d0
          A2=0.d0
          !CCCC  ----> A1
          A1(1,2)=1.d0
          A1(2,1)= GM/2.D0*RMOD2-VX*VX ; A1(2,2)=(3.D0-GAMA)*VX ; A1(2,3)=-GM*VY ; A1(2,4)=GM
          A1(3,1)= -VX*VY ; A1(3,2)=VY ; A1(3,3)=VX
          A1(4,1)=(GM*RMOD2-GAMA*ET)*VX ; A1(4,2)=GAMA*ET-GM/2.D0*RMOD2-GM*VX*VX ; A1(4,3)=-GM*VX*VY ; A1(4,4)=GAMA*VX

          !CCCC  ----> A2
          A2(1,3)=1.d0        
          A2(2,1)=-VX*VY ; A2(2,2)=VY ; A2(2,3)=VX
          A2(3,1)=GM/2.D0*RMOD2-VY*VY ; A2(3,2)=-GM*VX ; A2(3,3)=(3.D0-GAMA)*VY ; A2(3,4)=GM
          A2(4,1)=(GM*RMOD2-GAMA*ET)*VY ; A2(4,2)=-GM*VX*VY ; A2(4,3)=GAMA*ET-GM/2.D0*RMOD2-GM*VY*VY ; A2(4,4)=GAMA*VY 

          A(:,1:4)=a1
          A(:,5:8)=A2


!!$          print*,rnx
!!$          if (ielem.eq.1) then
!!$             do i=1,8
!!$                 write(*,'(12(f7.2))')( nn(i,j1),j1=1,12)
!!$              end do
!!$          end if

          !CCCC----------------------------------------
          !CCCC  ----> NO CALCULO LOS TERMINOS VISCOSOS CUANDO MU=0.0
          IF(FMU.GT.0.D0)THEN
             !CCCC  ----> TERMINOS VISCOSOS


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
             call dgemm('n','n',8,12,8,1.d0,RKMAT,8,nn,8,1.d0,kdN,8)
             call dgemm('t','n',12,12,8,1.d0,nn,8,kdN,12,1.d0,dNtkdN,12)

          END IF

          !CCCC-----------------------------------------------------------------------------------------------------------------
          !CCCC----> TERMINOS DE ESTABILIZACION Y CAPTURA DE CHOQUE                        
          ! ARR1=AR*TAU1
          ! ARR2=AR*TAU2
          ! ARR3=AR*TAU3

          !CCCC  ----> Captura de Choque

          CHOQ=alfa_mu*CTE!*ar 
          !print*,alfa_mu
          !CCCC----------------------------------------------------------------------------------------------------
          !CCCC---->************************************<----!CCCC 
          !CCCC---->  CALCULO DE MUa Y SUS COMPONENTES  <----!CCCC
          !CCCC---->************************************<----!CCCC


          !CCCC----------------------------------------------------------------------------------------------------
          !CCCC  ----> MULTIPLICO POR PARTES PARA SIMPLIFICAR EL ASUNTO


          !CCCC  ----> 'A'*dN
          AdN=0.d0
          AdNAux=0.d0
          dntat=0.d0
          dntataux=0.d0

          call dgemm('n','n',4,12,8,1.d0,A,4,nn,8,0.d0,AdN,4 )

          ! AdN1(1:4,:)=adn
          !aDN1(5:8,12)=ADNAUX
!!$          do i=1,4
!!$             do j1=1,12
!!$                do l=1,4
!!$                   AdN(i,j1)=AdN(i,j1)+ A1(i,l)*nnx(l,j1)
!!$                end do
!!$             end do
!!$          end do
          !  if (ielem.eq.6) then
          !     do i=1,4
          !        write(*,'(12(f15.8))')( adn(i,j1),j1=1,12)
          !     end do
          !     print*,"lala"
          !  end if

          !dNt*At
          call dgemm('t','t',12,4,8,1.d0,nn,8,A1,4,0.d0,dntat,12 )
          ! call dgemm('t','t',12,4,4,1.d0,nn(5:8,:),4,A2,4,0.d0,dntataux,12 )
          ! dntat1(12,1:4)=dntat
          !DNTAT1(12,5:8)=dntataux



          !CCCC  ----> dNt*'A't*'A'*dN = ('A'*dN)t*'A'*dN
          call mkl_set_num_threads(1)
          !Solestab=solucion para el termino de la estabilizacion (12*12)
          !dgemm= ver https://software.intel.com/es-es/node/520775#AE8380B9-CAC8-4C57-9AF3-2EAAC6ACFC1B
          call dgemm('n','n',12,12,4,1.d0,dntat,12,AdN,4,0.d0,Solestab,12 )

          !print*,rnx!
          !PARA MULTIPLICAR POR TAU A LAS FILAS
          do m=1,12
             do i=1,3

                Solestab(1+4*(i-1),m)=tau1*Solestab(1+4*(i-1),m)
                Solestab(2+4*(i-1):3+4*(i-1),m)=tau2*Solestab(2+4*(i-1):3+4*(i-1),m)
                Solestab(4+4*(i-1),m)=tau3*Solestab(4+4*(i-1),m)

             end do
          end do
!!$          do  l = 1, 12
!!$             write(*,'(12(f18.1))') (solestab(l,m), m=1,12)
!!$          end do
!!$

!!$          stop      
          ! print*,solestab
          !CCCC  ----> Nt*'A'*dN
          do l=1,4
             do m=1,3

                nt(l+(m-1)*4,l)=n(m,j)

             end do
          end do
          !    if(ielem.eq.6) then
          !       do  l = 1, 12
          !          write(*,'(12(f6.3))') (nt(l,m), m=1,4)
          !       end do
          !    end if
          !  print*,"lalal"

          call dgemm('n','n',12,12,4,1.d0,nt,12,AdN,4,0.d0,SolNtAdN,12)

          call dgemm('t','n',12,12,8,1.d0,nn,8,nn,8,0.d0,SoldNtdN,12 )

          !print*,soldntdn
          !CCCC  ----> Nt*N

          !call dgemm('n','t',12,12,4,1.d0,nt,12,nt,12,0.d0,SolNtN,12 )

          !SUMA DE TODOS LOS TERMINOS
          ! print*, lhs(1,1),solestab(1,1)
          !  result1= lhs(1,1)
          LHS=LHS+soldNtdN*CHOQ*ar+AR*(Solestab)+SolNtAdN*ar!+dntkdn*ru1)!
          !
          ! write(*,'(4(e15.8))')lhs(1,1),solntn(1,1)*AREA(IELEM)/3.D0,ar*solestab(1,1),result1+solntn(1,1)*AREA(IELEM)/3.D0+ar*solestab(1,1)
          ! do  l = 1, 12
          !   write(*,'(12(f9.2))') (solNtN(l,m), m=1,12)
          !  end do
          !pause



          !CALCULAR RHS QUE ES NT*N*UVIEJO+DELTAT*tau*dNt*Ai*TITA (VECTOR)

          !!TRANSPONE AdN y lo multiplica por el vector Theta_k (theta integrados en el elemento)
          !print*,size(adn),size(theta_k),size(solntn),size(uold)
          call dgemv('t', 4, 12, 1.d0, AdN, 4,theta_k , 1, 0.d0, dNAtheta, 1)
          !do  l = 1, 12
          !   write(*,'((f7.4))') dnatheta(l)
          !end do

          !print*,"lalal"
          !print*,theta_k
          !CALCULA NtN*Uviejo (da como resultado un vector)
          call dgemv('n', 12, 12, 1.d0, solNtN, 12,Uold1 , 1, 0.d0, NtNUold, 1)

          do i=1,3
             RHS(1+4*(i-1))= RHS(1+4*(i-1))+NtNUold(1+4*(i-1))*AREA(IELEM)/3.D0!+AR*dNAtheta(1+4*(i-1))*tau1
             RHS(2+4*(i-1))= RHS(2+4*(i-1))+NtNUold(2+4*(i-1))*AREA(IELEM)/3.D0!+AR*dNAtheta(2+4*(i-1))*tau2
             RHS(3+4*(i-1))= RHS(3+4*(i-1))+NtNUold(3+4*(i-1))*AREA(IELEM)/3.D0!+AR*dNAtheta(3+4*(i-1))*tau2
             RHS(4+4*(i-1))= RHS(4+4*(i-1))+NtNUold(4+4*(i-1))*AREA(IELEM)/3.D0!+AR*dNAtheta(4+4*(i-1))*tau3

          end do
          !print*,area(ielem)
          ! do  l = 1, 12
          !    write(*,'((f19.13))') rhs(l)
          ! end do
          ! print*,"lalal"
       END DO
       LHS=LHS+solNtN*area(ielem)


       !  do  l = 1, 12
       !    write(*,'(12(e15.8))') (lhs(l,m), m=1,12)
       ! end do
       !pause
       ! print*,ielem,"lala"

       call lrhsvector(lhs,rhs,ipoi1,ipoi2,ipoi3,ielem)
       !stop
    END DO
!!$      do i=1,28
!!$         print*,i
!!$         do p=1,4
!!$            print*,esup7((i-1)*4+p)
!!$         end do
!!$      end do
!!$stop
!!$    adns=0.d0
!!$    auxesup4=1
!!$    p=0
!!$    do i=2,npoin*4+1
!!$       cantcol=esup4(i)-auxesup4
!!$       do l=1,cantcol
!!$          p=p+1         
!!$          auxcol=esup3(p)
!!$          adns(i-1,auxcol)= esup6(p)
!!$       end do
!!$       auxesup4=esup4(i)
!!$    end do
!!$    ! call mkl_ddnscsr(job, rowa, cola, adns, 36, esup3, esup4, esup3, info)
!!$
!!$
!!$    do  l = 1, 4*npoin
!!$       write(101,'(72(e20.12))') (adns(l,m), m=1,4*npoin)
!!$    end do
!!$    do i=1,72
!!$       sum=0.d0
!!$       do l=1,72
!!$          sum=sum+adns(i,l)
!!$       end do
!!$       ! print*,sum
!!$    end do
!!$    ! l=0
!!$    do i=1,72
!!$       ! l=l+1
!!$       write(102,'(e20.12)')esup7(i)
!!$       !if (l.eq.4) print*,"lala"
!!$       !if (l.eq.4) l=0
!!$    end do



    !stop
  end subroutine ncalcRHS
!!$
  subroutine lrhsvector(lhs,rhs,ipoi1,ipoi2,ipoi3,ielem)
    use PointNeighbor, only: esup22,esup4,esup5,esup6,esup7
    !PARA PASAR DE LHS LOCAL A VECTOR DE LHS GLOBAL
    !ESUP6 GUARDA LOS VALORES DEL LHS EN FORMA DE VECTOR
    implicit none
    real*8, dimension(12,12):: lhs
    real*8, dimension(12):: rhs   
    integer*4, dimension(3):: nvector
    integer(4) orden,indices(9),indcolum1,indcolum,i,j,k,ielem
    integer(4) ipoi1,ipoi2,ipoi3

    nvector=(/ ipoi1,ipoi2,ipoi3 /)
    !print*,nvector
    !do ielem=1,3
    do i=1,9
       indices(i)=esup5((ielem-1)*9+i)
    end do
    ! print*,indices
    do i=1,3
       indcolum1=esup22(nvector(i))-esup22(1)
       ! print*,indcolum1
       indcolum=esup22(nvector(i)+1)-esup22(nvector(i))
       !print*,indcolum      
       do j=1,4
          do k=1,3
             orden=indcolum1*16+(j-1)*indcolum*4+(indices(k+(i-1)*3)-esup22(nvector(i)))*4
             esup6(orden+1) = esup6(orden+1) + lhs((i-1)*4+j,k*4-3)
             !print*,orden
             esup6(orden+2) = esup6(orden+2) + lhs((i-1)*4+j,k*4-2)
             esup6(orden+3) = esup6(orden+3) + lhs((i-1)*4+j,k*4-1)         
             esup6(orden+4) = esup6(orden+4) + lhs((i-1)*4+j,k*4)        
          end do
       end do
    end do
    !stop

    !ESUP 7 GUARDA EN UN VECTOR LOS VALORES DEL RHS
    do i=1,3
       do j=1,4
          esup7((nvector(i)-1)*4+j) = esup7((nvector(i)-1)*4+j) + rhs(j+(i-1)*4)
       end do
    end do
   ! print*,esup7
   ! stop
  end subroutine lrhsvector


end module calcRHS_mod

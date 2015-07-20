        !COMPILER-GENERATED INTERFACE MODULE: Mon Jun 29 12:32:11 2015
        MODULE ESTAB__genmod
          INTERFACE 
            SUBROUTINE ESTAB(U,T,GAMA,FR,FMU,DTMIN,RHOINF,TINF,UINF,VINF&
     &,GAMM)
              USE MESHDATA
              REAL(KIND=8) :: U(4,NPOIN)
              REAL(KIND=8) :: T(NPOIN)
              REAL(KIND=8) :: GAMA
              REAL(KIND=8) :: FR
              REAL(KIND=8) :: FMU
              REAL(KIND=8) :: DTMIN
              REAL(KIND=8) :: RHOINF
              REAL(KIND=8) :: TINF
              REAL(KIND=8) :: UINF
              REAL(KIND=8) :: VINF
              REAL(KIND=8) :: GAMM(NPOIN)
            END SUBROUTINE ESTAB
          END INTERFACE 
        END MODULE ESTAB__genmod

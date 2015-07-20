        !COMPILER-GENERATED INTERFACE MODULE: Mon Jun 29 12:32:11 2015
        MODULE DELTAT__genmod
          INTERFACE 
            SUBROUTINE DELTAT(DTMIN,DT)
              USE MESHDATA, ONLY :                                      &
     &          INPOEL,                                                 &
     &          NELEM,                                                  &
     &          AREA
              REAL(KIND=8) :: DTMIN
              REAL(KIND=8) :: DT(NELEM)
            END SUBROUTINE DELTAT
          END INTERFACE 
        END MODULE DELTAT__genmod

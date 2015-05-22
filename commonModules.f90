MODULE MVELOCIDADES
  REAL(8), DIMENSION(:), ALLOCATABLE:: VEL_X, VEL_Y, W_X, W_Y
END MODULE MVELOCIDADES

MODULE MVARIABGEN
  REAL(8), DIMENSION(:, :), ALLOCATABLE:: U, U1, U2, RHS, RHS1, RHS2, RHS3, UN
END MODULE MVARIABGEN

MODULE MVARIABLES
  REAL(8), DIMENSION(:), ALLOCATABLE:: P, T, RHO, E, RMACH
END MODULE MVARIABLES

MODULE MESTABILIZACION
  REAL(8) , DIMENSION(:), ALLOCATABLE:: SHOC, T_SUGN1, T_SUGN2, T_SUGN3
END MODULE MESTABILIZACION

MODULE TIMERS
  integer rate, start_t, end_t, sub_start, sub_end, sub_rate, m_start, m_end, m_rate
  real calcrhs_t, cuarto_t, output_t, masas_t, deriv_t, laplace_t, normales_t, forces_t, newmark_t, grad_t, transf_t, estab_t
  real spmv_t, residuo_t, fuente_t, total_t
END MODULE TIMERS


MODULE VARIMPLICIT
  use MeshData, only:  npoin
  integer*4, parameter :: im = 20 !size of krylov subspace
  REAL(8), dimension(:), allocatable:: alu,sol,rhsimplicit
  integer(4), dimension(:), allocatable:: jlu,ju
  real (8) vv(npoin+1,im+1)
END MODULE VARIMPLICIT

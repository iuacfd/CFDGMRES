module implicit
!!!NO INDEXARRRRRRRRR. SE MUEVE TODO LO QUE ESTA EN FORTRAN 77
contains



  subroutine intimpli
    use PointNeighbor, only: esup3,esup4,esup6,esup7
    use MeshData, only: nelem, npoin     
    use varimplicit

    integer(4) n,nz_num 
    real(8) eps,droptol
    integer(4) maxits,iout,lfil,iwk,ierr
    real(8), allocatable :: w(:)
    integer, allocatable :: jw(:)
    n = (npoin*4)
    allocate (w(n+1),jw(2*n))
    nz_num = (esup4(npoin*4+1)-esup4(1))
    eps=1.D-15
    maxits=100	
    iout=0
    lfil=n!poin
    droptol=1d-9
    iwk=n*100

    call ilut(n,esup6,esup3,esup4,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)

    if (ierr.ne.0) then
       print*, "PROBLEMA DE ILUT:", ierr
       ! stop
    end if

    call pgmres(n, im, esup7, sol, vv, eps, maxits, iout, esup6, esup3, esup4, alu, jlu, ju, ierr)

    if (ierr.ne.0)then
       print*, "PROBLEMA DE GMRES:", ierr
       ! stop
    end if

  end subroutine intimpli


!!$c-----------------------------------------------------------------------
  subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout, aa, ja, ia, alu, jlu, ju, ierr)
!!$c-----------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
    integer n, im, maxits, iout, ierr, ja(*), ia(n+1) , jlu(*), ju(n)
    real*8 vv(n,*), rhs(n), sol(n), aa(*), eps , alu(*)
!!$c----------------------------------------------------------------------*
!!$c     *
!!$c     *** ILUT - Preconditioned GMRES ***                  *
!!$c     *
!!$c----------------------------------------------------------------------*
!!$c     This is a simple version of the ILUT preconditioned GMRES algorithm. *
!!$c     The ILUT preconditioner uses a dual strategy for dropping elements   *
!!$c     instead  of the usual level of-fill-in approach. See details in ILUT *
!!$c     subroutine documentation. PGMRES uses the L and U matrices generated *
!!$c     from the subroutine ILUT to precondition the GMRES algorithm.        *
!!$c     The preconditioning is applied to the right. The stopping criterion  *
!!$c     utilized is based simply on reducing the residual norm by epsilon.   *
!!$c     This preconditioning is more reliable than ilu0 but requires more    *
!!$c     storage. It seems to be much less prone to difficulties related to   *
!!$c     strong nonsymmetries in the matrix. We recommend using a nonzero tol *
!!$c     (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
!!$c     lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
!!$c     more reliable the code is. Efficiency may also be much improved.     *
!!$c     Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
!!$c     Gaussian elimination without pivoting.                               *
!!$c     *
!!$c     ILU(0) and MILU(0) are also provided for comparison purposes         *
!!$c     USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
!!$c     then call pgmres.                                                    *
!!$c----------------------------------------------------------------------*
!!$c     Coded by Y. Saad - This version dated May, 7, 1990.                  *
!!$c----------------------------------------------------------------------*
!!$c     parameters                                                           *
!!$c-----------*
!!$c     on entry:                                                            *
!!$c==========*
!!$c     *
!!$c     n     == integer. The dimension of the matrix.                       *
!!$c     im    == size of krylov subspace:  should not exceed 50 in this      *
!!$c     version (can be reset by changing parameter command for     *
!!$c     kmax below)                                                 *
!!$c     rhs   == real vector of length n containing the right hand side.     *
!!$c     Destroyed on return.                                        *
!!$c     sol   == real vector of length n containing an initial guess to the  *
!!$c     solution on input. approximate solution on output           *
!!$c     eps   == tolerance for stopping criterion. process is stopped        *
!!$c     as soon as ( ||.|| is the euclidean norm):                  *
!!$c     || current residual||/||initial residual|| <= eps           *
!!$c     maxits== maximum number of iterations allowed                        *
!!$c     iout  == output unit number number for printing intermediate results *
!!$c     if (iout .le. 0) nothing is printed out.                    *
!!$c     *
!!$c     aa, ja,                                                              *
!!$c     ia    == the input matrix in compressed sparse row format:           *
!!$c     aa(1:nnz)  = nonzero elements of A stored row-wise in order *
!!$c     ja(1:nnz) = corresponding column indices.                   *
!!$c     ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
!!$c     here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
!!$c     *
!!$c     alu,jlu== A matrix stored in Modified Sparse Row format containing   *
!!$c     the L and U factors, as computed by subroutine ilut.       *
!!$c     *
!!$c     ju     == integer array of length n containing the pointers to       *
!!$c     the beginning of each row of U in alu, jlu as computed     *
!!$c     by subroutine ILUT.                                        *
!!$c     *
!!$c     on return:                                                           *
!!$c==========*
!!$c     sol   == contains an approximate solution (upon successful return).  *
!!$c     ierr  == integer. Error message with the following meaning.          *
!!$c     ierr = 0 --> successful return.                             *
!!$c     ierr = 1 --> convergence not achieved in itmax iterations.  *
!!$c     ierr =-1 --> the initial guess seems to be the exact        *
!!$c     solution (initial residual computed was zero)  *
!!$c     *
!!$c----------------------------------------------------------------------*
!!$c     *
!!$c     work arrays:                                                         *
!!$c=============*
!!$c     vv    == work array of length  n x (im+1) (used to store the Arnoli  *
!!$c     basis)                                                      *
!!$c----------------------------------------------------------------------*
!!$c     subroutines called :                                                 *
!!$c     amux   : SPARSKIT routine to do the matrix by vector multiplication  *
!!$c     delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
!!$c     lusol : combined forward and backward solves (Preconditioning ope.) *
!!$c     BLAS1  routines.                                                     *
!!$c----------------------------------------------------------------------*
    parameter (kmax=50)
    real*8 hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
!!$c-------------------------------------------------------------
!!$c     arnoldi size should not exceed kmax=50 in this version..
!!$c     to reset modify paramter kmax accordingly.
!!$c-------------------------------------------------------------
    data epsmac/1.d-16/
    n1 = n + 1
    its = 0
!!$c-------------------------------------------------------------
!!$c     outer loop starts here..
!!$c--------------compute initial residual vector --------------
    call amux (n, sol, vv, aa, ja, ia)
    do 21 j=1,n
       vv(j,1) = rhs(j) - vv(j,1)
21  continue
!!$c-------------------------------------------------------------
20  ro = dnrm2(n, vv, 1)
    if (iout .gt. 0 .and. its .eq. 0) write(iout, 199) its, ro
    if (ro .eq. 0.0d0) goto 999
    t = 1.0d0/ ro
    do 210 j=1, n
    vv(j,1) = vv(j,1)*t
210  continue
     if (its .eq. 0) eps1=eps*ro
!!$c     ** initialize 1-st term  of rhs of hessenberg system..
     rs(1) = ro
     i = 0
4    i=i+1
     its = its + 1
     i1 = i + 1
     call lusol (n, vv(1,i), rhs, alu, jlu, ju)
     call amux (n, rhs, vv(1,i1), aa, ja, ia)
!!$c-----------------------------------------
!!$c     modified gram - schmidt...
!!$c-----------------------------------------
     do 55 j=1, i
     t = ddot(n, vv(1,j),1,vv(1,i1),1)
     hh(j,i) = t
     call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
55   continue
     t = dnrm2(n, vv(1,i1), 1)
     hh(i1,i) = t
     if ( t .eq. 0.0d0) goto 58
     t = 1.0d0/t
     do 57  k=1,n
     vv(k,i1) = vv(k,i1)*t
57   continue
!!$c     
!!$c     done with modified gram schimd and arnoldi step..
!!$c     now  update factorization of hh
!!$c     
58   if (i .eq. 1) goto 121
!!$c--------perfrom previous transformations  on i-th column of h
     do 66 k=2,i
     k1 = k-1
     t = hh(k1,i)
     hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
     hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
66   continue
121  gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
!!$c     
!!$c     if gamma is zero then any small value will do...
!!$c     will affect only residual estimate
!!$c     
     if (gam .eq. 0.0d0) gam = epsmac
!!$c     
!!$c     get  next plane rotation
!!$c     
     c(i) = hh(i,i)/gam
     s(i) = hh(i1,i)/gam
     rs(i1) = -s(i)*rs(i)
     rs(i) =  c(i)*rs(i)
!!$c     
!!$c     detrermine residual norm and test for convergence-
!!$c     
     hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
     ro = abs(rs(i1))
131  format(1h ,2e14.4)
     if (iout .gt. 0)  write(*, 199) its, ro
     if (i .lt. im .and. (ro .gt. eps1))  goto 4
!!$c     
!!$c     now compute solution. first solve upper triangular system.
!!$c     
     rs(i) = rs(i)/hh(i,i)
     do 30 ii=2,i
     k=i-ii+1
     k1 = k+1
     t=rs(k)
     do 40 j=k1,i
     t = t-hh(k,j)*rs(j)
40   continue
     rs(k) = t/hh(k,k)
30   continue
!!$c     
!!$c     form linear combination of v(*,i)'s to get solution
!!$c     
     t = rs(1)
     do 15 k=1, n
     rhs(k) = vv(k,1)*t
15   continue
     do 16 j=2, i
     t = rs(j)
     do 161 k=1, n
     rhs(k) = rhs(k)+t*vv(k,j)
161  continue
16   continue
!!$c     
!!$c     call preconditioner.
!!$c     
  
     call lusol (n, rhs, rhs, alu, jlu, ju)
     do 17 k=1, n
     sol(k) = sol(k) + rhs(k) 
17   continue
!!$c     
!!$c     restart outer loop  when necessary
!!$c     
     !call setcondition
    
     if (ro .le. eps1) goto 990
     if (its .ge. maxits) goto 991
!!$c     
!!$c     else compute residual vector and continue..
!!$c     
     do 24 j=1,i
     jj = i1-j+1
     rs(jj-1) = -s(jj-1)*rs(jj)
     rs(jj) = c(jj-1)*rs(jj)
24   continue
     do 25  j=1,i1
     t = rs(j)
     if (j .eq. 1)  t = t-1.0d0
     call daxpy (n, t, vv(1,j), 1,  vv, 1)
25   continue
199  format('PGMRES its =', i4, ' res. norm =', d20.6)
!!$c     restart outer loop.
     goto 20
990  ierr = 0
     write(*, 199) its, ro
     return
991  ierr = 1
     return
999  continue
     ierr = -1
     return
!!$c-----------------end of pgmres ---------------------------------------
!!$c-----------------------------------------------------------------------
 end subroutine
!!$c-----------------------------------------------------------------------



!!$c----------------------------------------------------------------------c
    subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
!!$c-----------------------------------------------------------------------
    implicit none 
    integer n
    real*8 a(iwk),alu(iwk*15),w(n+1),droptol
    integer ja(iwk),ia(n+1),jlu(iwk*15),ju(n),jw(2*n),lfil,iwk,ierr
!!$c----------------------------------------------------------------------*
!!$c                      *** ILUT preconditioner ***                     *
!!$c      incomplete LU factorization with dual truncation mechanism      *
!!$c----------------------------------------------------------------------*
!!$c     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
!!$c----------------------------------------------------------------------*
!!$c PARAMETERS                                                           
!!$c-----------                                                           
!!$c
!!$c on entry:
!!$c========== 
!!$c n       = integer. The row dimension of the matrix A. The matrix 
!!$c
!!$c a,ja,ia = matrix stored in Compressed Sparse Row format.              
!!$c
!!$c lfil    = integer. The fill-in parameter. Each row of L and each row
!!$c           of U will have a maximum of lfil elements (excluding the 
!!$c           diagonal element). lfil must be .ge. 0.
!!$c           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!!$c           EARLIER VERSIONS. 
!!$c
!!$c droptol = real*8. Sets the threshold for dropping small terms in the
!!$c           factorization. See below for details on dropping strategy.
!!$c
!!$c  
!!$c iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!!$c           are not big enough to store the ILU factorizations, ilut
!!$c           will stop with an error message. 
!!$c
!!$c On return:
!!$c===========
!!$c
!!$c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!!$c           the L and U factors together. The diagonal (stored in
!!$c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!!$c           contains the i-th row of L (excluding the diagonal entry=1)
!!$c           followed by the i-th row of U.
!!$c
!!$c ju      = integer array of length n containing the pointers to
!!$c           the beginning of each row of U in the matrix alu,jlu.
!!$c
!!$c ierr    = integer. Error message with the following meaning.
!!$c           ierr  = 0    --> successful return.
!!$c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!!$c           ierr  = -1   --> Error. input matrix may be wrong.
!!$c                            (The elimination process has generated a
!!$c                            row in L or U whose length is .gt.  n.)
!!$c           ierr  = -2   --> The matrix L overflows the array al.
!!$c           ierr  = -3   --> The matrix U overflows the array alu.
!!$c           ierr  = -4   --> Illegal value for lfil.
!!$c           ierr  = -5   --> zero row encountered.
!!$c
!!$c work arrays:
!!$c=============
!!$c jw      = integer work array of length 2*n.
!!$c w       = real work array of length n+1.
!!$c  
!!$c----------------------------------------------------------------------
!!$c w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
!!$c jw(n+1:2n)  stores nonzero indicators
!!$c 
!!$c Notes:
!!$c ------
!!$c The diagonal elements of the input matrix must be  nonzero (at least
!!$c 'structurally'). 
!!$c
!!$c----------------------------------------------------------------------* 
!!$c---- Dual drop strategy works as follows.                             *
!!$c                                                                      *
!!$c     1) Theresholding in L and U as set by droptol. Any element whose *
!!$c        magnitude is less than some tolerance (relative to the abs    *
!!$c        value of diagonal element in u) is dropped.                   *
!!$c                                                                      *
!!$c     2) Keeping only the largest lfil elements in the i-th row of L   * 
!!$c        and the largest lfil elements in the i-th row of U (excluding *
!!$c        diagonal elements).                                           *
!!$c                                                                      *
!!$c Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
!!$c keeping  the largest  elements in  each row  of L  and U.   Taking   *
!!$c droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
!!$c (however, fill-in is then mpredictible).                             *
!!$c----------------------------------------------------------------------*
!!$c     locals
     integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
     real*8 tnorm, t, abs, s, fact 
     if (lfil .lt. 0) goto 998



!!$c-----------------------------------------------------------------------
!!$c     initialize ju0 (points to next element to be added to alu,jlu)
!!$c     and pointer array.
!!$c-----------------------------------------------------------------------

     ju0 = n+2
     jlu(1) = ju0

!!$c
!!$c     initialize nonzero indicator array. 
!!$c
    do 1 j=1,n
    jw(n+j)  = 0
1   continue
!!$c-----------------------------------------------------------------------
!!$c     beginning of main loop.
!!$c-----------------------------------------------------------------------
     do 500 ii = 1, n
     j1 = ia(ii)
     j2 = ia(ii+1) - 1
     tnorm = 0.0d0
     do 501 k=j1,j2
     tnorm = tnorm+abs(a(k))
501  continue
     if (tnorm .eq. 0.0) goto 999
     tnorm = tnorm/real(j2-j1+1)
!!$c     
!!$c     unpack L-part and U-part of row of A in arrays w 
!!$c     

     lenu = 1
     lenl = 0
     jw(ii) = ii
     w(ii) = 0.0
     jw(n+ii) = ii
!!$c
     do 170  j = j1, j2
     k = ja(j)
     t = a(j)
     if (k .lt. ii) then
     lenl = lenl+1
     jw(lenl) = k
     w(lenl) = t
     jw(n+k) = lenl
     else if (k .eq. ii) then
     w(ii) = t
     else
     lenu = lenu+1
     jpos = ii+lenu-1 
     jw(jpos) = k
     w(jpos) = t
     jw(n+k) = jpos
     endif
170  continue
     jj = 0
     len = 0 
!!$c     
!!$c     eliminate previous rows
!!$c     
150  jj = jj+1
     if (jj .gt. lenl) goto 160
!!$c-----------------------------------------------------------------------
!!$c     in order to do the elimination in the correct order we must select
!!$c     the smallest column index among jw(k), k=jj+1, ..., lenl.
!!$c-----------------------------------------------------------------------
     jrow = jw(jj)
     k = jj
!!$c     
!!$c     determine smallest column index
!!$c     
     do 151 j=jj+1,lenl
     if (jw(j) .lt. jrow) then
     jrow = jw(j)
     k = j
     endif
151  continue
!!$c
     if (k .ne. jj) then
!!$c     exchange in jw
     j = jw(jj)
     jw(jj) = jw(k)
     jw(k) = j
!!$c     exchange in jr
     jw(n+jrow) = jj
     jw(n+j) = k
!!$c     exchange in w
     s = w(jj)
     w(jj) = w(k)
     w(k) = s
                                                         end if
!!$c
!!$c     zero out element in row by setting jw(n+jrow) to zero.
!!$c     
     jw(n+jrow) = 0
!!$c
!!$c     get the multiplier for row to be eliminated (jrow).
!!$c     
     fact = w(jj)*alu(jrow)
     if (abs(fact) .le. droptol) goto 150
!!$c     
!!$c     combine current row and row jrow
!!$c
     do 203 k = ju(jrow), jlu(jrow+1)-1
     s = fact*alu(k)
     j = jlu(k)
     jpos = jw(n+j)
     if (j .ge. ii) then
!!$c     
!!$c     dealing with upper part.
!!$c     
     if (jpos .eq. 0) then
!!$c
!!$c     this is a fill-in element
!!$c     
     lenu = lenu+1
     if (lenu .gt. n) goto 995
     i = ii+lenu-1
     jw(i) = j
     jw(n+j) = i
     w(i) = - s
     else
!!$c
!!$c     this is not a fill-in element 
!!$c
     w(jpos) = w(jpos) - s
     endif
     else
!!$c     
!!$c     dealing  with lower part.
!!$c     
     if (jpos .eq. 0) then
!!$c
!!$c     this is a fill-in element
!!$c     
     lenl = lenl+1
     if (lenl .gt. n) goto 995
     jw(lenl) = j
     jw(n+j) = lenl
     w(lenl) = - s
     else
!!$c     
!!$c     this is not a fill-in element 
!!$c     
     w(jpos) = w(jpos) - s
     endif
     endif
203  continue
!!$c     
!!$c     store this pivot element -- (from left to right -- no danger of
!!$c     overlap with the working elements in L (pivots). 
!!$c     
     len = len+1 
     w(len) = fact
     jw(len)  = jrow
     goto 150
160  continue
!!$c     
!!$c     reset double-pointer to zero (U-part)
!!$c     
     do 308 k=1, lenu
     jw(n+jw(ii+k-1)) = 0
308  continue
!!$c     
!!$c     update L-matrix
!!$c     
     lenl = len 
     len = min0(lenl,lfil)
!!$c     
!!$c     sort by quick-split
!!$c
     call qsplit (w,jw,lenl,len)
!!$c
!!$c     store L-part
!!$c 
     do 204 k=1, len 
     if (ju0 .gt. iwk) goto 996
     alu(ju0) =  w(k)
     jlu(ju0) =  jw(k)
     ju0 = ju0+1
204  continue

!!$c     
!!$c     save pointer to beginning of row ii of U
!!$c     
     ju(ii) = ju0
!!$c
!!$c     update U-matrix -- first apply dropping strategy 
!!$c
     len = 0
     do k=1, lenu-1
     if (abs(w(ii+k)) .gt. droptol*tnorm) then 
     len = len+1
     w(ii+len) = w(ii+k) 
     jw(ii+len) = jw(ii+k) 
     endif
     enddo
     lenu = len+1
     len = min0(lenu,lfil)
!!$c
     call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
!!$c
!!$c     copy
!!$c 
     t = abs(w(ii))
     if (len + ju0 .gt. iwk) goto 997
     do 302 k=ii+1,ii+len-1 
     jlu(ju0) = jw(k)
     alu(ju0) = w(k)
     t = t + abs(w(k) )
     ju0 = ju0+1
302  continue
!!$c     
!!$c     store inverse of diagonal element of u
!!$c     
      if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
!!$c     
      alu(ii) = 1.0d0/ w(ii) 
!!$c     
!!$c     update pointer to beginning of next row of U.
!!$c     
      jlu(ii+1) = ju0
!!$c-----------------------------------------------------------------------
!!$c     end main loop
!!$c-----------------------------------------------------------------------
500   continue
      ierr = 0
      return
!!$c
!!$c     incomprehensible error. Matrix must be wrong.
!!$c     
995   ierr = -1
      return
!!$c     
!!$c     insufficient storage in L.
!!$c     
996   ierr = -2
      return
!!$c     
!!$c     insufficient storage in U.
!!$c     
997   ierr = -3
      return
!!$c     
!!$c     illegal lfil entered.
!!$c     
998   ierr = -4
      return
!!$c     
!!$c     zero row encountered
!!$c     
999   ierr = -5
      return
!!$c----------------end-of-ilut--------------------------------------------
!!$c-----------------------------------------------------------------------
      end subroutine
!!$c----------------------------------------------------------------------


      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
!!$c-----------------------------------------------------------------------
!!$c         A times a vector
!!$c----------------------------------------------------------------------- 
!!$c multiplies a matrix by a vector using the dot product form
!!$c Matrix A is stored in compressed sparse row storage.
!!$c
!!$c on entry:
!!$c----------
!!$c n     = row dimension of A
!!$c x     = real array of length equal to the column dimension of
!!$c         the A matrix.
!!$c a, ja,
!!$c    ia = input matrix in compressed sparse row format.
!!$c
!!$c on return:
!!$c-----------
!!$c y     = real array of length n, containing the product y=Ax
!!$c
!!$c-----------------------------------------------------------------------
!!$c local variables
!!$c
      real*8 t
      integer i, k
!!$c-----------------------------------------------------------------------
      do 100 i = 1,n
!!$c
!!$c     compute the inner product of row i with vector x
!!$c 
      t = 0.0d0
      do 99 k=ia(i), ia(i+1)-1 
      t = t + a(k)*x(ja(k))
99    continue
!!$c
!!$c     store result in y(i) 
!!$c
      y(i) = t
100   continue
!!$c
      return
!!$c---------end-of-amux---------------------------------------------------
!!$c-----------------------------------------------------------------------
      end subroutine
!!$c-----------------------------------------------------------------------

      subroutine lusol(n, y, x, alu, jlu, ju)
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)
!!$c-----------------------------------------------------------------------
!!$c
!!$c This routine solves the system (LU) x = y, 
!!$c given an LU decomposition of a matrix stored in (alu, jlu, ju) 
!!$c modified sparse row format 
!!$c
!!$c-----------------------------------------------------------------------
!!$c on entry:
!!$c n   = dimension of system 
!!$c y   = the right-hand-side vector
!!$c alu, jlu, ju 
!!$c     = the LU matrix as provided from the ILU routines. 
!!$c
!!$c on return
!!$c x   = solution of LU x = y.     
!!$c-----------------------------------------------------------------------
!!$c 
!!$c Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) 
!!$c       will solve the system with rhs x and overwrite the result on x . 
!!$c
!!$c-----------------------------------------------------------------------
!!$c local variables
!!$c
      integer i,k
!!$c
!!$c forward solve
!!$c

      do 40 i = 1, n
      x(i) = y(i)
      do 41 k=jlu(i),ju(i)-1
      x(i) = x(i) - alu(k)* x(jlu(k))
41    continue
40    continue
!!$c
!!$c     backward solve.
!!$c
      do 90 i = n, 1, -1
      do 91 k=ju(i),jlu(i+1)-1
      x(i) = x(i) - alu(k)*x(jlu(k))
91    continue
      x(i) = alu(i)*x(i)
90    continue
!!$c
      return
!!$c----------------end of lusol ------------------------------------------
!!$c-----------------------------------------------------------------------
      end subroutine
!!$c-----------------------------------------------------------------------


      subroutine qsplit(a,ind,n,ncut)
      real*8 a(n)
      integer ind(n), n, ncut
!!$c-----------------------------------------------------------------------
!!$c     does a quick-sort split of a real array.
!!$c     on input a(1:n). is a real array
!!$c     on output a(1:n) is permuted such that its elements satisfy:
!!$c
!!$c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
!!$c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
!!$c
!!$c     ind(1:n) is an integer array which permuted in the same way as a(*).
!!$c-----------------------------------------------------------------------
      real*8 tmp, abskey
      integer itmp, first, last
!!$c-----
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
!!$c
!!$c     outer loop -- while mid .ne. ncut do
!!$c
1     mid = first
      abskey = abs(a(mid))
      do 2 j=first+1, last
      if (abs(a(j)) .gt. abskey) then
      mid = mid+1
!!$c     interchange
      tmp = a(mid)
      itmp = ind(mid)
      a(mid) = a(j)
      ind(mid) = ind(j)
      a(j)  = tmp
      ind(j) = itmp
      endif
2     continue
!!$c
!!$c     interchange
!!$c
      tmp = a(mid)
      a(mid) = a(first)
      a(first)  = tmp
!!$c
      itmp = ind(mid)
      ind(mid) = ind(first)
      ind(first) = itmp
!!$c
!!$c     test for while loop
!!$c
      if (mid .eq. ncut) return
      if (mid .gt. ncut) then
      last = mid-1
      else
      first = mid+1
      endif
      goto 1
!!$c----------------end-of-qsplit------------------------------------------
!!$c-----------------------------------------------------------------------
      end subroutine

      subroutine daxpy(n,da,dx,incx,dy,incy)
!!$c
!!$c     constant times a vector plus a vector.
!!$c     uses unrolled loops for increments equal to one.
!!$c     jack dongarra, linpack, 3/11/78.
!!$c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
!!$c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!!$c
!!$c        code for unequal increments or equal increments
!!$c          not equal to 1
!!$c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
      dy(iy) = dy(iy) + da*dx(ix)
      ix = ix + incx
      iy = iy + incy
10    continue
      return
!!$c
!!$c        code for both increments equal to 1
!!$c
!!$c
!!$c        clean-up loop
!!$c
20    m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
      dy(i) = dy(i) + da*dx(i)
30    continue
      if( n .lt. 4 ) return
40    mp1 = m + 1
      do 50 i = mp1,n,4
      dy(i) = dy(i) + da*dx(i)
      dy(i + 1) = dy(i + 1) + da*dx(i + 1)
      dy(i + 2) = dy(i + 2) + da*dx(i + 2)
      dy(i + 3) = dy(i + 3) + da*dx(i + 3)
50    continue
      return
      end subroutine


                                                                                             
end module implicit

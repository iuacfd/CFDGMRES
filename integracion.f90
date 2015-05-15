PROGRAM INTEGRAL

  REAL(8) T(10000),Fx(10000),RINT,H
  real(8) fy(1000), rm(1000),pi,fint,wint,mint
  INTEGER I,J

  OPEN(1,FILE='DESPLAZAMIENTO',STATUS='old')
  OPEN(2,FILE='Resultado',STATUS='unknown')
  pi=dacos(-1.d0)

  DO I=1,29    
     READ(1,*)T(I),Fx(I),fy(I),rm(I)

  END DO
  CLOSE(1)

  RINT=0.D0
  j=1
  DO while (t(j).lt.1.d0) 
     H=T(J+1)-T(J)
     tmed=(T(J+1)+T(J))/2.d0
     RINT=RINT+H*(Fy(J+1)+Fy(J))/2.D0
     Fint=h*(fy(j+1)+fy(j))*pi*sin(pi*tmed)/2.d0
     mint=h*(rm(j+1)+rm(j))*pi**2.d0*sin(2.d0*pi*tmed)/6.d0
     wint=wint+fint+mint
     j=j+1
  END DO
  WRITE(2,*) RINT,wint

END PROGRAM INTEGRAL

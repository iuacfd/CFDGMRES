module PointNeighbor
  integer, dimension(:), allocatable :: esup1, esup2,esup3,esup4,esup5
  integer, dimension(:), allocatable :: psup1, psup2, lpoin
  real(8), dimension(:), allocatable :: esup6, esup7
contains
  subroutine getEsup(inpoel,nelem,npoin)
   ! use implicit
    implicit none
    integer istor
    integer nelem, npoin
    integer inpoel(3,nelem) 
    integer ielem, ipoi1, ipoin, i
    allocate(esup2(npoin+1),esup7(npoin*4))
    esup2 = 0

    !contar la cantidad de elementos vecinos de cada nodo (*histogram pattern*)
    do ielem=1,nelem
       do i=1,3
          ipoi1 = inpoel(i,ielem) + 1
          esup2(ipoi1) = esup2(ipoi1) + 1
       end do
    end do

    !reshuffle para obtener esup2 (*scan pattern*)
    do ipoin = 2,npoin+1
       esup2(ipoin) = esup2(ipoin) + esup2(ipoin-1)
    end do

    allocate(esup1(esup2(npoin+1)),esup4(npoin*4+1)) 
    call esup4aux
    allocate (esup3((esup4(npoin*4+1)-esup4(1))*4),esup5(9*nelem))
    allocate(esup6(size(esup3,1)))


    !allocate(esup2l((npoin+1)*4),esup1l(esup2(npoin+1)*4))
    !obtener el array esup1 de los elementos vecinos
    do ielem=1,nelem
       do i=1,3
          ipoin = inpoel(i,ielem)
          istor = esup2(ipoin) + 1 
          esup2(ipoin) = istor     
          esup1(istor) = ielem
       end do
    end do

    call esup3aux

    !restaurar esup2
    do ipoin=npoin+1,2,-1
       esup2(ipoin) = esup2(ipoin-1)
    end do
    esup2(1) = 0

    !do i=1,4
    !   esup2l(1+(n+1)*(i-1):(n+1)*i)=esup2
    !   esup1l(1+n*(i-1):n*i)=esup1
    !end do

    call esup5aux

  end subroutine getEsup

  subroutine getPsup(inpoel, nelem, npoin)
    implicit none
    integer istor
    integer nelem, npoin
    integer inpoel(3,nelem)
    integer ielem, jpoin, ipoin, iesup, i

    if(.not.allocated(esup1)) call getEsup(inpoel, nelem, npoin)

    allocate(psup2(npoin + 1))
    allocate(lpoin(npoin))
    lpoin = 0; psup2(1) = 0; istor = 0

    !calcular total de nodos vecinos
    do ipoin=1,npoin
       do iesup=esup2(ipoin)+1, esup2(ipoin+1)
          ielem=esup1(iesup)
          do i=1,3
             jpoin=inpoel(i,ielem)
             if(jpoin /= ipoin .and. lpoin(jpoin) /= ipoin) then
                istor = istor + 1
                lpoin(jpoin) = ipoin
             end if
          enddo
       end do
       psup2(ipoin+1) = istor
    end do

    allocate(psup1(istor))

    !obtener array psup1
    lpoin=0; istor=0
    do ipoin=1,npoin
       do iesup=esup2(ipoin)+1, esup2(ipoin+1)
          ielem=esup1(iesup)
          do i=1,3
             jpoin=inpoel(i,ielem)
             if(jpoin /= ipoin .and. lpoin(jpoin) /= ipoin) then
                istor = istor+1
                psup1(istor) = jpoin
                lpoin(jpoin) = ipoin
             end if
          enddo
       end do
    end do

    deallocate(lpoin)
  end subroutine getPsup



  subroutine esup4aux
    !se llama en PointNeighbor.f90
    use MeshData, only: npoin
   ! use PointNeighbor, only: esup2,esup4        
    integer(4) ipoin,factor
    !ESUP4 INDICE GLOBAL DE COLUMNAS
    esup2=esup2+1 
    esup4(1)=esup2(1)  
    do ipoin=1,npoin    
       factor=4*(esup2(ipoin+1)-esup2(ipoin))
       esup4(ipoin*4-2) = esup4(ipoin*4-3) + factor
       esup4(ipoin*4-1) = esup4(ipoin*4-2) + factor 
       esup4(ipoin*4)   = esup4(ipoin*4-1) + factor 
       esup4(ipoin*4+1) = esup4(ipoin*4)   + factor
    end do
  end subroutine esup4aux


  subroutine esup3aux
    !se llama en PointNeighbor.f90
    use MeshData, only: npoin
    !use PointNeighbor, only: esup1,esup2,esup3        
    integer(4) p,contador,contador1
    p=0
    contador =0
    do ipoin=1,npoin
       factor=(esup2(ipoin+1)-esup2(ipoin))
       do i=1,4
          do l= 1,factor
             contador=contador+1
             if (i.eq.1.and.l.eq.1) contador1=contador 
             if (l.eq.1) contador=contador1
             p=p+1
             esup3(p) = esup1(contador)*4-3
             p=p+1     
             esup3(p) = esup1(contador)*4-2
             p=p+1     
             esup3(p) = esup1(contador)*4-1
             p=p+1     
             esup3(p) = esup1(contador)*4
          end do
       end do
    end do
  end subroutine esup3aux

  subroutine esup5aux
!!$    !Esup5 guarda las posiciones de los ipoi por elemento teniendo en cuenta los esup2. Sirve para pasar de la matriz local 12*12 a la matriz global.
!!$    !se llama en PointNeighbor.f90
    !use PointNeighbor, only: esup2,esup5
    use MeshData, only: inpoel, nelem
    integer(4) ipoi(3)
    integer(4)p,l,m,i,ielem
    m=0
    do ielem=1,nelem

       ipoi(1)=inpoel(1,ielem); ipoi(2)=inpoel(2,ielem); ipoi(3)=inpoel(3,ielem)
       do i=1,3     !barre nodos de los elementos
          do p=1,3  !barre nodos dentro de un nodo de esup2
             do l=1,esup2(ipoi(i)+1)-esup2(ipoi(i))   !para buscar la posicion   
                if (ipoi(p).eq.esup1(esup2(ipoi(i))+l-1)) then
                   m=m+1
                   esup5(m)=esup2(ipoi(i))+l-1
                   exit
                end if
             end do
          end do
       end do
    end do
  end subroutine esup5aux
  
end module PointNeighbor

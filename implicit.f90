  module implicit

  contains

    subroutine CallAll
      use PointNeighbor
      use MeshData, only: inpoel, nelem, npoin
      implicit none
      call getEsup(inpoel,nelem,npoin)
      call getLocal
    end subroutine CallAll


    subroutine getLocal
      use PointNeighbor, only: esup2, esup1
      use MeshData, only: inpoel, nelem, npoin
      implicit none
      integer, dimension(:), allocatable :: lnum1
      integer ipoin, ielem, iesup, jpoin, i


      !Armar vector con las posiciones locales en un elemento del nodo i
      allocate(lnum1(esup2(npoin+1)))
      do ipoin=1,npoin
         do iesup=esup2(ipoin)+1, esup2(ipoin+1)
            ielem=esup1(iesup)
            do i=1,3
               jpoin=inpoel(i,ielem)
               if(jpoin==ipoin) then
                  lnum1(iesup)=i
               end if
            enddo
         end do
      end do
      
      do i=1, esup2(npoin+1)
         print*,esup1(i),lnum1(i)
      end do
      stop
    end subroutine getLocal

  end module implicit

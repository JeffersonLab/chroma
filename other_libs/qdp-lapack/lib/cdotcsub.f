c-----------------------------------------------------------------------
      subroutine cdotcsub( cdotcr, n, x, incx, y, incy ) 
      complex cdotcr, x(*), y(*)
      integer n, incx, incy
c-----------------------------------------------------------------------
c     interface for zdotc
c----------------------------------------------------------------------- 
      external cdotc
      complex cdotc
c
!        double complex z, a(100), b(100)
!        do i=1,n
!           a(i) = x(i)
!           b(i) = y(i)
!           print *, a(i),b(i)
!        enddo
!        print *, n, a(1)
!        z = zdotc( n, a, 1, b, 1 )
!
!        print *, z
!      print *, n, (x(i),i=1,n)
      cdotcr = cdotc( n, x, incx, y, incy )
!      print *, n, (x(i),i=1,n)
!      print *, cdotcr
c
      return
c-----end-of-cdotcsub---------------------------------------------------
c-----------------------------------------------------------------------
      end

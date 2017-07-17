c-----------------------------------------------------------------------
c CUBIC FUNCTION COEFFICIENTS y = ax^3 + bx^2 + cx + d
c-----------------------------------------------------------------------
      subroutine cubecoeff(x1,y1,m1,x2,y2,m2,a,b,c,d)
      implicit none
      real*8 x1, y1, m1, x2, y2, m2, a, b, c, d, num, den

      den = -(x1-x2)**3.d0

      num = 2.d0 * (y1-y2) - (m1+m2) * (x1-x2)
      a = num / den

      num = -3.d0 * (x1+x2) * (y1-y2) +
     8             m1 * (x1**2.d0 + x1*x2 - 2.d0*x2**2.d0) -
     8             m2 * (x2**2.d0 + x1*x2 - 2.d0*x1**2.d0)
      b = num / den

      num = 6.d0 * x1 * x2 * (y1-y2) +
     8            m1 * x2 * (x2**2.d0 + x1*x2 - 2.d0*x1**2.d0) -
     8            m2 * x1 * (x1**2.d0 + x1*x2 - 2.d0*x2**2.d0)
      c = num / den

      num = x1**2.d0 * x2**2.d0 * (m2-m1) + x1*x2*
     8 ( 3.d0 * (x1**2.d0*m2 - x2**2.d0*m1) + 6.d0 * (x1*y2-x2*y1) -
     8   2.d0 * (x1+x2) * (x1*m2 - x2*m1)) +
     8   2.d0 * (x1**3.d0*y2 - x2**3.d0*y1) -
     8   3.d0 * (x1+x2) * (x1**2.d0*y2 - x2**2.d0*y1)
      d = num / den

      return
      end

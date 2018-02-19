*
* $Id: lfit.F,v 1.1.1.1 1996/04/01 15:02:28 mclareni Exp $
*
* $Log: lfit.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:28  mclareni
* Mathlib gen
*
*
      SUBROUTINE LFIT(X,Y,L,KEY,A,B,E)
      implicit none
C
C     TO FIT A STRAIGHT LINE    Y=A*X+B    TO L POINTS WITH ERROR E
C     SEE MENZEL , FORMULAS OF PHYSICS P.116
C     POINTS WITH Y=0 ARE IGNOERD IF KEY=0
C     L IS NO. OF POINTS
C
      integer*4 L, J
      real*8 X(L),Y(L), KEY
      real*8 COUNT, A, B, E

      print *, 'L = ', L
      print *, 'KEY = ', KEY
      do j = 1, L
         print*, 'x, y = ', X(j), Y(j)
      enddo
      A = 0.
      B = 0.
      E = 0.
      return
      end

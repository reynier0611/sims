      implicit none
      include 'simulate.inc'
      character*255 header
      character*80 name
      real*8 qsq, ge, gm
      real*8 qmu, gep_in, gmp_in
      integer i

      read (5,*) header
      do i = 1,20 
         read(5,*) name, qmu, gep_in, gmp_in
         qsq = -qmu**2/hbarc**2
         call fofa_best_fit(qsq, ge, gm)
         print *, qmu, ge, gm
      enddo
      stop
      end

!-------------------------------------------------------------------------

      subroutine fofa_best_fit(qsquar,GE,GM)
      ! call fofa_bosted(qsquar,GE,GM)
      call fofa_JRA(qsquar,GE,GM)
      return
      end
      
!-------------------------------------------------------------------------
      subroutine fofa_JRA(qsquar,GE,GM)
      implicit none
      
      include 'simulate.inc'
      real*8 qsquar, GE, GM, mu_p, Q2, Q
      real*8 Q22, Q23, Q24, Q25, Q26
      
      Q2 = -qsquar*(hbarc**2.)*1.d-6
      Q  = sqrt(max(Q2,0.d0))
      
      Q22 = Q2**2 
      Q23 = Q2**3 
      Q24 = Q2**4 
      Q25 = Q2**5 
      Q26 = Q2**6 
      
      mu_p = 2.793
C     
C     -- This calculates the form factors Gep and Gmp using
C     John Arrington's fit to world data 
C     reference: Phys. Rev. C69, 022201, 2004
C     
      GM=mu_p/(1.D0
     1     + Q2 * (3.19)
     2     + Q22 * (1.355)
     3     + Q23 * (0.151)
     4     + Q24 * (-0.114E-01)
     5     + Q25 * (0.533E-03)
     6     + Q26 * (-0.900E-05)  )
      GE=1.D0/(1.D0
     1     + Q2 * (3.226)
     2     + Q22 * (1.508)
     3     + Q23 * (-0.3773)
     4     + Q24 * (0.611)
     5     + Q25 * (-0.1853)
     6     + Q26 * (0.1596E-01)  )
      return
      end
      
!-------------------------------------------------------------------------
      
      subroutine fofa_Bosted(qsquar,GE,GM)
      
*     csa 9/14/98 -- This calculates the form factors Gep and Gmp using
*     Peter Bosted's fit to world data (Phys. Rev. C 51, 409, Eqs. 4
*     and 5 or, alternatively, Eqs. 6)
      
      implicit none
      include 'simulate.inc'
      
      real*8  qsquar,GE,GM,mu_p
      real*8  Q,Q2,Q3,Q4,Q5,denom
      
      mu_p = 2.793
      
      Q2 = -qsquar*(hbarc**2.)*1.d-6
      Q  = sqrt(max(Q2,0.d0))
      
      Q3 = Q**3.
      Q4 = Q**4.
      Q5 = Q**5.
      
*     Use Eqs 4, 5:
      denom = 1. + 0.62*Q + 0.68*Q2 + 2.8*Q3 + 0.83*Q4
      GE = 1./denom
      denom = 1. + 0.35*Q + 2.44*Q2 + 0.5*Q3 + 1.04*Q4 + 0.34*Q5
      GM = mu_p/denom
      
*     OR Eqs 6:
*     denom = 1. + 0.14*Q + 3.01*Q2 + 0.02*Q3 + 1.20*Q4 + 0.32*Q5
*     GE = 1./denom
*     GM = mu_p/denom
      
      return
      end

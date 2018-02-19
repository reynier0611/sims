C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    Copyright (c) 2001 Institut des Sciences Nucleaires de Grenoble
C
C    Authors: J. Mougey, E. Voutier                 (Name@isn.in2p3.fr) 
C-----------------------------------------------------------------------
C
C    WARNING !!
C
C    These routines calculate the energy loss of electrons and hadrons
C    for elemental and compound materials, assuming that the  relevant
C    physics parameters are tabulated. In the contrary, the  tables at
C    the end of the section should be implemented using references:
C
C    M.J. Berger and S.M. Seltzer, National Bureau of Standards Report
C    82-2550A (1983).
C    R.M. Sternheimer, M.J. Berger and S.M. Seltzer,Atomic and Nuclear
C    Data Tables 30 (1984) 261.
C
C-----------------------------------------------------------------------
C
C    These routines are based on the written document  
C
C    ESPACE Energy Loss Corrections Revisited
C    J. Mougey, E89-044 Analysis Progress Report, November 2000.
C
C    which PostScript version can be downloaded from the WebSite
C
C    http://isnwww.in2p3.fr/hadrons/helium3/Anal/AnaPag.html
C
C-----------------------------------------------------------------------

      REAL*8 FUNCTION ELOSS_ELECTRON
     +                (beta,z_med,a_med,d_med,t_med)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C Energy loss of electrons taking into account the usual Bethe-Block    
C Stopping power + the Density Effect Corrections that are important
C at high electron energy. Additional  shell effects  can  be safely
C neglected for ultrarelativistic electrons.
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Passed variables
C
C            beta     electron velocity              (relative to light) 
C
C           z_med     effective charge of the medium
C           a_med     effective atomic mass          (AMU)
C           d_med     medium density                 (g/cm^3)
C           t_med     medium thickness               (g/cm^2)
C
C  z_med = sum(i)(w(i)*z(i))/sum(i)w(i)
C            w(i) is the abundacy of element i in the medium
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Function result
C
C  ELOSS_ELECTRON     Energy loss of electron        (MeV)
C
C-----------------------------------------------------------------------

C Passed variables

      REAL*8 z_med,a_med,d_med,t_med
      REAL*8 beta

C Local variables

      REAL*8 BETA2,BETA3,PLAS,EKIN,TAU,FTAU
      REAL*8 BETH,DENS,ESTP
      REAL*8 COEF,EMAS
      REAL*8 EXEN

*---- Constant factor corresponding to 2*pi*N_a*(r_e)^2*m_e*c^2
*     Its dimension is MeV.cm2/g          

      DATA COEF / 0.1535374581D+00 /

*---- Electron mass (MeV)

      DATA EMAS / 0.5109989020D+00 /


*---- Input variables consistency check

      IF((beta.LE.0.D+00) .OR.(beta.GE.1.D+00) .OR.(z_med.EQ.0.D+00).OR.
     +   (a_med.EQ.0.D+00).OR.(t_med.LE.0.D+00)) 
     +   THEN
         ELOSS_ELECTRON = 0.D+00
         GO TO 1
      ENDIF

      BETA2 = beta * beta
      BETA3 = 1.D+00 - BETA2

*---- Reduced Bethe-Block stopping power

      CALL EX_ENERG(z_med,d_med,EXEN)
      IF(EXEN.EQ.0.D+00) THEN
         ELOSS_ELECTRON = 0.D+00
         GO TO 1
      ENDIF
      EKIN = EMAS * ( (1.D+00/DSQRT(BETA3)) - 1.D+00 )
      TAU  = EKIN / EMAS
      FTAU = 1.D+00+(TAU*TAU/8.D+00)-((2.D+00*TAU+1.D+00)*DLOG(2.D+00))
      BETH = 2.D+00 * DLOG(1.D+06*EKIN/EXEN) + BETA3 * FTAU
      BETH = BETH + DLOG(1.D+00+(0.5D+00*TAU))

*---- Reduced density correction

      PLAS = 28.8084D+00 * DSQRT(d_med*z_med/a_med)
      DENS = 2.D+00*DLOG(PLAS)-DLOG(BETA3) - 2.D+00*DLOG(EXEN) - 1.D+00 
      IF(DENS.LT.0.D+00) DENS = 0.D+00

*---- Total electron stopping power

      ESTP = COEF * z_med * ( BETH - DENS ) / a_med / BETA2

*---- Electron energy loss

      ELOSS_ELECTRON = ESTP * t_med

    1 RETURN

      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      REAL*8 FUNCTION ELOSS_HADRON
     +                (z_inc,beta,z_med,a_med,d_med,t_med)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C Energy loss of hadrons taking into account the  usual  Bethe-Block    
C Stopping power + the Density Effect Corrections (important at high
C energy) + the Shell Corrections (important at small energy). 
C The approximation 2 * gamma * m_e / M << 1 is used for the  Bethe-
C Block formulae: for a 4 GeV/c pion, this leads to  an effect about
C 0.5-0.7 % on the total stopping power.
C 
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Passed variables
C
C           z_inc     hadron charge 
C            beta     hadron velocity                (relative to light) 
C
C           z_med     effective charge of the medium
C           a_med     effective atomic mass          (AMU)
C           d_med     medium density                 (g/cm^3)
C           t_med     medium thickness               (g/cm^2)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Function result
C
C  ELOSS_HADRON       Energy loss of hadron          (MeV)
C
C-----------------------------------------------------------------------

C Passed variables

c      INTEGER z_inc
      INTEGER z_inc

      REAL*8 z_med,a_med,d_med,t_med
      REAL*8 beta

C Local variables

      REAL*8 BETA2,GAMA2,GAMA,PLAS,ETA,ETA2,ETA4,ETA6
      REAL*8 BETH,DENS,SHE0,SHE1,SHEL,HSTP
      REAL*8 COEF,EMAS
      REAL*8 X,X0,X1,M
      REAL*8 EXEN,C,A

*---- Constant factor corresponding to 2*pi*N_a*(r_e)^2*m_e*c^2
*     Its dimension is MeV.cm2/g          

      DATA COEF / 0.3070749162D+00 /

*---- Electron mass (MeV)

      DATA EMAS / 0.5109989020D+00 /


*---- Input variables consistency check

      IF((z_inc.EQ.0.D+00).OR.(beta.LE.0.D+00) .OR.(beta.GE.1.D+00) .OR.
     +   (z_med.EQ.0.D+00).OR.(a_med.EQ.0.D+00).OR.(t_med.EQ.0.D+00)) 
     +   THEN
         ELOSS_HADRON = 0.D+00
         GO TO 1
      ENDIF

      BETA2 = beta * beta
      GAMA2 = 1.D+00 / (1.D+00 - BETA2)
      GAMA  = DSQRT( GAMA2 )
      ETA  = beta * GAMA 

*---- Reduced Bethe-Block stopping power

      CALL EX_ENERG(z_med,d_med,EXEN)
*-    Tabulated material check
      IF(EXEN.EQ.0.D+00) THEN
         ELOSS_HADRON = 0.D+00
         GO TO 1
      ENDIF
      BETH = DLOG(2.D+06*EMAS*BETA2*GAMA2/EXEN) - BETA2

*---- Reduced density correction

      PLAS = 28.8084D+00 * DSQRT(d_med*z_med/a_med)
      C = 2.D+00*DLOG(PLAS) - 2.D+00*DLOG(EXEN) - 1.D+00
      X = DLOG10(ETA)
      CALL HA_DENSI(z_med,d_med,-C,EXEN,X0,X1,M)
*-    Tabulated density consistency
      IF((X0+X1+M).EQ.0.D+00) THEN
         ELOSS_HADRON = 0.D+00
         GO TO 1
      ENDIF      
      A = -1.D+00 * ( C + DLOG(1.D+02)*X0 ) / ((X1-X0)**M)
      IF(X.LT.X0) THEN
         DENS = 0.D+00
      ELSEIF(X.LT.X1) THEN
        DENS = DLOG(1.D+02)*X + C + A*((X1-X)**M)
      ELSE
         DENS = DLOG(1.D+02)*X + C
      ENDIF
      DENS = 0.5D+00 * DENS

      ETA2 = 1.D+00 / (ETA*ETA) 
      ETA4 = ETA2 * ETA2 
      ETA6 = ETA4 * ETA2 

*---- Reduced shell correction

      SHE0 = 4.2237D-01*ETA2 + 3.040D-02*ETA4 - 3.80D-04*ETA6
      SHE1 = 3.8580D+00*ETA2 - 1.668D-01*ETA4 + 1.58D-03*ETA6
      SHEL = 1.D-06 * EXEN * EXEN * (SHE0 + 1.D-03*SHE1*EXEN) / z_med

*---- Total hadron stopping power
   
      HSTP = COEF * z_med * DFLOAT(z_inc) * DFLOAT(z_inc) / a_med
      HSTP = HSTP * ( BETH - DENS - SHEL ) / BETA2

*---- Electron energy loss

      ELOSS_HADRON = HSTP * t_med

    1 RETURN

      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      SUBROUTINE EX_ENERG(z_med,d_med,EXEN)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C Excitation energy of the elemental or compound material from
C 
C R.M. Sternheimer, M.J. Berger and S.M. Seltzer,Atomic and Nuclear Data
C Tables 30 (1984) 261.
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Passed_in variables
C
C           z_med     effective charge of the medium
C           d_med     medium density                 (g/cm^3)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Passed_out variable
C
C            EXEN     Excitation energy              (eV)
C
C-----------------------------------------------------------------------

C Passed_in variables

      REAL*8 z_med,d_med

C Passed_out variable

      REAL*8 EXEN

*---- Some Z=1 isotopes

      IF(DABS(z_med-1.D+00).LT.5.D-01) THEN
*-    Gazeous hydrogen
         IF(d_med.LT.1.D-02) THEN
            EXEN = 19.2D+00
*-    Liquid Hydrogen
         ELSEIF(d_med.LT.1.D-01) THEN
            EXEN = 21.8D+00
*---- Liquid Deuterium
         ELSEIF(d_med.GE.1.D-01) THEN
            EXEN = 21.8D+00
         ENDIF

*---- Some Z=2 isotopes

      ELSEIF(DABS(z_med-2.D+00).LT.1.D-01) THEN
*-    Gaseous/Liquid Helium
         EXEN = 41.8D+00

*---- Plastic scintillator (Polyvinyltolulene 2-CH[3]C[6]H[4]CH=CH[2])
*     
*     z_eff = 3.36842    a_eff =  6.21996    d_med = 1.03200

      ELSEIF(DABS(z_med-3.37D+00).LT.1.D-01) THEN
         EXEN = 64.7D+00

*---- Kapton (Polymide film C[22]H[10]N[2]O[5])
*     
*     z_eff = 5.02564    a_eff =  9.80345    d_med = 1.42000

      ELSEIF(DABS(z_med-5.03D+00).LT.1.D-01) THEN
         EXEN = 79.6D+00

*---- Some Z=6 isotopes

      ELSEIF(DABS(z_med-6.D+00).LT.1.D-01) THEN
         EXEN = 78.0D+00

*---- Air (dry, near sea level, 78% N2 + 22% O2)
*     
*     z_eff = 7.22000    a_eff = 14.46343    d_med = 1.20480E-03

      ELSEIF(DABS(z_med-7.22D+00).LT.1.D-01) THEN
         EXEN = 85.7D+00

*---- Aluminum

      ELSEIF(DABS(z_med-13.D+00).LT.1.D-01) THEN
         EXEN = 166.0D+00

*---- Titanium

      ELSEIF(DABS(z_med-22.D+00).LT.1.D-01) THEN
         EXEN = 233.0D+00

*---- Table overflow

      ELSE
c calculate the parameters according to Stermheimer
         IF (z_med < 13.) THEN
            EXEN = 12.d0 * z_med + 7.d0
         ELSE
            EXEN = 9.76d0*z_med + 58.8*z_med**(-0.19)
         ENDIF
      ENDIF

      RETURN

      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      SUBROUTINE HA_DENSI(z_med,d_med,CBAR,EXEN,X0,X1,M)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C Parameters of the density effect calculation of hadrons for elemental
C or compound materials from
C 
C R.M. Sternheimer, M.J. Berger and S.M. Seltzer,Atomic and Nuclear Data
C Tables 30 (1984) 261.
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Passed_in variables
C
C           z_med     effective charge of the medium
C           d_med     medium density                 (g/cm^3)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Passed_out variables
C
C              X0     Parameter of the density effect fit
C              X1     Parameter of the density effect fit
C               M     Parameter of the density effect fit
C
C-----------------------------------------------------------------------

C Passed_in variables

      REAL*8 z_med,d_med, CBAR, EXEN

C Passed_out variable

      REAL*8 X0,X1,M

*---- Some Z=1 isotopes

      IF(DABS(z_med-1.D+00).LT.1.D-01) THEN
*-    Gazeous hydrogen
         IF(d_med.LT.1.D-02) THEN
            X0 =  1.8639D+00
            X1 =  3.2718D+00
            M  =  5.7273D+00
*-    Liquid Hydrogen
         ELSEIF(d_med.LT.1.D-01) THEN
            X0 =  0.4759D+00
            X1 =  1.9215D+00
            M  =  5.6249D+00
*---- Liquid Deuterium
         ELSEIF(d_med.GE.1.D-01) THEN
            X0 =  0.4759D+00
            X1 =  1.9215D+00
            M  =  5.6249D+00
         ELSE
            call get_dens_par(EXEN, d_med, CBAR, X0, X1, M)
         ENDIF

*---- Some Z=2 isotopes

      ELSEIF(DABS(z_med-2.D+00).LT.1.D-01) THEN
*-    Gaseous/Liquid Helium
         X0 =   2.2017D+00
         X1 =   3.6122D+00
         M  =   5.8347D+00

*---- Plastic scintillator (Polyvinyltolulene 2-CH[3]C[6]H[4]CH=CH[2])
*     
*     z_eff = 3.36842    a_eff =  6.21996    d_med = 1.03200

      ELSEIF(DABS(z_med-3.37D+00).LT.1.D-01) THEN
         X0 =  0.1464D+00
         X1 =  2.4855D+00
         M  =  3.2393D+00

*---- Kapton (Polymide film C[22]H[10]N[2]O[5])
*     
*     z_eff = 5.02564    a_eff =  9.80345    d_med = 1.42000

      ELSEIF(DABS(z_med-5.03D+00).LT.1.D-01) THEN
         X0 =  0.1509D+00
         X1 =  2.5631D+00
         M  =  3.1921D+00

*---- Some Z=6 isotopes

      ELSEIF(DABS(z_med-6.D+00).LT.1.D-01) THEN
*-    Carbon (graphite, density 1.700 g/cm3)
         IF(d_med.LT.1.750D+00) THEN
            X0 =  0.0480D+00
            X1 =  2.5387D+00
            M  =  2.9532D+00
*-    Carbon (graphite, density 2.000 g/cm3)
         ELSEIF(d_med.LT.2.050D+00) THEN
            X0 = -0.0351D+00 
            X1 =  2.4860D+00
            M  =  3.0036D+00
*-    Carbon (graphite, density 2.265 g/cm3)
         ELSEIF(d_med.LT.2.270D+00) THEN
            X0 = -0.0178D+00 
            X1 =  2.3415D+00
            M  =  2.8697D+00
         ELSE
            call get_dens_par(EXEN, d_med, CBAR, X0, X1, M)
         ENDIF

*---- Air (dry, near sea level, 78% N2 + 22% O2)
*     
*     z_eff = 7.22000    a_eff = 14.46343    d_med = 1.20480E-03

      ELSEIF(DABS(z_med-7.22D+00).LT.1.D-01) THEN
         X0 =   1.7418D+00
         X1 =   4.2759D+00
         M  =   3.3994D+00

*---- Aluminum

      ELSEIF(DABS(z_med-13.D+00).LT.1.D-01) THEN
         X0 =  0.1708D+00
         X1 =  3.0127D+00
         M  =  3.6345D+00

*---- Titanium

      ELSEIF(DABS(z_med-22.D+00).LT.1.D-01) THEN
         X0 =  0.0957D+00
         X1 =  3.0386D+00
         M  =  3.0302D+00

*---- Table overflow

      ELSE
c calculate the parameters according to stermheimer 
c Phys.Rev. B Vol 3,3681 (1971)
         call get_dens_par(EXEN, d_med, CBAR, X0, X1, M) 
      ENDIF
      RETURN
 
      END

C-----------------------------------------------------------------------

      subroutine get_dens_par(exen, d_med, cbar, X0, X1, M)
c calculate the parameters according to stermheimer 
c Phys.Rev. B Vol 3,3681 (1971)
      implicit none
c input
      real*8 exen, d_med, cbar
c output
      real*8 X0, X1, M
      if (d_med .gt. 0.01) then
c     this is for solid or liquid
         M = 3.0
         if (EXEN .lt. 100.) then
            X1 = 2.0d0
            if (CBAR .lt. 3.681)then
               X0 = 0.2
            else
               X0 = 0.326*CBAR - 1.d0
            endif
         else
            X1 = 3.0d0
            if (CBAR .lt. 5.215)then
               X0 = 0.2
            else
               X0 = 0.326*CBAR - 1.5d0
            endif
         endif
      else
c     this is for gases
         M = 3.0
         if (CBAR .lt. 12.25)then
            X1 = 4.0
         else
            X1 = 5.0
         endif
         if (CBAR .lt. 10.0)then
            X0 = 1.6
         elseif (CBAR .lt. 10.5)then
            X0 = 1.7
         elseif (CBAR .lt. 11.0)then
            X0 = 1.8
         elseif (CBAR .lt. 11.5)then
            X0 = 1.9
         elseif (CBAR .lt. 12.25)then
            X0 = 2.0
         elseif (CBAR .lt. 13.804)then
            X0 = 2.0
         else
            X0 = 0.326*CBAR - 2.5
         endif
      endif
      return
      end
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

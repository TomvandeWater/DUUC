C=========================================================================
C DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoustic design
C and analysis tool for aircraft with propulsors in ducted fan
C configurations.
C 
C This software was developed under the auspices and sponsorship of the
C Tactical Technology Office (TTO) of the Defense Advanced Research
C Projects Agency (DARPA).
C 
C Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights Reserved
C
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version.
C 
C This program is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C 
C You should have received a copy of the GNU General Public License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
C
C Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (drela@mit.edu)
C Program Management: Brad Tousley, Paul Eremenko (eremenko@alum.mit.edu)
C
C=========================================================================


      SUBROUTINE DFDCINIT(LDEBUG)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCINIT
#endif
C--------------------------------------------------
C     DFDC interface routine - initialize DFDC analysis  
C     Input:  LDEBUG  logical debugging flag (T for debug output)
C--------------------------------------------------
      LOGICAL LDEBUG
C
C---- Set DEBUG flag and initialize
      CALL DFINIT(LDEBUG)
C
      RETURN
      END


      SUBROUTINE DFDCEXEC(LCONVG)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCEXEC
#endif
C--------------------------------------------------
C     DFDC interface routine - solve current case
C     Solve actuator disk or rotor problem using
C     iterative solver
C     Input:   None
C     Output:  LCONVG (logical) T if converged
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LCONVG
C
      RLX    = RLXSOLV
      WXEPS  = EPSSOLV
      ITRMAX = ITRMAXSOLV
C
C---- Solve current case (execute solver)
      RLXF = RLX
      IF(NROTOR.GT.0) THEN
       N = 1
       IF(IRTYPE(N).EQ.1) THEN
         WRITE(*,*) 'Using actuator disk...'
         IF(.NOT.LVMAV) THEN
           CALL ROTINITBGAM
C---- Put flow data into wake grid
           CALL SETGRDFLW
         ENDIF
         CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
cc         CALL CONVGTH(ITRMAX,RLXF,WXEPS)
       ELSEIF(IRTYPE(N).EQ.2) THEN
         WRITE(*,*) 'Using current blade...'
C---- check for uninitialized BGAM for blade
         DO IR = 1, NRC
           IF(BGAM(IR,N).NE.0.0) GO TO 32
         END DO
         CALL ROTINITBLD
C---- Put flow data into wake grid
         CALL SETGRDFLW
C
 32      CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
       ENDIF
C
       RLXF = RLX
       IF(LWRLX) THEN 
         CALL WAKERESET
         IF(IRTYPE(N).EQ.1) THEN
           CALL CONVGTH(ITRMAX,RLXF,WXEPS)
         ELSE
           CALL CONVGTHBG(ITRMAX,RLXF,WXEPS)
         ENDIF
       ENDIF
C
       ITYP = 1
       CALL TQCALC(ITYP)
       LCONVG = LCONV
      ENDIF
C
      RETURN
      END



      SUBROUTINE DFDCEXECTHR(TSPEC,ISPEC,LCONVG)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCEXECTHR
#endif
C--------------------------------------------------
C     DFDC interface routine - solve for specified thrust 
C     Input:  TSPEC - specified thrust (N)
C             ISPEC - 1 to specify rotor thrust
C                     2 to specify total thrust
C     Output:  LCONVG (logical) T if converged
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LCONVG
C
      RLX    = RLXSOLV
      WXEPS  = EPSSOLV
      ITRMAX = ITRMAXSOLV
C
      N = 1
      IF(IRTYPE(N).EQ.2) THEN
        WRITE(*,*) 'Driving PITCH to get specified thrust'
C----- initialize rotor from blade definition
cc        IF(.NOT.LCONV) CALL ROTINITBLD
C----- if rotor is defined converge operating point from rotor geometry
        RLXF = RLX
        CALL CONVGTHBGT(ITRMAX,RLXF,WXEPS,THRUST,ISPEC)
      ELSE
        WRITE(*,*) 'Driving BGAM to get specified thrust'
C----- initialize rotor BGAM from thrust spec (constant BGAM)
        IF(.NOT.LCONV) CALL ROTINITTHR(THRUST)
C---- Put flow data into wake grid
        CALL SETGRDFLW
C----- converge operating point from specified BGAM
        RLXF = RLX
        CALL CONVGTHBGT(ITRMAX,RLXF,WXEPS,THRUST,ISPEC)
cc        CALL CONVGTHT(ITRMAX,RLXF,WXEPS,THRUST,ISPEC)
      ENDIF
      ITYP = 1
      CALL TQCALC(ITYP)
      LCONVG = LCONV
C
      RETURN
      END



      SUBROUTINE DFDCCALCROTORFORCES
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCCALCROTORFORCES
#endif
C--------------------------------------------------
C     DFDC interface routine - calculate forces on rotor
C     Input:  None
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Call routine to calculate rotor forces
      ITYP = 1
      CALL TQCALC(ITYP)
C
      RETURN
      END



      SUBROUTINE DFDCDESBLADE(ND)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCDESBLADE
#endif
C--------------------------------------------------
C     DFDC interface routine - design blade with current B*GAM 
C                              circulation
C     Input:  None
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Design blade to current circulation
      CALL DESBLADE(ND)
C
      RETURN
      END


      SUBROUTINE DFDCWAKEALIGN
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCWAKEALIGN
#endif
C--------------------------------------------------
C     DFDC interface routine - realign current wakes to match flow
C     Input:  None
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Align wakes with flow
      CALL WAKERESET
C
      RETURN
      END



      SUBROUTINE DFDCLOAD(FNAME)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOAD
#endif
C--------------------------------------------------
C     DFDC interface routine - read case from file
C     Input:  FNAME
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAME
C
      WRITE(*,*) 'DFDCLOAD called with filename ',FNAME
      WRITE(*,*) '         length of FNAME    = ',LEN(FNAME)
C
C---- Read case from file
      CALL DFLOAD(FNAME)
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Case file loaded:'
        WRITE(LUNDBG,15) FNAME
      ENDIF
 15   FORMAT(1X,A)
C
      RETURN
      END



      SUBROUTINE DFDCSAVE(FNAME)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSAVE
#endif
C--------------------------------------------------
C     DFDC interface routine - save current case to file
C     Input:  FNAME
C--------------------------------------------------
      CHARACTER*128 FNAME
C
C---- Save case to file
      CALL DFSAVE(FNAME)
C
      RETURN
      END


      SUBROUTINE DFDCPRINTROTORFORCES
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCPRINTROTORFORCES
#endif
C--------------------------------------------------
C     DFDC interface routine - print current rotor/disk status
C     Input:  None
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
C---- Print rotor forces
      CALL ROTRPRT(6)
C
      RETURN
      END



      SUBROUTINE DFDCCPDUMP(FNAME1)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCCPDUMP
#endif
      INCLUDE 'DFDC.INC'
      CHARACTER*(*) FNAME1
C
      CHARACTER*80 FNAME, FILDEF
C
 1000 FORMAT(A)
C
      IF(FNAME1(1:1).NE.' ') THEN
       FNAME = FNAME1
      ELSE
C----- no argument... get it somehow
       IF(NPREFIX.GT.0) THEN
C------ offer default using existing prefix
        FILDEF = PREFIX(1:NPREFIX) // '.cp'
        WRITE(*,1100) FILDEF
 1100   FORMAT(/' Enter filename:  ', A)
        READ(*,1000) FNAME
        CALL STRIP(FNAME,NFN)
        IF(NFN.EQ.0) FNAME = FILDEF
       ELSE
C------ nothing available... just ask for filename
        CALL ASKS('Enter filename^',FNAME)
       ENDIF
      ENDIF
C
C
      LU = 19
      OPEN(LU,FILE=FNAME,STATUS='UNKNOWN')
      REWIND(LU)
C
      WRITE(LU,1000)
     & '#     x          y         Cp         Q  '
C         0.23451  0.23451  0.23451  0.23451
c             123456789012345678901234567890
C
ccc      CALL COMSET
C
ccc      BETA = SQRT(1.0 - MINF**2)
ccc      BFAC = 0.5*MINF**2 / (1.0 + BETA)
cc      BETA = 1.0
cc      BFAC = 0.0
C
      DO 10 IEL=1, NEL
        IF(ICFRST(IEL).LE.0) GO TO 10
C
        IC1 = ICFRST(IEL)
        IC2 = ICLAST(IEL)
        N = IC2 - IC1 + 1
C
        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 1      ) THEN
C-------- dump Cp on right side
          WRITE(LU,1000) ' '
          WRITE(LU,1000) ' '
          DO IC = IC1, IC2
            QMAG = SQRT(QCR(1,IC)**2 + QCR(2,IC)**2) 
            WRITE(LU,8500) XC(IC), YC(IC), CPR(IC), QMAG
          ENDDO
        ENDIF

        IF(ISPLOT(IEL).EQ. 0 .OR.
     &     ISPLOT(IEL).EQ. 2      ) THEN
C-------- dump Cp on left side
          WRITE(LU,1000) ' '
          WRITE(LU,1000) ' '
          DO IC = IC1, IC2
            QMAG = SQRT(QCL(1,IC)**2 + QCL(2,IC)**2) 
            WRITE(LU,8500) XC(IC), YC(IC), CPL(IC), QMAG
          ENDDO
        ENDIF
 10   CONTINUE
C
 8500 FORMAT(5(1X,F10.5))
C
      CLOSE(LU)
      RETURN
      END ! DFDCCPDUMP


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



      SUBROUTINE DFDCLOADNAME(NAMEIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADNAME
#endif
C--------------------------------------------------
C     DFDC interface routine - set case title
C     Input:  TITLIN - case title (80 characters)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER NAMEIN*80
C 
C---- Set case title
      NAME = NAMEIN
      CALL STRIP(NAME,NNAME)
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Case title loaded:'
        WRITE(LUNDBG,15) NAME
      ENDIF
 15   FORMAT(1X,A80)
C
      RETURN
      END


      SUBROUTINE DFDCLOADOPER(QINFIN,QREFIN,RPM1IN,RPM2IN,
     &                        RHOIN,VSOIN,RMUIN,ALTHIN,
     &                        XDWKLENIN,NWAKEIN,LWRLXIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADOPER
#endif
C--------------------------------------------------
C     DFDC interface routine - set case title
C     Input:  
C           QINF      freestream velocity
C           QREF      reference velocity
C           RPM1      RPM for disk #1   
C           RPM2      RPM for disk #2   
C           RHO       Rho   
C           VSO       VSound
C           RMU       rMU   
C           ALTH      altitude (in km)
C           XDWKLEN   X length of wake region downstream of duct and CB
C           NWAKE     # of wake points from TE to downstream boundary
C           LWRLX     wake relaxation flag (T to relax wake with iterations) 
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LWRLXIN
C 
C---- Set operational data
      QINF    = QINFIN   
      QREF    = QREFIN   
      RHO     = RHOIN    
      VSO     = VSOIN    
      RMU     = RMUIN    
      ALTH    = ALTHIN   
      XDWKLEN = XDWKLENIN
      NWAKE   = NWAKEIN  
C
      LWRLX   = .FALSE.
      IF(QREF.LE.0.0) THEN
         QREF = 1.0
      ENDIF
      OMEGA(1) = 2.0*PI*RPM1IN/60.0
      OMEGA(2) = 2.0*PI*RPM2IN/60.0
      IF(RHO.EQ.0.0 .OR. VSO.EQ.0.0 .OR. RMU.EQ.0.0) THEN
         DELTAT = 0.0
         CALL ATMO(ALTH,DELTAT,VSO,RHO,RMU) 
      ENDIF
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'OPER data loaded'
        WRITE(LUNDBG,*) ' Qinf   = ',QINFIN
        WRITE(LUNDBG,*) ' Qref   = ',QREFIN
        WRITE(LUNDBG,*) ' RPM1   = ',RPM1IN
        WRITE(LUNDBG,*) ' RPM2   = ',RPM2IN
        WRITE(LUNDBG,*) ' Rho    = ',RHOIN
        WRITE(LUNDBG,*) ' VSound = ',VSOIN
        WRITE(LUNDBG,*) ' rMU    = ',RMUIN
        WRITE(LUNDBG,*) ' Alt    = ',ALTHIN
        WRITE(LUNDBG,*) ' Xdwake = ',XDWKLENIN
        WRITE(LUNDBG,*) ' Nwake  = ',NWAKEIN
        WRITE(LUNDBG,*) ' LWRLX  = ',LWRLXIN
      ENDIF
C
      RETURN
      END




      SUBROUTINE DFDCLOADAEROALL(ND,NASECT,
     &                           XISECT,A0DEG,
     &                           CLMAX,CLMIN,
     &                           DCLDA,DCLDA_STALL,DCL_STALL,
     &                           CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                           CMCON,MCRIT,TOC,REREF,REXP)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADAEROALL
#endif
C--------------------------------------------------
C     DFDC interface routine - load aero data
C     Input: 
C           ND             disk # for aero sections
C           NASECT         # of aero sections specified
C           XISECT(.)      r/R for aero section
C           A0DEG(.)       A0 in degrees 
C           CLMAX(.)       CLmax for aero section
C           CLMIN(.)       CLmin for aero section
C           DCLDA(.)       dCL/dAlfa for aero section
C           DCLDA_STALL(.) dCL/dAlfa in stall regions
C           DCL_STALL(.)   deltaCL size of transition to stall
C           CDMIN(.)       minimum CD for aero section
C           CLDMIN(.)      CL at minimum CD for aero section
C           DCDCL2(.)      dCD/dCL^2 for aero section
C           DCDCL2S(.)     Secondary, annulus drag param dCD/dCL^2
C           CMCON(.)       constant CM for aero section
C           MCRIT(.)       critical Mach # for aero section
C           TOC(.)         thickness/chord ratio for aero section
C           REREF(.)       reference RN for aero section
C           REXP(.)        RN exponent for aero section
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XISECT(*),A0DEG(*),
     &          CLMAX(*),CLMIN(*),
     &          DCLDA(*),DCLDA_STALL(*),DCL_STALL(*),
     &          CDMIN(*),CLDMIN(*),DCDCL2(*),DCDCL2S(*),
     &          CMCON(*),MCRIT(*),TOC(*),REREF(*),REXP(*)
C 
      IF(ND.LT.1 .OR. ND.GT.NRX) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCLOADAEROALL ', ND
        RETURN
      ENDIF
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'DFDCLOADAEROALL'
      ENDIF
C---- Install AERO data
      DO N = 1, NASECT
        A0 = A0DEG(N)*DTR
        CALL PUTAERO(ND,N,XISECT(N),A0,CLMAX(N),CLMIN(N),
     &               DCLDA(N),DCLDA_STALL(N),DCL_STALL(N),
     &               CDMIN(N),CLDMIN(N),DCDCL2(N),DCDCL2S(N),
     &               CMCON(N),MCRIT(N),TOC(N),REREF(N),REXP(N))
        IF(LDBG) THEN
         WRITE(LUNDBG,1420) ND,N,XISECT(N)
         WRITE(LUNDBG,1500) A0DEG(N),DCDCL2(N)
         WRITE(LUNDBG,1510) DCLDA(N),DCDCL2S(N)
         WRITE(LUNDBG,1515) DCLDA_STALL(N),REREF(N)
         WRITE(LUNDBG,1520) CLMAX(N),REXP(N)
         WRITE(LUNDBG,1530) CLMIN(N),CMCON(N)
         WRITE(LUNDBG,1540) DCL_STALL(N),MCRIT(N)
         WRITE(LUNDBG,1550) CDMIN(N),TOC(N) 
         WRITE(LUNDBG,1560) CLDMIN(N)
        ENDIF
      END DO
      NAERO = NASECT
C---- reset aero section pointers for blade sections
      CALL SETIAERO
C
 1420 FORMAT(/' Disk# = ',I3,' Sect# = ',I3,' r/R = ',F10.4)
C
 1500 FORMAT( 1X,72('=')/
     &       ' 1) Zero-lift alpha (deg):',F7.2,6X,
     &       ' 9) dCd/dCl**2           :',F7.4)
 1510 FORMAT(' 2) dCl/dalpha           :',F7.3,6X,
     &       '10) Secondary dCd/dCl**2 :',F7.4)
 1515 FORMAT(' 3) dCl/dalpha@stall     :',F7.3,6X,
     &       '11) Reference Re number  :',F9.0)
 1520 FORMAT(' 4) Maximum Cl           :',F6.2,7X,
     &       '12) Re scaling exponent  :',F8.4)
 1530 FORMAT(' 5) Minimum Cl           :',F6.2,7X,
     &       '13) Cm                   :',F7.3)
 1540 FORMAT(' 6) Cl increment to stall:',F7.3,6X,
     &       '14) Mcrit                :',F7.3)
 1550 FORMAT(' 7) Minimum Cd           :',F7.4,6X,
     &       '15) Section t/c          :',F7.3)
 1560 FORMAT(' 8) Cl at minimum Cd     :',F6.3/
     &         1X,72('='))
C
      RETURN
      END


      SUBROUTINE DFDCLOADAERO1(ND,IASECT,
     &                         XISECT,A0DEG,
     &                         CLMAX,CLMIN,
     &                         DCLDA,DCLDA_STALL,DCL_STALL,
     &                         CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                         CMCON,MCRIT,TOC,REREF,REXP)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADAERO1
#endif
C--------------------------------------------------
C     DFDC interface routine - load a single aero section data
C     Input: 
C           ND          disk # for aero sections
C           IASECT      index of aero section specified
C           XISECT      r/R for aero section
C           A0DEG       A0 in degrees 
C           CLMAX       CLmax for aero section
C           CLMIN       CLmin for aero section
C           DCLDA       dCL/dAlfa for aero section
C           DCLDA_STALL dCL/dAlfa in stall regions
C           DCL_STALL   deltaCL size of transition to stall
C           CDMIN       minimum CD for aero section
C           CLDMIN      CL at minimum CD for aero section
C           DCDCL2      dCL/dAlfa^2 for aero section
C           DCDCL2S     Secondary, annulus drag param dCD/dCL^2
C           CMCON       constant CM for aero section
C           MCRIT       critical Mach # for aero section
C           TOC         thickness/chord ratio for aero section
C           REREF       reference RN for aero section
C           REXP        RN exponent for aero section
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(ND.LT.1 .OR. ND.GT.NRX) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCLOADAERO1 ', ND
        RETURN
      ENDIF
C
C---- If section # is input as 0 then reset # of aero sections
      IF(IASECT.LE.0) THEN
        NAERO = 0
        IF(LUNDBG) THEN
          WRITE(LUNDBG,*) 'DFDCLOADAERO1 clearing aero sections'
        ENDIF
        RETURN
      ENDIF
C 
C---- Install AERO data into aero section # IASECT
      A0 = A0DEG *DTR
      CALL PUTAERO(ND,IASECT,XISECT,A0,CLMAX,CLMIN,
     &             DCLDA,DCLDA_STALL,DCL_STALL,
     &             CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &             CMCON,MCRIT,TOC,REREF,REXP)
      NAERO = MAX(NAERO,IASECT)
C---- reset aero section pointers for blade sections
      CALL SETIAERO
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'DFDCLOADAERO1'
        A0DEG = A0/DTR
        WRITE(LUNDBG,1420) ND,IASECT,XISECT
        WRITE(LUNDBG,1500) A0DEG,DCDCL2
        WRITE(LUNDBG,1510) DCLDA,DCDCL2S
        WRITE(LUNDBG,1515) DCLDA_STALL,REREF
        WRITE(LUNDBG,1520) CLMAX,REXP
        WRITE(LUNDBG,1530) CLMIN,CMCON
        WRITE(LUNDBG,1540) DCL_STALL,MCRIT
        WRITE(LUNDBG,1550) CDMIN,TOC 
        WRITE(LUNDBG,1560) CLDMIN
      ENDIF
C
 1420 FORMAT(/' Disk# = ',I3,' Sect# = ',I3,' r/R = ',F10.4)
C
 1500 FORMAT( 1X,72('=')/
     &       ' 1) Zero-lift alpha (deg):',F7.2,6X,
     &       ' 9) dCd/dCl**2           :',F7.4)
 1510 FORMAT(' 2) dCl/dalpha           :',F7.3,6X,
     &       '10) Secondary dCd/dCl**2 :',F7.4)
 1515 FORMAT(' 3) dCl/dalpha@stall     :',F7.3,6X,
     &       '11) Reference Re number  :',F9.0)
 1520 FORMAT(' 4) Maximum Cl           :',F6.2,7X,
     &       '12) Re scaling exponent  :',F8.4)
 1530 FORMAT(' 5) Minimum Cl           :',F6.2,7X,
     &       '13) Cm                   :',F7.3)
 1540 FORMAT(' 6) Cl increment to stall:',F7.3,6X,
     &       '14) Mcrit                :',F7.3)
 1550 FORMAT(' 7) Minimum Cd           :',F7.4,6X,
     &       '15) Section t/c          :',F7.3)
 1560 FORMAT(' 8) Cl at minimum Cd     :',F6.3/
     &         1X,72('='))
C
      RETURN
      END



      SUBROUTINE DFDCLOADROTOR(ND,XDSKIN,NRBLDIN,NRPIN,
     &                         NRDIN, YRDIN, CHRDIN, BETAIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADROTOR
#endif
C--------------------------------------------------
C     DFDC interface routine - load rotor definition
C     Input:  
C           ND        disk # 
C           XDSKIN    axial coordinate of rotor disk
C           NRBLDIN   # of rotor blades
C           NRPIN     # of radial stations to use for calculations
C           NRDIN     # of rotor radial stations defined
C           YRDIN(.)  y (r) coordinates for actuator disk definition
C           CHRDIN(.) chord for rotor 
C           BETAIN(.) blade angle for rotor (deg)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION YRDIN(NRDIN), CHRDIN(NRDIN), BETAIN(NRDIN)
C
      IF(ND.LT.1 .OR. ND.GT.NRX) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCLOADROTOR ', ND
        RETURN
      ENDIF
C
C---- load actuator disk data
      XDISK(ND)  = XDSKIN
      NRBLD(ND)  = NRBLDIN
C
      DO IR = 1, NRDIN
        YRDEF(IR,ND)   = YRDIN(IR)
        CHRDEF(IR,ND)  = CHRDIN(IR)
        BETADEF(IR,ND) = BETAIN(IR)*DTR
      END DO   
      NRDEF(ND) = NRDIN
      IRTYPDEF(ND) = 2
      LROTOR = .TRUE.
      NROTOR = MAX(NROTOR,ND)
C---- reset aero sections associated with blade sections
      CALL SETIAERO
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Rotor/Stator blading for disk # ',ND,
     &                  ' at X =',XDISK(ND)
        WRITE(LUNDBG,*) 'Rotor/Stator blading defined for disk # ',ND
        WRITE(LUNDBG,*) 'Blades defined with #points = ',NRDEF(ND)
        WRITE(LUNDBG,*) '   R           CH           Beta'
        DO IR = 1, NRDEF(ND)
          WRITE(LUNDBG,15) YRDEF(IR,ND),CHRDEF(IR,ND),
     &                     BETADEF(IR,ND)/DTR,IAERO(IR,ND)
        END DO
      ENDIF
 15   FORMAT(3(1X,F12.6),2X,I4)
C
      RETURN
      END



      SUBROUTINE DFDCLOADACTDISK(ND,XDSKIN,NRPIN,
     &                           NADIN, YADIN, BGAMIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADACTDISK
#endif
C--------------------------------------------------
C     DFDC interface routine - load actuator disk definition
C     Input:  
C           ND        disk # 
C           XDSKIN    axial coordinate of actuator disk
C           NRPIN     # of radial stations to use for calculations
C           NADIN     # of actuator disk radial stations defined
C           YDIN(.)   y (r) coordinates for actuator disk definition
C           BGAMIN(.) circulation B*GAM for actuator disk
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION YADIN(NADIN), BGAMIN(NADIN)
C
      IF(ND.LT.1 .OR. ND.GT.NRX) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCLOADACTDISK ', ND
        RETURN
      ENDIF
C
C---- load actuator disk data
      XDISK(ND) = XDSKIN
C
      DO IR = 1, NADIN
        YRDEF(IR,ND)   = YADIN(IR)
        BGAMDEF(IR,ND) = BGAMIN(IR)
      END DO   
      NRDEF(ND) = NADIN
      LROTOR = .TRUE.
      IRTYPDEF(ND) = 1
      NROTOR = MAX(NROTOR,ND)
C
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'Rotor/Stator defined for disk # ',ND,
     &                  ' at X =',XDISK(ND)
        WRITE(LUNDBG,*) 'Actuator disk with #points = ',NRDEF(ND)
        WRITE(LUNDBG,*) '   R           B*GAM'
        DO IR = 1, NRDEF(ND)
          WRITE(LUNDBG,15) YRDEF(IR,ND),BGAMDEF(IR,ND)
        END DO
      ENDIF
 15   FORMAT(3(1X,F12.6),2X,I4)
C
      RETURN
      END



      SUBROUTINE DFDCLOADDRAGOBJ(NDDIN,
     &                           XDDIN,YDDIN,CDAIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADDRAGOBJ
#endif
C--------------------------------------------------
C     DFDC interface routine - load drag object definition
C     Input:  
C           NDDIN     # of points defining drag object
C           XDDIN(.)  X coordinates for drag object point
C           YDDIN(.)  y (r) coordinates for drag object point
C           CDAIN(.)  CDA for drag object point
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XDDIN(NDDIN),YDDIN(NDDIN),CDAIN(NDDIN)
C
C---- load drag object data
      N = NDOBJ + 1
      DO ID = 1, NDDIN
        XDDEF(ID,N)  = XDDIN(ID)
        YDDEF(ID,N)  = YDDIN(ID)
        CDADEF(ID,N) = CDAIN(ID)
      END DO   
      NDDEF(N) = NDDIN
      NDOBJ = NDOBJ + 1
      IF(LDBG) WRITE(LUNDBG,*) 'Drag object ',NDOBJ,NDDEF(NDOBJ)
C
C==========================================================
C---- Add drag object into buffer geometry as element 3 or higher
      IF(NDDEF(N).GT.1) THEN
C
        IF(NBEL.LT.2) THEN
          IEL = 3
         ELSE
          IEL = NBEL + 1
        ENDIF
C        
        IF(IBLAST(IEL-1).GE.1) THEN 
          IBNEXT = IBLAST(IEL-1) + 1
        ELSE
          IBNEXT = 1
        ENDIF
        IBFRST(IEL) = IBNEXT
        IB = IBNEXT - 1       
C--- Get drag object definition data
        DO ID = 1, NDDEF(N)
          IB = IB+1
          IF(IB.GT.IBX) STOP 'LOAD: Array overflow on IBX'
          XB(IB) = XDDEF(ID,N)
          YB(IB) = YDDEF(ID,N)
          IF(LDBG) WRITE(LUNDBG,*) IB,XB(IB),YB(IB)
        END DO 
C---- set buffer drag objects to NBTYPE=5  
        NBTYPE(IEL) = 5
        IBLAST(IEL) = IB
C---- set total number of elements and points in buffer geometry
        NB = IBLAST(IEL)
        IELDRGOBJ(NDOBJ) = IEL
        NBEL = IEL
        IF(LDBG) THEN
          WRITE(LUNDBG,*) 'Drag object element ',IEL
          WRITE(LUNDBG,*) '  IB1,IB2 ',IBFRST(IEL),IBLAST(IEL)
          WRITE(LUNDBG,*) '  NBTYPE  ',NBTYPE(IEL)
          WRITE(LUNDBG,*) '  NBEL,NB ',NBEL,NB
        ENDIF
C
      ENDIF
C
      RETURN
      END




      SUBROUTINE DFDCLOADGEOMETRY(NCBIN,XCBIN,YCBIN,
     &                            NDCTIN,XDCTIN,YDCTIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCLOADGEOMETRY
#endif
C--------------------------------------------------
C     DFDC interface routine - load duct geometry for CB and duct
C     Input:  
C           NCBIN     # of points in CB X,Y definitions
C           XCBIN(.)  X     coordinates for CB definition
C           YCBIN(.)  Y (r) coordinates for CB definition
C
C           NDCTIN    # of points in duct X,Y definitions
C           XDCTIN(.) X     coordinates for duct definition
C           YDCTIN(.) Y (r) coordinates for duct definition
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LCBIN, LDCTIN
      DIMENSION XCBIN(NCBIN),   YCBIN(NCBIN)
      DIMENSION XDCTIN(NDCTIN), YDCTIN(NDCTIN)
C
C---- Check for # of points to see which foils are specified
      LCBIN  = (NCBIN.GT.0)
      LDCTIN = (NDCTIN.GT.0)
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'DFDCLOADGEOMETRY '
        WRITE(LUNDBG,*) ' NCBIN  ',NCBIN
        WRITE(LUNDBG,*) ' NDCTIN ',NDCTIN
        WRITE(LUNDBG,*) ' NBEL,NB ',NBEL,NB
        DO IEL = 1, NBEL
          WRITE(LUNDBG,*) ' IEL, IB1,IB2 ',IEL,IBFRST(IEL),IBLAST(IEL)
        END DO
      ENDIF
C
      IF(LCBIN) THEN
         IF(LDBG) WRITE(LUNDBG,*) 'DFDCLOADGEOMETRY - load CB'
C
C---- Move element points up in buffer geometry
       IF(NBEL.GT.1) THEN
        IEL = 1
        NOLD = IBLAST(IEL) - IBFRST(IEL) + 1
        IF(NOLD.LT.0) NOLD = 0
        NMOV = NCBIN - NOLD
        IF(NMOV.NE.0) THEN
         IF(LDBG) WRITE(LUNDBG,*) 'CB moving buffer points NMOV=',NMOV
C---- Move element points by NMOV 
         DO IE = NBEL, 2, -1
           IB1 = IBFRST(IE)
           IB2 = IBLAST(IE)
           NBP = IB2-IB1+1
C---- Only move elements with points
           IF(NBP.GT.0) THEN 
cc            write(LUNDBG,*) 'before move IE,IB1,IB2 ',IE, IB1,IB2
            IF(NMOV.GT.0) THEN 
              DO IB=IB2,IB1,-1
                XB(IB+NMOV) = XB(IB)  
                YB(IB+NMOV) = YB(IB)  
              END DO
            ELSE
              DO IB=IB1,IB2
                XB(IB-NMOV) = XB(IB)  
                YB(IB-NMOV) = YB(IB)  
              END DO
            ENDIF
            IBFRST(IE) = IB1+NMOV
            IBLAST(IE) = IB2+NMOV
cc            write(LUNDBG,*) 'after move IE,IB1,IB2 ',IE,IB1+NMOV,IB2+NMOV
           ENDIF
         END DO
C---- reset total number of points in buffer geometry
         NB = IBLAST(NBEL)
        ENDIF
       ENDIF
C---- Process CB input points into buffer arrays as IEL=1
       IEL = 1
       IB1 = 1
       IB2 = IB1 + NCBIN-1
       IF(LDBG) WRITE(LUNDBG,*) 'CB Element IEL,IB1,IB2 ',IEL,IB1,IB2
       DO I=1, NCBIN
         IB = IB1 + I-1
         XB(IB) = XCBIN(I)  
         YB(IB) = YCBIN(I)  
         IF(LDBG) WRITE(LUNDBG,*) IB,XB(IB),YB(IB)
       END DO
       IBFRST(IEL) = IB1
       IBLAST(IEL) = IB2
       IF(NBEL.LE.0) THEN
         NBEL = 1
         NB   = IBLAST(IEL)
       ENDIF
       IF(LDBG) WRITE(LUNDBG,*) 
C---- set buffer airfoils to NBTYPE=0  
       NBTYPE(IEL) = 0
C
      ENDIF
C
C---- Now process the duct geometry inputs
      IF(LDCTIN) THEN
C
       IF(LDBG) WRITE(LUNDBG,*) 'DFDCLOADGEOMETRY - load Duct'
C---- Move element points up in buffer geometry
       IF(NBEL.GT.2) THEN
        IEL = 2
        NOLD = IBLAST(IEL) - IBFRST(IEL) + 1
        IF(NOLD.LT.0) NOLD = 0
        NMOV = NDCTIN - NOLD
        IF(NMOV.NE.0) THEN
         IF(LDBG) WRITE(LUNDBG,*) 'DCT moving buffer points NMOV=',NMOV
C---- Move element points by NMOV 
         DO IE = NBEL, 3, -1
           IB1 = IBFRST(IE)
           IB2 = IBLAST(IE)
           NBP = IB2-IB1+1
C---- Only move elements with points
           IF(NBP.GT.0) THEN 
cc            write(LUNDBG,*) 'before move IE,IB1,IB2 ',IE,IB1,IB2
            IF(NMOV.GT.0) THEN 
              DO IB=IB2,IB1,-1
                XB(IB+NMOV) = XB(IB)  
                YB(IB+NMOV) = YB(IB)  
              END DO
            ELSE
              DO IB=IB1,IB2
                XB(IB-NMOV) = XB(IB)  
                YB(IB-NMOV) = YB(IB)  
              END DO
            ENDIF
            IBFRST(IE) = IB1+NMOV
            IBLAST(IE) = IB2+NMOV
cc            write(LUNDBG,*) 'after move IE,IB1,IB2 ',IE,IB1+NMOV,IB2+NMOV
           ENDIF
         END DO
C---- reset total number of points in buffer geometry
         NB = IBLAST(NBEL)
        ENDIF
       ENDIF
C---- Process Duct input points into buffer arrays as IEL=2
       IEL = 2
       IB1 = IBLAST(IEL-1) + 1
       IB2 = IB1 + NDCTIN-1
       IF(LDBG) WRITE(LUNDBG,*) 'DCT Element IEL,IB1,IB2 ',IEL,IB1,IB2
       DO I=1, NDCTIN
         IB = IB1 + I-1
         XB(IB) = XDCTIN(I)  
         YB(IB) = YDCTIN(I)  
         IF(LDBG) WRITE(LUNDBG,*) IB,XB(IB),YB(IB)
       END DO
       IBFRST(IEL) = IB1
       IBLAST(IEL) = IB2
       IF(NBEL.LE.1) THEN
         NBEL = 2
C---- set total number of points in buffer geometry
         NB   = IBLAST(IEL)
       ENDIF
C---- set buffer airfoils to NBTYPE=0  
       NBTYPE(IEL) = 0
C
      ENDIF
      IF(LDBG) THEN
        WRITE(LUNDBG,*) 'After geometry loaded'
        WRITE(LUNDBG,*) ' NBEL,NB ',NBEL,NB
      ENDIF
C
      RETURN
      END


      SUBROUTINE DFDCGENGEOM(FNAMIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGENGEOM
#endif
C--------------------------------------------------
C     DFDC interface routine - generate paneled geometry
C     Input:  FNAMIN - file name for case (used for xxx.pan
C                      paneling parameter file
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER FNAMIN*(*)
C 
C---- process buffer geometry and set up paneled geometry and grid
      CALL GEPROC
      CALL GENGEOM(FNAMIN)
C
      RETURN
      END



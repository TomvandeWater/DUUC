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



      SUBROUTINE DFDCSETDEBUG(LDEBUG)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETDEBUG
#endif
C--------------------------------------------------
C     DFDC interface routine - set debug flag
C     Input:  LDEBUG  logical debugging flag (T for debug output)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LDEBUG
C
C---- Set DEBUG flag 
      LDBG = LDEBUG
C
      RETURN
      END



      SUBROUTINE DFDCSETXDISK(ND,XRAXIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETXDISK
#endif
C--------------------------------------------------
C     DFDC interface routine - set rotor axial X location
C     Input:  
C           ND -     index of disk (rotor or stator)
C           XRAXIN - rotor axial (X) location
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C 
      IF(ND.LT.1 .OR. ND.GT.NROTOR) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCSETXDISK ', ND
        RETURN
      ENDIF
C
C---- Set Rotor axial location
      XDISK(ND) = XRAXIN
C
C---- Adjust paneling on airfoils for rotor and wake grid
C     Re-initialize rotor setup
      CALL ADJPANL
CHHY remember to reset any existing rotor definition to new hub, tip ??
      CALL ROTORINIT
      CALL INIGRD
      CALL SETROTWAK
      CALL ROTPINIT
cc      CALL GETROTTHR
C----------------------------------------------------------------
C---- invalidate any existing solution
      LNCVP = .FALSE.
      LQAIC = .FALSE.
      LQGIC = .FALSE.
      LQCNT = .FALSE.
      LSYSP = .FALSE.
      LGSYS = .FALSE.
      LGAMU = .FALSE.
      LGAMA = .FALSE.
      LSIGP = .FALSE.
      LSIGM = .FALSE.
C
      NQSP = 0
      NSEG = 0
C
      RETURN
      END



      SUBROUTINE DFDCSETNRPTS(NRSTAIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETNRPTS
#endif
C--------------------------------------------------
C     DFDC interface routine - set # rotor points
C     Input:  NRSTAIN - # of points used on rotor
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Set # of rotor radial stations
      NRSTA = NRSTAIN
C
      CALL ROTORINIT
      CALL INIGRD
      CALL SETROTWAK
      CALL ROTPINIT
C
      RETURN
      END




      SUBROUTINE DFDCSETNBLD(ND,NBLDIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETNBLD
#endif
C--------------------------------------------------
C     DFDC interface routine - set blade count
C     Input:  ND - index of disk 
C             NBLDIN - # of rotor blades
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(ND.LT.1 .OR. ND.GT.NROTOR) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCSETNBLD ', ND
        RETURN
      ENDIF
C
C---- Set # of rotor blades on disk ND
      NRBLD(ND) = NBLDIN
C
      RETURN
      END


      SUBROUTINE DFDCSETCL(CL1,CL2)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETCL
#endif
C--------------------------------------------------
C     DFDC interface routine - set blade design CL
C     Input:  CL1 - CL at inboard  blade radial station
C             CL2 - CL at outboard blade radial station
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Set CL for rotor blade
      DO IR = 1, NRC
        CLDES(IR) = CL1 + FLOAT(IR-1)/FLOAT(NRC-1) * (CL2-CL1)
      END DO
C
      RETURN
      END


      SUBROUTINE DFDCSETPITCH(ND,DANG1,DANG2)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETPITCH
#endif
C--------------------------------------------------
C     DFDC interface routine - change blade pitch
C     Input:  ND - index of disk 
C             DANG1 - change in blade angle at  inboard station
C             DANG2 - change in blade angle at outboard station
C             (DANG1, DANG2 in degrees)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      IF(ND.LT.1 .OR. ND.GT.NROTOR) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCSETPITCH ', ND
        RETURN
      ENDIF
C
C---- Change blade angle by linear angle change (hub-tip)
      ANG1 = DANG1*DTR
      ANG2 = DANG2*DTR
      DO IR = 1, NRC
        BETAR(IR,ND) = BETAR(IR,ND) 
     &               + ANG1*DTR 
     &               + FLOAT(IR-1)/FLOAT(NRC-1) * (ANG2-ANG1)*DTR
      END DO
C
      RETURN
      END




      SUBROUTINE DFDCSETBGAM(ND,BGAM1,BGAM2)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETBGAM
#endif
C--------------------------------------------------
C     DFDC interface routine - change specified rotor circulation 
C     Input:  ND - index of disk 
C             BGAM1 - specified B*GAM at  inboard station
C             BGAM2 - specified B*GAM at outboard station
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      IF(ND.LT.1 .OR. ND.GT.NROTOR) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCSETBGAM ', ND
        RETURN
      ENDIF
C
C---- Set rotor circulation
      DO IR = 1, NRC
        BGAM(IR,ND) = BGAM1 + FLOAT(IR-1)/FLOAT(NRC-1) * (BGAM2-BGAM1)
      END DO
C----- initialize rotor operating point from B*GAMMA
      CALL ROTINITBGAM
C---- Reset rotor flag back to actuator disk
      IRTYPE(ND) = 1
C
      RETURN
      END



      SUBROUTINE DFDCSETRPM(ND,RPMIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETRPM
#endif
C--------------------------------------------------
C     DFDC interface routine - set rotor RPM
C     Input:  ND - index of disk 
C             RPMIN - rotor RPM
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
      IF(ND.LT.1 .OR. ND.GT.NROTOR) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCSETRPM ', ND
        RETURN
      ENDIF
C
C---- Set rotor RPM
      OMEGA(ND) = RPMIN*PI/30.
C
      RETURN
      END




      SUBROUTINE DFDCSETALT(ALTIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETALT
#endif
C--------------------------------------------------
C     DFDC interface routine - change fluid properties using 
C     standard atmosphere from altitude
C
C     Input:  ALTIN - altitude in km
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
C---- Set fluid properties using std atmosphere and altitude
      ALT = ALTIN
      CALL ATMO(ALT,VSO,RHO,RMU)
ccc      CALL FLOSHO(6, VSO, RHO, RMU)
C
      RETURN
      END



      SUBROUTINE DFDCSETFLUIDPROP(VSOIN,RHOIN,RMUIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETFLUIDPROP
#endif
C--------------------------------------------------
C     DFDC interface routine - change fluid properties
C     Input:  VSOIN - speed of sound (m/s)
C             RHOIN - density (kg/m^3)
C             RMU   - viscosity (kg/m-s)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
C---- Set fluid properties
      VSO = VSOIN
      RHO = RHOIN
      RMU = RMUIN
C
      RETURN
      END



      SUBROUTINE DFDCSETSOLVPARMS(RLXIN,EPSIN,ITRMAXIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETSOLVPARMS
#endif
C--------------------------------------------------
C     DFDC interface routine - change iterative solver parameters
C     Input:  RLXIN - relaxation factor
C             EPSIN - convergence criteria (for WX changes)
C             ITRMAXIN - maximum # of iterations
C
C     Suggested values:   RLXIN = 0.7
C                         EPSIN = 0.0002
C                         ITRMAXIN = 30
C     These values are installed by default if any of the 
C     inputs are 0.0 or 0
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
C---- Set solver parameters (0., 0., 0) gives defaults
      IF(RLXIN.LE.0.0 .OR.
     &   EPSIN.LE.0.0 .OR.
     &   ITRMAXIN.LE.0 ) THEN
        RLXSOLV = 0.7
        EPSSOLV = 0.0002
        ITRMAXSOLV = 30
      ELSE
        RLXSOLV = RLXIN
        EPSSOLV = EPSIN
        ITRMAXSOLV = ITRMAXIN
      ENDIF
C
      RETURN
      END




      SUBROUTINE DFDCSETQINF(QINFIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETQINF
#endif
C--------------------------------------------------
C     DFDC interface routine - change freestream velocity
C     Input:  QINFIN - freestream velocity
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      INCLUDE 'PLOT.INC'
C
C---- Set Qinf from specified value, limit to 0.0 min value
      QINF = MAX(0.0,QINFIN)
C
      RETURN
      END


      SUBROUTINE DFDCSETQREF(QREFIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETQREF
#endif
C--------------------------------------------------
C     DFDC interface routine - change freestream velocity
C     Input:  QREFIN - freestream velocity
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Set Qref from specified value, set to 1.0 if specified as 0.0
      QREF = QREFIN
      IF(QREF.LE.0.0) QREF = 1.0
      CALL QCPFOR
C
      RETURN
      END





      SUBROUTINE DFDCSETWRLX(LWRLXIN)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCSETWRLX
#endif
C--------------------------------------------------
C     DFDC interface routine - set wake relaxation flag
C     Input:  LWRLXIN - wake relaxation flag (T for relax wake in iteration)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LWRLXIN
C
C---- Set wake relaxation flag
      LWRLX = LWRLXIN
C
      RETURN
      END



C--------------------------------------------------------------------
C---Unimplemented commands
C     VGET command
C----- Enter slipstream velocity inflow
c        CALL GETVEL(FNAME)
C
C     PVEL command
C----- print slipstream velocities just downstream of rotor
c       CALL PRTVEL('Rotor velocities',.TRUE.,.TRUE.,.TRUE.)
C
C     SHOW command
C----- display duct, actuator disk or blade geometry, and drag objects
c        CALL FLOSHO(6, VSO, RHO, RMU)
c        CALL SHOWDUCT
c        CALL SHOWACTDSK
c        CALL SHOWBLADE
c        CALL SHOWDRAGOBJ
c
C     PRINT FORCES command
c        IF(NINPUT.EQ.0) THEN
C------- no element index specified... display all elements
c         IEL = 0
c        ELSE
C------- display only specified element
c         IEL = IINPUT(1)
c        ENDIF
C
c        CALL PRFORC(IEL)
C--------------------------------------------------------------------







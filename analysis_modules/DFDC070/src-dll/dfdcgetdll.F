cdC=========================================================================
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


      SUBROUTINE DFDCGETNAME(CNAME)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETNAME
#endif
C--------------------------------------------------
C     DFDC interface routine - get current case title
C     Input:  CNAME - case title (80 characters)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      CHARACTER CNAME*80
C 
C---- Get case title
      CNAME = NAME
C
      RETURN
      END


      SUBROUTINE DFDCGETFLOWCOND(VINF,VREF,
     &                           RHOINF,RMUINF,VSOINF)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETFLOWCOND
#endif
C--------------------------------------------------
C     DFDC interface routine - get current flow conditions
C     Input:      None
C     Output:  
C         VINF    freestream velocity
C         VREF    reference velocity
C         RHOINF  density
C         RMUINF  viscosity
C         VSOINF  speed of sound
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      VINF = QINF
      VREF = QREF
C
      RHOINF = RHO
      RMUINF = RMU
      VSOINF = VSO
C
      RETURN
      END


      SUBROUTINE DFDCGETRPM(ND,RRPM)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETRPM
#endif
C--------------------------------------------------
C     DFDC interface routine - get RPM on disk # ND
C     Input:     
C         ND     # of disk (rotor or stator)
C     Output:  
C         RRPM   disk RPM
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      RRPM = 60.0*OMEGA(ND)/PI/2.0
C
      RETURN
      END


      SUBROUTINE DFDCGETXDISK(ND,XDSK)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETXDISK
#endif
C--------------------------------------------------
C     DFDC interface routine - get axial location of disk # ND
C     Input:     
C         ND     # of disk (rotor or stator)
C     Output:  
C         XDSK   axial location of disk
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      XDSK = XDISK(ND)
C
      RETURN
      END



      SUBROUTINE DFDCGETNUMDISKS(NDSKS)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETNUMDISKS
#endif
C--------------------------------------------------
C     DFDC interface routine - get # of current rotor/stators defined 
C     Input:      None
C     Output:  
C         NDSKS       # of rotor/stator disks defined
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      NDSKS = NROTOR
C
      RETURN
      END


      SUBROUTINE DFDCGETROTORDEFSIZE(NRPT)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETROTORDEFSIZE
#endif
C--------------------------------------------------
C     DFDC interface routine - get # of current rotor 
C                              or actuator disk defining points
C     Input:      None
C     Output:  
C         NRPT       # of rotor points defined in radial station data
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      NRPT = NRC
C
      RETURN
      END

      SUBROUTINE DFDCGETROTORDEF(ND,
     &                           IRTYP,XDSK,RRTIP,
     &                           NRPT,
     &                           XI,CHI,BETAI,BGAMI)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETROTORDEF
#endif
C--------------------------------------------------
C     DFDC interface routine - get current rotor or actuator disk definition
C     Input:     
C         ND         # of disk (rotor or stator)
C     Output:  
C         IRTYP      flag for disk type (1 for act. disk, 2 for blade)
C         XDSK       axial location of disk
C         RRTIP      rotor tip radius
C
C         NRPT       # of rotor points defined in radial station data
C         XI(.)      r/R at rotor station
C         CHI(.)     c/R at rotor station (blade defined only)
C         BETAI(.)   blade angle (in deg) at rotor station (blade defined only)
C         BGAMI(.)   blade circulation B*GAM at rotor station
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LBLADEF
      DIMENSION XI(*),CHI(*),BETAI(*)
      DIMENSION BGAMI(*)
C
      IF(ND.LT.1 .OR. ND.GT.NRX) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCGETROTORDEF ', ND
        RETURN
      ENDIF
C
      RRTIP = RTIP(ND)
      IRTYP = IRTYPE(ND)
      XDSK = XDISK(ND)
C
      DO I=1, NRC
        XI(I)    = YRC(I,ND)/RRTIP
        CHI(I)   = CHR(I,ND)/RRTIP
        BETAI(I) = BETAR(I,ND)/DTR
        BGAMI(I) = BGAM(I,ND)
      END DO
      NRPT = NRC
C
      RETURN
      END

      SUBROUTINE DFDCGETROTORDATA(ND,
     &                            IRTYP,NBLD,XDSK,RRTIP,RRHUB,
     &                            VINF,VREF,RHOINF,RMUINF,VSOINF,
     &                            RRPM,ADV,
     &                            TDIM,TVDIM,TNACEL,
     &                            PDIM,PVDIM,QDIM,
     &                            EFFTOT,EFFIND,EIDEAL,
     &                            CT, CP, TC, PC,
     &                            CTH,CPH,FOMH,SIGMA,CTOS,
     &                            NRPT,
     &                            XI,CHI,BETAI,LSTALLI,
     &                            CLI,CDI,REI,MACHI,BGAMI)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETROTORDATA
#endif
C--------------------------------------------------
C     DFDC interface routine - get current rotor state
C     Input:      
C         ND      # of disk (rotor or stator)
C     Output:  
C         IRTYP   flag for disk type (1 for act. disk, 2 for blade)
C         NBLD    # of blades
C         XDSK    axial location for disk
C         RRTIP   rotor tip radius
C         RRHUB   rotor hub radius
C         VINF    freestream velocity
C         VREF    reference velocity
C         RHOINF  density
C         RMUINF  viscosity
C         VSOINF  speed of sound
C         RRPM    rotor RPM
C         ADV     rotor advance ratio (Vinf/(omega*Rtip)
C         TDIM    total thrust
C         TVDIM   viscous contribution to thrust
C         TNACEL  nacelle thrust (CB+duct)
C         PDIM    total power
C         PVDIM   viscous contribution to power
C         QDIM    total torque
C         EFFTOT  total efficiency
C         EFFIND  induced efficiency
C         EIDEAL  ideal efficiency
C         CT, CP  thrust, power coefficients based on rotational speed
C         TC, PC  thrust, power coefficients based on forward speed
C         CTH, CPH  thrust, power coefficients based on helicopter definitions
C         FOMH    figure of merit (only useful for hover)
C         SIGMA   blade solidity at 3/4 radius
C         CTOS    CTH/SIGMA  helicopter measure of blade loading
C
C         NRPT       # of rotor points defined in radial station data
C         XI(.)      r/R at rotor station
C         CHI(.)     c/R at rotor station
C         BETAI(.)   blade angle (in deg) at rotor station
C         LSTALLI(.) stall flag at rotor station (T for stalled section)
C         CLI(.)     blade CL at rotor station
C         CDI(.)     blade CD at rotor station
C         REI(.)     blade Reynolds number at rotor station
C         MACHI(.)   blade Mach number at rotor station
C         BGAMI(.)   blade circulation B*GAM at rotor station
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      LOGICAL LSTALLI
      DIMENSION XI(*),CHI(*),BETAI(*),LSTALLI(*)
      DIMENSION CLI(*),CDI(*),REI(*),MACHI(*),BGAMI(*)
C---- Local temporary storage
      DIMENSION WRK1(IRX)
C
      IF(ND.LT.1 .OR. ND.GT.NRX) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCGETROTORDATA ', ND
        RETURN
      ENDIF
C
      VINF = QINF
      VREF = QREF
      RHOINF = RHO
      RMUINF = RMU
      VSOINF = VSO
C
      IRTYP = IRTYPE(ND)
      NBLD  = NRBLD(ND)
      XDSK  = XDISK(ND)
      RRTIP = RTIP(ND)
      RRHUB = RHUB(ND)
      RRPM  = 60.0*OMEGA(ND)/PI/2.0
C
C---- advance ratio
      IF(OMEGA(ND).GT.0.0) THEN
        ADV = VINF/(OMEGA(ND)*RTIP(ND))
      ELSE
        ADV = 0.0
      ENDIF
C---- dimensional thrust, power, torque, rpm
      TNACEL = TDUCT
      TDIM   = TTOT + TNACEL
      QDIM   = QTOT
      PDIM   = PTOT
      TVDIM  = TVIS
      PVDIM  = PVIS
      TIDIM  = TDIM - TVIS
      PIDIM  = PDIM - PVIS
      IF(PDIM.NE.0.0) THEN
        EFFTOT = QINF*TDIM/PDIM
      ELSE
        EFFTOT = 0.0
      ENDIF
C---- blade solidity
      CALL SPLINE(CHR,WRK1,YRC,NRC)
      CH34  = SEVAL(0.75*RTIP(ND),CHR(1,ND),WRK1,YRC(1,ND),NRC)
      SIGMA = FLOAT(NRBLD(ND))*CH34/RTIP(ND)/PI
C---- standard coefficients based on forward speed
      IF(QINF.GT.0.0) THEN
        TC = TDIM/(0.5*RHO*QINF**2 * ATIP(ND))
        PC = PDIM/(0.5*RHO*QINF**3 * ATIP(ND))
      ELSE
        TC = 0.0
        PC = 0.0
      ENDIF
C---- standard coefficients based on rotational speed
      DIA = 2.0*RTIP(ND)
      EN = RRPM/60.0
      IF(EN.GT.0.0) THEN
        CT = TDIM/(RHO*EN**2*DIA**4)
        CP = PDIM/(RHO*EN**3*DIA**5)
      ELSE
        CT = 0.0
        CP = 0.0
      ENDIF
C---- induced efficiency (including nacelle thrust effect)
      IF(PIDIM.NE.0.0) THEN
        EFFIND = QINF*(TIDIM)/PIDIM
      ELSE
        EFFIND = 0.0
      ENDIF
C---- ideal (actuator disk) efficiency
      IF(TC.EQ.0) THEN
        EIDEAL = 0.0
      ELSE
        TCLIM = MAX( -1.0 , TC )
        EIDEAL = 2.0 / (1.0 + SQRT(TCLIM + 1.0))
      ENDIF
C---- define low advance ratio (hover?) data using helicopter nomenclature
      IF(OMEGA(ND).GT.0.0) THEN
       VTIP = OMEGA(ND)*RTIP(ND)
       CTH  = TDIM/(RHO*ADISK(ND)*VTIP**2)
       CPH  = PDIM/(RHO*ADISK(ND)*VTIP**3)
       CTOS = CTH / SIGMA
       IF(CTH.GE.0.0 .AND. CPH.NE.0.0) THEN
         FOMH = ABS(CTH)**1.5 / CPH / 2.0
       ELSE
         FOMH = 0.0
       ENDIF
      ENDIF
C
      DO I=1, NRC
        XI(I)    = YRC(I,ND)/RTIP(ND)
        CHI(I)   = CHR(I,ND)/RTIP(ND)
        BETAI(I) = BETAR(I,ND)/DTR
        REI(I)   = RER(I,ND)
        CLI(I)   = CLR(I,ND)
        CDI(I)   = CDR(I,ND)
        MACHI(I) = MACHR(I,ND)
        BGAMI(I) = BGAM(I,ND)
        LSTALLI(I) = LSTALLR(I,ND)
      END DO
      NRPT = NRC
C
      RETURN
      END


      SUBROUTINE DFDCGETFORCE1(IEL,ITYPE,
     &                         FX,FY,MZ,CXI,CYI,CMI,
     &                         VINF,VREF,RHOINF,QUREF)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETFORCE1
#endif
C--------------------------------------------------
C     DFDC interface routine - get forces on element IEL
C     Input:   IEL  - element index (0 for total forces)
C     Output:  IFORC  element force flag (0 for no forces, 1 for forces)
C              FX   - axial force (X direction)
C              FY   - radial force (Y direction)
C              MZ   - moment about element reference center (x,y)=(0.0,0.0)
C              CXI  - coefficient of axial force  (per unit QUEref)
C              CYI  - coefficient of radial force (per unit QUEref)
C              CMI  - moment coefficient about reference (per unit QUEref)
C
C              VINF - freestream velocity
C              VREF - reference velocity
C              RHOINF - density
C              QUREF- reference dynamic pressure (0.5*RHOINF*Vref^2)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Set dynamic pressure
      VINF = QINF
      VREF = QREF
      RHOINF = RHO
      QUREF = 0.5*RHOINF*VREF**2
C
C---- Exclude wake and vortex wake elements from force reporting
      FX  = 0.0
      FY  = 0.0
      MZ  = 0.0
      CXI = 0.0
      CYI = 0.0
      CMI = 0.0
      IFORC = 0
      IF(IEL.EQ.0 .OR. 
     &   NETYPE(IEL).EQ.0 .OR.
     &   NETYPE(IEL).EQ.5 .OR.
     &   NETYPE(IEL).EQ.6) THEN
C
        CXI = CX(IEL)
        CYI = CY(IEL)
        CMI = CM(IEL)
        FX = CXI*QUREF
        FY = CYI*QUREF
        MZ = CMI*QUREF
        IFORC = 1
      ENDIF
C
      RETURN
      END 


      SUBROUTINE DFDCGETFORCEALLSIZE(NELE1)
     &                         
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETFORCEALLSIZE
#endif
C--------------------------------------------------
C     DFDC interface routine - get # of elements with defined forces
C     Input:   None
C     Output:
C              NELE1      # elements + 1 (element #0 has total forces)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      NELE1 = NEL+1
C
      RETURN
      END 


      SUBROUTINE DFDCGETFORCEALL(VINF,VREF,RHOINF,QUREF,
     &                           NELE,IFORC,FX,FY,MZ,CXI,CYI,CMI)
     &                         
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETFORCEALL
#endif
C--------------------------------------------------
C     DFDC interface routine - get forces on elements (CB, duct, rotor 
C                              and drag elements)
C     Input:   None
C     Output:
C              VINF      freestream velocity
C              VREF      reference velocity
C              RHOINF    density
C              QUREF     reference dynamic pressure (0.5*RHOINF*Vref^2)
C
C              NELE      # elements (also, element 0 for total forces)
C              IFORC(.)  element force flag (0 for no forces, 1 for forces)
C              FX(.)     axial force (X direction)
C              FY(.)     radial force (Y direction)
C              MZ(.)     moment about element reference center (x,y)=(0.0,0.0)
C              CXI(.)    coefficient of axial force  (per unit QUEref)
C              CYI(.)    coefficient of radial force (per unit QUEref)
C              CMI(.)    moment coefficient about reference (per unit QUEref)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION IFORC(0:IRX),
     &          FX(0:IRX),FY(0:IRX),MZ(0:IRX),
     &          CXI(0:IRX),CYI(0:IRX),CMI(0:IRX)
C
C---- Set dynamic pressure
      VINF = QINF
      VREF = QREF
      RHOINF = RHO
      QUREF = 0.5*RHOINF*VREF**2
C
C---- Exclude wake and vortex wake elements from force reporting
      DO IE = 1, NEL
        FX(IE)  = 0.0
        FY(IE)  = 0.0
        MZ(IE)  = 0.0
        CXI(IE) = 0.0
        CYI(IE) = 0.0
        CMI(IE) = 0.0
        IFORC(IE) = 0
        IF(IE.EQ.0 .OR. 
     &     NETYPE(IE).EQ.0 .OR.
     &     NETYPE(IE).EQ.5 .OR.
     &     NETYPE(IE).EQ.6) THEN
C
          CXI(IE) = CX(IE)
          CYI(IE) = CY(IE)
          CMI(IE) = CM(IE)
          FX(IE) = CXI(IE)*QUREF
          FY(IE) = CYI(IE)*QUREF
          MZ(IE) = CMI(IE)*QUREF
          IFORC(IE) = 1
        ENDIF
      END DO
      NELE = NEL
C
      RETURN
      END 

 
      SUBROUTINE DFDCGETROTORVELSIZE(NRV)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETROTORVELSIZE
#endif
C--------------------------------------------------
C     DFDC interface routine - get velocities on current rotor 
C                              or actuator disk
C     Input:      None
C     Output:  
C         NRV       # of rotor stations for slipstream
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
      NRV = NRC
C
      RETURN
      END
 
      SUBROUTINE DFDCGETROTORVEL(ND,
     &                           VINF,
     &                           NRV,YRV,
     &                           VXAI,VRAI,VTAI,
     &                           VXA, VRA, VTA, VMA, VVA,
     &                           WXR, WRR, WTR, WMR, WWR)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETROTORVEL
#endif
C--------------------------------------------------
C     DFDC interface routine - get velocities on downstream side of
C                              rotor or actuator disk
C     Input:      
C         ND        # of disk (rotor or stator)
C     Output:  
C         VINF      freestream velocity
C         NRV       # of rotor stations for slipstream
C         YRV(.)    r coordinate of rotor stations
C---- induced velocity, absolute frame
C         VXAI(.)   axial  (X) induced velocity, absolute frame
C         VRAI(.)   radial (Y) induced velocity, absolute frame
C         VTAI(.)   tangential (theta) induced velocity, absolute frame
C---- total velocity, absolute frame
C         VXA(.)    axial  (X) velocity absolute frame
C         VRA(.)    radial (Y) velocity, absolute frame
C         VTA(.)    tangential (theta) velocity, absolute frame
C         VMA(.)    m-plane velocity, absolute frame
C         VVA(.)    velocity magnitude, absolute frame
C---- total velocity, relative frame
C         WXR(.)    axial  (X) velocity relative frame
C         WRR(.)    radial (Y) velocity, relative frame
C         WTR(.)    tangential (theta) velocity, relative frame
C         WMR(.)    m-plane velocity, relative frame
C         WWR(.)    velocity magnitude, relative frame
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION YRV(*),
     &          VXAI(*),VRAI(*),VTAI(*),
     &          VXA(*), VRA(*), VTA(*), VMA(*), VVA(*),
     &          WXR(*), WRR(*), WTR(*), WMR(*), WWR(*)
C
      IF(ND.LT.1 .OR. ND.GT.NRX) THEN
        WRITE(*,*) 'Disk index out of bounds, DFDCGETROTORVEL ', ND
        RETURN
      ENDIF
C
      VINF = QINF
C
      DO I = 1, NRC
        YRV(I) = YRC(I,ND)
C---- Absolute frame induced velocities
        VXAI(I) = VIND(1,I,ND)
        VRAI(I) = VIND(2,I,ND)
        VTAI(I) = VIND(3,I,ND)
C---- Absolute frame velocities, magnitude and flow angle
        VXA(I) = VABS(1,I,ND)
        VRA(I) = VABS(2,I,ND)
        VMA(I) = SQRT(VXA(I)**2 + VRA(I)**2)
        VTA(I) = VABS(3,I,ND)
        VVA(I) = SQRT(VMA(I)**2 + VTA(I)**2)
cc        ANGA(I) = ATAN2(VTA(I),VMA(I))
C---- Relative (rotor) frame velocities, magnitude and flow angle
        WXR(I) = VREL(1,I,ND)
        WRR(I) = VREL(2,I,ND)
        WMR(I) = SQRT(WXR(I)**2 + WRR(I)**2)
        WTR(I) = VREL(3,I,ND)
        WWR(I) = SQRT(WMR(I)**2 + WTR(I)**2)
cc        ANGR(I) = ATAN2(WMR(I),-WTR(I))
      END DO
      NRV = NRC
C
      RETURN
      END


      SUBROUTINE DFDCGETSURFFLOWSIZE(NCB, NDCT)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETSURFFLOWSIZE
#endif
C--------------------------------------------------
C     DFDC interface routine - get # stations for data on CB 
C                              and duct surfaces       
C     Input:      None
C     Output:  
C---- CB surface points, velocity, Cp 
C         NCB        # of CB panels (with velocities and Cp's)
C---- DUCT surface points, velocity, Cp 
C         NDCT       # of DUCT panels (with velocities and Cp's)
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
C
C---- Data for CB (element #1)
      IEL = 1
      IF(ICFRST(IEL).LE.0) THEN
        WRITE(*,*) 'No control points defined for element ',IEL
        STOP
      ENDIF
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      NCB = IC2 - IC1 + 1
C
C---- Data for duct (element #2)
      IEL = 2
      IF(ICFRST(IEL).LE.0) THEN
        WRITE(*,*) 'No control points defined for element ',IEL
        STOP
      ENDIF
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      NDCT = IC2 - IC1 + 1
C
      RETURN
      END

      SUBROUTINE DFDCGETSURFFLOW(VINF,VREF,RHOINF,
     &                           NCB, 
     &                           XCCB,  YCCB,
     &                           X1CB,  Y1CB,  X2CB,  Y2CB,
     &                           VXCB,  VYCB,  CPCB,
     &                           NDCT, 
     &                           XCDCT, YCDCT,
     &                           X1DCT, Y1DCT, X2DCT, Y2DCT,
     &                           VXDCT, VYDCT, CPDCT)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETSURFFLOW
#endif
C--------------------------------------------------
C     DFDC interface routine - get velocities on CB and duct surfaces       
C     Input:      None
C     Output:  
C         VINF       freestream velocity
C         VREF       reference velocity
C         RHOINF     density
C---- CB surface points, velocity, Cp 
C         NCB        # of CB panels (with velocities and Cp's)
C         XCCB(.)    X (axial)  coordinate of panel centerpoint
C         YCCB(.)    Y (radial) coordinate of panel centerpoint
C         X1CB(.)    X (axial)  coordinate of panel endpoint #1
C         Y1CB(.)    Y (radial) coordinate of panel endpoint #1
C         X2CB(.)    X (axial)  coordinate of panel endpoint #2
C         Y2CB(.)    Y (radial) coordinate of panel endpoint #2
C         VXCB(.)    axial  (X) velocity
C         VYCB(.)    radial (Y) velocity
C         CPCB(.)    Cp on CB (note Cp referenced to VREF, see below)
C---- DUCT surface points, velocity, Cp 
C         NDCT       # of DUCT panels (with velocities and Cp's)
C         XCDCT(.)   X (axial)  coordinate of panel centerpoint
C         YCDCT(.)   Y (radial) coordinate of panel centerpoint
C         X1DCT(.)   X (axial)  coordinate of panel endpoint #1
C         Y1DCT(.)   Y (radial) coordinate of panel endpoint #1
C         X2DCT(.)   X (axial)  coordinate of panel endpoint #2
C         Y2DCT(.)   Y (radial) coordinate of panel endpoint #2
C         VXDCT(.)   axial  (X) velocity
C         VYDCT(.)   radial (Y) velocity
C         CPDCT(.)   Cp on DCT (note Cp referenced to VREF, see below)
C
C     Note:
C         Cp referenced to VINF and VREF -> Cp=(Vinf^2-V^2)/Vref^2
C                         
C--------------------------------------------------
      INCLUDE 'DFDC.INC'
      DIMENSION XCCB(*),  YCCB(*),
     &          X1CB(*),  Y1CB(*),  X2CB(*),  Y2CB(*),
     &          VXCB(*),  VYCB(*),  CPCB(*)
      DIMENSION XCDCT(*), YCDCT(*),
     &          X1DCT(*), Y1DCT(*), X2DCT(*), Y2DCT(*),
     &          VXDCT(*), VYDCT(*), CPDCT(*)
C
C---- Set dynamic pressure
      VINF = QINF
      VREF = QREF
      RHOINF = RHO
C
C---- Dump data for CB (element #1)
      IEL = 1
      IF(ICFRST(IEL).LE.0) THEN
        WRITE(*,*) 'No control points defined for element ',IEL
        STOP
      ENDIF
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      NCB = IC2 - IC1 + 1
      DO IC = IC1, IC2
        I = IC - IC1 + 1
C---- get center point for panel
        XCCB(I) = XC(IC)
        YCCB(I) = YC(IC)
C---- get #1 and #2 corner points at ends of panel
        IP = IPCO(IC)
        X1CB(I) = XP(IP)
        Y1CB(I) = YP(IP)
        IP = IPCP(IC)
        X2CB(I) = XP(IP)
        Y2CB(I) = YP(IP)
C-------- use velocity and Cp on right side of panel
        IF(ISPLOT(IEL).EQ. 1) THEN
C---- use right side of panel for velocity, Cp
          VXCB(I) = QCR(1,IC)
          VYCB(I) = QCR(2,IC)
          CPCB(I) = CPR(IC)
        ELSEIF(ISPLOT(IEL).EQ. 2) THEN
C---- use left side of panel for velocity, Cp
          VXCB(I) = QCL(1,IC)
          VYCB(I) = QCL(2,IC)
          CPCB(I) = CPL(IC)
        ENDIF 
      ENDDO
C
C---- Dump data for duct (element #2)
      IEL = 2
      IF(ICFRST(IEL).LE.0) THEN
        WRITE(*,*) 'No control points defined for element ',IEL
        STOP
      ENDIF
      IC1 = ICFRST(IEL)
      IC2 = ICLAST(IEL)
      NDCT = IC2 - IC1 + 1
      DO IC = IC1, IC2
        I = IC - IC1 + 1
C---- get center point for panel
        XCDCT(I) = XC(IC)
        YCDCT(I) = YC(IC)
C---- get #1 and #2 corner points at ends of panel
        IP = IPCO(IC)
        X1DCT(I) = XP(IP)
        Y1DCT(I) = YP(IP)
        IP = IPCP(IC)
        X2DCT(I) = XP(IP)
        Y2DCT(I) = YP(IP)
C-------- use velocity and Cp on right side of panel
        IF(ISPLOT(IEL).EQ. 1) THEN
C---- use right side of panel for velocity, Cp
          VXDCT(I) = QCR(1,IC)
          VYDCT(I) = QCR(2,IC)
          CPDCT(I) = CPR(IC)
        ELSEIF(ISPLOT(IEL).EQ. 2) THEN
C---- use left side of panel for velocity, Cp
          VXDCT(I) = QCL(1,IC)
          VYDCT(I) = QCL(2,IC)
          CPDCT(I) = CPL(IC)
        ENDIF 
      ENDDO
C
      RETURN
      END


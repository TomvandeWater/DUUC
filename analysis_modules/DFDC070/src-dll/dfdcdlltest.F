      PROGRAM DFDCDLLDRIVER
C-------------------------------------------------------------------
C     DLL Development version 0.70
C
C     This version is based on DFDC version 0.70 10/05
C     This formulation assumes the vortex wakes strengths depend 
C     on a single "wake advance ratio" Wtheta/Wm for each rotor
C     station.   
C 
C-------------------------------------------------------------------
C----- DLL import declarations for DLL linkage
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCINIT
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOAD
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSAVE
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCEXEC
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCEXECTHR
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNUMDISKS
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETRPM
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETRPM
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETQINF
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETQREF
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETXDISK
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETXDISK
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETBGAM
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETCL
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETPITCH
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCDESBLADE
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGENGEOM
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETFLOWCOND
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADNAME
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCSETDEBUG
#endif
C
C-------------------------------------------------------------------
C
      IMPLICIT REAL (A-H,M,O-Z)
C
      CHARACTER*80 FNAME, LINE, ARGP1, ARGP2
      CHARACTER*24 PROMPT
      LOGICAL LCONVG, LDEBUG
C
C---- command input strings
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG
C
C---- local arrays for getting keyboard inputs
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C
 1000 FORMAT(A)
C
C---- print out a GNU license to start off
      CALL WRGNULIC
      LDEBUG = .FALSE.
      VERSION = 0.70
C
      TTOT = 0.0
      TNACEL = 0.0
      TSPEC = 0.0
      ISPEC = 1
      IDSK = 0
C
C---- initialize
      CALL DFDCINIT(LDEBUG)
C
C---- identify this version of DFDC
      WRITE(*,1010) VERSION
 1010 FORMAT(/' ========================='
     &       /'   DFDC DLL Version', F5.2             
     &       /' =========================')
C
C---- try to read input dataset from 1st Unix argument filename
      CALL GETARG0(1,ARGP1)
      FNAME = ARGP1
      IF(FNAME(1:1) .NE. ' ') CALL DFLOADDLL(FNAME)
C
C
C---- print top-level menu
      WRITE(*,1500)
C---- begin user-interaction loop
C.....................................................................
C
  500 CONTINUE
      PROMPT = ' DFDCDLL^'
C
      CALL ASKC(PROMPT,COMAND,COMARG)
      IF(COMAND.EQ.'    ') GO TO 500
C
C---- extract numbers (if any) from command argument string
      DO I=1, 20
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 0
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 0
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
C---------------------------------------------------------
      IF(    COMAND.EQ.'HELP'
     &  .OR. COMAND.EQ.'?   ') THEN
       WRITE(*,1500)
 1500 FORMAT(
     &  //' ========================================================='
     &  //'  QUIT     Exit program',
     &  //'  LOAD f   Load ducted fan case through DLL',
cc     &   /'  LOAF f   Read ducted fan case from file',
     &   /'  SAVE f   Write current ducted fan case to file',
     &  //'  E        Execute (solve) current case',
     &   /'  ETHR ri  Execute (solve) for specified thrust',
     &   /'  RPM  r   Set RPM',
     &   /'  QINF r   Set QINF',
     &   /'  QREF r   Set QREF',
     &   /'  BGAM rr  Set actuator disk circulation B*GAM',
     &   /'  CL   rr  Set blade CL for design',
     &   /'  DESI     Design blade to current CL and B*GAM',
     &   /'  PITC rr  Change blade pitch',
     &   /'  DFLO     Display flow condition',
     &   /'  DFOR     Display duct and CB forces',
     &   /'  DROT     Display rotor state and forces',
     &   /'  DVEL     Display rotor slipstream velocities',
     &   /'  DSRF     Display velocity, Cp on surfaces',
     &   /'  IDSK i   Set current default rotor/stator disk index',
     &   /'  XDSK r   Set rotor or actuator disk axial location',
     &   /'  NAME s   set case name',
     &   /'  NGEO     regenerate geometry', 
     &   /'  INIT     Reinitialize DFDC'
     &   /'  DEBU     toggle debugging flag')
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'QUIT') THEN
        STOP
C
C---------------------------------------------------------
C---- Load case from file using DLL
      ELSEIF(COMAND.EQ.'LOAD') THEN
        WRITE(*,*) 'COMARG ',COMARG
        CALL DFLOADDLL(COMARG)
        CALL DFDCGETNUMDISKS(NROTOR)
C
C---------------------------------------------------------
C---- Load case directly from file
      ELSEIF(COMAND.EQ.'LOAF') THEN
        CALL DFDCLOAD(COMARG)
        CALL DFDCGETNUMDISKS(NROTOR)
C
C---------------------------------------------------------
C---- Save current case to file
      ELSEIF(COMAND.EQ.'SAVE') THEN
        CALL DFDCSAVE(COMARG)
C
C---------------------------------------------------------
C---- Solve case for current circulation
      ELSEIF(COMAND.EQ.'E') THEN
        CALL DFDCEXEC(LCONVG)
        IF(.NOT.LCONVG) WRITE(*,*) '** Not converged **'
C
C---------------------------------------------------------
C---- Solve case for specified thrust (changing pitch)
      ELSEIF(COMAND.EQ.'ETHR') THEN
C---- Specify total thrust if one value specified
        IF(NINPUT.EQ.1) THEN
         TSPEC = RINPUT(1)
         ISPEC = 2
C---- Specify thrust and thrust flag (ISPEC=1 total, ISPEC=2 rotor)
        ELSEIF(NINPUT.EQ.2) THEN
         TSPEC = RINPUT(1)
         ISPEC = IINPUT(2)
        ELSE
         IF(ISPEC.EQ.2) THEN
           TSPEC = TTOT+TNACEL
         ELSE
           TSPEC = TTOT
         ENDIF
         CALL ASKR('Enter Tspec (thrust)^',TSPEC)
         CALL ASKI('Enter Ispec (1 rotor,2 total)^',ISPEC)
        ENDIF
        CALL DFDCEXECTHR(TSPEC,ISPEC,LCONVG)
        IF(.NOT.LCONVG) WRITE(*,*) '** Not converged **'
C
C---------------------------------------------------------
C---- Specify disk # for specific rotor
      ELSEIF(COMAND.EQ.'IDSK') THEN
        IF(NINPUT.EQ.1) THEN
         IDSK = IINPUT(1)
        ENDIF
        IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
          CALL ASKI('Enter disk #^',IDSK)
          IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
        ENDIF
C
C---------------------------------------------------------
C---- Change rotor RPM (or actuator disk RPM)
      ELSEIF(COMAND.EQ.'RPM'
     &  .OR. COMAND.EQ.'R   ') THEN
        IF(NINPUT.GE.1) THEN
         RPMIN = RINPUT(1)
        ELSE
         IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
           CALL ASKI('Enter disk #^',IDSK)
           IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
         ENDIF
         IF(IDSK.GT.0) CALL DFDCGETRPM(IDSK,RRPM)
         RPMIN = RRPM
         CALL ASKR('Enter rotor RPM^',RPMIN)
        ENDIF
        IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
          CALL ASKI('Enter disk #^',IDSK)
          IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
        ENDIF
        IF(IDSK.GT.0) CALL DFDCSETRPM(IDSK,RPMIN)
C
C---------------------------------------------------------
C---- Change freestream velocity
      ELSEIF(COMAND.EQ.'QINF') THEN
        IF(NINPUT.GE.1) THEN
         QINFIN = RINPUT(1)
        ELSE
         CALL DFDCGETFLOWCOND(VINF,VREF,RHOINF,RMUINF,VSOINF,RRPM)
         QINFIN = VINF
         CALL ASKR('Enter freestream velocity QINF^',QINFIN)
        ENDIF
        IF(NINPUT.GE.1) QINFIN = RINPUT(1)
        CALL DFDCSETQINF(QINFIN)
C
C---------------------------------------------------------
C---- Change reference velocity
      ELSEIF(COMAND.EQ.'QREF') THEN
        IF(NINPUT.GE.1) THEN
         QREFIN = RINPUT(1)
        ELSE
         CALL DFDCGETFLOWCOND(VINF,VREF,RHOINF,RMUINF,VSOINF,RRPM)
         QREFIN = VREF
         CALL ASKR('Enter reference velocity QREF^',QREFIN)
        ENDIF
        CALL DFDCSETQREF(QREFIN)
C
C---------------------------------------------------------
C---- Change rotor or actuator disk axial location
      ELSEIF(COMAND.EQ.'XDSK') THEN
        IF(NINPUT.GE.1) THEN
         XRAXIN = RINPUT(1)
        ELSE
         IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
           CALL ASKI('Enter disk #^',IDSK)
           IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
         ENDIF
         IF(IDSK.GT.0) CALL DFDCGETXDISK(IDSK,XDSK)
         XRAXIN = XDSK
         CALL ASKR('Enter disk X location^',XRAXIN)
        ENDIF
        IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
          CALL ASKI('Enter disk #^',IDSK)
          IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
        ENDIF
        IF(IDSK.GT.0) CALL DFDCSETXDISK(IDSK,XRAXIN)
C
C---------------------------------------------------------
C---- Set rotor circulation B*GAM
      ELSEIF(COMAND.EQ.'BGAM') THEN
C---- constant BGAM if one value specified
        IF(NINPUT.EQ.1) THEN
         BGAM1 = RINPUT(1)
         BGAM2 = RINPUT(1)
C---- linear variation of BGAM if two values specified
        ELSEIF(NINPUT.EQ.2) THEN
         BGAM1 = RINPUT(1)
         BGAM2 = RINPUT(2)
        ELSE
         LINE = ' '
 502     CALL ASKS('Enter B*GAMMA^',LINE)
         IF(LINE.NE.' ') THEN
           READ(LINE,*,ERR=502) BGAM1,BGAM2
         ENDIF
        ENDIF
        IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
          CALL ASKI('Enter disk #^',IDSK)
          IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
        ENDIF
        IF(IDSK.GT.0) CALL DFDCSETBGAM(IDSK,BGAM1,BGAM2)
C
C---------------------------------------------------------
C---- Change rotor CL for design
      ELSEIF(COMAND.EQ.'CL') THEN
C---- constant CL if one value specified
        IF(NINPUT.EQ.1) THEN
         CL1 = RINPUT(1)
         CL2 = RINPUT(1)
C---- linear variation of CL if two values specified
        ELSEIF(NINPUT.EQ.2) THEN
         CL1 = RINPUT(1)
         CL2 = RINPUT(2)
        ELSE
         LINE = ' '
 505     CALL ASKS('Enter CL (inboard,outboard)^',LINE)
         IF(LINE.NE.' ') THEN
           READ(LINE,*,ERR=505) CL1,CL2
         ENDIF
        ENDIF
        CALL DFDCSETCL(CL1,CL2)
C
C---------------------------------------------------------
C---- Change rotor blade pitch
      ELSEIF(COMAND.EQ.'PITC') THEN
C---- constant added pitch if one value specified
        IF(NINPUT.EQ.1) THEN
         ANG1 = RINPUT(1)
         ANG2 = RINPUT(1)
C---- linear variation of blade angle change if two values specified
        ELSEIF(NINPUT.EQ.2) THEN
         ANG1 = RINPUT(1)
         ANG2 = RINPUT(2)
        ELSE
         LINE = ' '
 510     CALL ASKS('Enter Pitch change (inboard,outboard)^',LINE)
         IF(LINE.NE.' ') THEN
           READ(LINE,*,ERR=510) ANG1,ANG2
         ENDIF
        ENDIF
        IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
          CALL ASKI('Enter disk #^',IDSK)
          IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
        ENDIF
        IF(IDSK.NE.0) CALL DFDCSETPITCH(IDSK,ANG1,ANG2)
C
C---------------------------------------------------------
C---- Design blade to current circulation and CL
      ELSEIF(COMAND.EQ.'DESI') THEN
        IF(IDSK.LE.0 .OR. IDSK.GT.NROTOR) THEN
          CALL ASKI('Enter disk #^',IDSK)
          IF(IDSK.LT.1 .OR. IDSK.GT.NROTOR) IDSK = 0
        ENDIF
        IF(IDSK.GT.0) CALL DFDCDESBLADE(IDSK)
C
C--------------------------------------------------------------------
C----- Display flow condition
      ELSEIF(COMAND.EQ.'DFLO') THEN
       CALL DISPFLOWCOND
C
C--------------------------------------------------------------------
C----- Display forces on duct and CB
      ELSEIF(COMAND.EQ.'DFOR') THEN
       CALL DISPFORCES
C
C--------------------------------------------------------------------
C----- Display forces and state on rotor
      ELSEIF(COMAND.EQ.'DRO') THEN
       CALL DISPROTORDATA1
C
C--------------------------------------------------------------------
C----- Display forces and state on rotor
      ELSEIF(COMAND.EQ.'DROT') THEN
       LU = 6
       CALL DISPROTORDATA2(LU)
C
C--------------------------------------------------------------------
C---- Display slipstream velocities just downstream of rotor
      ELSEIF(COMAND.EQ.'DVEL') THEN
       CALL DISPROTORVEL
C
C--------------------------------------------------------------------
C---- Display velocity and Cp on CB and duct surfaces
      ELSEIF(COMAND.EQ.'DSRF') THEN
       CALL DISPSURFFLOW
C
C---------------------------------------------------------
C---- Command to regenerate geometry from buffer (REMOVE LATER)
      ELSEIF(COMAND.EQ.'NGEO') THEN
        CALL DFDCGENGEOM(COMARG)
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'NAME') THEN
        CALL ASKS('Enter new case name^',LINE)
        CALL STRIP(LINE,NN)
        CALL DFDCLOADNAME(LINE)
C
C---------------------------------------------------------
C---- Reinitialize DFDC and solver
      ELSEIF(COMAND.EQ.'INIT') THEN
        CALL DFDCINIT(LDEBUG)
C
C---------------------------------------------------------
      ELSEIF(COMAND.EQ.'DEBU') THEN
        LDEBUG = .NOT.LDEBUG
        IF(LDEBUG) THEN
          WRITE(*,*) 'DEBUG output is enabled'
        ELSE
          WRITE(*,*) 'DEBUG output is disabled'
        ENDIF
        CALL DFDCSETDEBUG(LDEBUG)
C
C---------------------------------------------------------
C---------------------------------------------------------
C---- Debugging routines
c      ELSEIF(COMAND.EQ.'AERO') THEN
c       CALL AERO
C
c      ELSEIF(COMAND.EQ.'PANE') THEN
c        CALL PANGEN(.FALSE.)
C
c      ELSEIF(COMAND.EQ.'PPAR') THEN
c        CALL PANGEN(.TRUE. )
C
C---- initialize slipstream and circulation using blade CH,BETA geometry
c      ELSEIF(COMAND.EQ.'RINI') THEN
c       CALL ROTINITBLD
C---- initialize slipstream using B*GAM specified or current circulation
c      ELSEIF(COMAND.EQ.'RING') THEN
c       CALL ROTINITBGAM
C---- initialize slipstream using specified thrust
c      ELSEIF(COMAND.EQ.'RINT') THEN
c        IF(NINPUT.GE.1) THEN
c         THR = RINPUT(1)
c        ELSE
c         THR = THRUST
c         CALL ASKR('Enter rotor thrust^',THR)
c        ENDIF
c       CALL ROTINITTHR(THR)
C----- print slipstream velocities just downstream of rotor
c      ELSEIF(COMAND.EQ.'PVEL') THEN
c       CALL PRTVEL('Rotor velocities',.TRUE.,.TRUE.,.TRUE.)
C---------------------------------------------------------
C---------------------------------------------------------
C
      ELSE 
        WRITE(*,1050) COMAND
 1050   FORMAT(1X,A4,' command not recognized.  Type a "?" for list')
C
      ENDIF
C
C---- finished with last command... go back to menu prompt
      GO TO 500
      END


      SUBROUTINE DISPFLOWCOND
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNAME
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETFLOWCOND
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETRPM
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNUMDISKS
#endif
C-------------------------------------------------------------------
C     Display flow condition
C-------------------------------------------------------------------
      IMPLICIT REAL (A-H,M,O-Z)
      CHARACTER NAME*80
C
      CALL DFDCGETNAME(NAME)
      WRITE(*,*)
      WRITE(*,1010) 
     & '---------------------------------------------------------------'
      WRITE(*,1050) NAME(1:48), '  '
C
      CALL DFDCGETFLOWCOND(VINF,VREF,RHOINF,RMUINF,VSOINF)
C
      WRITE(*,1200) 'Qinf  =',VINF,  '     Rho  =',RHOINF 
      WRITE(*,1200) 'MUinf =',RMUINF,'     Vsou =',VSOINF 
C
      CALL DFDCGETNUMDISKS(NROTOR)
      DO ND = 1, NROTOR
        CALL DFDCGETRPM(ND,RRPM)
        WRITE(*,1250) 'Disk # ',ND,'  RPM   =',RRPM
      END DO
C
 1010 FORMAT(1X, A)
 1050 FORMAT(1X, A,A,A,A,A)
 1200 FORMAT(1X, 5(A, G13.5))
 1250 FORMAT(1X, (A, I3), 4(A, G13.5))
C
      RETURN
      END


      SUBROUTINE DISPFORCES
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNAME
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETFORCEALL
#endif
C-------------------------------------------------------------------
C     Display forces on duct and CB
C-------------------------------------------------------------------
      IMPLICIT REAL (A-H,M,O-Z)
      CHARACTER NAME*80
C
      PARAMETER (IRX=50)
      DIMENSION IFORC(0:IRX),
     &          FX(0:IRX),FY(0:IRX),MZ(0:IRX),
     &          CXI(0:IRX),CYI(0:IRX),CMI(0:IRX)
C
      CALL DFDCGETNAME(NAME)
      CALL DFDCGETFORCEALL(VINF,VREF,RHOINF,QUREF,
     &                     NEL,IFORC,FX,FY,MZ,CXI,CYI,CMI)
C
      WRITE(*,*)
      WRITE(*,1010) 
     & '---------------------------------------------------------------'
      WRITE(*,1050) NAME(1:48), '  '
      WRITE(*,1200) 'Qinf  =', VINF, '     Rho  =', RHOINF 
      WRITE(*,1200) 'Qref  =', VREF, '    q_ref =', QUREF
C
      DO IEL = 0, NEL
C---- Exclude elements without forces from reporting
        IF(IFORC(IEL).EQ.1) THEN
C
          IF(IEL.GT.0) THEN
           WRITE(*,1300) IEL
          ELSE
           WRITE(*,1310)
          ENDIF
C
          WRITE(*,1200) 'Fx   =', FX(IEL), '   Fy =', FY(IEL)
          WRITE(*,1200) 'MZ   =', MZ(IEL)
          WRITE(*,1200) 'Cx   =', CXI(IEL),'   Cy =', CYI(IEL)
          WRITE(*,1200) 'Cm   =', CMI(IEL)
C
        ENDIF
      END DO
C
 1010 FORMAT(1X, A)
 1050 FORMAT(1X, A,A,A,A,A)
 1100 FORMAT(1X, A, 5G13.5)
 1200 FORMAT(1X, 5(A, G13.5))
C
 1300 FORMAT(1X,
     & '-----------------------------------------------------'
     &/' Element', I3 )
 1310 FORMAT(1X,
     &/' Total forces', I3 )
C
      RETURN
      END


      SUBROUTINE DISPROTORDATA1
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETROTORDATA
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNUMDISKS
#endif
C-------------------------------------------------------------------
C     Display rotor data
C-------------------------------------------------------------------
cc      INCLUDE 'DFDC.INC'
      IMPLICIT REAL (A-H,M,O-Z)
      PARAMETER (IRX=50)
C
      LOGICAL LSTALLI
      DIMENSION XI(IRX),CHI(IRX),BETAI(IRX),LSTALLI(IRX)
      DIMENSION CLI(IRX),CDI(IRX),REI(IRX),MACHI(IRX),BGAMI(IRX)
C
      CALL DFDCGETNUMDISKS(NROTOR)
      DO ND = 1, NROTOR
       WRITE(*,*) ' '
C---- Get rotor forces and operating state
       CALL DFDCGETROTORDATA(ND,
     &                       IRTYP,NBLD,XDSK,RRTIP,RRHUB,
     &                       VINF,VREF,RHOINF,RMUINF,VSOINF,
     &                       RRPM,ADV,
     &                       TDIM,TVDIM,TNACEL,
     &                       PDIM,PVDIM,QDIM,
     &                       EFFTOT,EFFIND,EIDEAL,
     &                       CT, CP, TC, PC,
     &                       CTH,CPH,FOMH,SIGMA,CTOS,
     &                       NRPT,
     &                       XI,CHI,BETAI,LSTALLI,
     &                       CLI,CDI,REI,MACHI,BGAMI)
       IF(IRTYP.EQ.1) THEN
         WRITE(*,*) 'Actuator disk at X = ',XDSK
       ELSEIF(IRTYP.EQ.2) THEN
         WRITE(*,*) 'Bladed disk   at X = ',XDSK
       ENDIF
       WRITE(*,*) 'Vinf =',VINF 
       WRITE(*,*) 'RPM  =',RRPM,' #blades =',NBLD
       WRITE(*,*) 'Ttotal =',TDIM,' Tnacel =',TNACEL
       WRITE(*,*) 'Trotor =',TDIM-TNACEL,' Tvis =',TVDIM
       WRITE(*,*) 'Ptotal =',PDIM,' Pvis =',PVDIM
       WRITE(*,*) 'CT     =',CT , ' CP   =',CP 
       WRITE(*,*) 'CTh    =',CTH, ' CPh  =',CPH
       WRITE(*,*) 'FOM    =',FOMH
       WRITE(*,10) 
       DO IR = 1, NRPT
         WRITE(*,15) IR,CHI(IR),BETAI(IR),CLI(IR),CDI(IR),BGAMI(IR)
       ENDDO
      ENDDO
C
 10   FORMAT('  IR     Ch        Beta        CL         CD        BGAM')
 15   FORMAT(I4,6F11.5)
CIIII123456789011234567890112345678901123456789011234567890112345678901
C
      RETURN
      END


      SUBROUTINE DISPROTORDATA2(LU)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNAME
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETROTORDATA
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNUMDISKS
#endif
C-------------------------------------------------------------------
C     Display rotor data and operating state 
C-------------------------------------------------------------------
      IMPLICIT REAL (A-H,M,O-Z)
      PARAMETER (PI=3.1415926535897932384)
      PARAMETER (DTR=PI/180.)
      PARAMETER (IRX=50)
C
      LOGICAL LSTALLI
      CHARACTER SCHAR*1, NAME*80
      DIMENSION XI(IRX),CHI(IRX),BETAI(IRX),LSTALLI(IRX)
      DIMENSION CLI(IRX),CDI(IRX),REI(IRX),MACHI(IRX),BGAMI(IRX)
C
      WRITE(LU,1000)
      CALL DFDCGETNAME(NAME)
      WRITE(LU,1001) NAME
C
      CALL DFDCGETNUMDISKS(NROTOR)
      DO ND = 1, NROTOR
C---- Get rotor forces and operating state
       CALL DFDCGETROTORDATA(ND,
     &                       IRTYP,NBLD,XDSK,RRTIP,RRHUB,
     &                       VINF,VREF,RHOINF,RMUINF,VSOINF,
     &                       RRPM,ADV,
     &                       TDIM,TVDIM,TNACEL,
     &                       PDIM,PVDIM,QDIM,
     &                       EFFTOT,EFFIND,EIDEAL,
     &                       CT, CP, TC, PC,
     &                       CTH,CPH,FOMH,SIGMA,CTOS,
     &                       NRPT,
     &                       XI,CHI,BETAI,LSTALLI,
     &                       CLI,CDI,REI,MACHI,BGAMI)
C
       WRITE(LU,1005) ND, XDSK
       WRITE(LU,1023) 
       WRITE(LU,1010) VINF,RRPM,ADV,
     &                TDIM,PDIM,EFFTOT,
     &                TVDIM,PVDIM,EFFIND,
     &                TNACEL,QDIM,EIDEAL,
     &                NBLD,RRTIP,RRHUB,
     &                RHOINF,VSOINF,RMUINF
C
C---- low advance ratio (hover?) data
       IF(ADV.LT.0.1) THEN
         WRITE(LU,1116) CTH,CPH,FOMH,SIGMA,CTOS
       ENDIF
C
C---- coefficients based on rotational speed
       WRITE(LU,1015) CT, CP, ADV*PI
C---- coefficients based on forward speed
       WRITE(LU,1016) TC, PC, ADV
cc      IF(TERSE) RETURN

       IADD = 1
       IF(IRTYP.EQ.2) THEN
C---- Blade is defined, print blade data
C----- find maximum RE on blade
        REMAX = 0.0
        DO I=1, NRPT
          REMAX = MAX(REI(I),REMAX)
        END DO
        REEXP = 1.0
        IF(REMAX.GE.1.0E6) THEN
          REEXP = 6.0
        ELSEIF(REMAX.GE.1.0E3) THEN
          REEXP = 3.0
        ENDIF
        IF(REEXP.EQ.1.0) THEN
          WRITE(LU,1020) 
         ELSE
          WRITE(LU,1120) IFIX(REEXP)
        ENDIF
C
        DO I=1, NRPT, IADD
          BDEG = BETAI(I)
          XRE  = REI(I)/(10.0**REEXP)
          SCHAR = ' '
          IF(LSTALLI(I)) SCHAR = 's'
          WRITE(LU,1030) I,XI(I),CHI(I),BDEG,CLI(I),SCHAR,
     &                   CDI(I),XRE,MACHI(I),BGAMI(I)
cc     &                   EFFI,EFFP(I)
        END DO
C
       ELSEIF(IRTYP.EQ.1) THEN
C---- Print actuator disk datal
        WRITE(LU,1220) 
        DO I=1, NRPT, IADD
          WRITE(LU,1230)
     &      I,XI(I),MACHI(I),BGAMI(I)
        END DO
       ENDIF
C
      END DO
C
      RETURN
C....................................................................
C
 1000 FORMAT(/1X,75('='))
 1001 FORMAT(' Ducted Fan Rotor Solution:  ', A32)
 1005 FORMAT(' Disk # ',I3,' Xdisk = ',F11.4)
 1021 FORMAT(' (External slipstream present)',19X)
 1023 FORMAT(50X)
 1010 FORMAT( '  Vinf(m/s) :',F10.3,4X,'RPM        :',F11.3,3X,
     &   ' adv. ratio:',F9.4,
     &       /' Thrust(N)  :',G11.3,3X,'Power(W)   :',G11.3,3X,
     &   'Efficiency :',F9.4,
     &       /' Tvisc (N)  :',F11.4,3X,'Pvisc(W)   :',G11.3,3X,
     &   'Induced Eff:',F9.4,
     &       /' Tnacel(N)  :',F11.4,3X,'torQue(N-m):',G11.3,3X,
     &   'Ideal Eff  :',F9.4,
     &       /' no. blades :',I3,  11X,'radius(m)  :',F9.4, 5X,
     &   'hub rad.(m):',F9.4,
     &       /' rho(kg/m3) :',F10.5,4X,'Vsound(m/s):',F9.3, 5X,
     &   'mu(kg/m-s) :',E11.4,
     &      /1X,75('-'))
 1015 FORMAT('by(Rho,N,Dia)',
     &       '   Ct:', F11.5, '      Cp:', F11.5, '    J:', F11.5)
 1016 FORMAT('by(Rho,V,R)  ',
     &       '   Tc:', F11.5, '      Pc:', F11.5, '  adv:', F11.5)
 1116 FORMAT('Hover:       ',
     &       '  CTH:', F11.5, '     CPH:', F11.5,'  FOM:', F11.5,
     &      /'             ',
     &       'Sigma:', F11.5, ' CTH/Sig:', F11.5)
 1017 FORMAT(' Cpv:', F11.5, '    Cpi:', F11.5 )
C---- Rotor data
 1020 FORMAT('Rotor blade:'
     &      /'   i    r/R     c/R  beta(deg)      CL       Cd',
     &        '   RE      ',   '    Mach      B*Gam')
 1120 FORMAT('Rotor blade:'
     &      /'   i    r/R     c/R    beta(deg)    CL       Cd',
     &        '   REx10^',I1,  '    Mach      B*Gam')
 1030 FORMAT(2X,I2,F7.3,F8.4,F8.2,2X,F8.3,1X,A1,F8.4,1X,
     &       F7.2,4X,F7.3,1X,F10.3)
C---- Actuator disk data
 1220 FORMAT('Actuator disk:'
     &      /'   i    r/R      Mach      B*Gam')
 1230 FORMAT(2X,I2,F7.3,4X,F7.3,1X,F10.3)
 2000 FORMAT(/19X,'********** NOT CONVERGED **********'/)
C
C
c 1120 FORMAT(/'   i    r/R     c/R    beta(deg)',
c     & '    CL       Cd     RE        Mach        B*Gam')
c 1030 FORMAT(2X,I2,F7.3,8X,8X,8X,2X,1X,F8.4,1X,
c     &       F7.3,F10.3)
c   i    r/R     c/R    beta(deg)    CL       Cd     RE        Mach      B*Gam\
c   i    r/R     c/R    beta(deg)    CL       Cd     REx10^I   Mach      B*Gam')
cxxiiffffff7fffffff8fffffff8xxxfffffff8xSffffff8xffffff7xxxxffffff7xfffffffff9
cxx2i f7.3   f8.4    f8.2    3x  f8.3  x   f8.4 x  f7.2 4x    f7.3 x  f10.3
C
      END


      SUBROUTINE DISPROTORVEL
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNAME
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETFLOWCOND
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNUMDISKS
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETRPM
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETXDISK
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETROTORVEL
#endif
C----------------------------------------------------------
C     Print out velocities just downstream of rotor 
C     in either:
C        absolute frame induced velocities
C        absolute frame total velocities
C        blade relative frame total velocities
C     Also prints blade relative frame velocities on rotor blade
C----------------------------------------------------------
      IMPLICIT REAL (A-H,M,O-Z)
      PARAMETER (PI=3.1415926535897932384)
      PARAMETER (DTR=PI/180.)
      PARAMETER (IRX=50)
      CHARACTER NAME*80
C
      DIMENSION YRV(IRX),
     &          VXAI(IRX),VRAI(IRX),VTAI(IRX),
     &          VXA(IRX), VRA(IRX), VTA(IRX), VMA(IRX), VVA(IRX),
     &          WXR(IRX), WRR(IRX), WTR(IRX), WMR(IRX), WWR(IRX)
C
      CALL DFDCGETNAME(NAME)
      WRITE(*,20) NAME
      CALL DFDCGETFLOWCOND(VINF,VREF,RHOINF,RMUINF,VSOINF)
      WRITE(*,25) VINF,VREF
C
      CALL DFDCGETNUMDISKS(NROTOR)
      DO ND = 1, NROTOR
       CALL DFDCGETXDISK(ND,XDSK)
       CALL DFDCGETRPM(ND,RRPM)
       WRITE(*,27) ND,XDSK,RRPM
C
C---- Get slipstream velocity data
       CALL DFDCGETROTORVEL(ND,
     &                      VINF,
     &                      NRV,YRV,
     &                      VXAI,VRAI,VTAI,
     &                      VXA, VRA, VTA, VMA, VVA,
     &                      WXR, WRR, WTR, WMR, WWR)
C
C---- Absolute frame induced velocities
C---- Flow angles measured from axial direction
       WRITE(*,30)
       DO I = 1, NRV
         YY = YRV(I)
         VX = VXAI(I)
         VR = VRAI(I)
         VT = VTAI(I)
         VM = SQRT(VX**2 + VR**2)
         VV = SQRT(VM**2 + VT**2)
         ANG = ATAN2(VT,VM)
         WRITE(*,50) YY,VX,VR,VM,VT,VV,ANG/DTR
       END DO
C
C---- Absolute frame velocities
C---- Flow angles measured from axial direction
       WRITE(*,35)
       DO I = 1, NRV
         YY = YRV(I)
         VX = VXA(I)
         VR = VRA(I)
         VT = VTA(I)
         VM = SQRT(VX**2 + VR**2)
         VV = SQRT(VM**2 + VT**2)
         ANG = ATAN2(VT,VM)
         WRITE(*,50) YY,VX,VR,VM,VT,VV,ANG/DTR
       END DO
C
C---- Relative frame velocities downstream of rotor
C---- Blade relative velocity includes rotational speed
C---- Flow angles measured from plane of rotation
       WRITE(*,40)
       DO I = 1, NRV
         YY = YRV(I)
         VX = WXR(I)
         VR = WRR(I)
         VT = WTR(I)
         VM = SQRT(VX**2 + VR**2)
         VV = SQRT(VM**2 + VT**2)
         ANG = ATAN2(VM,-VT)
         WRITE(*,50) YY,VX,VR,VM,VT,VV,ANG/DTR
       END DO
C
C---- Relative frame velocities on blade lifting line
C---- Blade relative velocity includes rotational speed, 1/2 induced swirl
C---- Flow angles measured from plane of rotation
       WRITE(*,60)
       DO I = 1, NRV
         YY = YRV(I)
         VX = WXR(I)
         VR = WRR(I)
         VT = WTR(I) - 0.5*VTAI(I)
         VM = SQRT(VX**2 + VR**2)
         VV = SQRT(VM**2 + VT**2)
         ANG = ATAN2(VM,-VT)
         WRITE(*,50) YY,VX,VR,VM,VT,VV,ANG/DTR
       END DO
      END DO
C
 20   FORMAT(/,A)
 25   FORMAT(' QINF  =',F12.4,
     &      /' QREF  =',F12.4)
 27   FORMAT(/' Disk # ',I3,' at X =',F12.4,' RPM   =',F12.4)
 30   FORMAT(/'Induced vel, flow angles in absolute frame',
     &        ' (downstream of rotor)',
     &      /'     r          Vxi        Vri',
     &   '        Vmi        Vti        Vi     Swirl(deg)')
 35   FORMAT(/'Velocities, flow angles in absolute frame',
     &        ' (downstream of rotor)',
     &      /'     r          Vx         Vr',
     &   '         Vm         Vt         V     Swirl(deg)')
 40   FORMAT(/'Velocities, flow angles relative to rotor frame',
     &        ' (downstream of rotor)',
     &      /'     r          Wx         Wr',
     &   '         Wm         Wt         W       Phi(deg)')
 60   FORMAT(/'Velocities in rotor frame,',
     &        ' on rotor lifting line',
     &       /'flow angle from plane of rotation',
     &      /'     r          Wx         Wr',
     &   '         Wm         Wt         W       Phi(deg)')
C                  12345678901123456789011234567890112345678901
 50   FORMAT(7G11.4)
C
      RETURN
      END



      SUBROUTINE DISPSURFFLOW
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETNAME
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGETFLOWCOND
      !DEC$ ATTRIBUTES C,REFERENCE,DLLEXPORT :: DFDCGETSURFFLOW
#endif
C----------------------------------------------------------
C     Print out velocities and Cp's on CB and duct surfaces
C     Note:
C         Cp referenced to VINF and VREF -> Cp=(Vinf^2-V^2)/Vref^2
C----------------------------------------------------------
      IMPLICIT REAL (A-H,M,O-Z)
      PARAMETER (PI=3.1415926535897932384)
      PARAMETER (DTR=PI/180.)
      PARAMETER (IPX=1000)
      CHARACTER NAME*80
C
      DIMENSION XCCB(IPX),  YCCB(IPX),
     &          X1CB(IPX),  Y1CB(IPX),  X2CB(IPX),  Y2CB(IPX),
     &          VXCB(IPX),  VYCB(IPX),  CPCB(IPX)
      DIMENSION XCDCT(IPX), YCDCT(IPX),
     &          X1DCT(IPX), Y1DCT(IPX), X2DCT(IPX), Y2DCT(IPX),
     &          VXDCT(IPX), VYDCT(IPX), CPDCT(IPX)
C
      CALL DFDCGETNAME(NAME)
      WRITE(*,20) NAME
C
      CALL DFDCGETFLOWCOND(VINF,VREF,RHOINF,RMUINF,VSOINF,RRPM)
      WRITE(*,25) VINF,VREF,RRPM
C
 20   FORMAT(/,A)
 25   FORMAT(' QINF  =',F12.4,
     &      /' QREF  =',F12.4,
     &      /' RPM   =',F12.4)
C
      CALL DFDCGETSURFFLOW(VINF,VREF,RHOINF,
     &                     NCB, 
     &                     XCCB,  YCCB,
     &                     X1CB,  Y1CB,  X2CB,  Y2CB,
     &                     VXCB,  VYCB,  CPCB,
     &                     NDCT, 
     &                     XCDCT, YCDCT,
     &                     X1DCT, Y1DCT, X2DCT, Y2DCT,
     &                     VXDCT, VYDCT, CPDCT)
C
C---- Data on CB surface
      WRITE(*,30) NCB
      WRITE(*,40)
      DO I = 1, NCB
        WRITE(*,50) XCCB(I),YCCB(I),VXCB(I),VYCB(I),CPCB(I),
     &              X1CB(I),Y1CB(I),X2CB(I),Y2CB(I)
      END DO
C
C---- Data on DUCT surface
      WRITE(*,35) NDCT
      WRITE(*,40)
      DO I = 1, NDCT
        WRITE(*,50) XCDCT(I),YCDCT(I),VXDCT(I),VYDCT(I),CPDCT(I),
     &              X1DCT(I),Y1DCT(I),X2DCT(I),Y2DCT(I)
      END DO
C
 30   FORMAT(/'Flow data on CB   surface',
     &       /'NCB   =',I4)
 35   FORMAT(/'Flow data on DUCT surface',
     &       /'NDUCT =',I4)
 40   FORMAT('      Xc          Yc          Vx          Vy    ',
     &       '      Cp          X1          Y1          X2    ',
     &       '      Y2')
C             123456789012123456789012123456789012123456789012123456789012123456789012
 50   FORMAT(9(1X,G11.5))
C
      RETURN
      END


      

      SUBROUTINE DFLOADDLL(FNAME)
#ifdef _WIN32
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADNAME
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADOPER
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADAERO1
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADROTOR
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADACTDISK
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADDRAGOBJ
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCLOADGEOMETRY
      !DEC$ ATTRIBUTES C,REFERENCE,DLLIMPORT :: DFDCGENGEOM
#endif
C---------------------------------------------------------
C     Reads previously saved duct case data and geometry 
C     from file FNAME in DFDC version 0.5 format.  This
C     format separates data into groups by keywords.
C     Once read the geometry is processed from buffer to 
C     paneled geometry.
C---------------------------------------------------------
      IMPLICIT REAL (A-H,M,O-Z)
C
      PARAMETER (PI=3.1415926535897932384)
      PARAMETER (PI2I=0.5/PI)
      PARAMETER (DTR=PI/180.)
C
      CHARACTER*(*)   FNAME
      CHARACTER*128   ANAME, NAMEIN
      CHARACTER*128   IFILE, LINE, FILECB, FILEDCT
      LOGICAL LOPEN, LREPANEL
      DIMENSION RINPUT(10)
C
C---- local arrays for calling AREAD
      PARAMETER (ITX=300,NEX=2)
      DIMENSION XT(ITX,NEX), YT(ITX,NEX)
      DIMENSION NT(NEX)
C
      PARAMETER (IX=300)
      DIMENSION YRDIN(IX), CHRDIN(IX), BETAIN(IX)
      DIMENSION XDDIN(IX), YDDIN(IX),  CDAIN(IX)
      DIMENSION YADIN(IX), BGAMIN(IX)
C
      LOGICAL ERROR
C
      LU = 1
      CALL STRIP(FNAME,NF)
      IF(FNAME.NE.' ') THEN
        OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
        LOPEN = .TRUE.
       ELSE
        RETURN
      ENDIF
      ICNT = 0
      NROTOR = 0
C
C---- Read header line and version from DFDC input file
      CALL RDLINE(LU,LINE,ICNT)
      IF(LINE.EQ.'END' .OR. LINE.EQ.'ERR') GO TO 210
C
      IF(LINE(1:4).NE.'DFDC') THEN
        WRITE(*,*) 'Not DFDC input file...may be duct geometry file'
        CLOSE(LU)
cc        CALL LOADG(FNAME)
        RETURN
      ENDIF
      IV = INDEX(LINE,'Version') 
      READ(LINE(IV+8:LEN(LINE)),*) FILEVERS
      WRITE(*,1005) FILEVERS
C
 1000 FORMAT(A)
 1005 FORMAT(' Reading file from DFDC Version ',F5.2)
 1010 FORMAT(' File  ',A,' not found'/)
 1020 FORMAT(' File  ',A,' has incompatible format'/
     &       ' Loading not completed'/)
C
C
C--- Get Case title from line #2
      CALL RDLINE(LU,LINE,ICNT)
      IF(LINE.EQ.'END' .OR. LINE.EQ.'ERR') GO TO 210
      NAMEIN = LINE
      CALL STRIP(NAMEIN,NNAME)
C==== Use DLL to load 
      CALL DFDCLOADNAME(NAMEIN)
C
C-----------------------------------------------------
C--- Main loop to find keywords and start read actions 
C    for data groups
 10   CALL RDLINE(LU,LINE,ICNT)
      IF(LINE.EQ.'END') THEN
        GO TO 110
      ELSEIF(LINE.EQ.'ERR') THEN 
        GO TO 210
C
C
C---- OPER data
      ELSEIF(INDEX(LINE,'OPER').NE.0) THEN
        CALL RDLINE(LU,LINE,ICNT)
cc        READ(LINE,*,ERR=210) QINFIN,QREFIN,RPMIN
        NINPUT = 4
        CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
        IF(NINPUT.GE.4) THEN
          QINFIN = RINPUT(1)            
          QREFIN = RINPUT(2)            
          RPM1IN = RINPUT(3)            
          RPM2IN = RINPUT(4)            
        ELSEIF(NINPUT.GE.3) THEN
          QINFIN = RINPUT(1)            
          QREFIN = RINPUT(2)            
          RPM1IN = RINPUT(3)
          RPM2IN = 0.0            
        ELSEIF(NINPUT.GE.2) THEN
          QINFIN = RINPUT(1)            
          QREFIN = RINPUT(2)            
          RPM1IN = 0.0
          RPM2IN = 0.0            
        ENDIF
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) RHOIN,VSOIN,RMUIN,ALTHIN
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) XDWKLENIN,NWAKEIN
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) LWRLXIN
        CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDOPER').NE.0) THEN
          WRITE(*,*) 'OPER data read from file'
        ELSE
          WRITE(*,*) 'No ENDOPER to OPER section'
          STOP
        ENDIF
C==== Use DLL to load 
        CALL DFDCLOADOPER(QINFIN,QREFIN,RPM1IN,RPM2IN,
     &                    RHOIN,VSOIN,RMUIN,ALTHIN,
     &                    XDWKLENIN,NWAKEIN,LWRLXIN)
C
C
C---- AERO data
      ELSEIF(INDEX(LINE,'AERO').NE.0) THEN
        NR = NROTOR + 1
C---- reset # of aero sections defined
        N = 0
        CALL DFDCLOADAERO1(NR,N,XISECT,A0DEG,
     &                     CLMAX,CLMIN,
     &                     DCLDA,DCLDA_STALL,DCL_STALL,
     &                     CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                     CMCON,MCRIT,TOC,REREF,REXP)
C--- Read aero section definitions
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NAEROIN
        DO N = 1, NAEROIN
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) XISECT
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) A0DEG,DCLDA,CLMAX,CLMIN
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) DCLDA_STALL,DCL_STALL,CMCON,MCRIT
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,ERR=210) CDMIN,CLDMIN,DCDCL2
          CALL RDLINE(LU,LINE,ICNT)
cc          READ(LINE,*,ERR=210) REREF,REXP
          NINPUT = 4
          CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
          IF(NINPUT.GE.4) THEN
            REREF = RINPUT(1)            
            REXP  = RINPUT(2)            
            TOC   = RINPUT(3)            
            DCDCL2S = RINPUT(4)            
          ELSEIF(NINPUT.GE.3) THEN
            REREF = RINPUT(1)            

            REXP  = RINPUT(2)            
            TOC   = RINPUT(3)            
            DCDCL2S = 0.020
          ELSEIF(NINPUT.GE.2) THEN
            REREF = RINPUT(1)            
            REXP  = RINPUT(2)            
            TOC     = 0.1
            DCDCL2S = 0.020
          ENDIF
C==== Use DLL to load 
          CALL DFDCLOADAERO1(NR,N,
     &                       XISECT,A0DEG,
     &                       CLMAX,CLMIN,
     &                       DCLDA,DCLDA_STALL,DCL_STALL,
     &                       CDMIN,CLDMIN,DCDCL2,DCDCL2S,
     &                       CMCON,MCRIT,TOC,REREF,REXP)
        END DO
        CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDAERO').NE.0) THEN
          WRITE(*,*) 'AERO data read from file for # sections ',NAEROIN
        ELSE
          WRITE(*,*) 'No ENDAERO to AERO section'
          STOP
        ENDIF
C
C
C---- ROTOR data
      ELSEIF(INDEX(LINE,'ROTOR').NE.0) THEN
        NR = NROTOR + 1
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) XDSKIN, NRBLDIN, NRPIN
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NRDIN                
C--- Get rotor definition data
        DO IR = 1, NRDIN
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,END=210,ERR=210) YRDIN(IR),CHRDIN(IR),
     &                                 BETAIN(IR)
        END DO
C
 13     CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDROTOR').NE.0) THEN
          WRITE(*,*) 'ROTOR data read from file'
          IF(NRDIN.GE.2) THEN
           WRITE(*,*) 'Rotor blade defined with #points = ',NRDIN
          ELSE
           WRITE(*,*) 'Rotor blade defined by too few stations'
          ENDIF
        ELSE
          WRITE(*,*) 'No ENDROTOR found'
          STOP
        ENDIF
C==== Use DLL to load 
        NROTOR = NR
        CALL DFDCLOADROTOR(NR,
     &                     XDSKIN,NRBLDIN,NRPIN,
     &                     NRDIN, YRDIN, CHRDIN, BETAIN)
C
C
C---- ACTDISK data
      ELSEIF(INDEX(LINE,'ACTDISK').NE.0) THEN
        NR = NROTOR + 1
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) XDSKIN, NRPIN
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NADIN
C--- Get actuator disk definition data
        DO IR = 1, NADIN
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,END=210,ERR=210) YADIN(IR),BGAMIN(IR)
        END DO   
C
 18     CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDACTDISK').NE.0) THEN
          WRITE(*,*) 'ACTDISK data read from file'
          IF(NADIN.GE.2) THEN
           WRITE(*,*) 'Actuator disk defined with #points = ',NADIN
          ELSE
           WRITE(*,*) 'Actuator disk defined by too few stations'
          ENDIF
        ELSE
          WRITE(*,*) 'No ENDACTDISK found'
          STOP
        ENDIF
C==== Use DLL to load 
        NROTOR = NR
        CALL DFDCLOADACTDISK(NR,
     &                       XDSKIN,NRPIN,
     &                       NADIN, YADIN, BGAMIN)
C
C
C---- DRAGOBJ data
      ELSEIF(INDEX(LINE,'DRAGOBJ').NE.0) THEN
        CALL RDLINE(LU,LINE,ICNT)
        READ(LINE,*,ERR=210) NDDIN
C--- Get drag object definition data
        DO IR = 1, NDDIN
          CALL RDLINE(LU,LINE,ICNT)
          READ(LINE,*,END=210,ERR=210) XDDIN(IR),YDDIN(IR),CDAIN(IR)
        END DO   
C
 19     CALL RDLINE(LU,LINE,ICNT)
        IF(INDEX(LINE,'ENDDRAGOBJ').NE.0) THEN
          WRITE(*,*) 'DRAGOBJ data read from file'
          IF(NDDIN.GE.2) THEN
           WRITE(*,*) 'Drag object defined with #points = ',NDDIN
          ELSE
           WRITE(*,*) 'Drag object defined by too few stations'
          ENDIF
        ELSE
          WRITE(*,*) 'No ENDDRAGOBJ found'
          STOP
        ENDIF
C==== Use DLL to load 
        CALL DFDCLOADDRAGOBJ(NDDIN,
     &                       XDDIN,YDDIN,CDAIN)
C
C
C---- GEOM data (CB and duct coordinates in a single XFOIL multi-element file)
      ELSEIF(INDEX(LINE,'GEOM').NE.0) THEN
C---- Read the combined CB and duct airfoil file
        CALL AREADNR(LU, ITX,NEX, XT,YT,
     &               NT,NTEL,
     &               ANAME,ISPARS,IFTYPE)
        IF(IFTYPE.EQ.0 .OR. NTEL.NE.2) THEN
C----- read error occurred for two elements
         WRITE(*,*) 'Error reading GEOM data'
         STOP
        ENDIF
        WRITE(*,*) 'Coordinates read for ',NTEL,' foils'
C==== Use DLL to load 
        CALL DFDCLOADGEOMETRY(NT(1),XT(1,1),YT(1,1),
     &                        NT(2),XT(1,2),YT(1,2))
C
      ENDIF

      GO TO 10
C
C------------------------------------------------------------------------
 110  CONTINUE
C
C=======================================================================
C---- process buffer foil geometry and set up paneled geometry and grid
      CALL DFDCGENGEOM(FNAME)
C
      IF(LOPEN) CLOSE(LU)
      RETURN
C...............................................................
   98 CONTINUE
      WRITE(*,1050) FNAME(1:NF)
      RETURN
C
   99 CONTINUE
      WRITE(*,1100) FNAME(1:NF)
      IF(LOPEN) CLOSE(LU)
      RETURN
C
 210  CONTINUE
      WRITE(*,1150) ICNT,FNAME(1:NF)
      IF(LOPEN) CLOSE(LU)
      RETURN
C
 1050 FORMAT(/' File OPEN error:  ', A)
 1100 FORMAT(/' File READ error:  ', A)
 1150 FORMAT(/' File READ error on line ',I3,':  ', A)
C
      END ! DFLOADDLL





      SUBROUTINE WRGNULIC
C---- Print out a copy of the GNU license for DFDC
      WRITE(*,10)
 10   FORMAT(
     &/'==============================================================='
     &,'===============',
     &/' DFDC (Ducted Fan Design Code) is an aerodynamic and aeroacoust'
     &,'ic design',
     &/' and analysis tool for aircraft with propulsors in ducted fan'
     &,' configurations.',
     &/' ',
     &/' This software was developed under the auspices and sponsorship'
     &,' of the',
     &/' Tactical Technology Office (TTO) of the Defense Advanced Resea'
     &,'rch',
     &/' Projects Agency (DARPA).',
     &/' ',
     &/' Copyright (c) 2004, 2005, Booz Allen Hamilton Inc., All Rights'
     &,' Reserved'
     &/' ',
     &/' This program is free software; you can redistribute it and/or '
     &,'modify it',
     &/' under the terms of the GNU General Public License as published'
     &,' by the',
     &/' Free Software Foundation; either version 2 of the License, or '
     &,'(at your',
     &/' option) any later version.',
     &/' ',
     &/' This program is distributed in the hope that it will be useful'
     &,', but',
     &/' WITHOUT ANY WARRANTY; without even the implied warranty of'
     &/' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the '
     &,'GNU',
     &/' General Public License for more details.',
     &/' ',
     &/' You should have received a copy of the GNU General Public Lice'
     &,'nse along',
     &/' with this program; if not, write to the Free Software Foundati'
     &,'on, Inc.,',
     &/' 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA',
     &/' ',
     &/' Authors: Harold Youngren (guppy@maine.rr.com), Mark Drela (dre'
     &,'la@mit.edu)',
     &/' Program Management: Brad Tousley, Paul Eremenko (eremenko@alum'
     &,'.mit.edu)',
     &/'==============================================================='
     &,'===============')
      RETURN
      END


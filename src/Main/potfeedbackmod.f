      SUBROUTINE POTFEEDBACKMOD(TIMEHERE,RISF2,FMP,MP2,POTFB,FMPPOT)

      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      include 'tlist.h'

      REAL*8 TFB(5000),MCFB(5000),MATFB(1:3,5000)
      integer NFILES,TFBI,NPOT
      CHARACTER*200 FILEFB(2000)

      REAL*8 TSTAR2,SXFB,MSCL1,TSG,TRET
      CHARACTER*200 PATHALL
      COMMON/extpot/ TFB,MATFB,MCFB,FILEFB,TSTAR2,SXFB,MSCL1
      COMMON/extpot/ TSG,TRET,NFILES,TFBI,NPOT,PATHALL
      
      REAL TIMESN(10000),TCOL,MASSSN2(10000)!TSN,TCOL,MASSSN
      INTEGER ITSN,TITSN,MASSSN(10000)  !ISN,TOTSN,NAMESN
      COMMON/TMSN/ TIMESN,TCOL,MASSSN,MASSSN2,ITSN,TITSN,LECT
      
      REAL*8 TIMEHERE,RISF2,FMP,RHO_UNIF,MP2,POTFB,R_UNIF
      REAL*8 TIMEMYR,RISF,CTE,CTE2,G,MPDOT2,TDELAY2
      CHARACTER*200 FILE1,FILENOW
      INTEGER IFILE,JFILE,FMPPOT
      

      TIMEMYR=TIMEHERE*TSTAR2+TRET+TSG !NBODYTIME TO MYR AND VIRIALIZATION
      !MPDOT2=TSG/TSTAR2             !INITIAL TIME IN NB
      RISF=(RISF2**0.5)             !R=SQRT(R**2) STAR RADII
      !TDELAY2=-TRET/TSTAR2      !TIME OF CHANGE PLUMMER AFTER VIRIALIZATION
      !TCOL = 1000.
      IF(TIMEHERE.EQ.0.AND.NFILES.LT.10) THEN
         OPEN(UNIT=102,FILE="initcond.ini",STATUS='old') !OPENING FILE
         READ(102,*) PATHALL   !PATH TO POTENTIAL FILES
         READ(102,*) TSTAR2    !Scaling factor TIME
         READ(102,*) SXFB      !Scaling factor RADIUS CLUSTER INDEPENDENT OF EXTERNAL POTENTIAL
         READ(102,*) MSCL1     !Cluster Mass [Mo]
         READ(102,*) TSG       !TIME 2NDGEN STARTED (0 FOR 1STGEN) 
         READ(102,*) TRET      !VIRIALIZATION TIME, ALWAYS NEGATIVE
         READ(102,*) TCOL      !COLLAPSE TIME OR TEND FOR 2NDGEN
         CLOSE(102)
         TIMEMYR=TIMEHERE*TSTAR2+TRET+TSG
         !print*, PATHALL
         !print*,TSTAR2,SXFB,MSCL1,TSG,TRET
         TCOL = (TCOL-TRET)/TSTAR2
         TDELAY2=-TRET/TSTAR2   !TIME OF CHANGE PLUMMER AFTER VIRIALIZATION
         MPDOT2=TSG/TSTAR2
         TIMEMYR=TIMEHERE*TSTAR2+TRET+TSG
         FILE1="file_time_Mcloud_file.dat" !time(myr),M_cloud,file
         FILENOW=TRIM(PATHALL)//TRIM(FILE1) !COMPLETE PATH
         OPEN(UNIT=11,FILE=FILENOW,STATUS='old') !OPENING FILE
         NFILES=0                                !FOR COUNTING LINES
         DO IFILE=1,2000
            READ(11,*,END=11) TFB(IFILE),MCFB(IFILE),FILEFB(IFILE) !time(myr),M_cloud,file1
            NFILES=NFILES+1                                        !COUNTING LINES
         END DO
 11      CLOSE(11)
         OPEN(14,FILE="TableSN.dat",STATUS="OLD")                  !INFORMATION FOR SN REMOVING
         TITSN=0                                                   !
         ITSN=0                                                   !
         DO JFILE=1,10000
            READ(14,*,END=14) MASSSN(JFILE),TIMESN(JFILE),MASSSN2(JFILE) !SN MASS, TIME OF SN REMOVING
            MASSSN2(JFILE)=MASSSN2(JFILE)/MSCL1 !SCALING SN MASS
            TIMESN(JFILE)=(TIMESN(JFILE)-TRET)/TSTAR2  !(myr-->nb - ti(NB)) SCALING SN TIME TO NB
            !print*,TIMESN(JFILE),MASSSN(JFILE),"TableSN", MASSSN2(JFILE)   !FOR CHECKING DEBUGGING OFF
            IF(TIMESN(JFILE).GT.0.AND.ITSN.EQ.0) ITSN=JFILE ! CHECK
            IF(TIMESN(TITSN).GT.TCOL) EXIT
            TITSN=TITSN+1                                          !NUMBER OF POSSIBLE SN TO REMOVE
         END DO
 14      CLOSE(14)
         IF(TSG.GT.0) THEN                  !IF TIME SECOND GEN > 0                    
            DO IFILE=1,NFILES               !
               IF(TFB(IFILE).GT.TSG) THEN   !
                  TFBI=IFILE                !
                  EXIT                      !                       !CHECK!!!
               END IF                       
            END DO 
         ELSE
            DO IFILE=1,NFILES               !
               IF(TFB(IFILE).GT.TIMEMYR) THEN   !
                  TFBI=IFILE                !
                  EXIT                      !                       !CHECK!!!
               END IF                       
            END DO                          !IF TIME NO SECOND GEN, TIME FEED BACK = 0 
         END IF
         PRINT*,""                                          !INFO OF INFORMATION FOR SCALING EXTERNAL POTENTIAL
         PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         PRINT*,"INFO USED TO SCALE EXTERNAL POTENTIAL"
         PRINT*,"SCALING TIME",TSTAR2,"SCALING POSITION",SXFB !TIME SCALING, POSITION SCALING
         PRINT*,"TOTAL NUMBER OF SN",TITSN,"UNTIL",
     &        TIMESN(TITSN)*TSTAR2+TRET,"MYR" !,TCOL,mASSSN(TITSN)
 
         IF(TIMESN(TITSN).LT.TCOL) THEN
            PRINT*,"WARNING SN"
            PRINT*,"LAST SN OCURRING BEFORE END OF SIMULATION",
     &           MASSSN(TITSN)
            STOP
         END IF
         PRINT*,"TOTAL CLUSTER MASS",MSCL1,"[SOLARMASS]"            !CLUSTER MASS
         !PRINT*,"INITIAL RADIUS OF CLOUD",R_UNIF,"[PC]"             !RADIUS OF CLOUD (USELESS)
         !PRINT*,"DENSITY SPHERE",RHO_UNIF,"[SOLARMASS/PC**3]"       !RHO_0 OF CLOUD (USELESS)
         !&        ,"USELESS FOR EXPONENTIAL PROFILE"         
         PRINT*,"TIME FOR VIRIALIZATION",INT(-TRET),"[MYR]"         !TIME BEFORE OF EXTERNAL POTENTIAL EVOLUTION
         IF(TSG.GT.0) THEN
            PRINT*,"TIME OF SECOND GENERATION",TSG,"[MYR]" !INFORMATION OF SECOND GENERATION
         ELSE
            PRINT*,"NO SECOND GENERATION"
         END IF
         PRINT*,"TIME OF COLLAPSING",TCOL/TSTAR2,"MYR"
         PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         PRINT*,""
         MP2=MCFB(TFBI)/(MSCL1) 
         GOTO 101               !THIS IS DONE ONLY ONE TIME IN ALL THE VIRIALIZATION TIME
         
      END IF

      !IF(TIMEMYR.GT.TRET)THEN
      IF(TIMEMYR.GE.TFB(NFILES))THEN !IF TIME >= LASTSNAPSHOT 
         FMP = 0.               !EXTERNAL FORCES = 0
         POTFB = 0.             !EXTERNAL POTENTIAL = 0
         MCFB(TFBI)=0.
         IF(TFBI.EQ.(NFILES-1)) THEN
            PRINT*,"LAST WF FILE REACHED, ALL THE GAS IS GONE",TIMEMYR
            TFBI=NFILES
         END IF
         GOTO 104               !SKIP ALL THE REST (FASTER)
      END IF
      !END IF

      IF(TIMEMYR.GT.TRET)THEN                   !IF THE CLUSTER IS VIRIALIZED
         DO IFILE=TFBI,NFILES-1                 !CHECK IF GO TO NEXT SNAPSHOT
            IF(TIMEMYR.LT.TFB(TFBI+1)) GOTO 103 !IF TIME HAS NOT REACHED THE NEXT SNAPSHOT
            IF(TIMEMYR.GE.TFB(IFILE).AND.TIMEMYR.LT.TFB(IFILE+1))THEN !IF TIME HAS REACHED THE NEXT SNAPSHOT
               TFBI=IFILE                                            !USE NEXT SNAPSHOT
               MP2=MCFB(TFBI)/(MSCL1)                                !USE NEW INTEGRATED MASS
               GOTO 101                                              !STOP CHECKING
            END IF            
         END DO
      END IF

 101  IF(TIMEMYR.LE.TRET.AND.NPOT.GT.0) GOTO 103 !IF THE PROGRAM HAS ALREADY READ THE FIRST SNAPSHOT AND UNTIL GET VIRIAL STATE
      
      FILENOW=TRIM(PATHALL)//TRIM(FILEFB(TFBI))                          !PATH OF POTENTIAL FILE
      IF(TIMEMYR.EQ.TRET) PRINT*,"WORKING IN FOLDER ",PATHALL           !PRINTING WORKING PATH ONLY ONCE
      PRINT*,"READING WF FILE FOR FORCES AND POTENTIAL ",FILEFB(TFBI) !PRINTING WORKING POTENTIAL FILE
     &     ,TIMEMYR,"MYR",tfbi,mp2,MCFB(TFBI)
      OPEN(UNIT=12,FILE=FILENOW,STATUS='OLD')                         !OPENING POTENTIAL FILE
      NPOT=0                                                          !FOR COUNTING LINES POTENTIAL FILE
      G = 4.302*1E-3                                                  !GRAVITATIONAL CONSTANT 
      CTE=1/(MCFB(TFBI)*G)                                            !1/G*CLOUD_MASS FOR NBODY UNITS
      
      DO IFILE=1,5000
         READ(12,*,END=12) MATFB(1:3,IFILE)          !RADII,FORCES,POTENTIAL           
         MATFB(1,IFILE)=MATFB(1,IFILE)/(SXFB)        !radii /sx    --> Nbody Units 
         MATFB(2,IFILE)=MATFB(2,IFILE)*SXFB**2*CTE   !F*sx**2/(mass*G)    --> Nbody Units
         MATFB(3,IFILE)=MATFB(3,IFILE)*SXFB*CTE      !POT*sx/(mass*G)  --> Nbody Units
         NPOT=NPOT+1                                 !total lines file (5000)
      END DO
 12   CLOSE(12)
      
! 103  IF(TIMEHERE.LT.TDELAY2) THEN               !I THINK IS NOT USEFULL ANYMORE
!        MP2=MCFB(TFBI)/(MSCL1)
!      ELSE
!        MP2=MP
!      END IF

! FOR INITIAL UNIFORM POTENTIAL (USELESS FOR EXPONENTIAL PROFILE)
!      IF(TFBI.EQ.1.OR.TIMEMYR.LE.0) THEN 
!         FMP=4.*3.141592653589793*RHO_UNIF/(3.*MCFB(TFBI)) !remove M_cloud not scaled
!         POTFB=(3*R_UNIF**2-(RISF*SXFB)**2)/((2*R_UNIF**3)) !POTENTIAL WITHOUT M & G
!         NPOT=5000              !just to skip previous steps again
         !!!!!!IF(RISF.LT.MATFB(1,NPOT))
!         GOTO 104
!         print*,"paso"
!      END IF
      
 103  IF(RISF.GE.MATFB(1,NPOT).AND.FMPPOT.EQ.1) THEN !FORCES CALCULATION OUT OF SPHERE
         IFILE=NPOT-1
         FMP=((MATFB(2,IFILE+1)-MATFB(2,IFILE))/     !LINEAL FITTNG CHECK!!
     &        (MATFB(1,IFILE+1)-MATFB(1,IFILE)))*
     &        (RISF-MATFB(1,IFILE))+MATFB(2,IFILE)
         FMP=FMP/RISF                                !F/R --> M_CLOUD/R**3 *X(I)
         POTFB=0.                                    !POTENTIAL OUTISIDE = 0 !NOT NECESSARY FOR FMPOT = 1
         IF(FMP.LT.0) FMP=0.                         !FMP CAN NOT BE > 0
         GOTO 104                                    !SKIP THE REST
      END IF
      IF(RISF.GE.MATFB(1,NPOT).AND.FMPPOT.EQ.2) THEN        !POTENTIAL CALCULATION OUT OF SPHERE
         POTFB=0                                            !NOT INCLUDE BEFORE BECAUSE OF OPTIMIZATION, check later!
         FMP=0.                                             !NOT NECESSARY FOR FMPOT = 2
         GOTO 104
      END IF
      IF(RISF.LT.MATFB(1,1)) THEN           !FORCES AND POTENTIAL CALCULATION INSIDE BUBLE
         FMP=0.0                            !FORCES INSIDE BUBLE = 0
         POTFB=MATFB(3,1)                   !POTENTIAL = CTE !CHECK
         GOTO 104                           !SKIP THE REST
      END IF
      IF(FMPPOT.EQ.1) THEN                  !FORCES CALCULATIONS 
         DO IFILE=1,NPOT
            IF(RISF.GE.MATFB(1,IFILE).AND.RISF.LT.MATFB(1,IFILE+1))THEN !LOOKING FOR RADII
               FMP=((MATFB(2,IFILE+1)-MATFB(2,IFILE))/                  !LINEAL FITTNG 
     &              (MATFB(1,IFILE+1)-MATFB(1,IFILE)))*
     &              (RISF-MATFB(1,IFILE))+MATFB(2,IFILE)
               FMP=FMP/RISF                                             !F/R --> M_CLOUD/R**3 *X(I)
               POTFB=0.                                                 !NOT NECESSARY FOR FMPOT = 1
               GOTO 104                                                 !SKIP THE REST
            END IF
         END DO
      END IF
      IF(FMPPOT.EQ.2) THEN                                              !POTENTIAL CALCULATIONS
         DO IFILE=1,NPOT
            IF(RISF.GE.MATFB(1,IFILE).AND.RISF.LT.MATFB(1,IFILE+1))THEN !LOOKING FOR RADII
               POTFB=((MATFB(3,IFILE+1)-MATFB(3,IFILE))/                !LINEAL FITTNG 
     &              (MATFB(1,IFILE+1)-MATFB(1,IFILE)))*
     &              (RISF-MATFB(1,IFILE))+MATFB(3,IFILE)
               FMP=0                                                    !NOT NECESSARY FOR FMPOT = 2
               GOTO 104                                                 !SKIP THE REST
            END IF
         END DO
      END IF

 104  MP2=MCFB(TFBI)/(MSCL1)
      IF(MP2.NE.MP.AND.MP2.GE.0.AND.ABS(MP).LT.1000.) WRITE (6,*)
     &     "PLUMMER POTENTIAL NEW: MP=",MP2, "TIME=",TIMEMYR,"MYR"
     &     ,tfbi,MP
!     &     ,(TFB(NFILES)-TRET),TFB(NFILES),TFB(NFILES-1)
 !70   FORMAT (/,12X,'PLUMMER POTENTIAL NEW:    MP =',F7.3)
 !     IF(MP2.NE.MP.AND.MP2.GE.0) WRITE (6,*) "TIME",TIMEMYR,"MYR" 
      !print*,mp2,"caxa,tfbi
!print*,mp2,"caxa"

      RETURN
      
      END SUBROUTINE

#include <fintrf.h>

C
C Gateway function for pvaws2.F
C calls other function in file cpvaws2.F
C
        SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C
       IMPLICIT NONE
       MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

       INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C
       MWPOINTER MXGETPR
       MWPOINTER MXCREATEDOUBLEMATRIX

       real*8  mxGetScalar

C
C-----------------------------------------------------------------------
C

       MWPOINTER Y_IN_PR
       MWPOINTER THETA_IN_PR
       MWPOINTER III_IN_PR
       MWPOINTER JJJ_IN_PR

       MWPOINTER INVCOV_IN_PR
       MWPOINTER BI_IN_PR
       MWPOINTER WGHT_IN_PR

       MWPOINTER BI_OUT_PR
       MWPOINTER THNEW_OUT_PR

       MWPOINTER LWGHT_LOCAL_PR
       MWPOINTER SWJY_LOCAL_PR

       MWPOINTER LWGHT_LOCAL
       MWPOINTER SWJY_LOCAL

C-----------------------------------------------------------------------
C
C      LOCAL VARIABLES

       integer*4 N1_IN, N2_IN, N3_IN, PSIZE_IN,
     +           NV_IN, NVD_IN, NMASK_IN, NCORES_IN

       REAL*8 LAMBDA_IN, HAKT_IN, WGHT_IN(2)

       integer*4 IH1,IH2,IH3,DWL1,DWL2,DWL3

       MWSIZE NMASK, NV, NVD, ONE, DWL, NCORES

       real*8 w(2)
       mwsize numrel
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
       IF (NRHS .NE. 17) THEN
         CALL MEXERRMSGTXT('VPAWS requires 17 input arguments')
       ELSEIF (NLHS .GT. 2) THEN
         CALL MEXERRMSGTXT('VAWS requires 2 output argument')
       ENDIF

       CALL MEXPRINTF(' Start Iteration of patchbased VAWS ... ')

       NV_IN = INT(MXGETSCALAR(PRHS(8)))
       NVD_IN = INT(MXGETSCALAR(PRHS(9)))
       N1_IN = INT(MXGETSCALAR(PRHS(10)))
       N2_IN = INT(MXGETSCALAR(PRHS(11)))
       N3_IN = INT(MXGETSCALAR(PRHS(12)))
       NMASK_IN = INT(MXGETSCALAR(PRHS(13)))
       PSIZE_IN = INT(MXGETSCALAR(PRHS(14)))
       HAKT_IN = MXGETSCALAR(PRHS(15))
       LAMBDA_IN = MXGETSCALAR(PRHS(16))
       NCORES_IN = INT(MXGETSCALAR(PRHS(17)))

C       getting all the other INPUT pointers

       Y_IN_PR = MXGETPR(PRHS(1))
       III_IN_PR = MXGETPR(PRHS(2))
       JJJ_IN_PR = MXGETPR(PRHS(3))
       THETA_IN_PR = MXGETPR(PRHS(4))
       INVCOV_IN_PR = MXGETPR(PRHS(5))
       BI_IN_PR = MXGETPR(PRHS(6))
       WGHT_IN_PR = MXGETPR(PRHS(7))

       numrel=2
       call MXCOPYPTRTOREAL8(WGHT_IN_PR,w,numrel)
       ih3=FLOOR(HAKT_IN/w(2))
       ih2=FLOOR(HAKT_IN/w(1))
       ih1=FLOOR(HAKT_IN)

       DWL1 = 2*IH1+1
       DWL2 = 2*IH2+1
       DWL3 = 2*IH3+1
       DWL = DWL1*DWL2*DWL3
C CREATE MATRICES FOR RETURN ARGUMENT
       ONE=1
       NMASK=NMASK_IN
       NV=NV_IN
       NVD=NVD_IN
       NCORES=NCORES_IN

C     BI
      PLHS(1) = MXCREATEDOUBLEMATRIX(ONE, NMASK, 0)

C     THNEW
      PLHS(2) = MXCREATEDOUBLEMATRIX(NV, NMASK, 0)
C
C ASSIGN POINTERS TO THE OUTPUT PARAMETERS
C
      BI_OUT_PR = MXGETPR(PLHS(1))
      THNEW_OUT_PR = MXGETPR(PLHS(2))

C
      LWGHT_LOCAL = mxcreatedoublematrix(ONE,DWL,0)
      SWJY_LOCAL = mxcreatedoublematrix(NV,NCORES,0)

      LWGHT_LOCAL_PR = mxgetpr(LWGHT_LOCAL)
      SWJY_LOCAL_PR = mxgetpr(SWJY_LOCAL)
C
C DO THE ACTUAL COMPUTATIONS ALSO USING VAL
C

        CALL CPVAWS2(%val(Y_IN_PR),
     +               NV_IN,
     +               NVD_IN,
     +               N1_IN,
     +               N2_IN,
     +               N3_IN,
     +               NMASK_IN,
     +               %val(III_IN_PR),
     +               %val(JJJ_IN_PR),
     +               HAKT_IN,
     +               LAMBDA_IN,
     +               %val(THETA_IN_PR),
     +               %val(INVCOV_IN_PR),
     +               %val(BI_IN_PR),
     +               %val(BI_OUT_PR),
     +               %val(THNEW_OUT_PR),
     +               NCORES_IN,
     +               %val(LWGHT_LOCAL_PR),
     +               %val(WGHT_IN_PR),
     +               %val(SWJY_LOCAL_PR),
     +               PSIZE_IN)

       call mxDestroyArray(SWJY_LOCAL)
       call mxDestroyArray(LWGHT_LOCAL)
      CALL MEXPRINTF('iteration completed: time to return \n')

      RETURN
      END

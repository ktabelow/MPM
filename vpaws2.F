#include <fintrf.h>

C
C Gateway function for pvaws2.F
C calls other function in file cpvaws2.F
C
        SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C
C
       IMPLICIT NONE
       MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

       INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C
C
       MWPOINTER MXGETPR
       MWPOINTER MXCREATEDOUBLEMATRIX
       MWPOINTER MXCREATENUMERICMATRIX
       MWPOINTER MXCREATENUMERICARRAY

       mwPointer mxGetNumberOfElements
       mwPointer mxGetDimensions
       real*8  mxGetScalar
       integer*4 mxIsNumeric, mxClassIDFromClassName
       character*8 double

C-----------------------------------------------------------------------
C

       MWSIZE MXGETM, MXGETN
       mwSize mxGetNumberOfDimensions
C
C-----------------------------------------------------------------------
C
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
       MWPOINTER THPATCH_LOCAL_PR
       MWPOINTER INVCOVP_LOCAL_PR
       MWPOINTER BIIPATCH_LOCAL_PR
       MWPOINTER SWJY_LOCAL_PR

C-----------------------------------------------------------------------
C
C      LOCAL VARIABLES



       INTEGER*4 N1_IN, N2_IN, N3_IN, NP1_IN, NP2_IN, NP3_IN,
     +           NV_IN, NVD_IN, NMASK_IN, NCORES_IN


       REAL*8 LAMBDA_IN, HAKT_IN, WGHT_IN(2)


       INTEGER*4 IH

       MWSIZE NMASK, NV, NVD, ONE, DWL, PSIZE, NCORES, PNCORES

       real*8 x
       integer*4 nx,i
       mwsize numrel
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
       IF (NRHS .NE. 19) THEN
         CALL MEXERRMSGTXT('VPAWS requires 19 input arguments')
       ELSEIF (NLHS .GT. 2) THEN
         CALL MEXERRMSGTXT('VAWS requires 2 output argument')
       ENDIF


C       getting the values of n1, n2, n3, nV, NCORES, NMASK
C       HAKT, lambda
C       using get scalar
C       need to cast since it returns a double
        CALL MEXPRINTF(' BEFORE getting SCALARS \n')

       NV_IN = INT(MXGETSCALAR(PRHS(8)))
       NVD_IN = INT(MXGETSCALAR(PRHS(9)))
       N1_IN = INT(MXGETSCALAR(PRHS(10)))
       N2_IN = INT(MXGETSCALAR(PRHS(11)))
       N3_IN = INT(MXGETSCALAR(PRHS(12)))
       NMASK_IN = INT(MXGETSCALAR(PRHS(13)))
       NP1_IN = INT(MXGETSCALAR(PRHS(14)))
       NP2_IN = INT(MXGETSCALAR(PRHS(15)))
       NP3_IN = INT(MXGETSCALAR(PRHS(16)))
       HAKT_IN = MXGETSCALAR(PRHS(17))
       LAMBDA_IN = MXGETSCALAR(PRHS(18))
       NCORES_IN = INT(MXGETSCALAR(PRHS(19)))

       CALL MEXPRINTF(' AFTER getting SCALARS \n')

C       getting all the other INPUT pointers


       Y_IN_PR = MXGETPR(PRHS(1))
       III_IN_PR = MXGETPR(PRHS(2))
       JJJ_IN_PR = MXGETPR(PRHS(3))
       THETA_IN_PR = MXGETPR(PRHS(4))
       INVCOV_IN_PR = MXGETPR(PRHS(5))
       BI_IN_PR = MXGETPR(PRHS(6))
       WGHT_IN_PR = MXGETPR(PRHS(7))

       IH = INT(HAKT_IN)
       DWL = 2*IH+1
       DWL = DWL*DWL*DWL
       CALL MEXPRINTF(' AFTER getting OTHER input \n')

       CALL MEXPRINTF(' BEFORE CREATING OUTPUT ')
C CREATE MATRICES FOR RETURN ARGUMENT
       ONE=1
       NMASK=NMASK_IN
       NV=NV_IN
       NVD=NVD_IN
       NCORES=NCORES_IN
       PSIZE=NP1_IN*NP2_IN*NP3_IN
       PNCORES=NCORES_IN*NP1_IN*NP2_IN*NP3_IN
C
C
C     BI
      PLHS(1) = MXCREATEDOUBLEMATRIX(ONE, NMASK, 0)
      CALL MEXPRINTF(' BI done \n')

C     THNEW
      PLHS(2) = MXCREATEDOUBLEMATRIX(NV, NMASK, 0)
      CALL MEXPRINTF(' THNEW done \n')

      CALL MEXPRINTF(' AFTER PLACE OUTPUT ')
C
C ASSIGN POINTERS TO THE OUTPUT PARAMETERS
C
      BI_OUT_PR = MXGETPR(PLHS(1))
      THNEW_OUT_PR = MXGETPR(PLHS(2))

C
      LWGHT_LOCAL_PR = mxcreatedoublematrix(ONE,DWL,0)
      SWJY_LOCAL_PR = mxcreatedoublematrix(NV,NCORES,0)
      THPATCH_LOCAL_PR = mxcreatedoublematrix(NV,PNCORES,0)
      INVCOVP_LOCAL_PR = mxcreatedoublematrix(NVD,PNCORES,0)
      BIIPATCH_LOCAL_PR = mxcreatedoublematrix(ONE,PNCORES,0)

C      LWGHT_LOCAL_PR = mxgetpr(LWGHT_LOCAL)
C      SWJY_LOCAL_PR = mxgetpr(SWJY_LOCAL)
C      THPATCH_LOCAL_PR = mxgetpr(THPATCH_LOCAL)
C      INVCOVP_LOCAL_PR = mxgetpr(INVCOVP_LOCAL)
C      BIIPATCH_LOCAL_PR = mxgetpr(BIIPATCH_LOCAL)

C
C DO THE ACTUAL COMPUTATIONS ALSO USING VAL
C
       numrel=1
C       call MXCOPYPTRTOREAL8(Y_IN_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(THETA_IN_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(INVCOV_IN_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(BI_IN_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(WGHT_IN_PR,x,numrel)
C       call MXCOPYPTRTOINTEGER4(III_IN_PR,nx,numrel)
C       call MXCOPYPTRTOINTEGER4(JJJ_IN_PR,nx,numrel)

C       call MXCOPYPTRTOREAL8(LWGHT_LOCAL_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(SWJY_LOCAL_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(THPATCH_LOCAL_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(INVCOVP_LOCAL_PR,x,numrel)
C       call MXCOPYPTRTOREAL8(BIIPATCH_LOCAL_PR,x,numrel)

       CALL MEXPRINTF('BEFORE CALL \n')


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
     +               NP1_IN,
     +               NP2_IN,
     +               NP3_IN,
     +               %val(THPATCH_LOCAL_PR),
     +               %val(INVCOVP_LOCAL_PR),
     +               %val(BIIPATCH_LOCAL_PR))

       CALL MEXPRINTF(' call END \n')
       i=0
       DO WHILE(i.lt.1e8)
         i=i+1
       END DO
       call mxDestroyArray(SWJY_LOCAL_PR)
       i=0
       DO WHILE(i.lt.1e8)
         i=i+1
       END DO
       CALL MEXPRINTF(' destroyed 4\n')
       call mxDestroyArray(BIIPATCH_LOCAL_PR)
       i=0
       DO WHILE(i.lt.1e8)
         i=i+1
       END DO
       CALL MEXPRINTF(' destroyed 1 \n')

       call mxDestroyArray(INVCOVP_LOCAL_PR)
       i=0
       DO WHILE(i.lt.1e8)
         i=i+1
       END DO
       CALL MEXPRINTF(' destroyed 2 \n')
       call mxDestroyArray(THPATCH_LOCAL_PR)
       i=0
       DO WHILE(i.lt.1e8)
         i=i+1
       END DO
       CALL MEXPRINTF(' destroyed 3 \n')
       call mxDestroyArray(LWGHT_LOCAL_PR)
       i=0
       DO WHILE(i.lt.1e8)
         i=i+1
       END DO

      CALL MEXPRINTF(' time to return ')

      RETURN
      END

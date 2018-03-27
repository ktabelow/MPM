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
       MWPOINTER LWGHT_IN_PR
       MWPOINTER WGHT_IN_PR
       MWPOINTER SWJY_IN_PR
       MWPOINTER HAKT_IN_PR, HAKT_OUT_PR
       MWPOINTER BI_IN_PR, BI_OUT_PR
       MWPOINTER THNEW_IN_PR, THNEW_OUT_PR
       MWPOINTER IHPATCH_IN_PR
       MWPOINTER INVCOVP_IN_PR
       MWPOINTER BIIPATCH_IN_PR
       MWPOINTER THPATCH_IN_PR

C-----------------------------------------------------------------------
C
C      LOCAL VARIABLES



       INTEGER*4 N1_IN, N2_IN, N3_IN, NP1_IN, NP2_IN, NP3_IN,
     +           NV_IN, NVD_IN, NMASK_IN, NCORES_IN

C       INTEGER*4 III_IN(NMASK_IN), JJJ_IN(NMASK_IN)

C       REAL*8 BI_IN(NMASK_IN)
C       REAL*8 Y_IN(NV_IN,NMASK_IN), THETA_IN(NV_IN,NMASK_IN),
C     +        INVCOV_IN(NVD_IN,NMASK_IN)

       REAL*8 LAMBDA_IN, HAKT_IN, WGHT_IN(2)

       REAL*8, DIMENSION(:), ALLOCATABLE :: BI_OUT
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: THNEW_OUT
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: LWGHT_LOCAL
       REAL*8, DIMENSION(:), ALLOCATABLE :: SWJY_LOCAL
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: THPATCH_LOCAL
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: INVCOVP_LOCAL
       REAL*8, DIMENSION(:), ALLOCATABLE :: BIIPATCH_LOCAL

       INTEGER*4 DWL, IH
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
       IF (NRHS .NE. 19) THEN
         CALL MEXERRMSGTXT('VPAWS requires 19 input arguments')
       ELSEIF (NLHS .GT. 2) THEN
         CALL MEXERRMSGTXT('VAWS requires 2 output argument')
       ENDIF

C      CHECK IF THE IMPUTS ARE NUMERIC

C       CALL MEXPRINTF(' BEFORE Checking numerics ')
       IF(mxIsNumeric(prhs(1)) .eq. 0) then
         CALL MEXERRMSGTXT('Y HAS TO BE A NUMERIC MATRIX')
       ENDIF

       IF(mxIsNumeric(prhs(2)) .eq. 0) then
         CALL MEXERRMSGTXT('III HAS TO BE A NUMERIC MATRIX')
       ENDIF

       IF(mxIsNumeric(prhs(3)) .eq. 0) then
         CALL MEXERRMSGTXT('JJJ HAS TO BE A NUMERIC MATRIX')
       ENDIF

       IF(mxIsNumeric(prhs(4)) .eq. 0) then
         CALL MEXERRMSGTXT('THETA HAS TO BE A NUMERIC MATRIX')
       ENDIF

       IF(mxIsNumeric(prhs(5)) .eq. 0) then
         CALL MEXERRMSGTXT('INVCOV HAS TO BE A NUMERIC MATRIX')
       ENDIF

       IF(mxIsNumeric(prhs(6)) .eq. 0) then
         CALL MEXERRMSGTXT('BI HAS TO BE A NUMERIC MATRIX')
       ENDIF

       IF(mxIsNumeric(prhs(7)) .eq. 0) then
         CALL MEXERRMSGTXT('WGHT HAS TO BE A NUMERIC MATRIX')
       ENDIF

       IF(mxIsNumeric(prhs(8)) .eq. 0) then
         CALL MEXERRMSGTXT('NV HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(9)) .eq. 0) then
         CALL MEXERRMSGTXT('NVD HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(10)) .eq. 0) then
         CALL MEXERRMSGTXT('N1 HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(11)) .eq. 0) then
         CALL MEXERRMSGTXT('N2 HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(12)) .eq. 0) then
         CALL MEXERRMSGTXT('N3 HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(13)) .eq. 0) then
         CALL MEXERRMSGTXT('NMASK HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(14)) .eq. 0) then
         CALL MEXERRMSGTXT('NP1 HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(15)) .eq. 0) then
         CALL MEXERRMSGTXT('NP2 HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(16)) .eq. 0) then
         CALL MEXERRMSGTXT('NP3 HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(17)) .eq. 0) then
         CALL MEXERRMSGTXT('HAKT HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(18)) .eq. 0) then
         CALL MEXERRMSGTXT('LAMBDA HAS TO BE A NUMBER')
       ENDIF

       IF(mxIsNumeric(prhs(19)) .eq. 0) then
         CALL MEXERRMSGTXT('NCORES HAS TO BE A NUMBER')
       ENDIF

        CALL MEXPRINTF(' after checking numerics \n')

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
       ALLOCATE(BI_OUT(NMASK_IN))
       ALLOCATE(THNEW_OUT(NV_IN,NMASK_IN))
       ALLOCATE(LWGHT_LOCAL(DWL,DWL))
       ALLOCATE(SWJY_LOCAL(NV_IN))
       ALLOCATE(THPATCH_LOCAL(NV_IN,NP1_IN*NP2_IN*NP3_IN))
       ALLOCATE(INVCOVP_LOCAL(NVD_IN,NP1_IN*NP2_IN*NP3_IN))
       ALLOCATE(BIIPATCH_LOCAL(NP1_IN*NP2_IN*NP3_IN))

C       HAKT_IN_PR = MXGETPR(PRHS(17))
C       LAMBDA_IN_PR = MXGETPR(PRHS(18))
       CALL MEXPRINTF(' AFTER getting OTHER input \n')

       CALL MEXPRINTF(' BEFORE CREATING OUTPUT ')
C CREATE MATRICES FOR RETURN ARGUMENT
C
C
C     BI
      PLHS(1) = MXCREATEDOUBLEMATRIX(1, NMASK_IN, 0)
      CALL MEXPRINTF(' BI done \n')

C     THNEW
      PLHS(2) = MXCREATEDOUBLEMATRIX(NV_IN, NMASK_IN, 0)
      CALL MEXPRINTF(' THNEW done \n')

      CALL MEXPRINTF(' AFTER PLACE OUTPUT ')
C
C ASSIGN POINTERS TO THE OUTPUT PARAMETERS
C
      BI_OUT_PR = MXGETPR(PLHS(1))
      THNEW_OUT_PR = MXGETPR(PLHS(2))

C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES


C       CALL MXCOPYPTRTOREAL8(BI_IN_PR, BI_INOUT,
C     +                       MXGETNUMBEROFELEMENTS(PRHS(14)))C
C       CALL MXCOPYPTRTOREAL8(THNEW_IN_PR, THNEW_INOUT,
C     +                       MXGETNUMBEROFELEMENTS(PRHS(15)))
C       CALL MXCOPYPTRTOREAL8(HAKT_IN_PR, HAKT_INOUT,
C     +                       MXGETNUMBEROFELEMENTS(PRHS(10)))
C       CALL MXCOPYPTRTOREAL8(Y_IN_PR, Y_IN,
C     +                       MXGETNUMBEROFELEMENTS(PRHS(1)))



C       CALL MXCOPYPTRTOINTEGER4(III_IN_PR, III_IN,
C     +                          MXGETNUMBEROFELEMENTS(PRHS(8)))

C       CALL MXCOPYPTRTOINTEGER4(JJJ_IN_PR, JJJ_IN,
C     +                          MXGETNUMBEROFELEMENTS(PRHS(9)))


      CALL MEXPRINTF('BEFORE CALL \n')
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
     +               BI_OUT,
     +               THNEW_OUT,
     +               NCORES_IN,
     +               LWGHT_LOCAL,
     +               %val(WGHT_IN_PR),
     +               SWJY_LOCAL,
     +               NP1_IN,
     +               NP2_IN,
     +               NP3_IN,
     +               THPATCH_LOCAL,
     +               INVCOVP_LOCAL,
     +               BIIPATCH_LOCAL)

       CALL MEXPRINTF(' call END \n')

       CALL MXCOPYREAL8TOPTR(BI_OUT, BI_OUT_PR,
     +                       MXGETNUMBEROFELEMENTS(PRHS(14)))
       CALL MXCOPYREAL8TOPTR(THNEW_OUT, THNEW_OUT_PR,
     +                       MXGETNUMBEROFELEMENTS(PRHS(15)))

       DEALLOCATE(BI_OUT)
       DEALLOCATE(THNEW_OUT)
       DEALLOCATE(LWGHT_LOCAL)
       DEALLOCATE(SWJY_LOCAL)
       DEALLOCATE(THPATCH_LOCAL)
       DEALLOCATE(INVCOVP_LOCAL)
       DEALLOCATE(BIIPATCH_LOCAL)

      CALL MEXPRINTF(' time to return ')

      RETURN
      END

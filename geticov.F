#include <fintrf.h>

C
C Gateway function geticov
C calls other function in file esticov.F
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
       MWPOINTER  MXGETPR
       MWPOINTER MXCREATEDOUBLEMATRIX

       real*8  mxGetScalar

C
C-----------------------------------------------------------------------
C
C
       MWPOINTER THETA_IN_PR
       MWPOINTER TE_IN_PR
       MWPOINTER IND_IN_PR
       MWPOINTER SIGMAI_IN_PR

       MWPOINTER ICOV_OUT_PR

C-----------------------------------------------------------------------
C
C      LOCAL VARIABLES

      integer*4 N_IN, ND_IN, NV_IN

      MWSIZE N, NV, NVD


C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
       IF (NRHS .NE. 7) THEN
         CALL MEXERRMSGTXT('VAWS requires 16 input arguments')
       ELSEIF (NLHS .GT. 1) THEN
         CALL MEXERRMSGTXT('VAWS requires 2 output argument')
       ENDIF

       CALL MEXPRINTF(' Start Iteration of VAWS ... ')

       NV_IN = INT(MXGETSCALAR(PRHS(5)))
       N_IN = INT(MXGETSCALAR(PRHS(6)))
       ND_IN = INT(MXGETSCALAR(PRHS(7)))

C       getting all the other INPUT pointers

       THETA_IN_PR = MXGETPR(PRHS(1))
       TE_IN_PR = MXGETPR(PRHS(2))
       IND_IN_PR = MXGETPR(PRHS(3))
       SIGMAI_IN_PR = MXGETPR(PRHS(4))

C CREATE MATRICES FOR RETURN ARGUMENT
       NVD=NV_IN*(NV_IN+1)/2
       N=N_IN
       NV=NV_IN

C     ICOV
       PLHS(1) = MXCREATEDOUBLEMATRIX(NVD, N, 0)

C
      ICOV_OUT_PR = MXGETPR(PLHS(1))
C
      IF(NV_IN.eq.3) THEN
C
C   ESTATICS model without MT
C
        CALL INVCOV3(%val(THETA_IN_PR),
     +            %val(TE_IN_PR),
     +            %val(IND_IN_PR),
     +            %val(SIGMAI_IN_PR),
     +            ND_IN,
     +            N_IN,
     +            %val(ICOV_OUT_PR))
      ELSE
C
C   ESTATICS model without MT
C
          CALL INVCOV4(%val(THETA_IN_PR),
     +            %val(TE_IN_PR),
     +            %val(IND_IN_PR),
     +            %val(SIGMAI_IN_PR),
     +            ND_IN,
     +            N_IN,
     +            %val(ICOV_OUT_PR))
      END IF
C


      CALL MEXPRINTF('inverse covariance matrix created \n')

      RETURN
      END

#include <fintrf.h>

C
C Gateway function for cvaws2.F 
C calls other function in file cvaws2.F
C
C 
C
C 
C 
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
       
       MWPOINTER Y_IN_PR, THETA_IN_PR
       MWPOINTER III_IN_PR, JJJ_IN_PR

       MWPOINTER SI2_IN_PR, LWGHT_IN_PR
       MWPOINTER WGHT_IN_PR, SWIJ_IN_PR
       MWPOINTER HAKT_IN_PR, HAKT_OUT_PR
       MWPOINTER BI_IN_PR, BI_OUT_PR
       MWPOINTER THNEW_IN_PR, THNEW_OUT_PR

C      POINTERS TO M (rows) AND N (columns) of INPUTS

       MWPOINTER M_N1, N_N1, M_N2, N_N2, M_N3, N_N3,
     +  M_NMASK, N_NMASK,
     + M_HAKT, N_HAKT, M_NV, N_NV, M_NCORES, N_NCORES,
     + M_LAMBDA, N_LAMBDA,
     + M_WGHT, N_WGHT

       
C-----------------------------------------------------------------------
C
C      LOCAL VARIABLES



       INTEGER*4 N1_IN, N2_IN, N3_IN, NV_IN, NCORES_IN, NMASK_IN
       
        INTEGER*4, DIMENSION(:), ALLOCATABLE   :: III_IN, JJJ_IN
        
        REAL*8, DIMENSION(:), ALLOCATABLE   :: BI_INOUT
       
        REAL*8, DIMENSION(:,:), ALLOCATABLE   :: THNEW_INOUT, Y_IN
       
       REAL*8 LAMBDA_IN, HAKT_INOUT, HAKT_IN, WGHT_IN(2)
       
C      SIZE OF THE INPUT ARGUMENTS (MXN)      

       MWSIZE DIM_N1, DIM_N2, DIM_N3, DIM_NV, DIM_NCORES, 
     + NVSIZE,NSIZE, ONE, DIM_LAMBDA, DIM_HAKT, DIM_WGHT,
     + DIM_NMASK, NMASKSIZE
       

       
C       MWSIZE PRODN, NARRAY(3) 

C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
       IF (NRHS .NE. 18) THEN
         CALL MEXERRMSGTXT('VAWS requires 18 input arguments')
       ELSEIF (NLHS .GT. 3) THEN
         CALL MEXERRMSGTXT('VAWS requires 3 output argument')
       ENDIF
      
C      CHECK IF THE IMPUTS ARE NUMERIC

C       CALL MEXPRINTF(' BEFORE Checking numerics ')
       IF(mxIsNumeric(prhs(1)) .eq. 0) then 
           CALL MEXERRMSGTXT('S HAS TO BE A NUMERIC MATRIX')
        ENDIF

        IF(mxIsNumeric(prhs(2)) .eq. 0) then 
           CALL MEXERRMSGTXT('NV HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(3)) .eq. 0) then 
           CALL MEXERRMSGTXT('N1 HAS TO BE A NUMBER')
        ENDIF
       
        IF(mxIsNumeric(prhs(4)) .eq. 0) then 
           CALL MEXERRMSGTXT('N2 HAS TO BE A NUMBER')
        ENDIF
  
        IF(mxIsNumeric(prhs(5)) .eq. 0) then 
           CALL MEXERRMSGTXT('N3 HAS TO BE A NUMBER')
        ENDIF

         IF(mxIsNumeric(prhs(6)) .eq. 0) then 
           CALL MEXERRMSGTXT('NMASK HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(7)) .eq. 0) then 
           CALL MEXERRMSGTXT('III HAS TO BE A NUMERIC ARRAY')
        ENDIF

        IF(mxIsNumeric(prhs(8)) .eq. 0) then 
           CALL MEXERRMSGTXT('JJJ HAS TO BE A NUMERIC ARRAY')
        ENDIF

        IF(mxIsNumeric(prhs(9)) .eq. 0) then 
           CALL MEXERRMSGTXT('HAKT HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(10)) .eq. 0) then 
           CALL MEXERRMSGTXT('LAMBDA HAS TO BE A NUMBER')
        ENDIF
       
        IF(mxIsNumeric(prhs(11)) .eq. 0) then 
           CALL MEXERRMSGTXT('THETA HAS TO BE A NUMERIC MATRIX')
        ENDIF
  
        IF(mxIsNumeric(prhs(12)) .eq. 0) then 
           CALL MEXERRMSGTXT('SI2 HAS TO BE A NUMERIC MATRIX')
        ENDIF

        IF(mxIsNumeric(prhs(13)) .eq. 0) then 
           CALL MEXERRMSGTXT('BI HAS TO BE A NUMERIC ARRAY')
        ENDIF
       
        IF(mxIsNumeric(prhs(14)) .eq. 0) then 
           CALL MEXERRMSGTXT('THNEW HAS TO BE A NUMERIC MATRIX')
        ENDIF
  
        IF(mxIsNumeric(prhs(15)) .eq. 0) then 
           CALL MEXERRMSGTXT('NCORES HAS TO BE A NUMBER')
        ENDIF

         IF(mxIsNumeric(prhs(16)) .eq. 0) then 
           CALL MEXERRMSGTXT('LWGHT HAS TO BE A NUMERIC ARRAY')
        ENDIF

         IF(mxIsNumeric(prhs(17)) .eq. 0) then 
           CALL MEXERRMSGTXT('WGHT HAS TO BE A NUMERIC ARRAY')
        ENDIF

         IF(mxIsNumeric(prhs(18)) .eq. 0) then 
           CALL MEXERRMSGTXT('SWIJ HAS TO BE A NUMERIC MATRIX')
        ENDIF

        
C        CALL MEXPRINTF(' after checking numerics \n')
C
C CHECK THE DIMENSIONS OF INPUT VALUES.  
C N1,N2,N3,NV,NMASK NCORES, LAMBDA, HAKT HAVE TO BE 1
C WGHT HAS TO BE 2
C 
C
       M_N1 = MXGETM(PRHS(3))
       N_N1 = MXGETN(PRHS(3))
       M_N2 = MXGETM(PRHS(4))
       N_N2 = MXGETN(PRHS(4))
       M_N3 = MXGETM(PRHS(5))
       N_N3 = MXGETN(PRHS(5))
       
       M_NMASK = MXGETM(PRHS(6))
       N_NMASK = MXGETN(PRHS(6))
       

       M_WGHT = MXGETM(PRHS(17))
       N_WGHT = MXGETN(PRHS(17))
    
       M_NV = MXGETM(PRHS(2))
       N_NV = MXGETN(PRHS(2))
       M_HAKT = MXGETM(PRHS(9))
       N_HAKT = MXGETN(PRHS(9))
       M_LAMBDA = MXGETM(PRHS(10))
       N_LAMBDA = MXGETN(PRHS(10))
       M_NCORES = MXGETM(PRHS(15))
       N_NCORES = MXGETN(PRHS(15))
      
       

      
C       NARRAY = MXGETDIMENSIONS(PRHS(1))
C       PRODN = MXGETNUMBEROFELEMENTS(PRHS(3))
       
       
       DIM_N1 = M_N1*N_N1
       DIM_N2 = M_N2*N_N2
       DIM_N3 = M_N3*N_N3       
     
       DIM_NMASK = M_NMASK*N_NMASK     

       DIM_WGHT = M_WGHT*N_WGHT
       DIM_NV = M_NV*N_NV
       DIM_HAKT = M_HAKT*N_HAKT
       DIM_LAMBDA = M_LAMBDA*N_LAMBDA
       DIM_NCORES = M_NCORES*N_NCORES
      
      
C       CALL MEXPRINTF(' checking dimensions \n')
C
      
      IF ((DIM_N1 .NE. 1).OR.(DIM_N2 .NE. 1).OR.(DIM_N3 .NE. 1)) THEN
        CALL MEXERRMSGTXT('N1, N2, N3 have to be scalar')
      ENDIF

      IF ((DIM_NV .NE. 1).OR.(DIM_NMASK .NE. 1)) THEN
        CALL MEXERRMSGTXT('NV and NMASK  have to be a SCALAR')
      ENDIF
      IF (DIM_WGHT .NE. 2) THEN
        CALL MEXERRMSGTXT('WGHT has to be a vector of lenght 2')
      ENDIF
      
      IF (DIM_NCORES .NE. 1) THEN
        CALL MEXERRMSGTXT('NCORES HAS to be scalar')
      ENDIF
    
      IF ((DIM_HAKT .NE. 1).OR.(DIM_LAMBDA .NE. 1)) THEN
        CALL MEXERRMSGTXT('HAKT and LAMBDA have to be scalar')
      ENDIF

     

      
C       getting the values of n1, n2, n3, nV, NCORES, NMASK 
C       HAKT, lambda
C       using get scalar
C       need to cast since it returns a double
C        CALL MEXPRINTF(' BEFORE getting SCALARS \n')


        N1_IN = INT(MXGETSCALAR(PRHS(3)))
       N2_IN = INT(MXGETSCALAR(PRHS(4)))
       N3_IN = INT(MXGETSCALAR(PRHS(5)))
       
       NMASK_IN = INT(MXGETSCALAR(PRHS(6)))

        HAKT_IN = MXGETSCALAR(PRHS(9))
         NV_IN = INT(MXGETSCALAR(PRHS(2)))
       NCORES_IN = INT(MXGETSCALAR(PRHS(15)))
       LAMBDA_IN = MXGETSCALAR(PRHS(10))
        
       
      
    
C       CALL MEXPRINTF(' AFTER getting SCALARS \n')

       ALLOCATE(BI_INOUT(NMASK_IN))
       ALLOCATE(THNEW_INOUT(NV_IN,NMASK_IN)) 
C       ALLOCATE(Y_IN(NV_IN,NMASK_IN))
       ALLOCATE(III_IN(NMASK_IN))
       ALLOCATE(JJJ_IN(N1_IN*N2_IN*N3_IN))


       ONE=1
       NVSIZE=NV_IN
       NSIZE = N1_IN*N2_IN*N3_IN
       NMASKSIZE = NMASK_IN

C       getting all the other INPUT pointers
       

       Y_IN_PR = MXGETPR(PRHS(1))
       III_IN_PR = MXGETPR(PRHS(7))
        JJJ_IN_PR = MXGETPR(PRHS(8))
       
       
       THETA_IN_PR = MXGETPR(PRHS(11))
       SI2_IN_PR = MXGETPR(PRHS(12))
       
       BI_IN_PR = MXGETPR(PRHS(13))
       THNEW_IN_PR = MXGETPR(PRHS(14))
       HAKT_IN_PR = MXGETPR(PRHS(9))
       
       LWGHT_IN_PR = MXGETPR(PRHS(16))
       WGHT_IN_PR= MXGETPR(PRHS(17))
       SWIJ_IN_PR= MXGETPR(PRHS(18))
      
     
C       CALL MEXPRINTF(' AFTER getting OTHER input \n')

C       CALL MEXPRINTF(' BEFORE CREATING OUTPUT ')
C CREATE MATRICES FOR RETURN ARGUMENT
C
C    
C     BI
      PLHS(1) = MXCREATEDOUBLEMATRIX(ONE,
     + MXGETNUMBEROFELEMENTS(PRHS(13)),0)
C      CALL MEXPRINTF(' BI done \n')

C     THNEW 
      PLHS(2) = MXCREATEDOUBLEMATRIX(NVSIZE,
     + NMASKSIZE,0)
C      CALL MEXPRINTF(' THNEW done \n') 

C    HAKT 
      PLHS(3) = MXCREATEDOUBLEMATRIX(ONE,
     + ONE,0)
C        CALL MEXPRINTF(' HAKT done \n')


C      CALL MEXPRINTF(' AFTER PLACE OUTPUT ')
C
C ASSIGN POINTERS TO THE OUTPUT PARAMETERS
C
      BI_OUT_PR = MXGETPR(PLHS(1))
      THNEW_OUT_PR = MXGETPR(PLHS(2))
      HAKT_OUT_PR = MXGETPR(PLHS(3))
     
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
     


       CALL MXCOPYPTRTOREAL8(BI_IN_PR, BI_INOUT, 
     + MXGETNUMBEROFELEMENTS(PRHS(13)))
       CALL MXCOPYPTRTOREAL8(THNEW_IN_PR, THNEW_INOUT, 
     + MXGETNUMBEROFELEMENTS(PRHS(14)))
       CALL MXCOPYPTRTOREAL8(HAKT_IN_PR, HAKT_INOUT, 
     + MXGETNUMBEROFELEMENTS(PRHS(9)))
C       CALL MXCOPYPTRTOREAL8(Y_IN_PR, Y_IN, 
C     + MXGETNUMBEROFELEMENTS(PRHS(1)))



       CALL MXCOPYPTRTOINTEGER4(III_IN_PR, III_IN, 
     + MXGETNUMBEROFELEMENTS(PRHS(7)))

       CALL MXCOPYPTRTOINTEGER4(JJJ_IN_PR, JJJ_IN, 
     + MXGETNUMBEROFELEMENTS(PRHS(8)))
       
       
C      CALL MEXPRINTF('BEFORE CALL \n')
C
C DO THE ACTUAL COMPUTATIONS ALSO USING VAL
C



        CALL CVAWS2(%val(Y_IN_PR),
     + NV_IN,N1_IN,N2_IN,N3_IN,NMASK_IN, 
     + III_IN, JJJ_IN, HAKT_INOUT,     
     + LAMBDA_IN,%val(THETA_IN_PR),
     + %val(SI2_IN_PR), BI_INOUT, THNEW_INOUT,   
     + NCORES_IN,%val(LWGHT_IN_PR),
     + %val(WGHT_IN_PR),%val(SWIJ_IN_PR))

C         CALL MEXPRINTF(' call END \n')
    
       CALL MXCOPYREAL8TOPTR(BI_INOUT, BI_OUT_PR, 
     + MXGETNUMBEROFELEMENTS(PRHS(13)))
       CALL MXCOPYREAL8TOPTR(THNEW_INOUT, THNEW_OUT_PR, 
     + MXGETNUMBEROFELEMENTS(PRHS(14)))
       CALL MXCOPYREAL8TOPTR(HAKT_INOUT, HAKT_OUT_PR, ONE)

      DEALLOCATE(BI_INOUT)     
      DEALLOCATE(THNEW_INOUT)
C      DEALLOCATE(Y_IN)
      DEALLOCATE(III_IN)
      DEALLOCATE(JJJ_IN) 
C        CALL MEXPRINTF(' time to return ')

      RETURN
      END





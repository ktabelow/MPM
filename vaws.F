#include <fintrf.h>

C
C Gateway function for vaws.F 
C awslchi2 calls other function in file cvaws.F
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
       
       MWPOINTER Y_IN_PR, MASK_IN_PR
       MWPOINTER THETA_IN_PR

       MWPOINTER SI2_IN_PR, LWGHT_IN_PR
       MWPOINTER WGHT_IN_PR, SWIJ_IN_PR
       MWPOINTER HAKT_IN_PR, HAKT_OUT_PR
       MWPOINTER BI_IN_PR, BI_OUT_PR
       MWPOINTER THNEW_IN_PR, THNEW_OUT_PR

C      POINTERS TO M (rows) AND N (columns) of INPUTS

       MWPOINTER M_N1, N_N1, M_N2, N_N2, M_N3, N_N3,
     + M_HAKT, N_HAKT, M_NV, N_NV, M_NCORES, N_NCORES,
     + M_LAMBDA, N_LAMBDA,
     + M_WGHT, N_WGHT

       
C-----------------------------------------------------------------------
C
C      LOCAL VARIABLES



       INTEGER*4 N1_IN, N2_IN, N3_IN, NV_IN, NCORES_IN
       
        INTEGER*8, DIMENSION(:), ALLOCATABLE   :: MASK_IN
        
        REAL*8, DIMENSION(:), ALLOCATABLE   :: BI_INOUT
       
        REAL*8, DIMENSION(:,:), ALLOCATABLE   :: THNEW_INOUT, Y_IN
       
       REAL*8 LAMBDA_IN, HAKT_INOUT, HAKT_IN, WGHT_IN(6)
       
C      SIZE OF THE INPUT ARGUMENTS (MXN)      

       MWSIZE DIM_N1, DIM_N2, DIM_N3, DIM_NV, DIM_NCORES, 
     + NVSIZE,NSIZE, ONE, DIM_LAMBDA, DIM_HAKT, DIM_WGHT
       

       
C       MWSIZE PRODN, NARRAY(3) 

C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
       IF (NRHS .NE. 16) THEN
         CALL MEXERRMSGTXT('VAWS requires 16 input arguments')
       ELSEIF (NLHS .GT. 3) THEN
         CALL MEXERRMSGTXT('VAWS requires 3 output argument')
       ENDIF
      
C      CHECK IF THE IMPUTS ARE NUMERIC

C       CALL MEXPRINTF(' BEFORE Checking numerics ')
       IF(mxIsNumeric(prhs(1)) .eq. 0) then 
           CALL MEXERRMSGTXT('S HAS TO BE A NUMERIC MATRIX')
        ENDIF

        IF(mxIsNumeric(prhs(2)) .eq. 0) then 
           CALL MEXERRMSGTXT('MASK HAS TO BE A MATRIX')
        ENDIF

        IF(mxIsNumeric(prhs(3)) .eq. 0) then 
           CALL MEXERRMSGTXT('NV HAS TO BE A NUMBER')
        ENDIF
       
        IF(mxIsNumeric(prhs(4)) .eq. 0) then 
           CALL MEXERRMSGTXT('N1 HAS TO BE A NUMBER')
        ENDIF
  
        IF(mxIsNumeric(prhs(5)) .eq. 0) then 
           CALL MEXERRMSGTXT('N2 HAS TO BE A NUMBER')
        ENDIF

         IF(mxIsNumeric(prhs(6)) .eq. 0) then 
           CALL MEXERRMSGTXT('N3 HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(7)) .eq. 0) then 
           CALL MEXERRMSGTXT('HAKT HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(8)) .eq. 0) then 
           CALL MEXERRMSGTXT('LAMBDA HAS TO BE A NUMBER')
        ENDIF

        IF(mxIsNumeric(prhs(9)) .eq. 0) then 
           CALL MEXERRMSGTXT('THETA HAS TO BE A NUMERIC MATRIX')
        ENDIF

        IF(mxIsNumeric(prhs(10)) .eq. 0) then 
           CALL MEXERRMSGTXT('SI2 HAS TO BE A 4X4XN MATRIX')
        ENDIF
       
        IF(mxIsNumeric(prhs(11)) .eq. 0) then 
           CALL MEXERRMSGTXT('BI HAS TO BE A NUMERIC VECTOR')
        ENDIF
  
        IF(mxIsNumeric(prhs(12)) .eq. 0) then 
           CALL MEXERRMSGTXT('THNEW HAS TO BE A NUMERIC MATRIX')
        ENDIF

        IF(mxIsNumeric(prhs(13)) .eq. 0) then 
           CALL MEXERRMSGTXT('NCORES HAS TO BE A NUMBER')
        ENDIF
       
        IF(mxIsNumeric(prhs(14)) .eq. 0) then 
           CALL MEXERRMSGTXT('LWGHT HAS TO BE A NUMERIC MATRIX')
        ENDIF
  
        IF(mxIsNumeric(prhs(15)) .eq. 0) then 
           CALL MEXERRMSGTXT('WGHT HAS TO BE A NUMERIC VECTOR')
        ENDIF

         IF(mxIsNumeric(prhs(16)) .eq. 0) then 
           CALL MEXERRMSGTXT('SWIJ HAS TO BE A NUMERIC MATRIX')
        ENDIF

        
C        CALL MEXPRINTF(' after checking numerics \n')
C
C CHECK THE DIMENSIONS OF INPUT VALUES.  
C N1,N2,N3,NV, NCORES, LAMBDA, HAKT HAVE TO BE 1
C WGHT HAS TO BE 2
C 
C
       M_N1 = MXGETM(PRHS(4))
       N_N1 = MXGETN(PRHS(4))
       M_N2 = MXGETM(PRHS(5))
       N_N2 = MXGETN(PRHS(5))
       M_N3 = MXGETM(PRHS(6))
       N_N3 = MXGETN(PRHS(6))
       
       

       M_WGHT = MXGETM(PRHS(15))
       N_WGHT = MXGETN(PRHS(15))
    
       M_NV = MXGETM(PRHS(3))
       N_NV = MXGETN(PRHS(3))
       M_HAKT = MXGETM(PRHS(7))
       N_HAKT = MXGETN(PRHS(7))
       M_LAMBDA = MXGETM(PRHS(8))
       N_LAMBDA = MXGETN(PRHS(8))
       M_NCORES = MXGETM(PRHS(13))
       N_NCORES = MXGETN(PRHS(13))
      
       

      
C       NARRAY = MXGETDIMENSIONS(PRHS(1))
C       PRODN = MXGETNUMBEROFELEMENTS(PRHS(3))
       
       
       DIM_N1 = M_N1*N_N1
       DIM_N2 = M_N2*N_N2
       DIM_N3 = M_N3*N_N3       
     

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

      IF (DIM_NV .NE. 1) THEN
        CALL MEXERRMSGTXT('NV  has to be a SCALAR')
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

     

      
C       getting the values of n1, n2, n3, nV, NCORES 
C       HAKT, lambda
C       using get scalar
C       need to cast since it returns a double
C        CALL MEXPRINTF(' BEFORE getting SCALARS \n')


        N1_IN = INT(MXGETSCALAR(PRHS(4)))
       N2_IN = INT(MXGETSCALAR(PRHS(5)))
       N3_IN = INT(MXGETSCALAR(PRHS(6)))
       
        HAKT_IN = MXGETSCALAR(PRHS(7))
         NV_IN = INT(MXGETSCALAR(PRHS(3)))
       NCORES_IN = INT(MXGETSCALAR(PRHS(13)))
       LAMBDA_IN = MXGETSCALAR(PRHS(8))
        
       
      
    
C       CALL MEXPRINTF(' AFTER getting SCALARS \n')

       ALLOCATE(BI_INOUT(N1_IN*N2_IN*N3_IN))
       ALLOCATE(THNEW_INOUT(NV_IN,N1_IN*N2_IN*N3_IN)) 
       ALLOCATE(Y_IN(NV_IN,N1_IN*N2_IN*N3_IN))
       ALLOCATE(MASK_IN(N1_IN*N2_IN*N3_IN))

C       THREE=mxgetnumberofdimensions(PRHS(23))
       ONE=1
       NVSIZE=NV_IN
       NSIZE = N1_IN*N2_IN*N3_IN

C       getting all the other INPUT pointers
       

       Y_IN_PR = MXGETPR(PRHS(1))
       
       MASK_IN_PR = MXGETPR(PRHS(2))
       
       THETA_IN_PR = MXGETPR(PRHS(9))
       SI2_IN_PR = MXGETPR(PRHS(10))
       
       BI_IN_PR = MXGETPR(PRHS(11))
       THNEW_IN_PR = MXGETPR(PRHS(12))
       HAKT_IN_PR = MXGETPR(PRHS(7))
       
       LWGHT_IN_PR = MXGETPR(PRHS(14))
       WGHT_IN_PR= MXGETPR(PRHS(15))
       SWIJ_IN_PR= MXGETPR(PRHS(16))
      
     
C       CALL MEXPRINTF(' AFTER getting OTHER input \n')

C       CALL MEXPRINTF(' BEFORE CREATING OUTPUT ')
C CREATE MATRICES FOR RETURN ARGUMENT
C
C    
C     BI
      PLHS(1) = MXCREATEDOUBLEMATRIX(ONE,
     + MXGETNUMBEROFELEMENTS(PRHS(11)),0)
C      CALL MEXPRINTF(' BI done \n')

C     THNEW 
      PLHS(2) = MXCREATEDOUBLEMATRIX(NVSIZE,
     + NSIZE,0)
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
     + MXGETNUMBEROFELEMENTS(PRHS(11)))
       CALL MXCOPYPTRTOREAL8(THNEW_IN_PR, THNEW_INOUT, 
     + MXGETNUMBEROFELEMENTS(PRHS(12)))
       CALL MXCOPYPTRTOREAL8(HAKT_IN_PR, HAKT_INOUT, 
     + MXGETNUMBEROFELEMENTS(PRHS(7)))
       CALL MXCOPYPTRTOREAL8(Y_IN_PR, Y_IN, 
     + MXGETNUMBEROFELEMENTS(PRHS(1)))
       CALL MXCOPYPTRTOINTEGER8(MASK_IN_PR, MASK_IN, 
     + MXGETNUMBEROFELEMENTS(PRHS(2)))
       
       
C      CALL MEXPRINTF('BEFORE CALL \n')
C
C DO THE ACTUAL COMPUTATIONS ALSO USING VAL
C

C        CALL CVAWS2(%val(Y_IN_PR),MASK_IN,
C     + NV_IN,N1_IN,N2_IN,N3_IN,HAKT_INOUT,     
C     + LAMBDA_IN,%val(THETA_IN_PR),
C     + %val(SI2_IN_PR), %val(BI_IN_PR), %val(THNEW_IN_PR),   
C     + NCORES_IN,%val(LWGHT_IN_PR),
C     + %val(WGHT_IN_PR),%val(SWIJ_IN_PR))

C at the moment y and mask can be called both by value 
C %val(Y_IN_PR)
C or by reference Y_IN
        CALL CVAWS(%val(Y_IN_PR),MASK_IN,
     + NV_IN,N1_IN,N2_IN,N3_IN,HAKT_INOUT,     
     + LAMBDA_IN,%val(THETA_IN_PR),
     + %val(SI2_IN_PR), BI_INOUT, THNEW_INOUT,   
     + NCORES_IN,%val(LWGHT_IN_PR),
     + %val(WGHT_IN_PR),%val(SWIJ_IN_PR))

C         CALL MEXPRINTF(' call END \n')
    
       CALL MXCOPYREAL8TOPTR(BI_INOUT, BI_OUT_PR, 
     + MXGETNUMBEROFELEMENTS(PRHS(11)))
       CALL MXCOPYREAL8TOPTR(THNEW_INOUT, THNEW_OUT_PR, 
     + MXGETNUMBEROFELEMENTS(PRHS(12)))
       CALL MXCOPYREAL8TOPTR(HAKT_INOUT, HAKT_OUT_PR, ONE)

      DEALLOCATE(BI_INOUT)     
      DEALLOCATE(THNEW_INOUT)
      DEALLOCATE(Y_IN)
      DEALLOCATE(MASK_IN)
C        CALL MEXPRINTF(' time to return ')

      RETURN
      END





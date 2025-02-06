MODULE TRACER_MOD
  IMPLICIT NONE

  PUBLIC

  REAL*8, ALLOCATABLE            :: TCVV      (:)
  REAL*8, ALLOCATABLE            :: HR_STT_ADJ(:, :, :, :)
 
  CONTAINS

  SUBROUTINE INIT_TRACER()

    USE CONTROL_VAR_MOD,      ONLY : SPE_MASS, RUN_DIR, SPE_NUMS
    USE FORCING_MOD,          ONLY : INIT_FIELD

    include  "CMN_SIZE"

    !local var
    INTEGER              :: as, I

    allocate(HR_STT_ADJ(IIPAR, JJPAR, LLPAR, SPE_NUMS), stat=as)
    HR_STT_ADJ = 0d0

    allocate(TCVV(SPE_NUMS), stat=as)
    TCVV = 0d0

    DO I = 1, SPE_NUMS
      TCVV(I) = 28.97 / SPE_MASS(I)
    ENDDO
    IF ( RUN_DIR ) THEN
      HR_STT_ADJ(:, :, 1, 1) = INIT_FIELD
    ENDIF
    print *, 'NOW is check HR', sum(HR_STT_ADJ)

  END SUBROUTINE

  SUBROUTINE SAVE_TRACER()
    USE CONTROL_VAR_MOD,      ONLY : OUT_PATH, SPE_NAME
    USE TIME_MOD,             ONLY : TIME_OBJ
    USE NETCDF_IO_MOD,        ONLY : SAVE_NC
    USE PBL_MIX_MOD,          ONLY : GET_PBL_TOP_L
    USE PRESSURE_MOD,         ONLY : GET_PEDGE


    INCLUDE  "CMN_SIZE"

    !local var
    CHARACTER(LEN=255)      :: FILE_PATH
    REAL*8                  :: TMP_DATA(IIPAR, JJPAR)
    INTEGER                 :: I, J, L, TOP
    REAL*8                  :: TOTPRES, DELTPRES
    
    write(FILE_PATH, "('STT_ADJ.', I4.4, I2.2, I2.2, I2.2, I2.2,'.nc')")   &
    TIME_OBJ%YEAR, TIME_OBJ%MONTH, TIME_OBJ%DAY, TIME_OBJ%HOUR, TIME_OBJ%MINUTE

    FILE_PATH = TRIM(OUT_PATH) // TRIM(FILE_PATH)
    TMP_DATA = 0d0

    print *, "Now is right"
    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J, L, TOP, TOTPRES, DELTPRES )
    DO J = 1, JJPAR
      DO I = 1, IIPAR
        !TOP = FLOOR( GET_PBL_TOP_L( I, J ) )
        !TOP = MAX(TOP, 1)
        !TOTPRES = GET_PEDGE(I,J,1) - GET_PEDGE(I,J,TOP+1)
        !DO L = 1, TOP
          !DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
          !TMP_DATA(I, J) = TMP_DATA(I, J) +   &
          !HR_STT_ADJ(I, J, L, 1) * ( DELTPRES / TOTPRES )
        !ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    print *, "over"
    !CALL SAVE_NC(FILE_PATH, TRIM(SPE_NAME), HR_STT_ADJ(:, :, :, 1), IIPAR, JJPAR, LLPAR)
    !print *, SUM(TMP_DATA)
    CALL SAVE_NC(FILE_PATH, TRIM(SPE_NAME(1)), HR_STT_ADJ(:, :, 1, 1), IIPAR, JJPAR)

  END SUBROUTINE

END MODULE

module control_var_mod

  USE READER_MOD

  IMPLICIT NONE
  PUBLIC

  !module var
  logical            :: RUN_DIR, DIS_PLAN, HR_CHECK_MASS, LR_CHECK_MASS
  logical            :: Flux_C, FILL_NEG
  integer            :: y_s, m_s, d_s, h_s, min_s
  integer            :: y_e, m_e, d_e, h_e, min_e
  integer            :: IORD, JORD, KORD, SPE_NUMS

  integer            :: HR_STEP, LR_STEP, CON_STEP
  character(len=255) :: HR_MET_PATH, LR_MET_PATH, OUT_PATH
  character(len=255) :: HR_EMS_PATH, LR_EMS_PATH
  character(len=255) :: HR_RESTART_PATH, LR_RESTART_PATH
  character(len=255) :: HR_MASS_PATH, LR_MASS_PATH
  character(len=10), allocatable  :: SPE_NAME(:)
  integer, allocatable            :: SPE_MASS(:)


  CONTAINS

  SUBROUTINE INIT_CONTROL_VAR

    USE ERROR_MOD,      ONLY : STOP_ERROR

    !LOCAL VAR
    CHARACTER(LEN=100)       :: PAR_FILE_PATH
    INTEGER                  :: IOS, IU_GHSM
    CHARACTER(LEN=255)        :: LINE

    PAR_FILE_PATH = '/work/ese-zhengzy/model_z/GHSM/Parameter'
    OPEN( IU_GHSM, FILE=TRIM( PAR_FILE_PATH ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL STOP_ERROR( 'INIT_CONTROL_VAR-OPEN' )

    DO

      LINE = READ_ONE_LINE( IU_GHSM, 'INIT_CONTROL' )
      print *, TRIM(LINE)
      IF      ( INDEX( LINE, 'TIME SET'       ) > 0 ) THEN
            CALL READ_TIME(IU_GHSM)
      ELSE IF ( INDEX( LINE, 'PATH'           ) > 0 ) THEN
            CALL READ_PATH(IU_GHSM)
      ELSE IF ( INDEX( LINE, 'Corrective plan') > 0 ) THEN
            CALL READ_COR_PLAN(IU_GHSM)
      ELSE IF ( INDEX( LINE, 'Check_MASS'     ) > 0 ) THEN
            CALL READ_CHECK_MASS(IU_GHSM)
      ELSE IF ( INDEX( LINE, 'TRANSPORT MENU' ) > 0 ) THEN
            CALL READ_TRANSPORT(IU_GHSM)
      ELSE IF ( INDEX( LINE, 'SPECIE MENU'    ) > 0 ) THEN
            CALL READ_SPECIE(IU_GHSM)
      ELSE IF ( INDEX( LINE, 'END OF FILE'    ) > 0 ) THEN
            EXIT
      END IF

    END DO

    CLOSE( IU_GHSM )
    
  END SUBROUTINE

  SUBROUTINE READ_TIME( FILE_ID )
    INTEGER,      INTENT(IN)  :: FILE_ID

    !LOCAL VAR
    INTEGER                   :: N
    CHARACTER(LEN=255)        :: SUBSTRS(10)

    !Running direction, T represents forward direction,
    !F represents reverse direction

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_TIME:1' )
    READ( SUBSTRS(1:N), * ) RUN_DIR 

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 5, 'READ_TIME:2' )
    READ( SUBSTRS(1:N), * ) y_s, m_s, d_s, h_s, min_s

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 5, 'READ_TIME:3' )
    READ( SUBSTRS(1:N), * ) y_e, m_e, d_e, h_e, min_e

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_TIME:4' )
    READ( SUBSTRS(1:N), * ) HR_STEP

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_TIME:5' )
    READ( SUBSTRS(1:N), * ) LR_STEP

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_TIME:6' )
    READ( SUBSTRS(1:N), * ) CON_STEP
  END SUBROUTINE

  SUBROUTINE READ_PATH( FILE_ID )
    INTEGER,      INTENT(IN)  :: FILE_ID

    !LOCAL VAR
    INTEGER                   :: N
    CHARACTER(LEN=255)        :: SUBSTRS(10)

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_PATH:1' )
    READ( SUBSTRS(1:N), '(a)' ) HR_MET_PATH

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_PATH:2' )
    READ( SUBSTRS(1:N), '(a)' ) LR_MET_PATH

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_PATH:3' )
    READ( SUBSTRS(1:N), '(a)' ) OUT_PATH

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_PATH:4' )
    READ( SUBSTRS(1:N), '(a)' ) HR_EMS_PATH

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_PATH:5' )
    READ( SUBSTRS(1:N), '(a)' ) LR_EMS_PATH

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_PATH:5' )
    READ( SUBSTRS(1:N), '(a)' ) HR_RESTART_PATH

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_PATH:5' )
    READ( SUBSTRS(1:N), '(a)' ) LR_RESTART_PATH

  END SUBROUTINE

  SUBROUTINE READ_COR_PLAN( FILE_ID )
    INTEGER,      INTENT(IN)  :: FILE_ID

    !LOCAL VAR
    INTEGER                   :: N
    CHARACTER(LEN=255)        :: SUBSTRS(10)

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_COR_PLAN:1' )
    READ( SUBSTRS(1:N), * ) DIS_PLAN
    !there has an check process

  END SUBROUTINE

  SUBROUTINE READ_CHECK_MASS( FILE_ID )
    INTEGER,      INTENT(IN)  :: FILE_ID

    !LOCAL VAR
    INTEGER                   :: N
    CHARACTER(LEN=255)        :: SUBSTRS(10)

    !Running direction, T represents forward direction,
    !F represents reverse direction
    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_CHECK_MASS:1' )
    READ( SUBSTRS(1:N), * ) HR_CHECK_MASS

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_CHECK_MASS:2' )
    READ( SUBSTRS(1:N), * ) LR_CHECK_MASS

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_CHECK_MASS:3' )
    READ( SUBSTRS(1:N), '(a)' ) HR_MASS_PATH

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_CHECK_MASS:4' )
    READ( SUBSTRS(1:N), '(a)' ) LR_MASS_PATH

  END SUBROUTINE

  SUBROUTINE READ_TRANSPORT( FILE_ID )
    INTEGER,      INTENT(IN)  :: FILE_ID

    !LOCAL VAR
    INTEGER                   :: N
    CHARACTER(LEN=255)        :: SUBSTRS(10)

    !Running direction, T represents forward direction,
    !F represents reverse direction
    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_TRANSPORT:1' )
    READ( SUBSTRS(1:N), * ) Flux_C 

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_TRANSPORT:2' )
    READ( SUBSTRS(1:N), * ) FILL_NEG

    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 3, 'READ_TRANSPORT:3' )
    READ( SUBSTRS(1:N), * ) IORD, JORD, KORD

  END SUBROUTINE

  SUBROUTINE READ_SPECIE( FILE_ID )
    INTEGER,      INTENT(IN)  :: FILE_ID

    !LOCAL VAR
    INTEGER                   :: N, as, I, SPE_COUNT
    CHARACTER(LEN=255)        :: SUBSTRS(10)

    !Running direction, T represents forward direction,
    !F represents reverse direction
    SPE_COUNT = 0
    CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 1, 'READ_SPECIE:1' )
    READ( SUBSTRS(1:N), * ) SPE_NUMS

    allocate(SPE_NAME(SPE_NUMS), stat=as)

    allocate(SPE_MASS(SPE_NUMS), stat=as)
    SPE_MASS = 0

    DO I = 1, SPE_NUMS

      CALL SPLIT_ONE_LINE( FILE_ID, SUBSTRS, N, 3, 'READ_SPECIE:2' )
      READ( SUBSTRS(1:N), * ) SPE_COUNT, SPE_NAME(I), SPE_MASS(I)

    END DO

    print *, SPE_NAME
  END SUBROUTINE

end module

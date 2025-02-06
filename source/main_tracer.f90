PROGRAM GHSM

  USE CONTROL_VAR_MOD,  ONLY : INIT_CONTROL_VAR, CON_STEP, RUN_DIR
  USE TIME_MOD
  USE FORCING_MOD
  USE TOOL_MOD
  USE NETCDF_IO_MOD

  IMPLICIT NONE

  INCLUDE "CMN_SIZE"
  INCLUDE "CMN_GCTM"

  !local var
  INTEGER                   :: I, J, EW_INDEX, NS_INDEX
  INTEGER                   :: MASK_VAL(IIPAR, JJPAR)
  character(len=100)        :: file_path
  REAL*8                    :: EW(IIPAR, JJPAR)
  REAL*8                    :: NS(IIPAR, JJPAR)
  REAL*8                    :: HR_STT_ADJ(IIPAR, JJPAR)
  REAL*8                    :: SUM_VAL

  EW = 0d0
  NS = 0d0
  HR_STT_ADJ = 0d0

  CALL INIT_CONTROL_VAR()

  CALL INIT_TIME()

  CALL INIT_FORCING()

  DO WHILE( TIME_OBJ%RUN_FLAG )

    IF ( MOD(TIME_OBJ%MINUTE, 10) == 0 ) THEN
      !$OMP PARALLEL DO          &
      !$OMP DEFAULT( SHARED )    &
      !$OMP PRIVATE( I, J )
        DO J = 1, JJPAR
          DO I = 1, IIPAR
            HR_STT_ADJ(I, J) = HR_STT_ADJ(I, J) + hr_forcing(I, J) * MASK_VAL(I, J)
          END DO
        END DO
      !$OMP END PARALLEL DO
    END IF

    CALL RENEW_TIME()
    CALL OUTPUT_TIME()

    WRITE(file_path, "('/data/ese-zhengzy/gc/transport/transport_ratio.EW',  &
                            I4.4, I2.2, I2.2, '.', I2.2, I2.2)")             &
    TIME_OBJ%year, TIME_OBJ%month, TIME_OBJ%day, TIME_OBJ%hour, TIME_OBJ%MINUTE

    file_path = trim(file_path) // ".nc"

    CALL READ_NC( FILE_PATH, "EW", EW, 1152, 721 )

    WRITE(file_path, "('/data/ese-zhengzy/gc/transport/transport_ratio.NS',  &
                            I4.4, I2.2, I2.2, '.', I2.2, I2.2)")             &
    TIME_OBJ%year, TIME_OBJ%month, TIME_OBJ%day, TIME_OBJ%hour, TIME_OBJ%MINUTE

    file_path = trim(file_path) // ".nc"
    CALL READ_NC( FILE_PATH, "NS", NS, 1152, 721 )

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J )
      DO I = 1, IIPAR
        DO J = 9, JJPAR - 8
          IF ( HR_STT_ADJ(I, J) > 0 ) THEN
            IF (EW(I, J) > 0) THEN 
              EW_INDEX = I + 1 
            ELSE 
              EW_INDEX = I - 1
            ENDIF
            IF (NS(I, J) > 0) THEN 
              NS_INDEX = J + 1 
            ELSE 
              NS_INDEX = J - 1
            ENDIF
            IF (EW_INDEX > IIPAR) EW_INDEX = 1
            HR_STT_ADJ(I, J) = HR_STT_ADJ(I, J) * (1 - EW(I, J) - NS(I, J))
            HR_STT_ADJ(I, NS_INDEX) = HR_STT_ADJ(I, NS_INDEX) + NS(I, J) * HR_STT_ADJ(I, J)
            HR_STT_ADJ(EW_INDEX, J) = HR_STT_ADJ(EW_INDEX, J) + EW(I, J) * HR_STT_ADJ(I, J)
          ENDIF
        ENDDO
      ENDDO
    !$OMP END PARALLEL DO

  END DO

  file_path = trim(file_path) // ".nc"
  SUM_VAL = SUM(HR_STT_ADJ)

  !$OMP PARALLEL DO          &
  !$OMP DEFAULT( SHARED )    &
  !$OMP PRIVATE( I, J )
    DO I = 1, IIPAR
      DO J = 1, JJPAR
        HR_STT_ADJ(I, J) = HR_STT_ADJ(I, J) / SUM_VAL * 100
      ENDDO
    ENDDO
  !$OMP END PARALLEL DO
END PROGRAM

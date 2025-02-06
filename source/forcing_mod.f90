MODULE FORCING_MOD
  !Declare controlled var
  USE CONTROL_VAR_MOD,    ONLY : HR_EMS_PATH, LR_EMS_PATH
  USE CONTROL_VAR_MOD,    ONLY : HR_RESTART_PATH, RUN_DIR

  IMPLICIT NONE
  PUBLIC

  real*8, allocatable            :: hr_con     (:, :, :)
  real*8, allocatable            :: sum_mat    (:, :, :)
  real*8, allocatable            :: obs_data   (:, :, :)
  real*8, allocatable            :: hr_emis    (:, :, :)
  real*8, allocatable            :: lr_emis    (:, :, :)
  real*8, allocatable            :: hr_forcing (:, :)
  real*8, allocatable            :: lr_forcing (:, :, :)
  real*8, allocatable            :: init_field (:, :)
  real*8, allocatable            :: lr_con     (:, :, :)
  real*8, allocatable            :: lr_con_tran(:, :, :)
  real*8, allocatable            :: lr_25      (:, :)
  real*8, allocatable            :: hr_25      (:, :)

  PROCEDURE(), POINTER           :: add_forcing_hr
  PROCEDURE(), POINTER           :: add_forcing_lr

  CONTAINS

  SUBROUTINE INIT_FORCING()

    USE  NETCDF_IO_MOD,    ONLY : READ_NC

    INCLUDE  "CMN_SIZE"

    integer                :: as

    allocate(hr_con(IIPAR, JJPAR, LLPAR), stat=as)
    hr_con = 0d0

    allocate(lr_con(72, 46, 72), stat=as)
    lr_con = 0d0

    allocate(lr_con_tran(72, 46, 72), stat=as)
    lr_con_tran = 0d0

    allocate(obs_data(744, 72, 46), stat=as)
    obs_data = 0d0

    allocate(lr_25(144, 91), stat=as)
    lr_25 = 0d0

    allocate(hr_25(1152, 721), stat=as)
    hr_25 = 0d0

    IF ( RUN_DIR ) THEN

      allocate(hr_emis(IIPAR, JJPAR, 12), stat=as)
      hr_emis = 0d0

      allocate(init_field(IIPAR, JJPAR), stat=as)
      init_field = 0d0

      allocate(lr_emis(72, 46, 12), stat=as)
      lr_emis = 0d0

      call read_nc('/data/ese-zhengzy/gc/SO2.nc', 'SO2', hr_emis, IIPAR, JJPAR, 12)
      call read_nc('/data/ese-zhengzy/restart_SO2.nc', 'SO2', init_field, IIPAR, JJPAR)
      CALL READ_NC('/data/ese-zhengzy/china_so2_5.nc', 'SO2', obs_data, 744, 72, 46)
      print *, 'init_field'
      print *, sum(init_field)
      CALL RENEW_LR_STT()
    ELSE
      allocate(hr_forcing(IIPAR, JJPAR), stat=as)
      hr_forcing = 0d0

      call read_nc('/data/ese-hejl/gc/PKU/pm_grad_025.nc', 'gradient', hr_forcing, IIPAR, JJPAR)
    END IF

  END SUBROUTINE

  subroutine add_forcing_f(stt_adj, timestep, month)

    USE AIR_STATE_MOD,     ONLY : AD
    USE GRID_MOD,          ONLY : GET_AREA_M2

    INCLUDE  "CMN_SIZE"

    real*8, intent(inout)    :: stt_adj(IIPAR, JJPAR, LLPAR)
    integer, intent(in)      :: timestep, month
    integer                  :: I, J
    real*8                   :: AREA_M2

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J, AREA_M2 )
      DO J = 1, JJPAR
        AREA_M2 = GET_AREA_M2( J )
        DO I = 1, IIPAR
          !stt_adj(I, J, 1) = stt_adj(I, J, 1) + hr_emis(I, J, month)  &
          !                   * timestep * AREA_M2 / AD(I, J, 1)
          stt_adj(I, J, 1) = stt_adj(I, J, 1) + hr_forcing(I, J)
          !                   * timestep * AREA_M2 / AD(I, J, 1)
        END DO
      END DO
    !$OMP END PARALLEL DO

  end subroutine

  subroutine add_forcing_b(stt_adj)

    USE AIR_STATE_MOD,     ONLY : AD
    USE GRID_MOD,          ONLY : GET_AREA_M2

    INCLUDE  "CMN_SIZE"

    real*8, intent(inout)    :: stt_adj(IIPAR, JJPAR, LLPAR)
    integer                  :: I, J
    real*8                   :: AREA_M2

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J, AREA_M2 )
      DO J = 1, JJPAR
        AREA_M2 = GET_AREA_M2( J )
        DO I = 1, IIPAR
          stt_adj(I, J, 1) = stt_adj(I, J, 1) + hr_forcing(I, J)  &
                              * AREA_M2 / AD(I, J, 1)
        END DO
      END DO
    !$OMP END PARALLEL DO

  end subroutine

  subroutine add_forcing_45_f(lr_stt_adj)
    real*8, intent(inout)     :: lr_stt_adj(72, 46, 72)
    integer                   :: I, J

    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED )   &
    !$OMP PRIVATE( I, J )
      DO I = 1, 72
        DO J = 1, 46
          lr_stt_adj(I, J, 1) = lr_stt_adj(I, J, 1) + lr_emis(I, J, 2)
        END DO
      END DO
    !$OMP END PARALLEL DO

  end subroutine


  SUBROUTINE RENEW_LR_STT()

    USE TIME_MOD,        ONLY : TIME_OBJ, GET_NEXT_HOUR
    USE NETCDF_IO_MOD,   ONLY : READ_NC
    USE CONTROL_VAR_MOD, ONLY : SPE_NAME

    character(len=100)                :: file_path
    integer                           :: year, month, day, hour

    CALL GET_NEXT_HOUR( year, month, day, hour )

    WRITE(file_path, "('/data/ese-hejl/gc/new_forward/model.result.',   &
                        I4.4, I2.2, I2.2, '.', I2.2, I2.2)")            &
    year, month, day, hour

    file_path = trim(file_path) // "00.nc"

    CALL READ_NC( FILE_PATH, SPE_NAME(1), lr_con, 72, 46, 72 )

    IF ( .TRUE. ) THEN
      CALL REPLACE_OBS( LR_CON, OBS_DATA(day * 24 + hour, :, :) )
    END IF
  END SUBROUTINE

  SUBROUTINE RENEW_LR_25()

    USE TIME_MOD,        ONLY : TIME_OBJ
    USE NETCDF_IO_MOD,   ONLY : READ_NC
    USE CONTROL_VAR_MOD, ONLY : SPE_NAME
    USE ALGORITHM_MOD,   ONLY : double_line_interp
    

    character(len=100)                :: file_path

    WRITE(file_path, "('/data/ese-hejl/gc/new_forward/ourstep.ems.adj.', &
                        I4.4, I2.2, I2.2, '.', I2.2, I2.2, I2.2)")       &
    TIME_OBJ%year, TIME_OBJ%month, TIME_OBJ%day, TIME_OBJ%hour

    file_path = trim(file_path) // "00.nc"

    CALL READ_NC( FILE_PATH, SPE_NAME(1), lr_25, 144, 91 )

    !print *, '899999'
    !print *, sum(lr_25)

    CALL double_line_interp(lr_25, hr_25)
    !print *, '6666666'

  END SUBROUTINE


  SUBROUTINE REPLACE_OBS( LR_DATA, OBS_DATA )

    USE AIR_STATE_MOD,   ONLY : AIRVOL, AD

    REAL*8, INTENT(INOUT)    :: LR_DATA(72, 46, 72)
    REAL*8, INTENT(IN)       :: OBS_DATA(72, 46)

    INTEGER                  :: I, J

    DO I = 1,72
      DO J = 1,46
        IF ( NOT(ISNAN(OBS_DATA(I, J))) ) THEN
          LR_DATA(I, J, 1) = OBS_DATA(I, J) * AIRVOL(I, J, 1)  &
                             / AD(I, J, 1) * 1d-9
        END IF
      END DO
    END DO
  END SUBROUTINE

END MODULE

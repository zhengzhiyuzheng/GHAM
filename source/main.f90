PROGRAM GHSM

  USE PJC_PFIX_MOD,     ONLY : DO_PJC_PFIX
  USE PRESSURE_MOD,     ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
  USE PRESSURE_MOD,     ONLY : GET_AP, GET_BP
  USE PRESSURE_MOD,     ONLY : INIT_PRESSURE
  USE TPCORE_FVDAS_MOD, ONLY : TPCORE_FVDAS, Init_Tpcore
  USE GRID_MOD,         ONLY : COMPUTE_GRID, GET_YMID_R
  USE GRID_MOD,         ONLY : GET_AREA_M2
  USE OMP_LIB
  USE CONTROL_VAR_MOD,  ONLY : INIT_CONTROL_VAR, CON_STEP, RUN_DIR
  USE PBL_MIX_MOD,      ONLY : DO_PBL_MIX, DO_PBL_MIX_ADJ
  USE CONVECTION_MOD,   ONLY : NFCLDMX, NFCLDMX_ADJ
  USE TIME_MOD
  USE ALGORITHM_MOD,    ONLY : INIT_ALG, distribute, machine_predict
  USE ALGORITHM_MOD,    ONLY : interp_normal_surface
  USE FORCING_MOD
  USE AIR_STATE_MOD
  USE TRACER_MOD
  USE TOOL_MOD
  USE MACHINE_LEARN

  IMPLICIT NONE

  INCLUDE "CMN_SIZE"
  INCLUDE "CMN_GCTM"

  !local var
  INTEGER                   :: I, J, L, L2, N, N_DYN, K, as
  INTEGER                   :: JFIRST, JLAST, val_count
  REAL*8                    :: A_DIFF, D_DYN, TR_DIFF
  REAL*8                    :: SUMADA, decay_val
  REAL*8                    :: YMID_R(JJPAR)
  REAL*8                    :: Ap(LLPAR + 1), Bp(LLPAR + 1)
  REAL*8                    :: A_M2(JJPAR)
  REAL*8                    :: start_time, end_time, elapsed_time
  REAL*8                    :: TRACER_MASS(IIPAR, JJPAR, LLPAR)

  TRACER_MASS = 0d0
  D_DYN = 300

  !OPEN(UNIT=77, FILE='./mass.txt', STATUS='UNKNOWN', ACTION='WRITE', FORM='FORMATTED')

  !CALL OMP_SET_NUM_THREADS(5)

  !Basic program initialization
  CALL COMPUTE_GRID()
  CALL INIT_PRESSURE()

  !Custom program initialization
  print *, 'INIT_CONTROL_VAR'

  CALL INIT_CONTROL_VAR()
  print *, 'INIT_TIME'

  CALL INIT_TIME()
  print *, 'INIT_ALG'

  CALL INIT_ALG
  print *, TIME_OBJ%RUN_FLAG
  print *, 'INIT_AIR_STATE'

  CALL INIT_AIR_STATE()

  !Set surface pressure
  print *, 'INIT_FORCING'
  !CALL SET_FLOATING_PRESSURE( PSC )

  CALL INIT_FORCING()

  DO L = 1, LLPAR+1
    K = ( LLPAR + 1 ) - L + 1
    Ap(L) = GET_AP(K)
    Bp(L) = GET_BP(K)
  ENDDO

  DO J = 1, JJPAR
    A_M2(J) = GET_AREA_M2( J )
  ENDDO

  DO J = 1,JJPAR
    YMID_R(J) = GET_YMID_R(J)
  ENDDO

  print *, 'Init_Tpcore'
  call Init_Tpcore(IIPAR, JJPAR, LLPAR, JFIRST, JLAST, 0, 0, D_DYN, Re, YMID_R)
  print *, 'INIT_TRACER'
  CALL INIT_TRACER()
  print *, 'END INIT'
  print *, TIME_OBJ%RUN_FLAG

  CALL INIT_XG()

  DO WHILE( TIME_OBJ%RUN_FLAG )

    CALL OMP_SET_NUM_THREADS(5)
    start_time = omp_get_wtime()

    !CALL OUTPUT_TIME()
    IF ( MOD(TIME_OBJ%MINUTE, 10) == 0 ) THEN
      CALL ADD_FORCING_F( HR_STT_ADJ, 5, TIME_OBJ%MONTH )
    END IF
    PRINT *, 'HR_STT_ADJ1 IS ', SUM(HR_STT_ADJ)
    CALL RENEW_TIME()
    CALL OUTPUT_TIME()
    CALL RENEW_AIR_STATE()

    TRACER_MASS = HR_STT_ADJ(:,:,:,1)

    CALL CONVERT_UNITS( 2,  1, TCVV, AD, HR_STT_ADJ )

    IF ( NOT(RUN_DIR) ) THEN
      !print *, "Now is NFCLDMX_ADJ"
      !PRINT *, 'BEFORE IS ', SUM(HR_STT_ADJ(:, :, 1, 1))
      CALL NFCLDMX_ADJ(SPECIE, TCVV, CMFMC(:, :, 2:LLPAR+1), DTRAIN, HR_STT_ADJ, AD, BXHEIGHT)
      !PRINT *, 'AFTER IS ', SUM(HR_STT_ADJ(:, :, 1, 1))
      CALL DO_PBL_MIX_ADJ(HR_STT_ADJ, PBL, TCVV, AD, BXHEIGHT, CON_STEP)
    ENDIF

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED )     &
    !$OMP PRIVATE( I, J )
    DO J = 1, JJPAR
      DO I = 1, IIPAR

       ! True surface pressure at midpoint of dynamic timestep [hPa]
       P_TP1(I,J) = GET_PEDGE(I,J,1)

       ! True surface pressure at end of dynamic timestep [hPa]
       P_TP2(I,J) = PSC(I,J)

      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    XMASS = 0d0
    YMASS = 0d0

    CALL DO_PJC_PFIX( D_DYN, P_TP1, P_TP2,                        &
                        -u_wind,   -v_wind,       XMASS,     YMASS )

    P_TEMP_NEW = P_TP2

    IF (.FALSE.) THEN
      CALL SAVE_TRANSPORT_RATIO(XMASS(:, :, 1), YMASS(:, :, 1), P_TP1(:, :), AP, BP)
    END IF

    !PRINT *, 'u_wind IS ', sum(u_wind)
    !PRINT *, 'PS IS ', sum(P_TEMP_NEW)
    !PRINT *, 'DIFF ', sum(abs(P_TP1 - P_TP2))
    !PRINT *, 'SUM MASS ', sum(XMASS)

    UTMP(:,:,1:LLPAR) = -u_wind(:,:,LLPAR:1:-1)
    VTMP(:,:,1:LLPAR) = -v_wind(:,:,LLPAR:1:-1)

    CALL TPCORE_FVDAS( D_DYN,    Re,        IIPAR,    JJPAR,   &
                       LLPAR,    JFIRST,    JLAST,    0,       &
                       0,        SPECIE,         Ap,       Bp, &
                       UTMP,     VTMP,      P_TP1,    P_TP2,   &
                       P_TEMP,   HR_STT_ADJ(:,:,LLPAR:1:-1,:), &
                       3,        3,       7,          0,       &
                       XMASS(:,:,LLPAR:1:-1),                  &
                       YMASS(:,:,LLPAR:1:-1),                  &
                       MASSFLEW,                               &
                       MASSFLNS,                               &
                       MASSFLUP,    A_M2,  0,      0,     0,   &
                       .FALSE. )
    CALL SET_FLOATING_PRESSURE( P_TEMP_NEW )

    CALL CONVERT_UNITS( 1,  1, TCVV, AD, HR_STT_ADJ )

    PRINT *, 'TRANSPORT DIFF IS ', SUM(HR_STT_ADJ(:, :, :, 1) - TRACER_MASS)


    IF ( RUN_DIR ) THEN
      CALL DO_PBL_MIX(HR_STT_ADJ, PBL, TCVV, AD, BXHEIGHT, CON_STEP)
      print *, "Now is NFCLDMX"
      CALL NFCLDMX(SPECIE, TCVV, CMFMC, DTRAIN, HR_STT_ADJ, CON_STEP, AD)
    ENDIF


    PRINT *, 'HR_STT_ADJ4 IS ', SUM(HR_STT_ADJ)

    !IF ( MOD(TIME_OBJ%MINUTE, 10) == 0 ) THEN
    !  CALL ADD_FORCING_F( HR_STT_ADJ, 5, TIME_OBJ%MONTH )
    !END IF

    !CALL distribute( HR_STT_ADJ(:,:,:,1), lr_con )

    IF ( TIME_OBJ%MINUTE == 55 .OR. MOD(TIME_OBJ%MINUTE, 10) == 0 ) THEN
      !CALL RENEW_LR_25()
      !CALL interp_normal_surface(HR_STT_ADJ(:, :, :, 1), HR_25)

      !print *, 'Now is start'
      !CALL OMP_SET_NUM_THREADS(20)
      !CALL machine_predict(HR_STT_ADJ(:, :, :, 1))
      !CALL distribute( HR_STT_ADJ(:,:,:,1), lr_con )
      !CALL RENEW_LR_STT()
      
      IF ( TIME_OBJ%MINUTE == 0 ) THEN 
        !CALL SAVE_TRACER()
      ENDIF
      !stop
    END IF

    !CALL CHECK_MASS(HR_STT_ADJ(:, :, :, 1), AD, TRACER_MASS)
    !WRITE(77, *) SUM(TRACER_MASS)

    PRINT *, 'HR_STT_ADJ IS ', SUM(HR_STT_ADJ)

    end_time = omp_get_wtime()
    elapsed_time = end_time - start_time
    PRINT *, "Elapsed time (in seconds):", elapsed_time


  END DO

END PROGRAM

MODULE TOOL_MOD

  !This is a tool module that does not rely on any other modules. The
  !implementation of functions inside and the passing of external 
  !parameters must be displayed

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE CHECK_MASS(STT_ADJ, AD, RES)

    INCLUDE  "CMN_SIZE"
    REAL*8, intent(in)     :: STT_ADJ(:, :, :)
    REAL*8, intent(in)     :: AD(:, :, :)
    REAL*8, intent(inout)  :: RES(:, :, :)

    !local var
    INTEGER           :: I, J, L

    !$OMP PARALLEL DO NUM_THREADS(NUM_CORE)      &
    !$OMP DEFAULT( SHARED )                      &
    !$OMP PRIVATE( I, J, L )
      DO L = 1, LLPAR
        DO I = 1, IIPAR
          DO J = 1, JJPAR
            RES(I, J, L) = STT_ADJ(I, J, L) * AD(I, J, L)
          END DO
        END DO
      END DO
    !$OMP END PARALLEL DO
    
  END SUBROUTINE

  SUBROUTINE CONVERT_UNITS( IFLAG, N_TRACERS, TCVV, AD, STT )
!
      include "CMN_SIZE"     ! Size parameters
!
      INTEGER, INTENT(IN)    :: IFLAG

      ! Number of tracers
      INTEGER, INTENT(IN)    :: N_TRACERS

      REAL*8,  INTENT(IN)    :: TCVV(N_TRACERS)

      REAL*8,  INTENT(IN)    :: AD(IIPAR,JJPAR,LLPAR)
!
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
!
      INTEGER :: I, J, L, N

      SELECT CASE ( IFLAG )

         CASE ( 1 )

           !$OMP PARALLEL DO             &
           !$OMP DEFAULT( SHARED )       &
           !$OMP PRIVATE( I, J, L, N )   
            DO N = 1, N_TRACERS
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT(I,J,L,N) = STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
           !$OMP END PARALLEL DO

         CASE ( 2 )

            !$OMP PARALLEL DO           &
            !$OMP DEFAULT( SHARED )     &
            !$OMP PRIVATE( I, J, L, N ) 
            DO N = 1, N_TRACERS
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT(I,J,L,N) = STT(I,J,L,N) * AD(I,J,L) / TCVV(N)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO

         CASE ( 5 )
            !$OMP PARALLEL DO           &
            !$OMP DEFAULT( SHARED )     &
            !$OMP PRIVATE( I, J, L, N ) 
            DO N = 1, N_TRACERS
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT(I,J,L,N) = STT(I,J,L,N)/AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO

         CASE ( 6 )

            !$OMP PARALLEL DO            &
            !$OMP DEFAULT( SHARED )      &
            !$OMP PRIVATE( I, J, L, N )  
            DO N = 1, N_TRACERS
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT(I,J,L,N) = STT(I,J,L,N) * AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO

      END SELECT

      END SUBROUTINE CONVERT_UNITS


      SUBROUTINE SAVE_TRANSPORT_RATIO(X_MASS, Y_MASS, P, AP, BP)

        USE TIME_MOD,        ONLY : TIME_OBJ
        USE NETCDF_IO_MOD,   ONLY : SAVE_NC
        include "CMN_SIZE"     ! Size parameters

        REAL*8,  INTENT(IN)    :: X_MASS(IIPAR,JJPAR)
        REAL*8,  INTENT(IN)    :: Y_MASS(IIPAR,JJPAR)

        REAL*8,  INTENT(IN)    :: P(IIPAR,JJPAR)
        REAL*8,  INTENT(IN)    :: AP(LLPAR), BP(LLPAR)

        !local var
        REAL*8                 :: EW(IIPAR,JJPAR)
        REAL*8                 :: NS(IIPAR,JJPAR)

        INTEGER                :: I, J
        character(len=100)     :: file_path

        WRITE(file_path, "('/data/ese-zhengzy/gc/transport/transport_ratio.EW', &
                            I4.4, I2.2, I2.2, '.', I2.2, I2.2)")             &
        TIME_OBJ%year, TIME_OBJ%month, TIME_OBJ%day, TIME_OBJ%hour, TIME_OBJ%MINUTE

        file_path = trim(file_path) // ".nc"
        print *, maxval(X_MASS)
        !stop

        !$OMP PARALLEL DO           &
        !$OMP DEFAULT( SHARED )     &
        !$OMP PRIVATE( I, J )
          DO I = 1, 1152
            DO J = 1, 721
              EW(I, J) = X_MASS(I, J) / (AP(72) + (BP(72) - BP(71)) * P(I, J))
              NS(I, J) = Y_MASS(I, J) / (AP(72) + (BP(72) - BP(71)) * P(I, J))
            END DO
          END DO
        !$OMP END PARALLEL DO

        CALL SAVE_NC( FILE_PATH, "EW", EW, 1152, 721 )

        WRITE(file_path, "('/data/ese-zhengzy/gc/transport/transport_ratio.NS', &
                            I4.4, I2.2, I2.2, '.', I2.2, I2.2)")             &
        TIME_OBJ%year, TIME_OBJ%month, TIME_OBJ%day, TIME_OBJ%hour, TIME_OBJ%MINUTE

        file_path = trim(file_path) // ".nc"
        CALL SAVE_NC( FILE_PATH, "NS", NS, 1152, 721 )

      END SUBROUTINE

END MODULE

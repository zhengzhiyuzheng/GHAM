MODULE ALGORITHM_MOD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_alg
  PUBLIC :: distribute
  PUBLIC :: interp_normal_surface
  PUBLIC :: double_line_interp
  PUBLIC :: machine_predict

  TYPE, PUBLIC :: INTP
    PROCEDURE(), NOPASS, POINTER :: D2
    PROCEDURE(), NOPASS, POINTER :: D3
  END TYPE INTP

  TYPE(INTP), PUBLIC             :: INTERP_OBJ

  REAL*8                         :: RATIO(1152, 721)

  CONTAINS

  subroutine init_alg()

    INTERP_OBJ%D2 => linear_interp_2d
    INTERP_OBJ%D3 => linear_interp_3d

  end subroutine

  subroutine linear_interp_3d(data1, data2, res_con, I_SIZE, J_SIZE, L_SIZE)

    USE TIME_MOD,    ONLY : TIME_OBJ

    real*8, intent(in)    :: data1(I_SIZE, J_SIZE, L_SIZE)
    real*8, intent(in)    :: data2(I_SIZE, J_SIZE, L_SIZE)
    real*8, intent(inout) :: res_con(I_SIZE, J_SIZE, L_SIZE)
    integer, intent(in)   :: I_SIZE, J_SIZE, L_SIZE

    integer     :: I, J, L

    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED )           &
    !$OMP PRIVATE( I, J, L )
      do L = 1, L_SIZE
        do J = 1, J_SIZE
           do I = 1, I_SIZE
             res_con(I, J, L) = data1(I, J, L)    &
              + (data2(I, J, L) - data1(I, J, L)) &
              * TIME_OBJ%PRE_RATE
           end do
        end do
      end do
    !$OMP END PARALLEL DO
  end subroutine

  subroutine linear_interp_2d(data1, data2, res_con, I_SIZE, J_SIZE)

    USE TIME_MOD,    ONLY : TIME_OBJ

    real*8, intent(in)    :: data1(I_SIZE, J_SIZE)
    real*8, intent(in)    :: data2(I_SIZE, J_SIZE)
    real*8, intent(inout) :: res_con(I_SIZE, J_SIZE)
    integer, intent(in)   :: I_SIZE, J_SIZE

    integer     :: I, J

    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED )           &
    !$OMP PRIVATE( I, J )
        do J = 1, J_SIZE
           do I = 1, I_SIZE
             res_con(I, J) = data1(I, J)   &
              + (data2(I, J) - data1(I, J)) &
              * TIME_OBJ%PRE_RATE
           end do
        end do
    !$OMP END PARALLEL DO
  end subroutine

  subroutine machine_predict(hr_data)
    
    USE TIME_MOD,      ONLY : TIME_OBJ
    USE AIR_STATE_MOD, ONLY : INPUT_DATA, GET_PAR 
    USE MACHINE_LEARN
    use iso_c_binding

    include "CMN_SIZE"

    real*8, intent(inout)    :: hr_data(IIPAR, JJPAR, LLPAR)
    REAL(C_FLOAT)            :: dot_data(12), ans
    integer                  :: I, J, point_index

    IF ( TIME_OBJ%MINUTE == 55 ) THEN
    RATIO = 0d0
    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED )   &
    !$OMP PRIVATE( I, J, dot_data, ans)
      do I = 1, IIPAR
        do J = 1, JJPAR
          dot_data = get_par(input_data, i, j)
          ans = 0d0
          call predict(booster_ptr, dot_data, ans)
          !print *, 'this is ans', ans
          !print *, 'Now is stio'
          !stop
          !hr_data(i, j, 1) = hr_data(i, j, 1) * ans
          RATIO(i, j) = ans ** (0.2)
        enddo
      enddo
    !$OMP END PARALLEL DO
    ELSE
    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED )   &
    !$OMP PRIVATE( I, J)
      do I = 1, IIPAR
        do J = 1, JJPAR
          hr_data(i, j, 1) = hr_data(i, j, 1) * RATIO(I, J)
        enddo
      enddo
    !$OMP END PARALLEL DO
    ENDIF
  end subroutine

  subroutine distribute( hr_data, lr_data )

    USE TIME_MOD,    ONLY : TIME_OBJ

    include "CMN_SIZE"

    real*8, intent(inout)    :: hr_data(IIPAR, JJPAR, LLPAR)
    real*8, intent(in)       :: lr_data(72, 46, 72)

    !local var
    real*8        :: num_count(72, 46, 72)
    real*8        :: sum_mat  (72, 46, 72)
    integer       :: I, J, L, x, y

    !init
    num_count = 0d0
    sum_mat   = 0d0

    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED )   &
    !$OMP PRIVATE( I, J, L, x, y )
      do L = 1, 1
        do J = 1, JJPAR
          do I = 1, IIPAR
            if (I < 8 .or. I >= 1144) then
              x = 1
            else
              x = (I * 0.3125 - 2.5) / 5 + 2
            end if
            if (J < 8) then
              y = 1
            else
              y = (J * 0.25 - 2) / 4 + 2
            end if
            sum_mat(x, y, L) = sum_mat(x, y, L) + abs(hr_data(I, J, L))
            num_count(x, y, L) = num_count(x, y, L) + 1
          end do
        end do
      end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J, L, x, y )
      do L = 1, 1
        do J = 1, JJPAR
          do I = 1, IIPAR
            if (I < 8 .or. I >= 1144) then
              x = 1
            else
              x = (I * 0.3125 - 2.5) / 5 + 2
            end if
            if (J < 8) then
              y = 1
            else
              y = (J * 0.25 - 2) / 4 + 2
            end if
            if (sum_mat(x, y, L) /= 0) then
              hr_data(I, J, L) = abs(hr_data(I, J, L)) / sum_mat(x, y, L)  &
              * lr_data(x, y, L) * num_count(x, y, L)
             !* (sum_mat(x, y, L) + (lr_data(x, y, L) - sum_mat(x, y, L)) *  &
             !DBLE(TIME_OBJ%MINUTE + TIME_OBJ%TIME_STEP) / 60 ) * num_count(x, y, L)
            else
              hr_data(I, J, L) = lr_data(x, y, L)
            end if
          end do
        end do
      end do
    !$OMP END PARALLEL DO
  end subroutine

  subroutine interp_normal_surface(hr_data, hr_25)

    include "CMN_SIZE"

    real*8, intent(in)       :: hr_25(:, :)
    real*8, intent(inout)    :: hr_data(IIPAR, JJPAR, LLPAR)

    !local var
    integer                  :: I, J
    real*8                   :: total_sen, total_all


    !hr_data(:, :, 1) = hr_25
    !return
    total_sen = 0.9 * sum(hr_data(:, :, 1))

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J )
      DO I = 1, IIPAR
        DO J = 1, JJPAR
          if (hr_data(I, J, 1) /= 0 .AND. hr_25(i, j) /= 0) then
            hr_data(I, J, 1) = sqrt(abs(hr_data(I, J, 1) * hr_25(I, J)))  &
                               * hr_25(I, J) / abs(hr_25(I, J))
          end if
        END DO
      END DO
    !$OMP END PARALLEL DO

    !total_sen = sum(hr_25)
    total_all = sum(hr_data(:, :, 1))
    total_all = max(total_all, 1.0)

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J )
      DO I = 1, IIPAR
        DO J = 1, JJPAR
          hr_data(I, J, 1) = total_sen * hr_data(I, J, 1) / total_all
        END DO
      END DO
    !$OMP END PARALLEL DO

  end subroutine

  real*8 function deg2rad(degree)
      real*8, intent(in) :: degree
      real*8, parameter :: pi = 3.141592653589793
      deg2rad = degree * pi / 180.0
  end function deg2rad

  function cal_distance(lon1_in, lat1_in, lon2_in, lat2_in) result(distance)
    real*8, intent(in)  :: lon1_in, lat1_in, lon2_in, lat2_in
    real*8              :: lon1, lat1, lon2, lat2
    real*8              :: distance, earth_r, dlat, dlon, c, a

    lon1 = deg2rad(lon1_in)
    lat1 = deg2rad(lat1_in)
    lon2 = deg2rad(lon2_in)
    lat2 = deg2rad(lat2_in)

    earth_r = 6317

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))

    distance = earth_r * c

  end function


  subroutine double_line_interp(in_data, out_data)

    include "CMN_SIZE"

    real*8, intent(in)  :: in_data(:, :)
    real*8, intent(out) :: out_data(:, :)

    !local var
    integer             :: i, j, lon_index_left, lat_index_low, lon_index_right, lat_index_up
    real*8              :: distance, now_lon, now_lat
    real*8              :: lon_left, lon_right, lat_low, lat_up
    real*8              :: up_val, low_val


    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( PRIVATE )   &
    !$OMP SHARED( IN_DATA, OUT_DATA )
    do i = 1, IIPAR
      do j = 9, JJPAR - 8
        !print *, i, j
        now_lon = 0.3125 * (i - 1)
        now_lat = 0.25 * (j - 1)

        lon_index_left = 1 + floor(now_lon / 2.5)
        lat_index_low = 1 + floor(now_lat / 2)

        lon_index_right = 1 + lon_index_left
        lat_index_up = 1 + lat_index_low

        if (lon_index_right > 144) then
          lon_index_right = 1
        end if
        
        lon_left = (lon_index_left - 1) * 2.5
        lon_right = (lon_index_right - 1) * 2.5
        lat_low = (lat_index_low - 1) * 2
        lat_up = (lat_index_up - 1) * 2

        if ( lon_right == 0 ) lon_right = 360

        !print *, lon_index_right, lat_index_up
        !if (now_lon > lon_right .or. now_lon < lon_left) then
        !  print *, now_lon, lon_left, lon_right
        !  print *, 'lon is wrong'
        !  stop
        !end if

        !if (now_lat > lat_up .or. now_lat < lat_low) then
        !  print *, now_lat, lat_low, lat_up
        !  print *, 'lat is wrong'
        !  stop
        !end if
        !print *, 'Now is ', now_lon, now_lat, lon_right, lat_up
        
        up_val = (now_lon - lon_left) / (lon_right - lon_left) *  &
                  in_data(lon_index_right, lat_index_up) +         &
                 (lon_right - now_lon) / (lon_right - lon_left) * &
                  in_data(lon_index_left, lat_index_up)

        low_val = (now_lon - lon_left) / (lon_right - lon_left) *  &
                  in_data(lon_index_right, lat_index_low) +         &
                 (lon_right - now_lon) / (lon_right - lon_left) * &
                  in_data(lon_index_left, lat_index_low)

        out_data(i, j) = (now_lat - lat_low) / (lat_up - lat_low) *  up_val &
                         + (lat_up - now_lat) / (lat_up - lat_low) * low_val

      end do
    end do
    !$OMP END PARALLEL DO
    !stop
  end subroutine

END MODULE

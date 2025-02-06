MODULE NETCDF_IO_MOD

  IMPLICIT NONE
  PUBLIC

  INTERFACE READ_NC
    module procedure read_4d
    module procedure read_3d
    module procedure read_2d
  END INTERFACE

  INTERFACE READ_NC_VL
    module procedure read_4d_vl
    module procedure read_3d_vl
    module procedure read_2d_vl
  END INTERFACE

  INTERFACE SAVE_NC
    module procedure save_4d
    module procedure save_3d
    module procedure save_2d
  END INTERFACE

  CONTAINS

  SUBROUTINE READ_4D( FILE_PATH, VAR_NAME, DATA_CON, I, J, L, K )

    USE NETCDF

    character(len=*), intent(in)     :: file_path
    character(len=*), intent(in)     :: var_name
    integer,          intent(in)     :: I, J, L, K
    REAL*8,           intent(inout)  :: DATA_CON(I, J, L, K)

    integer                          :: kid, var_id, res

    print *, file_path

    res = NF90_OPEN(TRIM(file_path),nf90_nowrite,kid)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error opening NetCDF file"
      STOP
    END IF

    print *, var_name
    res = nf90_inq_varid(kid, TRIM(var_name), var_id)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error getting variable ID"
      STOP
    END IF

    res = nf90_get_var(kid, var_id, data_con)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error reading variable data"
      STOP
    END IF

    res = nf90_close(kid)
    
  END SUBROUTINE

  SUBROUTINE READ_4D_VL( FILE_PATH, VAR_NAME, DATA_CON, START_LOC, START_COUNT )

    USE NETCDF

    character(len=*), intent(in)     :: file_path
    character(len=*), intent(in)     :: var_name
    integer,          intent(in)     :: START_LOC(:), START_COUNT(:)
    REAL*8,           intent(inout)  :: DATA_CON(:, :, :, :)

    integer                          :: kid, var_id, res

    print *, file_path

    res = NF90_OPEN(TRIM(file_path),nf90_nowrite,kid)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error opening NetCDF file"
      STOP
    END IF

    print *, var_name
    res = nf90_inq_varid(kid, TRIM(var_name), var_id)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error getting variable ID"
      STOP
    END IF

    res = nf90_get_var(kid, var_id, data_con, start=START_LOC, count=START_COUNT)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error reading variable data"
      STOP
    END IF
    
    res = nf90_close(kid)

  END SUBROUTINE

  SUBROUTINE READ_3D_VL( FILE_PATH, VAR_NAME, DATA_CON, START_LOC, START_COUNT )

    USE NETCDF

    character(len=*), intent(in)     :: file_path
    character(len=*), intent(in)     :: var_name
    integer,          intent(in)     :: START_LOC(:), START_COUNT(:)
    REAL*8,           intent(inout)  :: DATA_CON(:, :, :)

    integer                          :: kid, var_id, res

    print *, file_path
    print *, "READ_3D_VL"

    res = NF90_OPEN(TRIM(file_path),nf90_nowrite,kid)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error opening NetCDF file"
      STOP
    END IF

    print *, var_name
    res = nf90_inq_varid(kid, TRIM(var_name), var_id)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error getting variable ID"
      STOP
    END IF

    res = nf90_get_var(kid, var_id, data_con, start=START_LOC, count=START_COUNT)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error reading variable data"
      STOP
    END IF

    res = nf90_close(kid)

  END SUBROUTINE

  SUBROUTINE READ_3D( FILE_PATH, VAR_NAME, DATA_CON, I, J, L )

    USE NETCDF

    character(len=*), intent(in)     :: file_path
    character(len=*), intent(in)     :: var_name
    integer,          intent(in)     :: I, J, L
    REAL*8,           intent(inout)  :: DATA_CON(I, J, L)

    integer                          :: kid, var_id, res

    print *, file_path

    res = NF90_OPEN(TRIM(file_path),nf90_nowrite,kid)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error opening NetCDF file"
      STOP
    END IF

    print *, var_name
    res = nf90_inq_varid(kid, TRIM(var_name), var_id)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error getting variable ID"
      STOP
    END IF

    res = nf90_get_var(kid, var_id, data_con)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error reading variable data"
      STOP
    END IF

    res = nf90_close(kid)

  END SUBROUTINE

  SUBROUTINE READ_2D( FILE_PATH, VAR_NAME, DATA_CON, I, J )

    USE NETCDF

    character(len=*), intent(in)     :: file_path
    character(len=*), intent(in)     :: var_name
    integer,          intent(in)     :: I, J
    REAL*8,           intent(inout)  :: DATA_CON(I, J)

    integer                          :: kid, var_id, res

    print *, file_path

    res = NF90_OPEN(TRIM(file_path),nf90_nowrite,kid)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error opening NetCDF file"
      STOP
    END IF

    print *, var_name
    res = nf90_inq_varid(kid, TRIM(var_name), var_id)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error getting variable ID"
      STOP
    END IF
    res = nf90_get_var(kid, var_id, data_con)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error reading variable data"
      STOP
    END IF
    res = nf90_close(kid)

  END SUBROUTINE

  SUBROUTINE READ_2D_VL( FILE_PATH, VAR_NAME, DATA_CON, START_LOC, START_COUNT)

    USE NETCDF

    character(len=*), intent(in)     :: file_path
    character(len=*), intent(in)     :: var_name
    integer,          intent(in)     :: START_LOC(:), START_COUNT(:)
    REAL*8,           intent(inout)  :: DATA_CON(:, :)

    integer                          :: kid, var_id, res

    print *, file_path
    print *, "READ_2D_VL"

    res = NF90_OPEN(TRIM(file_path),nf90_nowrite,kid)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error opening NetCDF file"
      STOP
    END IF

    res = nf90_inq_varid(kid, TRIM(var_name), var_id)
    IF (res /= nf90_noerr) THEN
      PRINT *, "Error getting variable ID"
      STOP
    END IF

    res = nf90_get_var(kid, var_id, data_con, start=START_LOC, count=START_COUNT)

    IF (res /= nf90_noerr) THEN
      PRINT *, "Error reading variable data"
      STOP
    END IF

    res = nf90_close(kid)

  END SUBROUTINE


  subroutine save_4d(file_name, var_name, data_con, I, J, L, K)
    use netcdf
    character(len=*)     :: file_name
    character(len=*)     :: var_name
    integer              :: I, J, L, K
    real*8               :: data_con(I, J, L, K)
    !local var
    integer :: ncid, varid, dimid_x, dimid_y, dimid_z, dimid_s, dimids(4)

    call check(nf90_create(trim(file_name), NF90_CLOBBER, ncid))

    call check(nf90_def_dim(ncid, "lon", I, dimid_x))
    call check(nf90_def_dim(ncid, "lat", J, dimid_y))
    call check(nf90_def_dim(ncid, "alt", L, dimid_z))
    call check(nf90_def_dim(ncid, "spe", K, dimid_s))

    dimids = (/dimid_x, dimid_y, dimid_z, dimid_s/)

    call check(nf90_def_var(ncid, trim(var_name), NF90_REAL, dimids, varid))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid, data_con))

    call check(nf90_close(ncid))
  end subroutine

  subroutine save_3d(file_name, var_name, data_con, I, J, L)
    use netcdf
    character(len=*)     :: file_name
    character(len=*)     :: var_name
    integer              :: I, J, L
    real*8               :: data_con(I, J, L)
    !local var
    integer :: ncid, varid, dimid_x, dimid_y, dimid_z, dimids(3)
      
    call check(nf90_create(trim(file_name), NF90_CLOBBER, ncid))
  
    call check(nf90_def_dim(ncid, "lon", I, dimid_x))
    call check(nf90_def_dim(ncid, "lat", J, dimid_y))
    call check(nf90_def_dim(ncid, "alt", L, dimid_z))

    dimids = (/dimid_x, dimid_y, dimid_z/)

    call check(nf90_def_var(ncid, trim(var_name), NF90_REAL, dimids, varid))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid, data_con))

    call check(nf90_close(ncid))
  end subroutine

  subroutine save_2d(file_name, var_name, data_con, I, J)
    use netcdf
    character(len=*)     :: file_name
    character(len=*)     :: var_name
    integer              :: I, J
    real*8               :: data_con(I, J)
    !local var
    integer :: ncid, varid, dimid_x, dimid_y, dimids(2)

    call check(nf90_create(trim(file_name), NF90_CLOBBER, ncid))

    call check(nf90_def_dim(ncid, "lon", I, dimid_x))
    call check(nf90_def_dim(ncid, "lat", J, dimid_y))

    dimids = (/dimid_x, dimid_y/)

    call check(nf90_def_var(ncid, trim(var_name), NF90_REAL, dimids, varid))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid, data_con))

    call check(nf90_close(ncid))
  end subroutine

  subroutine check(status)
    use netcdf
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      write(*, *) "Error in NETCDF SAVE", status
      stop
    end if
  end subroutine check

END MODULE


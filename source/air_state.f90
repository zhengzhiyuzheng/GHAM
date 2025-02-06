MODULE AIR_STATE_MOD

  USE CONTROL_VAR_MOD,      ONLY : SPE_NAME, SPE_MASS
  USE CONTROL_VAR_MOD,      ONLY : HR_MET_PATH, RUN_DIR
  USE TIME_MOD,             ONLY : TIME_OBJ
  use iso_c_binding
  IMPLICIT NONE
  PUBLIC

  interface READ_GEOS_FP
    module procedure READ_GEOS_FP_4D
    module procedure READ_GEOS_FP_3D
  end interface


  type point_array
    real*8, pointer              :: val  (:, :)
  endtype


  logical                        :: first
  logical                        :: XG_FIRST = .TRUE.
  real*8, allocatable            :: u_wind     (:, :, :)
  real*8, allocatable            :: v_wind     (:, :, :)
  real*8, allocatable            :: psc        (:, :)
  real*8, pointer                :: u_wind1    (:, :, :)
  real*8, pointer                :: v_wind1    (:, :, :)
  real*8, pointer                :: u_wind2    (:, :, :)
  real*8, pointer                :: v_wind2    (:, :, :)
  real*8, pointer                :: psc1       (:, :)
  real*8, pointer                :: psc2       (:, :)
  real*8, pointer                :: PBL        (:, :)
  real*8, allocatable, target    :: org_u_wind1(:, :, :, :)
  real*8, allocatable, target    :: org_v_wind1(:, :, :, :)
  real*8, allocatable, target    :: org_u_wind2(:, :, :, :)
  real*8, allocatable, target    :: org_v_wind2(:, :, :, :)
  real*8, allocatable, target    :: org_p1     (:, :, :)
  real*8, allocatable, target    :: org_p2     (:, :, :)

  real*8, allocatable            :: P_TP1      (:, :)
  real*8, allocatable            :: P_TP2      (:, :)
  real*8, allocatable            :: P_TEMP     (:, :)
  real*8, allocatable            :: P_TEMP_NEW (:, :)
  real*8, allocatable            :: AD_A       (:, :, :)
  real*8, allocatable            :: AD_B       (:, :, :)
  real*8, allocatable            :: UTMP       (:, :, :)
  real*8, allocatable            :: VTMP       (:, :, :)
  real*8, allocatable            :: XMASS      (:, :, :)
  real*8, allocatable            :: YMASS      (:, :, :)
  real*8, allocatable            :: MASS       (:, :, :)
  real*8, allocatable            :: MASS_45    (:, :, :)
  real*8, allocatable            :: MASSFLEW   (:, :, :, :)
  real*8, allocatable            :: MASSFLNS   (:, :, :, :)
  real*8, allocatable            :: MASSFLUP   (:, :, :, :)
  real*8, allocatable            :: MONTH_DATA (:, :, :)
  real*8, allocatable            :: TEM_DATA   (:, :, :)
  real*8, allocatable            :: T       (:, :, :)
  real*8, pointer                :: T1      (:, :, :)
  real*8, pointer                :: T2      (:, :, :)
  real*8, allocatable, target    :: org_t1  (:, :, :, :)
  real*8, allocatable, target    :: org_t2  (:, :, :, :)

  real*8, allocatable            :: BXHEIGHT(:, :, :)
  real*8, allocatable            :: AD      (:, :, :)
  real*8, allocatable            :: AIRVOL  (:, :, :)
  real*8, allocatable, target    :: ORG_PBL (:, :, :)

  real*8, allocatable, target    :: ORG_CMFMC (:, :, :, :)
  real*8, allocatable, target    :: ORG_DTRAIN(:, :, :, :)

  real*8, pointer                :: CMFMC     (:, :, :)
  real*8, pointer                :: DTRAIN    (:, :, :)


  !Now define XGboost input var
  type(point_array)              :: INPUT_DATA(12)
  !real*8, pointer                 :: INPUT_DATA(12, :, :)

  real*8, allocatable, target    :: XG_U10M(:, :) !A1
  real*8, allocatable, target    :: XG_PBL(:, :)
  !real*8, allocatable, target    :: XG_LWI(:, :) !A1
  real*8, allocatable, target    :: XG_USTAR(:, :) !A1
  !real*8, allocatable, target    :: XG_UWND(:, :)
  real*8, allocatable, target    :: XG_ALBD(:, :) !A1
  real*8, allocatable, target    :: XG_Z0(:, :) !A1
  real*8, allocatable, target    :: XG_QL(:, :) !A3cld
  real*8, allocatable, target    :: XG_SEAICE80(:, :) !A1
  real*8, allocatable, target    :: XG_RH(:, :) !A3dyn
  real*8, allocatable, target    :: XG_TSKIN(:, :) !A1
  real*8, allocatable, target    :: XG_BXHEIGHT(:, :)
  real*8, allocatable, target    :: XG_AD      (:, :)
  real*8, allocatable, target    :: XG_AIRVOL  (:, :)

  CONTAINS

  SUBROUTINE INIT_XG_VAR()

    include  "CMN_SIZE"

    !local var
    INTEGER              :: as


    allocate(XG_U10M(IIPAR, JJPAR), stat=as)
    XG_U10M = 0d0

    allocate(XG_PBL(IIPAR, JJPAR), stat=as)
    XG_PBL = 0d0

    !allocate(XG_LWI(IIPAR, JJPAR), stat=as)
    !XG_LWI = 0d0

    allocate(XG_USTAR(IIPAR, JJPAR), stat=as)
    XG_USTAR = 0d0

    !allocate(XG_UWND(IIPAR, JJPAR), stat=as)
    !XG_UWND = 0d0

    allocate(XG_ALBD(IIPAR, JJPAR), stat=as)
    XG_ALBD = 0d0

    allocate(XG_Z0(IIPAR, JJPAR), stat=as)
    XG_Z0 = 0d0

    allocate(XG_QL(IIPAR, JJPAR), stat=as)
    XG_QL = 0d0

    allocate(XG_SEAICE80(IIPAR, JJPAR), stat=as)
    XG_SEAICE80 = 0d0

    allocate(XG_RH(IIPAR, JJPAR), stat=as)
    XG_RH = 0d0

    allocate(XG_TSKIN(IIPAR, JJPAR), stat=as)
    XG_TSKIN = 0d0

    allocate(XG_BXHEIGHT(IIPAR, JJPAR), stat=as)
    XG_BXHEIGHT = 0d0

    allocate(XG_AD(IIPAR, JJPAR), stat=as)
    XG_AD = 0d0

    allocate(XG_AIRVOL(IIPAR, JJPAR), stat=as)
    XG_AIRVOL = 0d0

    call RENEW_XG_VAR()

    INPUT_DATA(1 )%VAL => XG_U10M
    INPUT_DATA(2 )%VAL => XG_PBL
    INPUT_DATA(3 )%VAL => XG_USTAR !A1
    INPUT_DATA(4 )%VAL => XG_ALBD !A1
    INPUT_DATA(5 )%VAL => XG_Z0 !A1
    INPUT_DATA(6 )%VAL => XG_QL !A3cld
    INPUT_DATA(7 )%VAL => XG_SEAICE80 !A1
    INPUT_DATA(8 )%VAL => XG_RH !A3dyn
    INPUT_DATA(9 )%VAL => XG_TSKIN !A1
    INPUT_DATA(10)%VAL => XG_BXHEIGHT
    INPUT_DATA(11)%VAL => XG_AD
    INPUT_DATA(12)%VAL => XG_AIRVOL

    !INPUT_DATA(1, :, :) => XG_U10M
    !INPUT_DATA(2, :, :) => XG_PBL
    !INPUT_DATA(3, :, :) => XG_USTAR !A1
    !INPUT_DATA(4, :, :) => XG_ALBD !A1
    !INPUT_DATA(5, :, :) => XG_Z0 !A1
    !INPUT_DATA(6, :, :) => XG_QL !A3cld
    !INPUT_DATA(7, :, :) => XG_SEAICE80 !A1
    !INPUT_DATA(8, :, :) => XG_RH !A3dyn
    !INPUT_DATA(9, :, :) => XG_TSKIN !A1
    !INPUT_DATA(10, :, :) => XG_BXHEIGHT
    !INPUT_DATA(11, :, :) => XG_AD
    !INPUT_DATA(12, :, :) => XG_AIRVOL

  END SUBROUTINE

  FUNCTION GET_PAR(INPUT_PO, I, J) RESULT(VAL)

  TYPE(point_array), INTENT(IN)    :: INPUT_PO(12)

  REAL(C_FLOAT)                    :: VAL(12)

  INTEGER, INTENT(IN)              :: I, J
  INTEGER                          :: K

  VAL = 0d0
  DO K = 1,12
    VAL(K) = INPUT_PO(K)%VAL(I, J)
  ENDDO

  END FUNCTION

  SUBROUTINE RENEW_XG_VAR()

    USE NETCDF_IO_MOD,   ONLY : READ_NC_VL

    include  "CMN_SIZE"

    character(len=100)               :: file_path
    integer                          :: YEAR, MONTH, DAY
    integer                          :: LOC_START_3(3), LOC_COUNT_3(3)
    integer                          :: LOC_START_4(4), LOC_COUNT_4(4)
    integer                          :: today_minute
  
    IF (TIME_OBJ%IS_NEW_HOUR .or. XG_FIRST) THEN

      !A1 var
      YEAR  = TIME_OBJ%YEAR
      MONTH = TIME_OBJ%MONTH
      DAY   = TIME_OBJ%DAY

      LOC_START_3 = (/1, 1, TIME_OBJ%HOUR+1/)
      LOC_COUNT_3 = (/IIPAR, JJPAR, 1/)
     
      WRITE(file_path, "('GEOSFP.', I4.4, I2.2, I2.2, '.')") &
            YEAR, MONTH, DAY

      file_path = trim(HR_MET_PATH) // trim(file_path) //    &
                  trim("A1") // ".025x03125.nc"

      print *, "Now is read" 
      CALL READ_NC_VL( FILE_PATH, "U10M", XG_U10M, LOC_START_3, LOC_COUNT_3)
      CALL READ_NC_VL( FILE_PATH, "USTAR", XG_USTAR, LOC_START_3, LOC_COUNT_3)
      CALL READ_NC_VL( FILE_PATH, "ALBEDO", XG_ALBD, LOC_START_3, LOC_COUNT_3)
      CALL READ_NC_VL( FILE_PATH, "Z0M", XG_Z0, LOC_START_3, LOC_COUNT_3)
      CALL READ_NC_VL( FILE_PATH, "TS", XG_TSKIN, LOC_START_3, LOC_COUNT_3)
      CALL READ_NC_VL( FILE_PATH, "SEAICE80", XG_SEAICE80, LOC_START_3, LOC_COUNT_3)

    END IF

    IF ( TIME_OBJ%IS_RENEW_A3 .or. XG_FIRST ) THEN
      !A3 VAR
      !today_minute = TIME_OBJ%HOUR * 60 + TIME_OBJ%MINUTE
      !IF (  ) THEN
      !ELSE
      YEAR  = TIME_OBJ%YEAR
      MONTH = TIME_OBJ%MONTH
      DAY   = TIME_OBJ%DAY
      !END IF

      LOC_START_4 = (/1, 1, 1, TIME_OBJ%RIGHT_INDEX/)
      LOC_COUNT_4 = (/IIPAR, JJPAR, 1, 1/)

      WRITE(file_path, "('GEOSFP.', I4.4, I2.2, I2.2, '.')") &
            YEAR, MONTH, DAY

      file_path = trim(HR_MET_PATH) // trim(file_path) //    &
                  trim("A3cld") // ".025x03125.nc"
     
      CALL READ_NC_VL( FILE_PATH, "QL", XG_QL, LOC_START_4, LOC_COUNT_4)

      !dyn
      WRITE(file_path, "('GEOSFP.', I4.4, I2.2, I2.2, '.')") &
            YEAR, MONTH, DAY

      file_path = trim(HR_MET_PATH) // trim(file_path) //    &
                  trim("A3dyn") // ".025x03125.nc"
      CALL READ_NC_VL( FILE_PATH, "RH", XG_RH, LOC_START_4, LOC_COUNT_4)

    END IF

    XG_BXHEIGHT = BXHEIGHT(:, :, 1)
    XG_PBL      = PBL(:, :)
    XG_AD       = 16 * 14 * AD(:, :, 1)
    XG_AIRVOL   = 16 * 14 * AIRVOL(:, :, 1)
    
    IF ( XG_FIRST ) XG_FIRST = .FALSE.
  END SUBROUTINE


  SUBROUTINE READ_GEOS_FP_3D(data_con, var_name, produce, FLAG)

    USE NETCDF_IO_MOD,   ONLY : READ_NC

    include  "CMN_SIZE"
    
    character(len=*), intent(in)     :: var_name
    character(len=*), intent(in)     :: produce
    real*8,           intent(inout)  :: data_con(:, :, :)
    INTEGER,intent(in)               :: FLAG

    !local var
    character(len=100)               :: file_path
    integer                          :: YEAR, MONTH, DAY

    IF ( FLAG == -1 ) THEN
      YEAR  = TIME_OBJ%PRE_YEAR
      MONTH = TIME_OBJ%PRE_MONTH
      DAY   = TIME_OBJ%PRE_DAY
    ELSE IF ( FLAG == 1 ) THEN
      YEAR  = TIME_OBJ%NEXT_YEAR
      MONTH = TIME_OBJ%NEXT_MONTH
      DAY   = TIME_OBJ%NEXT_DAY
    ELSE
      YEAR  = TIME_OBJ%YEAR
      MONTH = TIME_OBJ%MONTH
      DAY   = TIME_OBJ%DAY
    END IF

    WRITE(file_path, "('GEOSFP.', I4.4, I2.2, I2.2, '.')") &
    YEAR, MONTH, DAY

    file_path = trim(HR_MET_PATH) // trim(file_path) //    &
                trim(produce) // ".025x03125.nc"

    IF ( TRIM(var_name) == TRIM("PBLH") ) THEN
      call READ_NC(trim(file_path), trim(var_name), data_con, IIPAR, JJPAR, 24)
      RETURN
    END IF

    call READ_NC(trim(file_path), trim(var_name), data_con, IIPAR, JJPAR, 8)

  END SUBROUTINE

  SUBROUTINE READ_GEOS_FP_4D(data_con, var_name, produce, FLAG)

    USE NETCDF_IO_MOD,   ONLY : READ_NC

    include  "CMN_SIZE"

    character(len=*), intent(in)     :: var_name
    character(len=*), intent(in)     :: produce
    real*8,           intent(inout)  :: data_con(:, :, :, :)
    INTEGER,intent(in)               :: FLAG

    !local var
    character(len=100)               :: file_path
    integer                          :: YEAR, MONTH, DAY

    IF ( FLAG == -1 ) THEN
      YEAR  = TIME_OBJ%PRE_YEAR
      MONTH = TIME_OBJ%PRE_MONTH
      DAY   = TIME_OBJ%PRE_DAY
    ELSE IF ( FLAG == 1 ) THEN
      YEAR  = TIME_OBJ%NEXT_YEAR
      MONTH = TIME_OBJ%NEXT_MONTH
      DAY   = TIME_OBJ%NEXT_DAY
    ELSE
      YEAR  = TIME_OBJ%YEAR
      MONTH = TIME_OBJ%MONTH
      DAY   = TIME_OBJ%DAY
    END IF

    WRITE(file_path, "('GEOSFP.', I4.4, I2.2, I2.2, '.')") &
    YEAR, MONTH, DAY

    file_path = trim(HR_MET_PATH) // trim(file_path) //    &
                trim(produce) // ".025x03125.nc"

    IF ( TRIM("CMFMC") == TRIM(VAR_NAME) ) THEN
      call READ_NC(trim(file_path), trim(var_name), data_con, IIPAR, JJPAR, LLPAR+1, 8)
      RETURN
    END IF

    call READ_NC(trim(file_path), trim(var_name), data_con, IIPAR, JJPAR, LLPAR, 8)

  END SUBROUTINE

  SUBROUTINE INIT_AIR_STATE()

    !warn this subroutine should be call before init_time
    include  "CMN_SIZE"

    !local var
    INTEGER              :: as


    allocate(u_wind(IIPAR, JJPAR, LLPAR), stat=as)
    u_wind = 0d0

    allocate(v_wind(IIPAR, JJPAR, LLPAR), stat=as)
    v_wind = 0d0

    allocate(psc(IIPAR, JJPAR), stat=as)
    psc = 0d0

    allocate(org_u_wind1(IIPAR, JJPAR, LLPAR, 8), stat=as)
    org_u_wind1 = 0d0

    allocate(org_v_wind1(IIPAR, JJPAR, LLPAR, 8), stat=as)
    org_v_wind1 = 0d0

    allocate(org_u_wind2(IIPAR, JJPAR, LLPAR, 8), stat=as)
    org_u_wind2 = 0d0

    allocate(org_v_wind2(IIPAR, JJPAR, LLPAR, 8), stat=as)
    org_v_wind2 = 0d0

    allocate(org_p1(IIPAR, JJPAR, 8), stat=as)
    org_p1 = 0d0

    allocate(org_p2(IIPAR, JJPAR, 8), stat=as)
    org_p2 = 0d0

    allocate(P_TP1(IIPAR, JJPAR), stat=as)
    P_TP1 = 0d0

    allocate(P_TP2(IIPAR, JJPAR), stat=as)
    P_TP2 = 0d0

    allocate(P_TEMP(IIPAR, JJPAR), stat=as)
    P_TEMP = 0d0

    allocate(P_TEMP_NEW(IIPAR, JJPAR), stat=as)
    P_TEMP_NEW = 0d0

    allocate(AD_A(IIPAR, JJPAR, LLPAR), stat=as)
    AD_A = 0d0

    allocate(AD_B(IIPAR, JJPAR, LLPAR), stat=as)
    AD_B = 0d0

    allocate(UTMP(IIPAR, JJPAR, LLPAR), stat=as)
    UTMP = 0d0

    allocate(VTMP(IIPAR, JJPAR, LLPAR), stat=as)
    VTMP = 0d0

    allocate(XMASS(IIPAR, JJPAR, LLPAR), stat=as)
    XMASS = 0d0

    allocate(YMASS(IIPAR, JJPAR, LLPAR), stat=as)
    YMASS = 0d0

    allocate(MASS(IIPAR, JJPAR, LLPAR), stat=as)
    MASS = 0d0

    allocate(MASS_45(IIPAR, JJPAR, LLPAR), stat=as)
    MASS_45 = 0d0

    allocate(MONTH_DATA(IIPAR, JJPAR, LLPAR), stat=as)
    MONTH_DATA = 0d0

    allocate(TEM_DATA(IIPAR, JJPAR, LLPAR), stat=as)
    TEM_DATA = 0d0

    allocate(MASSFLEW(IIPAR, JJPAR, LLPAR, 1), stat=as)
    MASSFLEW = 0d0

    allocate(MASSFLNS(IIPAR, JJPAR, LLPAR, 1), stat=as)
    MASSFLNS = 0d0

    allocate(MASSFLUP(IIPAR, JJPAR, LLPAR, 1), stat=as)
    MASSFLUP = 0d0

    allocate(T(IIPAR, JJPAR, LLPAR), stat=as)
    T = 0d0

    allocate(org_t1(IIPAR, JJPAR, LLPAR, 8), stat=as)
    org_t1 = 0d0

    allocate(org_t2(IIPAR, JJPAR, LLPAR, 8), stat=as)
    org_t2 = 0d0

    allocate(AD(IIPAR, JJPAR, LLPAR), stat=as)
    AD = 0d0

    allocate(BXHEIGHT(IIPAR, JJPAR, LLPAR), stat=as)
    BXHEIGHT = 0d0

    allocate(ORG_PBL(IIPAR, JJPAR, 24), stat=as)
    ORG_PBL = 0d0

    allocate(AIRVOL(IIPAR, JJPAR, LLPAR), stat=as)
    AIRVOL = 0d0

    allocate(ORG_CMFMC(IIPAR, JJPAR, LLPAR+1, 8), stat=as)
    ORG_CMFMC = 0d0

    allocate(ORG_DTRAIN(IIPAR, JJPAR, LLPAR+1, 8), stat=as)
    ORG_DTRAIN = 0d0

    first = .TRUE.

    CALL READ_GEOS_FP( ORG_PBL, 'PBLH',  "A1", 0 )
    ! READ TODAY MET
    IF ( TIME_OBJ%IS_READ_LEFT ) THEN
      CALL READ_GEOS_FP( org_u_wind2, 'U',  "A3dyn", 0 )
      CALL READ_GEOS_FP( org_v_wind2, 'V',  "A3dyn", 0 )
      CALL READ_GEOS_FP( org_p2, 'PS', "I3"   , 0 )
      CALL READ_GEOS_FP( org_t2, 'T',  "I3"   , 0 )
      CALL READ_GEOS_FP( org_u_wind1, 'U',  "A3dyn", -1 )
      CALL READ_GEOS_FP( org_v_wind1, 'V',  "A3dyn", -1 )
      CALL READ_GEOS_FP( org_p1, 'PS', "I3"   , -1 )
      CALL READ_GEOS_FP( org_t1, 'T',  "I3"   , -1 )
    ELSE IF ( TIME_OBJ%IS_READ_RIGHT ) THEN
      CALL READ_GEOS_FP( org_u_wind1, 'U',  "A3dyn", 0 )
      CALL READ_GEOS_FP( org_v_wind1, 'V',  "A3dyn", 0 )
      CALL READ_GEOS_FP( org_p1, 'PS', "I3"   , 0 )
      CALL READ_GEOS_FP( org_t1, 'T',  "I3"   , 0 )
      CALL READ_GEOS_FP( org_u_wind2, 'U',  "A3dyn", 1 )
      CALL READ_GEOS_FP( org_v_wind2, 'V',  "A3dyn", 1 )
      CALL READ_GEOS_FP( org_p2, 'PS', "I3"   , 1 )
      CALL READ_GEOS_FP( org_t2, 'T',  "I3"   , 1 )
    ELSE 
      IF ( RUN_DIR ) THEN
        CALL READ_GEOS_FP( org_u_wind2, 'U',  "A3dyn", 0 )
        CALL READ_GEOS_FP( org_v_wind2, 'V',  "A3dyn", 0 )
        CALL READ_GEOS_FP( org_p2, 'PS', "I3"   , 0 )
        CALL READ_GEOS_FP( org_t2, 'T',  "I3"   , 0 )
      ELSE
        CALL READ_GEOS_FP( org_u_wind1, 'U',  "A3dyn", 0 )
        CALL READ_GEOS_FP( org_v_wind1, 'V',  "A3dyn", 0 )
        CALL READ_GEOS_FP( org_p1, 'PS', "I3"   , 0 )
        CALL READ_GEOS_FP( org_t1, 'T',  "I3"   , 0 )
      ENDIF
    END IF

    IF ( TIME_OBJ%IS_READ_COULD ) THEN
      CALL READ_GEOS_FP( ORG_DTRAIN, 'DTRAIN',   "A3dyn", -1 )
      CALL READ_GEOS_FP( ORG_CMFMC ,  'CMFMC',  "A3mstE", -1 )
      TIME_OBJ%IS_READ_COULD = .FALSE.
    ELSE
      CALL READ_GEOS_FP( ORG_DTRAIN, 'DTRAIN',   "A3dyn", 0 )
      CALL READ_GEOS_FP( ORG_CMFMC ,  'CMFMC',  "A3mstE", 0 )
    END IF
    CALL RENEW_AIR_STATE()
    CALL INIT_XG_VAR()
  END SUBROUTINE

  SUBROUTINE RENEW_AIR_STATE()

    USE PRESSURE_MOD, ONLY : SET_FLOATING_PRESSURE

    IF ( RUN_DIR ) THEN

      IF ( TIME_OBJ%IS_NEW_DAY ) THEN
        CALL READ_GEOS_FP( ORG_PBL, 'PBLH',  "A1", 0 )
      END IF

      IF ( TIME_OBJ%IS_READ_COULD ) THEN
        CALL READ_GEOS_FP( ORG_DTRAIN, 'DTRAIN',   "A3dyn", 0 )
        CALL READ_GEOS_FP( ORG_CMFMC ,  'CMFMC',  "A3mstE", 0 )
      END IF

      IF ( TIME_OBJ%IS_READ_RIGHT ) THEN
        org_u_wind1 = org_u_wind2
        org_u_wind1 = org_u_wind2
        org_p1      = org_p2
        org_t1      = org_t2
        CALL READ_GEOS_FP( org_u_wind2, 'U',  "A3dyn", 1 )
        CALL READ_GEOS_FP( org_v_wind2, 'V',  "A3dyn", 1 )
        CALL READ_GEOS_FP( org_p2, 'PS', "I3"   , 1 )
        CALL READ_GEOS_FP( org_t2, 'T',  "I3"   , 1 )
      END IF
      IF ( TIME_OBJ%LEFT_INDEX < TIME_OBJ%RIGHT_INDEX) THEN
        u_wind1 => org_u_wind2(:,:,:,TIME_OBJ%LEFT_INDEX)
        v_wind1 => org_v_wind2(:,:,:,TIME_OBJ%LEFT_INDEX)
        psc1    => org_p2     (:,:,  TIME_OBJ%LEFT_INDEX)
        t1      => org_t2     (:,:,:,TIME_OBJ%LEFT_INDEX)
      ELSE
        u_wind1 => org_u_wind1(:,:,:,TIME_OBJ%LEFT_INDEX)
        v_wind1 => org_v_wind1(:,:,:,TIME_OBJ%LEFT_INDEX)
        psc1    => org_p1     (:,:,  TIME_OBJ%LEFT_INDEX)
        t1      => org_t1     (:,:,:,TIME_OBJ%LEFT_INDEX)
      END IF
      u_wind2 => org_u_wind2(:,:,:,TIME_OBJ%RIGHT_INDEX)
      v_wind2 => org_v_wind2(:,:,:,TIME_OBJ%RIGHT_INDEX)
      psc2    => org_p2     (:,:,  TIME_OBJ%RIGHT_INDEX)
      t2      => org_t2     (:,:,:,TIME_OBJ%RIGHT_INDEX)
      PBL     => ORG_PBL    (:,:,  TIME_OBJ%HOUR + 1)
    ELSE

      IF ( TIME_OBJ%IS_NEW_DAY ) THEN
        CALL READ_GEOS_FP( ORG_PBL, 'PBLH',  "A1", 0 )
      END IF

      IF ( TIME_OBJ%IS_READ_COULD ) THEN
        CALL READ_GEOS_FP( ORG_DTRAIN, 'DTRAIN',   "A3dyn", -1 )
        CALL READ_GEOS_FP( ORG_CMFMC ,  'CMFMC',  "A3mstE", -1 )
      END IF

      IF ( TIME_OBJ%IS_READ_LEFT ) THEN
        org_u_wind2 = org_u_wind1
        org_u_wind2 = org_u_wind1
        org_p2      = org_p1
        org_t2      = org_t1
        CALL READ_GEOS_FP( org_u_wind1, 'U',  "A3dyn", -1 )
        CALL READ_GEOS_FP( org_v_wind1, 'V',  "A3dyn", -1 )
        CALL READ_GEOS_FP( org_p1, 'PS', "I3"   , -1 )
        CALL READ_GEOS_FP( org_t1, 'T',  "I3"   , -1 )
      END IF

      IF ( TIME_OBJ%LEFT_INDEX < TIME_OBJ%RIGHT_INDEX) THEN
        u_wind2 => org_u_wind1(:,:,:,TIME_OBJ%RIGHT_INDEX)
        v_wind2 => org_v_wind1(:,:,:,TIME_OBJ%RIGHT_INDEX)
        psc2    => org_p1     (:,:,  TIME_OBJ%RIGHT_INDEX)
        t2      => org_t1     (:,:,:,TIME_OBJ%RIGHT_INDEX)
      ELSE
        u_wind2 => org_u_wind2(:,:,:,TIME_OBJ%RIGHT_INDEX)
        v_wind2 => org_v_wind2(:,:,:,TIME_OBJ%RIGHT_INDEX)
        psc2    => org_p2     (:,:,  TIME_OBJ%RIGHT_INDEX)
        t2      => org_t2     (:,:,:,TIME_OBJ%RIGHT_INDEX)
      END IF

      u_wind1 => org_u_wind1(:,:,:,TIME_OBJ%LEFT_INDEX)
      v_wind1 => org_v_wind1(:,:,:,TIME_OBJ%LEFT_INDEX)
      psc1    => org_p1     (:,:,  TIME_OBJ%LEFT_INDEX)
      t1      => org_t1     (:,:,:,TIME_OBJ%LEFT_INDEX)
      PBL     => ORG_PBL    (:,:,  TIME_OBJ%HOUR + 1)
    END IF

    DTRAIN => ORG_DTRAIN(:, :, :, TIME_OBJ%LEFT_INDEX)
    CMFMC  => ORG_CMFMC (:, :, :, TIME_OBJ%LEFT_INDEX)
    CALL INTERP_AIR_MET()
    IF ( FIRST ) THEN
      CALL SET_FLOATING_PRESSURE( PSC )
      FIRST = .FALSE.
    END IF
    CALL CAL_AIR_STATE()
    if ( .not. XG_FIRST) CALL RENEW_XG_VAR()

  END SUBROUTINE

  SUBROUTINE CAL_AIR_STATE()

    USE GRID_MOD,     ONLY : GET_AREA_M2
    USE PRESSURE_MOD, ONLY : GET_BP, GET_PEDGE

    INCLUDE "CMN_SIZE"
    INCLUDE "CMN_GCTM"

    !LOCAL VAR
    INTEGER              :: I, J, L
    REAL*8               :: P1, P2, D_VAL, AREA_M2

    !$OMP PARALLEL DO               &
    !$OMP DEFAULT( SHARED )         &
    !$OMP PRIVATE( I, J, L, AREA_M2, P1, P2, D_VAL )
      DO J = 1, JJPAR
        AREA_M2 = GET_AREA_M2( J )
        DO I = 1, IIPAR
          DO L = 1, LLPAR
            P1    = GET_PEDGE(I,J,L)
            P2    = GET_PEDGE(I,J,L+1)
            D_VAL = P1 - P2
            AD(I, J, L) = D_VAL * G0_100 * AREA_M2
            BXHEIGHT(I,J,L) = Rdg0 * T(I,J,L) * LOG( P1 / P2 )
            AIRVOL(I,J,L) = BXHEIGHT(I,J,L) * AREA_M2
          END DO
        END DO
      END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE

  SUBROUTINE INTERP_AIR_MET()

    USE ALGORITHM_MOD,        ONLY : INTERP_OBJ

    INCLUDE "CMN_SIZE"

    CALL INTERP_OBJ%D3(u_wind1, u_wind2, u_wind, IIPAR, JJPAR, LLPAR)
    CALL INTERP_OBJ%D3(v_wind1, v_wind2, v_wind, IIPAR, JJPAR, LLPAR)
    CALL INTERP_OBJ%D2(psc1,    psc2,    psc,    IIPAR, JJPAR)
    CALL INTERP_OBJ%D3(T1,      T2,      T,      IIPAR, JJPAR, LLPAR)

  END SUBROUTINE
END MODULE

! $Id: grid_mod.f,v 1.1 2011/02/23 00:08:48 daven Exp $
      MODULE GRID_MOD
!
!******************************************************************************
!  Module GRID_MOD contains variables and routines which are used to specify
!  the parameters of a GEOS-CHEM horizontal grid. (bmy, 3/11/03, 4/20/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) IS_NESTED  (LOGICAL) : =T if we are using a nested-grid
!  (2 ) IIGLOB     (INTEGER) : Global longitude extent [# of boxes]
!  (3 ) JJGLOB     (INTEGER) : Global latitude  extent [# of boxes]
!  (4 ) I0         (INTEGER) : Nested-grid offset in longitude (X) dimension
!  (5 ) J0         (INTEGER) : Nested-grid offset in latitude  (Y) dimension
!  (6 ) XMID_G     (REAL*8 ) : GLOBAL array of grid-box lon centers [degrees]
!  (7 ) XEDGE_G    (REAL*8 ) : GLOBAL array of grid-box lon edges   [degrees]
!  (8 ) YMID_G     (REAL*8 ) : GLOBAL array of grid-box lat centers [degrees]
!  (9 ) YEDGE_G    (REAL*8 ) : GLOBAL array of grid-box lat edges   [degrees]
!  (10) YMID_R_G   (REAL*8 ) : GLOBAL array of grid-box lat centers [radians]
!  (11) YEDGE_R_G  (REAL*8 ) : GLOBAL array of grid-box lat edges   [radians]
!  (12) AREA_M2_G  (REAL*8 ) : GLOBAL array of grid-box surface areas [m2]
!  (13) AREA_CM2_G (REAL*8 ) : GLOBAL array of grid-box surface areas [cm2]
!  (14) XMID       (REAL*8 ) : WINDOW array of grid-box lon centers [degrees]
!  (15) XEDGE      (REAL*8 ) : WINDOW array of grid-box lon edges   [degrees]
!  (16) YMID       (REAL*8 ) : WINDOW array of grid-box lat centers [degrees]
!  (17) YEDGE      (REAL*8 ) : WINDOW array of grid-box lat edges   [degrees]
!  (18) YMID_R     (REAL*8 ) : WINDOW array of grid-box lat centers [radians]
!  (19) YEDGE_R    (REAL*8 ) : WINDOW array of grid-box lat edges   [radians]
!  (20) AREA_M2    (REAL*8 ) : WINDOW array of grid-box surface areas [m2]
!  (21) AREA_CM2   (REAL*8 ) : WINDOW array of grid-box surface areas [cm2]
!
!  Module Routines:
!  ============================================================================
!  (1 ) COMPUTE_GRID      : Computes all lon, lat, surface area quantities
!  (2 ) SET_XOFFSET       : Initializes nested grid longitude (X) offset
!  (3 ) SET_YOFFSET       : Initializes nested grid latitude  (Y) offset
!  (4 ) GET_XOFFSET       : Returns     nested grid longitude (X) offset
!  (5 ) GET_XOFFSET       : Returns     nested grid latitude  (Y) offset
!  (6 ) GET_XMID          : Returns grid box center  longitude [degrees]
!  (7 ) GET_XEDGE         : Returns grid box W. edge longitude [degrees]
!  (8 ) GET_YMID          : Returns grid box center  latitude  [degrees]
!  (9 ) GET_YEDGE         : Returns grid box S. edge latitude  [degrees]
!  (10) GET_YMID_R        : Returns grid box center  latitude  [radians]
!  (11) GET_YEDGE_R       : Returns grid box S. edge latitude  [radians]
!  (12) GET_AREA_M2       : Returns grid box surface area      [m2]
!  (13) GET_AREA_CM2      : Returns grid box surface area      [cm2]
!  (14) GET_BOUNDING_BOX  : Returns indices for a region denoted by lats/lons
!  (15) ITS_A_NESTED_GRID : Returns T for nested grid simulations; F otherwise
!  (16) INIT_GRID         : Allocates and zeroes all module arrays
!  (17) CLEANUP_GRID      : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by grid_mod.f:
!  ============================================================================
!  (1 ) error_mod.f       : Module containing I/O error and NaN check routines
!
!  NOTES:
!  (1 ) Fixed typos in "grid_mod.f" (bmy, 4/28/03)
!  (2 ) Added routine GET_BOUNDING_BOX.  Now define 1x125 grid. (bmy, 12/1/04)
!  (3 ) Modified for GCAP 4x5 horizontal grid (swu, bmy, 5/24/05)
!  (4 ) Added comments re: surface area derivation (bmy, 4/20/06)
!  (5 ) Modifications for GEOS-5 nested grids (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "grid_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC ::  COMPUTE_GRID
      PUBLIC ::  SET_XOFFSET
      PUBLIC ::  SET_YOFFSET
      PUBLIC ::  GET_XOFFSET
      PUBLIC ::  GET_YOFFSET
      PUBLIC ::  GET_XMID
      PUBLIC ::  GET_XEDGE
      PUBLIC ::  GET_XEDGE_G
      PUBLIC ::  GET_YMID
      PUBLIC ::  GET_YEDGE
      PUBLIC ::  GET_YSIN
      PUBLIC ::  GET_YEDGE_G
      PUBLIC ::  GET_YMID_R
      PUBLIC ::  GET_YMID_R_W
      PUBLIC ::  GET_YEDGE_R
      PUBLIC ::  GET_AREA_M2
      PUBLIC ::  GET_AREA_CM2
      PUBLIC ::  GET_BOUNDING_BOX
      ! (dkh, 10/07/09)
      PUBLIC ::  GET_IJ
      PUBLIC ::  GET_IJ_GLOBAL
      PUBLIC ::  ITS_A_NESTED_GRID
      PUBLIC ::  CLEANUP_GRID
      PUBLIC ::  GET_IIIPAR
      PUBLIC ::  GET_JJJPAR
      PUBLIC ::  DX_COV, DY_COV
      PUBLIC ::  YMID, XMID

      !==================================================================
      ! MODULE VARIABLES
      !==================================================================
      LOGICAL              :: IS_NESTED
      INTEGER              :: I0,           J0
      INTEGER              :: IIGLOB,       JJGLOB
      REAL*8,  ALLOCATABLE :: XMID_G(:),    XEDGE_G(:)
      REAL*8,  ALLOCATABLE :: YMID_G(:),    YEDGE_G(:)
      REAL*8,  ALLOCATABLE :: YMID_R_G(:),  YEDGE_R_G(:)
      REAL*8,  ALLOCATABLE :: AREA_M2_G(:), AREA_CM2_G(:)
      REAL*8,  ALLOCATABLE :: YSIN(:)
      REAL*8,  ALLOCATABLE :: XMID(:),      XEDGE(:)
      REAL*8,  ALLOCATABLE :: YMID(:),      YEDGE(:)
      REAL*8,  ALLOCATABLE :: YMID_R(:),    YEDGE_R(:)
      REAL*8,  ALLOCATABLE :: YMID_R_W(:),  YEDGE_R_W(:)
      REAL*8,  ALLOCATABLE :: AREA_M2(:),   AREA_CM2(:)
      REAL*8,  ALLOCATABLE :: DX_COV(:),   DY_COV(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_GRID
!
!******************************************************************************
!  Subroutine COMPUTE_GRID initializes the longitude, latitude and surface
!  area arrays. (bmy, 3/11/03, 4/20/06)
!
!  NOTES:
!  (1 ) Added fancy output (bmy, 4/26/04)
!  (2 ) Suppress some output lines (bmy, 7/20/04)
!  (3 ) Now also support 1 x 1.25 grid (bmy, 12/1/04)
!  (4 ) Now modified for GCAP 4x5 horizontal grid (swu, bmy, 5/24/05)
!  (5 ) Added comments re: surface area derivation (bmy, 4/20/06)
!  (6 ) Compute YMID, YEDGE for 0.5x0.666 nested grids (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Physical constants
#     include "define.h"

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I,    J
      REAL*8              :: FMID, FEDGE

      !=================================================================
      ! COMPUTE_GRID begins here!
      !=================================================================

      ! Allocate variables on first call
      IF ( FIRST ) THEN
         CALL INIT_GRID
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Compute latitude centers & edges (algorithm from old "input.f")
      !=================================================================
      FMID  = 0.5d0 * DBLE( JJGLOB + 1 )
      FEDGE = 0.5d0 * DBLE( JJGLOB + 2 )

      DO J = 1, JJGLOB
         YMID_G(J)  = DJSIZE * ( DBLE(J) - FMID  )
         YEDGE_G(J) = DJSIZE * ( DBLE(J) - FEDGE )
      ENDDO

#if   defined( GRID4x5 ) && defined( GCAP )
      ! Overwrite YMID at poles for GCAP 4 x 5 grid (swu, bmy, 5/24/05)
      YMID_G(1)      = -88.d0
      YMID_G(JJGLOB) = +88.d0

#elif defined( GRID4x5 )
      ! Overwrite YMID at poles for 4 x 5 grid
      YMID_G(1)      = -89.d0
      YMID_G(JJGLOB) = +89.d0

#elif defined( GRID2x25 )
      ! Overwrite YMID at poles for 2 x 2.5 grid
      YMID_G(1)      = -89.5d0
      YMID_G(JJGLOB) = +89.5d0

#elif defined ( GRID1x1 ) || defined( GRID1x125 )
      ! Overwrite YMID at poles for 1 x 1 and 1 x 1.25 grids
      YMID_G(1)      = -89.75d0
      YMID_G(JJGLOB) = +89.75d0

#elif defined ( GRID05x0666 )
      ! Overwrite YMID at poles for 0.5 x 0.666 grids
      YMID_G(1)      = -89.875d0
      YMID_G(JJGLOB) = +89.875d0

#elif defined ( GRID025x03125 )
      ! Overwrite YMID at poles for 0.25 x 0.3125 grids
      ! (lzh,02/01/2015)
      YMID_G(1)      = -89.9375d0
      YMID_G(JJGLOB) = +89.9375d0

#endif

      ! Overwrite YEDGE at poles
      YEDGE_G(1)        = -90d0
      YEDGE_G(JJGLOB+1) = +90d0

      ! Compute latitude center/edges in radians
      DO J = 1, JJGLOB
         YMID_R_G(J)  = ( PI / 180d0 ) * YMID_G(J)
         YEDGE_R_G(J) = ( PI / 180d0 ) * YEDGE_G(J)
      ENDDO
      ! Overwrite RLATV at N. pole
      YEDGE_R_G(JJGLOB+1) = PI / 2d0

      !=================================================================
      ! Compute longitude centers & edges (algorithm from old "input.f")
      !=================================================================
      XMID_G(1)  = -180d0
      XEDGE_G(1) = XMID_G(1) - ( DISIZE / 2d0 )

      DO I = 1, IIGLOB-1
         XMID_G(I+1)  = XMID_G(I)  + DISIZE
      ENDDO

      DO I = 1, IIGLOB
         XEDGE_G(I+1) = XEDGE_G(I) + DISIZE
      ENDDO
      !=================================================================
      ! Compute grid box surface areas (algorithm from old "input.f")
      !
      ! The surface area of a grid box is derived as follows:
      !
      !    Area = dx * dy
      !
      ! Where:
      !
      !    dx is the arc length of the box in longitude
      !    dy is the arc length of the box in latitude
      !
      ! Which are computed as:
      !
      !    dx = r * delta-longitude
       !       = ( Re * cos[ YMID[J] ] ) * ( 2 * PI / IIIPAR )
      !
      !    dy = r * delta-latitude
      !       = Re * ( YEDGE[J+1] - YEDGE[J] )
      !
      ! Where:
      !
      !    Re         is the radius of the earth
      !    YMID[J]    is the latitude at the center of box J
      !    YEDGE[J+1] is the latitude at the N. Edge of box J
      !    YEDGE[J]   is the latitude at the S. Edge of box J
      !
      ! So, the surface area is thus:
      !
       !    Area = ( Re * cos( YMID[J] ) * ( 2 * PI / IIIPAR ) *
      !             Re * ( YEDGE[J+1] - YEDGE[J] )
      !
      !    2*PI*Re^2    {                                            }
      ! = ----------- * { cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] ) }
       !     IIIPAR      {                                            }
      !
      ! And, by using the trigonometric identity:
      !
      !    d sin(x) = cos x * dx
      !
      ! The following term:
      !
      !    cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] )
      !
      ! May also be written as a difference of sines:
      !
      !    sin( YEDGE[J+1] ) - sin( YEDGE[J] )
      !
      ! So the final formula for surface area of a grid box is:
      !
      !            2*PI*Re^2    {                                     }
      !    Area = ----------- * { sin( YEDGE[J+1] ) - sin( YEDGE[J] ) }
       !              IIIPAR     {                                     }
      !
      !
      ! NOTES:
      ! (1) The formula with sines is more numerically stable, and will
      !      yield identical global total surface areas for all grids.
      ! (2) The units are determined by the radius of the earth Re.
      !      if you use Re [m], then surface area will be in [m2], or
      !      if you use Re [cm], then surface area will be in [cm2], etc.
      ! (3) The grid box surface areas only depend on latitude, as they
      !      are symmetric in longitude.  To compute the global surface
      !      area, multiply the surface area arrays below by the number
       !      of longitudes (e.g. IIIPAR).
      ! (4) At present, assumes that GEOS-Chem will work on a
      !      Cartesian grid.
      !
      ! (bmy, 4/20/06, 2/24/12)
      !=================================================================
      DO J = 1, JJGLOB

         ! Grid box surface areas [m2]
         AREA_M2_G(J) = 2d0 * PI * Re * Re / DBLE( IIGLOB ) *
     &                  ( SIN( YEDGE_R_G(J+1) ) - SIN( YEDGE_R_G(J) ) )

         ! Grid box surface areas [cm2]
         AREA_CM2_G(J) = AREA_M2_G(J) * 1d4

      ENDDO
      !=================================================================
      ! Save to local size arrays so that we can index for all grids
      !=================================================================
      ! XMID
      DO I = 1, IIPAR
         XMID(I) = XMID_G(I+I0)
      ENDDO
      ! XEDGE
      DO I = 1, IIPAR+1
         XEDGE(I) = XEDGE_G(I+I0)
      ENDDO
      ! YMID, YMID_R, AREA_M2, AREA_CM2
      DO J = 1, JJPAR
         YMID(J)     = YMID_G(J+J0)
         YMID_R(J)   = YMID_R_G(J+J0)
         AREA_M2(J)  = AREA_M2_G(J+J0)
         AREA_CM2(J) = AREA_CM2_G(J+J0)
      ENDDO
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )
      ! This only needs to be done for GEOS-5 nested grid (dan, bmy, 11/18/08)
      DO J = 0, JJPAR+1
         YMID_R_W(J) = YMID_R_G(J+J0)
      ENDDO
#endif

      ! YEDGE, YEDGE_R
      DO J = 1, JJPAR+1
         YEDGE(J)   = YEDGE_G(J+J0)
         YEDGE_R(J) = YEDGE_R_G(J+J0)
         YSIN(J) = SIN ( YEDGE_R_G(J+J0) )
      ENDDO

      !=================================================================
      !  Added By : Asandu, KSingh, Mjardak   04/01/09
      !  Desc: Grid box (I,J) has size DX_COV[J] x DY_COV[J] in Km
      !        (same for all I = 1, IIGLOB)
      !=================================================================
      !replace JJGLOB by JJPAR (fp)
      DO J = 1, JJPAR
         DX_COV(J)    =( Re * cos( YMID(J)/180.d0*PI ) ) *
     &        ( 2 * PI / IIGLOB ) * 1.d-3
      END DO

      DO J = 1, JJPAR
          DY_COV(J)    = Re*(YEDGE(J+1)-YEDGE(J) )/180.d0*PI * 1.d-3
      END DO

      !=================================================================
      ! Echo info to stdout
      !=================================================================
      WRITE( 6, '(''Nested-Grid X-offset [boxes]:'', i4 )' ) I0
      WRITE( 6, '(''Nested-Grid Y-offset [boxes]:'', i4 )' ) J0
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( XMID(I),  I=1,IIPAR )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( XEDGE(I), I=1,IIPAR+1 )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( YMID(J),  J=1,JJPAR )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( YEDGE(J), J=1,JJPAR+1 )

      !=================================================================
      ! Deallocate global arrays -- we don't need these anymore
      !=================================================================
      IF ( ALLOCATED( XMID_G     ) ) DEALLOCATE( XMID_G     )
!      IF ( ALLOCATED( XEDGE_G    ) ) DEALLOCATE( XEDGE_G    )
      IF ( ALLOCATED( YMID_G     ) ) DEALLOCATE( YMID_G     )
!      IF ( ALLOCATED( YEDGE_G    ) ) DEALLOCATE( YEDGE_G    )
      IF ( ALLOCATED( YMID_R_G   ) ) DEALLOCATE( YMID_R_G   )
      IF ( ALLOCATED( YEDGE_R_G  ) ) DEALLOCATE( YEDGE_R_G  )
      IF ( ALLOCATED( AREA_M2_G  ) ) DEALLOCATE( AREA_M2_G  )
      IF ( ALLOCATED( AREA_CM2_G ) ) DEALLOCATE( AREA_CM2_G )

      END SUBROUTINE COMPUTE_GRID

!------------------------------------------------------------------------------

      SUBROUTINE SET_XOFFSET( X_OFFSET )
!
!******************************************************************************
!  Function SET_XOFFSET initializes the nested-grid latitude offset
!  variable I0. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) X_OFFSET (INTEGER) : Nested grid longitude offset (# of grid boxes)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: X_OFFSET

      !=================================================================
      ! SET_XOFFSET begins here!
      !=================================================================
      I0 = X_OFFSET

      END SUBROUTINE SET_XOFFSET

!------------------------------------------------------------------------------

      SUBROUTINE SET_YOFFSET( Y_OFFSET )
!
!******************************************************************************
!  Function SET_YOFFSET initializes the nested-grid latitude offset
!  variable J0. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) Y_OFFSET (INTEGER) : Nested grid latitude offset (# of grid boxes)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: Y_OFFSET

      !=================================================================
      ! SET_XOFFSET begins here!
      !=================================================================
      J0 = Y_OFFSET

      END SUBROUTINE SET_YOFFSET

!------------------------------------------------------------------------------
!          Prasad Kasibhatla - Duke University                                !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_iiipar
!
! !DESCRIPTION: Function GET\_IIIPAR returns IIIPAR.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_IIIPAR() RESULT( I )
!
! !RETURN VALUE:
!
      INTEGER              :: I  ! IIIPAR
!
! !REVISION HISTORY:
!  01 Mar 2012 - P. Kasibhatla
!EOP
!------------------------------------------------------------------------------
!BOC
      I =  IIGLOB

      END FUNCTION GET_IIIPAR
!EOC
!------------------------------------------------------------------------------
!          Prasad Kasibhatla - Duke University                                !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_jjjpar
!
! !DESCRIPTION: Function GET\_JJJPAR returns JJJPAR.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_JJJPAR() RESULT( J )
!
! !RETURN VALUE:
!
      INTEGER              :: J  ! JJJPAR
!
! !REVISION HISTORY:
!  01 Mar 2012 - P. Kasibhatla
!EOP
!------------------------------------------------------------------------------
!BOC
      J =  JJGLOB

      END FUNCTION GET_JJJPAR
!EOC
!------------------------------------------------------------------------------

      FUNCTION GET_XOFFSET( GLOBAL ) RESULT( X_OFFSET )
!
!******************************************************************************
!  Function GET_XOFFSET returns the nested-grid longitude offset to the
!  calling program. (bmy, 3/11/03)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: GLOBAL

      ! Function value
      INTEGER :: X_OFFSET

      !=================================================================
      ! GET_XOFFSET begins here!
      !=================================================================
      IF ( PRESENT( GLOBAL ) ) THEN

         ! If GLOBAL is passed, then return the actual window offset.
         ! This is necessary for certain instances (e.g. diagnostics)
         X_OFFSET = I0

      ELSE

         ! Otherwise, if we have a nested grid, then all of the met
         ! fields have been cut down to size already.  Return 0.
         X_OFFSET = 0

      ENDIF

      END FUNCTION GET_XOFFSET

!------------------------------------------------------------------------------

      FUNCTION GET_YOFFSET( GLOBAL ) RESULT( Y_OFFSET )
!
!******************************************************************************
!  Function GET_YOFFSET returns the nested-grid latitude offset to the
!  calling program. (bmy, 3/11/03)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: GLOBAL

      ! Function value
      INTEGER :: Y_OFFSET

      !=================================================================
      ! GET_XOFFSET begins here!
      !=================================================================
      IF ( PRESENT( GLOBAL ) ) THEN

         ! If GLOBAL is passed, then return the actual window offset.
         ! This is necessary for certain instances (e.g. diagnostics)
         Y_OFFSET = J0

      ELSE

         ! Otherwise, if we have a nested grid, then all of the met
         ! fields have been cut down to size already.  Return 0.
         Y_OFFSET = 0

      ENDIF

      END FUNCTION GET_YOFFSET

!------------------------------------------------------------------------------

      FUNCTION GET_XMID( I ) RESULT( X )
!
!******************************************************************************
!  Function GET_XMID returns the longitude in degrees at the center of a
!  GEOS-CHEM grid box.  Works for nested-grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) I (INTEGER) : GEOS-CHEM grid-box index for longitude
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I

      ! Function value
      REAL*8              :: X

      !=================================================================
      ! GET_XMID begins here!
      !=================================================================
      X = XMID( I )

      END FUNCTION GET_XMID

!------------------------------------------------------------------------------

      FUNCTION GET_XEDGE( I ) RESULT( X )
!
!******************************************************************************
!  Function GET_XEDGE returns the longitude in degrees at the western edge of
!  a GEOS-CHEM grid box.  Works for nested-grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) I (INTEGER) : GEOS-CHEM grid-box index for longitude
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I

      ! Function value
      REAL*8              :: X

      !=================================================================
      ! GET_XEDGE begins here!
      !=================================================================
      X = XEDGE( I )

      END FUNCTION GET_XEDGE

!------------------------------------------------------------------------------
!          Prasad Kasibhatla - Duke University                                !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xedge_g
!
! !DESCRIPTION: Function GET\_XEDGE\_G returns the longitude in degrees at the
!  western edge of a GEOS-Chem grid box
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_XEDGE_G( I ) RESULT( X )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I  ! Longitude index
!
! !RETURN VALUE:
!
      REAL*4              :: X  ! Corresponding lon value @ W edge of grid box
!
! !REVISION HISTORY:
!  01 Mar 2012 - P. Kasibhatla
!EOP
!------------------------------------------------------------------------------
!BOC
      X = XEDGE_G( I )

      END FUNCTION GET_XEDGE_G
!EOC
!------------------------------------------------------------------------------

      FUNCTION GET_YMID( J ) RESULT( Y )
!
!******************************************************************************
!  Function GET_YMID returns the latitude in degrees at the center of
!  a GEOS-CHEM grid box.  Works for nested grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) I (INTEGER) : GEOS-CHEM grid-box index for latitude
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: J

      ! Function value
      REAL*8              :: Y

      !=================================================================
      ! GET_YMID begins here!
      !=================================================================
      Y = YMID( J )

      END FUNCTION GET_YMID

!------------------------------------------------------------------------------

      FUNCTION GET_YEDGE( J ) RESULT( Y )
!
!******************************************************************************
!  Function GET_EDGE returns the latitude in degrees at the southern edge of
!  a GEOS-CHEM grid box.  Works for nested grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) J (INTEGER) : GEOS-CHEM grid-box index for latitude
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: J

      ! Function value
      REAL*8              :: Y

      !=================================================================
      ! GET_YEDGE begins here!
      !=================================================================
      Y = YEDGE( J )

      END FUNCTION GET_YEDGE

!------------------------------------------------------------------------------
!          Prasad Kasibhatla - Duke University                                !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_yedge_g
!
! !DESCRIPTION: Function GET\_YEDGE\_G returns the latitude in degrees at the
!  southern edge of a GEOS-Chem grid box.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YEDGE_G( J ) RESULT( Y )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*4              :: Y  ! Latitude value @ S edge of grid box [degrees]
!
! !REVISION HISTORY:
!  01 Mar 2012 - P. Kasibhatla
!EOP
!------------------------------------------------------------------------------
!BOC
      Y = YEDGE_G( J )

      END FUNCTION GET_YEDGE_G
!EOC
!------------------------------------------------------------------------------

      FUNCTION GET_YMID_R( J ) RESULT( Y )
!
!******************************************************************************
!  Function GET_YMID_R returns the latitude in radians at the center of
!  a GEOS-CHEM grid box.  Works for nested grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) J (INTEGER) : GEOS-CHEM grid-box index for latitude
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: J

      ! Function value
      REAL*8              :: Y

      !=================================================================
      ! GET_YMID_R begins here!
      !=================================================================
      Y = YMID_R( J )

      END FUNCTION GET_YMID_R

!------------------------------------------------------------------------------

      FUNCTION GET_YMID_R_W( J ) RESULT( Y )
!
!******************************************************************************
!  Function GET_YMID_R_W returns the latitude in radians at the center of
!  a GEOS-CHEM grid box for the GEOS-5 nested grid. (dan, bmy, 11/6/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1) J (INTEGER) : GEOS-CHEM grid-box index for latitude
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: J

      ! Function value
      REAL*8              :: Y

      !=================================================================
      ! dan for tpcore window GET_YMID_R_W begins here!
      !=================================================================
      Y = YMID_R_W( J )

      ! Return to calling program
      END FUNCTION GET_YMID_R_W

!------------------------------------------------------------------------------

      FUNCTION GET_YEDGE_R( J ) RESULT( Y )
!
!******************************************************************************
!  Function GET_YEDGE_R returns the latitude in radians at the southern edge
!  of a GEOS-CHEM grid box.  Works for nested grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) J (INTEGER) : GEOS-CHEM grid-box index for latitude
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: J

      ! Function value
      REAL*8              :: Y

      !=================================================================
      ! GET_YEDGE_R begins here!
      !=================================================================
      Y = YEDGE_R( J )

      END FUNCTION GET_YEDGE_R

!------------------------------------------------------------------------------
      FUNCTION Get_ySin( I, J, L ) RESULT( Y )


      INTEGER, INTENT(IN) :: I  ! Longitude index
      INTEGER, INTENT(IN) :: J  ! Latitude index
      INTEGER, INTENT(IN) :: L  ! Level index


      REAL*8              :: Y  ! Sine of Latitude value @ S edge of grid box

      Y = YSIN(J)

      END FUNCTION Get_ySin

!------------------------------------------------------------------------------

      FUNCTION GET_AREA_M2( J ) RESULT( A )
!
!******************************************************************************
!  Function GET_AREA_M2 returns the surface area [m2] of a GEOS-CHEM
!  grid box.  Works for nested grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J (INTEGER) : GEOS-CHEM grid-box index for latitude
!
!  NOTES:
!  (1 ) Surface area is only a function of latitude, since all grid boxes are
!        symmetrical in longitude. (bmy, 3/11/03)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: J

      ! Function value
      REAL*8              :: A

      !=================================================================
      ! GET_AREA_M2 begins here!
      !=================================================================
      A = AREA_M2( J )

      END FUNCTION GET_AREA_M2

!------------------------------------------------------------------------------

      FUNCTION GET_AREA_CM2( J ) RESULT( A )
!
!******************************************************************************
!  Function GET_AREA_CM2 returns the surface area [cm2] of a GEOS-CHEM
!  grid box.  Works for nested grids too. (bmy, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) J (INTEGER) : GEOS-CHEM grid-box index for latitude
!
!  NOTES:
!  (1 ) Surface area is only a function of latitude, since all grid boxes are
!        symmetrical in longitude. (bmy, 3/11/03)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: J

      ! Function value
      REAL*8              :: A

      !=================================================================
      ! GET_AREA_CM2 begins here!
      !=================================================================
      A = AREA_CM2( J )

      END FUNCTION GET_AREA_CM2
!-----------------------------------------------------------------------------

      SUBROUTINE GET_BOUNDING_BOX( COORDS, INDICES )

!******************************************************************************
!  Subroutine GET_BOUNDING_BOX returns the indices which specify the lower
!  left (LL) and upper right (UR) corners of a rectangular region, given the
!  corresponding longitude and latitude values. (bmy, 12/1/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) COORDS  (REAL*8)  : Contains (/LON_LL, LAT_LL, LON_UR, LAT_UR/) values
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) INDICES (INTEGER) : Contains indices (/I_LL, J_LL, I_UR, J_UR/)
!
!  NOTES:
!******************************************************************************
!

#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*8,  INTENT(IN)  :: COORDS(4)
      INTEGER, INTENT(OUT) :: INDICES(4)

      ! Local variables
      INTEGER              :: I, J
      CHARACTER(LEN=255)   :: LOCATION

      !=================================================================
      ! GET_BOUNDING_BOX begins here!
      !=================================================================

      ! Location
      LOCATION = 'GET_BOUNDING_BOX (grid_mod.f)'

      ! Initialize
      INDICES(:) = 0

      !=================================================================
      ! Longitude search
      !=================================================================
      DO I = 1, IIPAR

         ! Locate index corresponding to the lower-left longitude
         IF ( COORDS(1) >  XEDGE(I  )   .and.
     &        COORDS(1) <= XEDGE(I+1) ) INDICES(1) = I

         ! Locate index corresponding to upper-right longitude
         IF ( COORDS(3) >  XEDGE(I  )   .and.
     &        COORDS(3) <= XEDGE(I+1) ) INDICES(3) = I

      ENDDO

      ! Error check lower-left longitude
      IF ( INDICES(1) == 0 ) THEN
         print *, 'Invalid lower-left lon index!'
      ENDIF

      ! Error check upper-right longitude
      IF ( INDICES(3) == 0 ) THEN
         print *, 'Invalid upper-right lon index!'
      ENDIF

      !=================================================================
      ! Latitude search
      !=================================================================
      DO J = 1, JJPAR

         ! Locate index corresponding to the lower-left latitude
         IF ( COORDS(2) >  YEDGE(J  )   .and.
     &        COORDS(2) <= YEDGE(J+1) ) INDICES(2) = J

         ! Locate index corresponding to the upper-right latitude
         IF ( COORDS(4) >  YEDGE(J  )   .and.
     &        COORDS(4) <= YEDGE(J+1) ) INDICES(4) = J

      ENDDO

      ! Error check lower-left longitude
      IF ( INDICES(2) == 0 ) THEN
         print *, 'Invalid lower-left lat index!'
      ENDIF

      ! Error check upper-right longitude
      IF ( INDICES(4) == 0 ) THEN
         print *, 'Invalid upper-right lat index!'
      ENDIF

      END SUBROUTINE GET_BOUNDING_BOX

!------------------------------------------------------------------------------

      FUNCTION ITS_A_NESTED_GRID() RESULT( IT_IS_NESTED )
!
!******************************************************************************
!  Subroutine ITS_A_NESTED_GRID returns TRUE if we are using a nested-grid
!  (i.e. a subset of a global grid) or FALSE otherwise. (bmy, 3/11/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: IT_IS_NESTED

      !=================================================================
      ! ITS_A_NESTED_GRID begins here!
      !=================================================================
      IT_IS_NESTED = IS_NESTED

      END FUNCTION ITS_A_NESTED_GRID

!------------------------------------------------------------------------------

      FUNCTION GET_IJ( LON, LAT ) RESULT ( IIJJ )
!
!******************************************************************************
!  Subroutine GET_IJ returns I and J index for a LON, LAT coord. (dkh, 11/16/06)
!  Reformatted (zhe, dkh, 02/04/11).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LON (REAL*8) : Longitude                          [degrees]
!  (2 ) LAT (REAL*8) : Latitude                           [degrees]
!
!  Function result
!  ============================================================================
!  (1 ) IIJJ(1) (INTEGER) : Long index                    [none]
!  (2 ) IIJJ(2) (INTEGER) : Lati index                    [none]
!
!  NOTES:
!  (1 ) Updated to support nested domains and made much simpler  (zhe, 1/19/11)
!
!******************************************************************************


#     include "CMN_SIZE"     ! DISIZE, DJSIZE

      ! Arguments
      REAL*4  :: LAT, LON

      ! Return
      INTEGER ::  IIJJ(2)

      ! Local variables
      REAL*8  :: TLON, TLAT
      REAL*8  :: I0,   J0

      !=================================================================
      ! GET_IJ begins here!
      !=================================================================
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      TLON = INT( ( LON + 180d0 ) / DISIZE + 1.5d0 )
      TLAT = INT( ( LAT +  90d0 ) / DJSIZE + 1.5d0 )

#if defined( NESTED_CH ) || defined( NESTED_NA )
      TLON = TLON - I0
      TLAT = TLAT - J0
      IF ( TLAT < 1 .or. TLAT > JJPAR ) THEN
         print *, 'Beyond the nested window GET_IJ'
      ENDIF
#else
      IF ( TLON > IIPAR ) TLON = TLON - IIPAR

      ! Check for impossible values
      IF ( TLON > IIPAR .or. TLAT > JJPAR .or.
     &     TLON < 1     .or. TLAT < 1          ) THEN
         print *, 'Error finding grid box GET_IJ'
      ENDIF

#endif

      IIJJ(1) = TLON
      IIJJ(2) = TLAT

      END FUNCTION GET_IJ

!------------------------------------------------------------------------------

      FUNCTION GET_IJ_GLOBAL( LON, LAT ) RESULT ( IIJJ )
!
!******************************************************************************
!  Subroutine GET_IJ returns I and J index for a LON, LAT coord. (dkh, 11/16/06)
!  Reformatted (zhe, dkh, 02/04/11).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LON (REAL*8) : Longitude                          [degrees]
!  (2 ) LAT (REAL*8) : Latitude                           [degrees]
!
!  Function result
!  ============================================================================
!  (1 ) IIJJ(1) (INTEGER) : Long index                    [none]
!  (2 ) IIJJ(2) (INTEGER) : Lati index                    [none]
!
!  NOTES:
!  (1 ) Updated to support nested domains and made much simpler  (zhe, 1/19/11)
!
!******************************************************************************
!


#     include "CMN_SIZE"     ! DISIZE, DJSIZE

      ! Arguments
      REAL*4  :: LAT, LON

      ! Return
      INTEGER ::  IIJJ(2)

      ! Local variables
      REAL*8  :: TLON, TLAT

      !=================================================================
      ! GET_IJ_GLOBAL begins here!
      !=================================================================

      TLON = INT( ( LON + 180d0 ) / DISIZE + 1.5d0 )
      TLAT = INT( ( LAT +  90d0 ) / DJSIZE + 1.5d0 )

      IF ( TLON > IIPAR_L ) TLON = TLON - IIPAR_L

      ! Check for impossible values
      IF ( TLON > IIPAR_L .or. TLAT > JJPAR_L .or.
     &     TLON < 1     .or. TLAT < 1          ) THEN
         print *, 'Error finding grid box GET_IJ_GLOBAL'
      ENDIF

      IIJJ(1) = TLON
      IIJJ(2) = TLAT

      ! Return to calling program
      END FUNCTION GET_IJ_GLOBAL

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GRID
!
!******************************************************************************
!  Subroutine INIT_GRID initializes variables and allocates module arrays.
!  (bmy, 3/11/03, 11/6/08)
!
!  NOTES:
!  (1 ) Fixed typos that caused AREA_CM2_G and AREA_CM2 to be initialized
!        before they were allocated. (bmy, 4/28/03)
!  (2 ) Now define IIGLOB & JJGLOB for 1 x 1.25 grid (bmy, 12/1/04)
!  (3 ) Modified for GCAP 4x5 horizontal grid (swu, bmy, 5/24/05)
!  (4 ) Modifications for 0.5 x 0.666 nested grids (dan, bmy, 11/6/08)
!******************************************************************************

#     include "CMN_SIZE"

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_GRID begins here!
      !=================================================================

      ! Define global sizes for grid.  We need to redefine these here
      ! since for the nested grid, we set IGLOB=IIPAR and JGLOB=JJPAR

      !=================================================================
      ! Allocate global-sized arrays (e.g. use IIGLOB, JJGLOB)
      !=================================================================

#if   defined( GRID1x1 )
      IIGLOB = 360
      JJGLOB = 181
#elif defined( GRID05x0666 )
      IIGLOB = 540
      JJGLOB = 361
#elif defined( GRID025x03125 )
      IIGLOB = 1152    ! (lzh,02/01/2015)
      JJGLOB = 721
#elif defined( GRID1x125 )
      IIGLOB = 288
      JJGLOB = 181
#elif defined( GRID2x25 )
      IIGLOB = 144
      JJGLOB = 91
#elif defined( GRID4x5 ) && defined( GCAP )
      IIGLOB = 72
      JJGLOB = 45
#elif defined( GRID4x5 )
      IIGLOB = 72
      JJGLOB = 46
#endif

      ALLOCATE( XMID_G( IIGLOB ), STAT=AS )
      XMID_G = 0

      ALLOCATE( XEDGE_G( IIGLOB+1 ), STAT=AS )
      XEDGE_G = 0d0

      ALLOCATE( YMID_G( JJGLOB ), STAT=AS )
      YMID_G = 0d0

      ALLOCATE( YEDGE_G( JJGLOB+1 ), STAT=AS )
      YEDGE_G = 0d0

      ALLOCATE( YMID_R_G( JJGLOB ), STAT=AS )
      YMID_R_G = 0d0

      ALLOCATE( YEDGE_R_G( JJGLOB+1 ), STAT=AS )
      YEDGE_R_G = 0d0

      ALLOCATE( AREA_M2_G( JJGLOB ), STAT=AS )
      AREA_M2_G = 0d0

      ALLOCATE( AREA_CM2_G( JJGLOB ), STAT=AS )
      AREA_CM2_G = 0d0

      !=================================================================
      ! Allocate window-sized arrays (e.g. use IIPAR, JJPAR)
      !=================================================================
      ALLOCATE( XMID( IIPAR ), STAT=AS )
      XMID = 0

      ALLOCATE( XEDGE( IIPAR+1 ), STAT=AS )
      XEDGE = 0d0

      ALLOCATE( YMID( JJPAR ), STAT=AS )
      YMID = 0d0

      ALLOCATE( YEDGE( JJPAR+1 ), STAT=AS )
      YEDGE = 0d0

      ALLOCATE( YSIN( JJPAR+1 ), STAT=AS )
      YSIN = 0d0

      ALLOCATE( YMID_R( 1:JJPAR ), STAT=AS )
      YMID_R = 0d0

      ALLOCATE( YMID_R_W( 0:JJPAR+1 ), STAT=AS )
      YMID_R_W = 0d0

      ALLOCATE( YEDGE_R( JJPAR+1 ), STAT=AS )
      YEDGE_R = 0d0

      ALLOCATE( AREA_M2( JJPAR ), STAT=AS )
      AREA_M2 = 0d0

      ALLOCATE( AREA_CM2( JJPAR ), STAT=AS )
      AREA_CM2 = 0d0

      ALLOCATE( DX_COV(JJPAR), STAT=AS )
      DX_COV= 0d0

      ALLOCATE( DY_COV(JJPAR), STAT=AS )
      DY_COV= 0d0
      !=================================================================
      ! Also test for 1x1 nested-grid (smaller than global size)
      !=================================================================
      IF ( IIPAR == IIGLOB .and. JJPAR == JJGLOB ) THEN
         IS_NESTED = .FALSE.
      ELSE
         IS_NESTED = .TRUE.
      ENDIF

      END SUBROUTINE INIT_GRID

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GRID
!
!******************************************************************************
!  Subroutine CLEANUP_GRID deallocates all module arrays
!  (bmy, 3/11/03, 11/6/08)
!
!  NOTES:
!  (1 ) Now also deallocate YMID_R_W (dan, bmy, 11/6/08)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_GRID begins here!
      !=================================================================

      ! Deallocate window arrays
      IF ( ALLOCATED( XMID       ) ) DEALLOCATE( XMID       )
      IF ( ALLOCATED( XEDGE      ) ) DEALLOCATE( XEDGE      )
      IF ( ALLOCATED( YMID       ) ) DEALLOCATE( YMID       )
      IF ( ALLOCATED( YEDGE      ) ) DEALLOCATE( YEDGE      )
      IF ( ALLOCATED( YMID_R     ) ) DEALLOCATE( YMID_R     )
      IF ( ALLOCATED( YMID_R_W   ) ) DEALLOCATE( YMID_R_W   )
      IF ( ALLOCATED( YEDGE_R    ) ) DEALLOCATE( YEDGE_R    )
      IF ( ALLOCATED( AREA_M2    ) ) DEALLOCATE( AREA_M2    )
      IF ( ALLOCATED( AREA_CM2   ) ) DEALLOCATE( AREA_CM2   )
      IF ( ALLOCATED( YSIN       ) ) DEALLOCATE( YSIN       )

      END SUBROUTINE CLEANUP_GRID

!------------------------------------------------------------------------------

      END MODULE GRID_MOD

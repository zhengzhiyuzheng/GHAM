MODULE MACHINE_LEARN
  use iso_c_binding
  implicit none

  interface

    subroutine predict(booster, in_data, res) bind(C, name='predict')
      use iso_c_binding
      type(C_PTR), intent(in)             :: booster
      real(C_FLOAT), intent(in)           :: in_data(*)
      real(C_FLOAT), intent(inout)        :: res
    end subroutine predict

    subroutine load_model(model_path, booster) bind(C, name="load_model")
      use iso_c_binding
      character(kind=C_CHAR, len=1), intent(in), dimension(*)   :: model_path(*)
      type(C_PTR), intent(out)                                  :: booster
    end subroutine

    subroutine predict_tr(model, input_size, in_data, res) bind(C, name='predict_tr')
      use iso_c_binding
      type(C_PTR), intent(in), value      :: model
      integer(C_INT), intent(in), value   :: input_size
      real(C_FLOAT), intent(in)           :: in_data(*)
      real(C_FLOAT), intent(inout)        :: res(*)
    end subroutine predict_tr

    subroutine load_model_tr(model_path, model) bind(C, name="load_model_tr")
      use iso_c_binding
      character(kind=C_CHAR, len=1), intent(in), dimension(*)   :: model_path(*)
      type(C_PTR), intent(out)                                  :: model
    end subroutine

  end interface

  CHARACTER(kind=c_char, len=255) :: model_path
  TYPE(C_PTR)                   :: booster_ptr

  contains

  SUBROUTINE INIT_XG()

    model_path = "/work/ese-zhengzy/model_z/GHSM/XGboost/BC_model.json"//c_null_char
    call load_model(trim(model_path), booster_ptr)

  END SUBROUTINE

END MODULE

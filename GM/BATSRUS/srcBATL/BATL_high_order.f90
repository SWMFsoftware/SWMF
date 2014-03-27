!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_high_order
  implicit none

  save

contains

  !======================================================================
  real function correct_face_value(FaceValue, CellValue_I)

    ! FaceValue is at cell face. CellValue_I are cell centered.
    ! Return 6th order approximation

    real, intent(in):: CellValue_I(4), FaceValue
    real:: Der2, Der4
    real, parameter:: c1over6 = 1./6, c1over180 = 1./180
    !--------------------------------------------------------------------
    Der2 = c1over6*(CellValue_I(2) - 2*FaceValue + CellValue_I(3))
    Der4 = c1over180*(16*FaceValue - &
         9*(CellValue_I(2) + CellValue_I(3)) + &
         CellValue_I(1) + CellValue_I(4))
    correct_face_value = FaceValue - Der2 + Der4

  end function correct_face_value
  !======================================================================

  real function calc_center_first_derivate(CellValue_I, DxIn)

    ! Calculate df/dx at x=x_i with f(k), where k = i-3,i-2 ... i+3.
    ! Directly combine CellValue_I to get df/dx can save some 
    ! computation compare to the approach in this subroutine.

    real, intent(in):: CellValue_I(7)
    real, optional, intent(in):: DxIn
    real:: Dx
    real:: FaceL, FaceR ! FaceL = F_{i-1/2}, FaceR = F_{i+1/2}
    real:: CorrectedFaceL, CorrectedFaceR
    real, parameter:: c3over256 = 3./256, c25over256 = 25./256, &
         c150over256 = 150./256
    integer, parameter:: i = 4
    ! ----------------------------------------------------------------------

    Dx = 1.0
    if(present(DxIn)) Dx = DxIn

    FaceL = calc_face_value(CellValue_I(1:6))
    FaceR = calc_face_value(CellValue_I(2:7))

    CorrectedFaceL = correct_face_value(FaceL, CellValue_I(i-2:i+1))
    CorrectedFaceR = correct_face_value(FaceR, CellValue_I(i-1:i+2)) 

    calc_center_first_derivate = &
         (CorrectedFaceR - CorrectedFaceL)/Dx
  end function calc_center_first_derivate
  !======================================================================

  real function calc_face_value(CellValue_I)
    ! Calculate f_{i+1/2} with f(k), where k = i-2,i-1 ... i+3

    real, intent(in):: CellValue_I(6)
    real, parameter:: c3over256 = 3./256, c25over256 = 25./256, &
         c150over256 = 150./256
    !----------------------------------------------------------------------

    calc_face_value = c3over256*(CellValue_I(1) + CellValue_I(6)) - &
         c25over256*(CellValue_I(2) + CellValue_I(5)) + &
         c150over256*(CellValue_I(3) + CellValue_I(4))
  end function calc_face_value

  !===========================================================================
end module BATL_high_order

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_realkinds

  use , intrinsic :: iso_fortran_env
  use , intrinsic :: ieee_arithmetic

#define __SYSTEM_NAN_32__            -4194304_int32
#define __SYSTEM_INF_32__          2139095040_int32
#define __SYSTEM_NAN_64__   -2251799813685248_int64
#define __SYSTEM_INF_64__ 9218868437227405312_int64

  implicit none

  public

  integer , parameter :: rk4  = real32
  integer , parameter :: rk8  = real64
  integer , parameter :: rk16 = real128
  integer , parameter :: rkx = rk8
  real(rk8) , parameter :: nan = transfer(__SYSTEM_NAN_64__, 1._real64)
  real(rk8) , parameter :: inf = transfer(__SYSTEM_INF_64__, 1._real64)

  interface is_nan
    module procedure is_nan_single
    module procedure is_nan_double
  end interface

  interface is_inf
    module procedure is_inf_single
    module procedure is_inf_double
  end interface

  contains

  logical elemental function is_nan_double(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_nan_double = (ieee_class(x) == ieee_quiet_nan .or. &
                     ieee_class(x) == ieee_signaling_nan)
  end function is_nan_double

  logical elemental function is_inf_double(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_inf_double = (ieee_class(x) == ieee_negative_inf .or. &
                     ieee_class(x) == ieee_positive_inf)
  end function is_inf_double

  logical elemental function is_nan_single(x)
    implicit none
    real(rk4) , intent(in) :: x
    is_nan_single = (ieee_class(x) == ieee_quiet_nan .or. &
                     ieee_class(x) == ieee_signaling_nan)
  end function is_nan_single

  logical elemental function is_inf_single(x)
    implicit none
    real(rk4) , intent(in) :: x
    is_inf_single = (ieee_class(x) == ieee_negative_inf .or. &
                     ieee_class(x) == ieee_positive_inf)
  end function is_inf_single

end module mod_realkinds

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

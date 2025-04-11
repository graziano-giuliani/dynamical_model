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

module mod_regcm_types
  use mod_realkinds
  use mod_intkinds

  implicit none

  public

  type model_area
    logical :: bandflag
    logical :: crmflag
    logical :: has_bdy
    logical :: has_bdyleft , has_bdyright , has_bdytop , has_bdybottom
    logical :: has_bdytopleft , has_bdytopright
    logical :: has_bdybottomleft , has_bdybottomright
    integer(ik4) , dimension(2) :: location
    integer(ik4) :: left , right , top , bottom
    integer(ik4) :: topleft , topright , bottomleft , bottomright
    integer(ik4) :: ibt1 , ibt2 , ibt3 , ibt4 , ibt6
    integer(ik4) :: ibb1 , ibb2 , ibb3 , ibb4 , ibb6
    integer(ik4) :: jbl1 , jbl2 , jbl3 , jbl4 , jbl6
    integer(ik4) :: jbr1 , jbr2 , jbr3 , jbr4 , jbr6
  end type model_area

  type domain
    real(rkx) , pointer , dimension(:,:) :: ht => null( )
    real(rkx) , pointer , dimension(:,:) :: lndcat => null( )
    real(rkx) , pointer , dimension(:,:) :: lndtex => null( )
    real(rkx) , pointer , dimension(:,:) :: mask => null( )
    real(rkx) , pointer , dimension(:,:) :: area => null( )
    real(rkx) , pointer , dimension(:,:) :: dlat => null( )
    real(rkx) , pointer , dimension(:,:) :: dlon => null( )
    real(rkx) , pointer , dimension(:,:) :: ulat => null( )
    real(rkx) , pointer , dimension(:,:) :: ulon => null( )
    real(rkx) , pointer , dimension(:,:) :: vlat => null( )
    real(rkx) , pointer , dimension(:,:) :: vlon => null( )
    real(rkx) , pointer , dimension(:,:) :: xlat => null( )
    real(rkx) , pointer , dimension(:,:) :: xlon => null( )
    real(rkx) , pointer , dimension(:,:) :: msfu => null( )
    real(rkx) , pointer , dimension(:,:) :: msfv => null( )
    real(rkx) , pointer , dimension(:,:) :: hx => null( )
    real(rkx) , pointer , dimension(:,:) :: hy => null( )
    real(rkx) , pointer , dimension(:,:) :: msfx => null( )
    real(rkx) , pointer , dimension(:,:) :: msfd => null( )
    real(rkx) , pointer , dimension(:,:) :: coriol => null( )
    real(rkx) , pointer , dimension(:,:) :: coriou => null( )
    real(rkx) , pointer , dimension(:,:) :: coriov => null( )
    real(rkx) , pointer , dimension(:,:) :: ef => null( )
    real(rkx) , pointer , dimension(:,:) :: ddx => null( )
    real(rkx) , pointer , dimension(:,:) :: ddy => null( )
    real(rkx) , pointer , dimension(:,:) :: ex => null( )
    real(rkx) , pointer , dimension(:,:) :: crx => null( )
    real(rkx) , pointer , dimension(:,:) :: cry => null( )
    real(rkx) , pointer , dimension(:,:) :: dmdy => null( )
    real(rkx) , pointer , dimension(:,:) :: dmdx => null( )
    real(rkx) , pointer , dimension(:,:) :: xmsf => null( )
    real(rkx) , pointer , dimension(:,:) :: dmsf => null( )
    integer(ik4) , pointer , dimension(:,:) :: ldmsk => null( )
  end type domain

  type atmosphere
    real(rkx) , pointer , dimension(:,:,:) :: u => null( )
    real(rkx) , pointer , dimension(:,:,:) :: v => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ux => null( )
    real(rkx) , pointer , dimension(:,:,:) :: vx => null( )
    real(rkx) , pointer , dimension(:,:,:) :: w => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pai => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:) :: p => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rho => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tvirt => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tetav => null( )
    real(rkx) , pointer , dimension(:,:,:) :: qs => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qx => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zeta => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zetaf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dz => null( )
    real(rkx) , pointer , dimension(:,:,:) :: fmz => null( )
    real(rkx) , pointer , dimension(:,:,:) :: fmzf => null( )
  end type atmosphere

  type surfstate
    real(rkx) , pointer , dimension(:,:) :: ps => null( )
  end type surfstate

end module mod_regcm_types

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

module nlte_aliarms

!
! provides calculation of non-LTE heating rates by ALI-ARMS non-LTE code
!
  use ppgrid,             only: pcols, pver, pverp
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use cam_abortutils,     only: endrun
  use cam_logfile,        only: iulog
  use spmd_utils,         only: masterproc

  implicit none
  private
  save

! Public interfaces
  public &
     nlte_aliarms_init, &
     nlte_aliarms_calc

contains

!-----------------------------------------------------------------
  subroutine nlte_aliarms_init
!-----------------------------------------------------------------
!
!
!-----------------------------------------------------------------

  use cam_history,  only: addfld

  implicit none


  if (masterproc) then
    write(iulog,*) 'init: ALI-ARMS non-LTE code'
  end if

  call addfld ('ALIARMS_Q',(/ 'lev' /), 'A','K/s','Non-LTE LW CO2 heating')

  end subroutine nlte_aliarms_init
  
!-----------------------------------------------------------------
  subroutine nlte_aliarms_calc (lchnk,ncol,pmid,t,xo2,xo,xn2,xco2,cool)
!-----------------------------------------------------------------
!
!
!-----------------------------------------------------------------

  use cam_history,  only: outfld

  implicit none

! Input variables
  integer, intent(in) :: ncol                          ! number of atmospheric columns
  integer, intent(in) :: lchnk                         ! chunk identifier

  real(r8), intent(in) :: pmid(pcols,pver)             ! model pressure at mid-point
  real(r8), intent(in) :: t(pcols,pver)                ! Neutral temperature (K)
  real(r8), intent(in) :: xco2(pcols,pver)             ! CO2 profile
  real(r8), intent(in) :: xn2(pcols,pver)              ! N2 profile
  real(r8), intent(in) :: xo(pcols,pver)               ! O profile
  real(r8), intent(in) :: xo2(pcols,pver)              ! O2 profile

! Output variables
  real(r8), intent(out) :: cool(pcols,pver)            ! CO2 NLTE cooling rate

! local variables

  real(r8), dimension(pver) :: p, tn, zkm
  real(r8), dimension(pver) :: co2_vmr, o_vmr, n2_vmr, o2_vmr
  real(r8), dimension(pver) :: ali_cool
  
  integer :: icol
  
  cool(:,:) = 0.0_r8
 
  do icol=1,ncol
  
      p = pmid(icol,:)*1.0e-5_r8 ! conver pmid in Pa to bars
      zkm = -7.0_r8*log(p/1.013_r8)
      tn = t(icol,:)
      
      co2_vmr = xco2(icol,:)
      o_vmr = xo(icol,:)
      n2_vmr = xn2(icol,:)
      o2_vmr = xo2(icol,:)
      
      call ali(zkm, p, tn, co2_vmr, o_vmr, n2_vmr, o2_vmr, ali_cool, pver)
      
      cool(icol,:) = ali_cool(:) 
  
  enddo
  
  call outfld ('ALIARMS_Q', cool, pcols, lchnk)
  
  end subroutine nlte_aliarms_calc
  
  end module nlte_aliarms


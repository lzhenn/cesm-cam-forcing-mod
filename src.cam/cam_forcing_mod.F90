module cam_forcing_mod
!-----------------------------------------------------------------------
! Purpose:
!
! CESM CAM FORCING MODULE (CCFM) is an open source modification module for the
! CESM users using external forcing files to conduct sensitive experiments in
! the CAM workflow. This flexible architecture enables you to deploy the module
! to a series of different versions of CESM, although the source code was
! developed based on CESM1.2.2.
!
!---------------------------Code history-------------------------------
!
! 2019-04-30  Zhenning Li,    Creation of module
! 2019-05-05  Zhenning Li,    Add cam_forc_init()
!
!-----------------------------------------------------------------------

use spmd_utils,    only: masterproc
use cam_logfile,   only: iulog
use abortutils,    only: endrun
use shr_kind_mod,  only: r8 => shr_kind_r8


implicit none
private
save

public :: &
    cam_forc_readnl,             &! read namelist from file
    cam_forc_init                 ! initiate the CAM forcing module

! Private module data
character(len=16), parameter :: unset_str = 'UNSET'

character(len=256)  :: forcing_file_loc     = unset_str
character(len=32)   :: forcing_process      = unset_str  
character(len=32)   :: forcing_field1       = unset_str  
character(len=32)   :: forcing_field2       = unset_str  
character(len=32)   :: forcing_field3       = unset_str  
character(len=32)   :: forcing_method      = unset_str  

integer :: forcing_relax_time   = 86400
!======================================================================= 
contains
!======================================================================= 


subroutine cam_forc_readnl(nlfile)
!-----------------------------------------------------------------------
!     
! Purpose:     
! Read and check the cam_forcing needed namelist variables
!     
!---------------------------Code history-------------------------------
!
! 2019-04-30  Zhenning Li,    Creation of module
!
!-----------------------------------------------------------------------

    use namelist_utils,  only: find_group_name
    use units,              only: getunit, freeunit
    use mpishorthand
    
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cam_forc_readnl'

    namelist /cam_forcing_mod/ forcing_file_loc, forcing_process, forcing_method, &
    forcing_field1, forcing_field2, forcing_field3, forcing_relax_time
    
   !-----------------------------------------------------------------------------

    if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(nlfile), status='old' )
        call find_group_name(unitn, 'cam_forcing_mod', status=ierr)
        if (ierr == 0) then
            read(unitn, cam_forcing_mod, iostat=ierr)
            if (ierr /= 0) then
                call endrun(subname // ':: ERROR reading namelist')
            end if
        end if
        close(unitn)
        call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(forcing_file_loc, len(forcing_file_loc) , mpichar, 0, mpicom)
    call mpibcast(forcing_process,  len(forcing_process)  , mpichar, 0, mpicom)
    call mpibcast(forcing_method,   len(forcing_method)   , mpichar, 0, mpicom)
    call mpibcast(forcing_field1,   len(forcing_field1)   , mpichar, 0, mpicom)
    call mpibcast(forcing_field2,   len(forcing_field2)   , mpichar, 0, mpicom)
    call mpibcast(forcing_field3,   len(forcing_field3)   , mpichar, 0, mpicom)
    call mpibcast(forcing_relax_time,                   1 ,  mpiint, 0, mpicom)
#endif

    ! Error checking:
    if (.not. (forcing_process == 'convect_deep' .or. forcing_process == 'convect_shallow' .or. forcing_process == 'vertical_diffusion' .or. forcing_process == 'micro_phy')) then
        write(iulog,*)'forcing_mod: illegal value of forcing_process:',forcing_process
        call endrun('forcing_mod: illegal value of forcing_process')
    endif
    
end subroutine cam_forc_readnl 

subroutine cam_forc_init()
!-----------------------------------------------------------------------
!     
! Purpose:     
! Initiate the cam_forcing with input file and parameters
!     
!---------------------------Code history-------------------------------
!
! 2019-05-05  Zhenning Li, Creation of the subroutine
!
!-----------------------------------------------------------------------

    write(iulog,*) forcing_file_loc, ' forcing_process:',forcing_process
end subroutine cam_forc_init

end module cam_forcing_mod



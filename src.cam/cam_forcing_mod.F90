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

use spmd_utils,     only: masterproc
use cam_logfile,    only: iulog
use abortutils,     only: endrun
use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pcols, begchunk, endchunk, pver
use infnan,         only: nan, assignment(=)
use physconst,      only: cpair
implicit none
private
save

public :: &
    cam_forc_readnl,             &! read namelist from file
    cam_forc_reg,                &! register the CAM forcing module with variables
    cam_forc_init,               &! initiate the CAM forcing module
    cam_forc_exe,                &! exceute the CAM forcing module
    cam_forc_debug                ! output the debug info for CAM forcing module

public :: forc_shf                    ! public the 2-D sensible heat flux forcing
public :: forc_lhf                    ! public the 2-D latent heat flux forcing

public :: forc_utend                    ! public the 3-D utend forcing 
public :: forc_vtend                    ! public the 3-D vtend forcing 
public :: forc_ttend                    ! public the 3-D ttend forcing 

real(r8), allocatable:: forc_shf(:,:,:)    ! state var
real(r8), allocatable:: forc_lhf(:,:,:)    ! state var

real(r8), allocatable:: forc_utend(:,:,:,:)    ! state var
real(r8), allocatable:: forc_vtend(:,:,:,:)    ! state var
real(r8), allocatable:: forc_ttend(:,:,:,:)    ! state var

! Private module data
character(len=16), parameter :: unset_str = 'UNSET'

character(len=256)  :: forcing_file_loc     = unset_str
character(len=32)   :: forcing_process      = unset_str  
character(len=32)   :: forcing_method       = unset_str  
character(len=32)   :: forcing_time_res     = unset_str  

logical :: forcing_utend_flag   = .false.
logical :: forcing_vtend_flag   = .false.
logical :: forcing_ttend_flag   = .false.
logical :: forcing_shf_flag     = .false.
logical :: forcing_lhf_flag     = .false.

integer :: forcing_relax_time   = 86400
integer :: forcing_timeframe    = 12
integer :: forcing_start_day    = 1

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
    forcing_utend_flag, forcing_vtend_flag, forcing_ttend_flag, forcing_relax_time,&
    forcing_shf_flag, forcing_lhf_flag, forcing_time_res, forcing_start_day,&
    forcing_timeframe
    
    !-----------------------------------------------------------------------------
    
    ! Try to read namelist variables
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
    
    forcing_process=trim(forcing_process)
    forcing_method=trim(forcing_method)
    forcing_time_res=trim(forcing_time_res)

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(forcing_file_loc, len(forcing_file_loc) , mpichar, 0, mpicom)
    call mpibcast(forcing_process,  len(forcing_process)  , mpichar, 0, mpicom)
    call mpibcast(forcing_method,   len(forcing_method)   , mpichar, 0, mpicom)
    call mpibcast(forcing_time_res, len(forcing_time_res) , mpichar, 0, mpicom)
    call mpibcast(forcing_relax_time,                   1 ,  mpiint, 0, mpicom)
    call mpibcast(forcing_timeframe,                    1 ,  mpiint, 0, mpicom)
    call mpibcast(forcing_start_day,                    1 ,  mpiint, 0, mpicom)
    call mpibcast(forcing_utend_flag,                   1 ,  mpilog, 0, mpicom)
    call mpibcast(forcing_vtend_flag,                   1 ,  mpilog, 0, mpicom)
    call mpibcast(forcing_ttend_flag,                   1 ,  mpilog, 0, mpicom)
    call mpibcast(forcing_shf_flag,                     1 ,  mpilog, 0, mpicom)
    call mpibcast(forcing_lhf_flag,                     1 ,  mpilog, 0, mpicom)
#endif

    ! Error checking:
    if (.not. (forcing_process == 'convect_deep' .or. forcing_process == 'convect_shallow' .or. forcing_process == 'vert_diff_flx' &
    .or. forcing_process == 'vert_diff_tend' .or. forcing_process == 'micro_phy')) then
        write(iulog,*)'forcing_mod: illegal value of forcing_process:',forcing_process
        call endrun('forcing_mod: illegal value of forcing_process')
    endif
    
end subroutine cam_forc_readnl 

subroutine cam_forc_reg()
!-----------------------------------------------------------------------
!     
! Purpose:     
! Register the cam_forcing with three 2-D and three 3-D
!     
!---------------------------Code history-------------------------------
!
! 2019-05-08  Zhenning Li, Creation of the subroutine
!
!-----------------------------------------------------------------------
     
    write(iulog,*) 'CAM_FORC_REG: Register varibales.'

    ! register the forcing variables
    if (forcing_ttend_flag == .true.) then
        allocate (forc_ttend(pcols,pver,begchunk:endchunk,forcing_timeframe))
        forc_ttend=nan
    end if
    if (forcing_utend_flag == .true.) then
        allocate (forc_utend(pcols,pver,begchunk:endchunk,forcing_timeframe))
        forc_utend=nan
    end if
    if (forcing_vtend_flag == .true.) then
        allocate (forc_vtend(pcols,pver,begchunk:endchunk,forcing_timeframe))
        forc_vtend=nan
    end if
    if (forcing_shf_flag == .true.) then
        write(iulog,*) 'CAM_FORC_REG: Register SHFLX'
        allocate (forc_shf(pcols,begchunk:endchunk,forcing_timeframe))
        forc_shf=nan
    end if
    if (forcing_lhf_flag == .true.) then
        allocate (forc_lhf(pcols,begchunk:endchunk,forcing_timeframe))
        forc_lhf=nan
    end if
end subroutine cam_forc_reg

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
    use ncdio_atm,          only: infld
    use pio,                only: file_desc_t, pio_closefile
    use cam_pio_utils,      only: cam_pio_openfile

    type(file_desc_t), pointer ::  fh_forc   ! forcing file handle
    logical :: found=.false.
    integer :: nm   ! time frame controller 
    
    write(iulog,*) 'CAM_FORC_INIT: Reading forcing file:', forcing_file_loc
    write(iulog,*) 'CAM_FORC_INIT: Forcing process:', forcing_process
    write(iulog,*) 'CAM_FORC_INIT: Forcing time resolution:', forcing_time_res
    write(iulog,*) 'CAM_FORC_INIT: Forcing method:', forcing_method
    write(iulog,*) 'CAM_FORC_INIT: Forcing frames:', forcing_timeframe


    allocate(fh_forc)
    call cam_pio_openfile(fh_forc, forcing_file_loc, 0)
   
    ! Read the forcing file
    if (forcing_shf_flag) then
        do nm=1, forcing_timeframe
            write(iulog,*) 'CAM_FORC_INIT: Reading Attemptation', nm
            call infld('SHFLXF', fh_forc, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                forc_shf(:,:,nm), found, grid_map='PHYS', timelevel=nm)
        end do
    end if
    if (forcing_lhf_flag) then
        do nm=1, forcing_timeframe
            call infld('LHFLXF', fh_forc, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
                forc_lhf(:,:,nm), found, grid_map='PHYS', timelevel=nm)
        end do
    end if   
    if (forcing_utend_flag) then
        do nm=1, forcing_timeframe
            call infld('UTEND', fh_forc, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                forc_utend(:,:,:,nm), found, grid_map='PHYS', timelevel=nm)
        end do
    end if
    if (forcing_vtend_flag) then
        do nm=1, forcing_timeframe
            call infld('VTEND', fh_forc, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                forc_vtend(:,:,:,nm), found, grid_map='PHYS', timelevel=nm)
        end do
    end if
    if (forcing_ttend_flag) then
        do nm=1, forcing_timeframe
            call infld('TTEND', fh_forc, 'lon', 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                forc_ttend(:,:,:,nm), found, grid_map='PHYS', timelevel=nm)
        end do
        ! Convert T forcing to dry static energy
        forc_ttend=forc_ttend*cpair
    end if
    
    call pio_closefile(fh_forc)
    deallocate(fh_forc)
    nullify(fh_forc)

end subroutine cam_forc_init

subroutine cam_forc_debug(phy_proc_now, ptend, state, lchnk, exe_flag, cam_in)
!-----------------------------------------------------------------------
!     
! Purpose:     
! Output the debug info for cam_forcing module 
!     
!---------------------------Code history-------------------------------
!
! 2019-05-10  Zhenning Li, Creation of the subroutine
!
!-----------------------------------------------------------------------
    use physics_types,    only: physics_state, physics_ptend
    use camsrfexch,       only: cam_in_t
    use time_manager,       only: get_curr_calday

    type(physics_state), intent(inout)  :: state        ! state variable for nudging
    type(physics_ptend)                 :: ptend        ! indivdual parameterization tendencies
    type(cam_in_t), optional            :: cam_in       ! cam_in for flux operation
    
    character(len=*)                   :: phy_proc_now ! which physical process now in
    character(len=16)                   :: exe_state    ! whether after or before the CAM_FORC_EXE CALL
    integer                             :: lchnk        ! current chunk # the model is operated on 
    integer                             :: icol         ! column loop indicator
    logical                             :: exe_flag     ! has the forcing module already executed?
    
    real(r8)              :: calday   ! current model day
    real(r8)              :: stat_lat, stat_lon
    
    calday = get_curr_calday()
    
    if (exe_flag) then
        exe_state='*AFTER EXE*'
    else
        exe_state='*BEFORE EXE*'
    end if    

    if (phy_proc_now == forcing_process) then
        ! only execute in the first step of one certain calendar day and over
        ! (EQ, 120E)
        if ((calday-floor(calday))<0.01) then 
            do icol=1, pcols ! cols in chunk
                stat_lat=state%lat(icol)*180.0/3.1415926
                stat_lon=state%lon(icol)*180.0/3.1415926
                if ((abs(stat_lon-120.0)<1.5).and. &
                (abs(stat_lat-0.0)<1.0)) then
                    if ((forcing_utend_flag)) then
                        write(iulog,"(A16,A20,F8.2,A5,F6.2,A5,F6.2,A15,F6.2,A15)") exe_state, "CAM_FORC_DEBUG: Calday:",calday,&
                        "Lat:",stat_lat,"Lon:",stat_lon," ptend%u(lv13)=", ptend%u(icol,13)*86400,"m/s/day"
                    end if
                    if ((forcing_vtend_flag)) then
                        write(iulog,"(A16, A20,F8.2,A5,F6.2,A5,F6.2,A15,F6.2,A15)") exe_state, "CAM_FORC_DEBUG: Calday:",calday,&
                        "Lat:",stat_lat,"Lon:",stat_lon," ptend%v(lv13)=", ptend%v(icol,13)*86400,"m/s/day"
                    end if
                    if ((forcing_ttend_flag)) then
                        write(iulog,"(A16, A20,F8.2,A5,F6.2,A5,F6.2,A15,F6.2,A15)") exe_state, "CAM_FORC_DEBUG: Calday:",calday,&
                        "Lat:",stat_lat,"Lon:",stat_lon," ptend%t(lv13)=", ptend%s(icol,13)/cpair*86400,"K/day"
                    end if
                    
                    if ((forcing_shf_flag)) then
                        write(iulog,"(A16, A20,F8.2,A5,F6.2,A5,F6.2,A15,F6.2,A15)") exe_state, "CAM_FORC_DEBUG: Calday:",calday,&
                        "Lat:",stat_lat,"Lon:",stat_lon," cam_in%shf=", cam_in%shf(icol),"W/m^2"
                    end if       
                    if ((forcing_lhf_flag)) then
                        write(iulog,"(A16, A20,F8.2,A5,F6.2,A5,F6.2,A15,F6.2,A15)") exe_state, "CAM_FORC_DEBUG: Calday:",calday,&
                        "Lat:",stat_lat,"Lon:",stat_lon," cam_in%lhf=", cam_in%lhf(icol),"W/m^2"
                    end if      
                end if
            end do
       end if
    end if

end subroutine cam_forc_debug


subroutine cam_forc_exe(phy_proc_now, ptend, state, lchnk, cam_in)
!-----------------------------------------------------------------------
!     
! Purpose:     
! Execute the cam_forcing regarding different forcing process
! and forcing method.
!     
!---------------------------Code history-------------------------------
!
! 2019-05-08  Zhenning Li, Creation of the subroutine
!
!-----------------------------------------------------------------------
    use physics_types,   only: physics_state, physics_ptend
    use camsrfexch,       only: cam_in_t

    type(physics_state), intent(inout)  :: state        ! state variable for nudging
    type(physics_ptend)                 :: ptend        ! indivdual parameterization tendencies
    type(cam_in_t), optional            :: cam_in       ! cam_in for flux operation
    
    character(len=*)                   :: phy_proc_now
    integer     :: lchnk      ! current chunk # the model is operated on 
    
    if (phy_proc_now == forcing_process) then
        if ((forcing_utend_flag)) then
            call var_forc_exe3d(ptend%u, forc_utend(:,:,lchnk,:), state%u)
        end if
        if ((forcing_vtend_flag)) then
            call var_forc_exe3d(ptend%v, forc_vtend(:,:,lchnk,:), state%v)
        end if
        if ((forcing_ttend_flag)) then
            call var_forc_exe3d(ptend%s, forc_ttend(:,:,lchnk,:), state%t)
        end if
        
        if ((forcing_shf_flag)) then
            call var_forc_exe2d(cam_in%shf, forc_shf(:,lchnk,:))
        end if       
        if ((forcing_lhf_flag)) then
            call var_forc_exe2d(cam_in%lhf, forc_lhf(:,lchnk,:))
        end if      
    end if

end subroutine cam_forc_exe

subroutine var_forc_exe3d(iner_tend, ext_forc, iner_state)
!-----------------------------------------------------------------------
!     
! Purpose:     
! Execute the 3-D variable cam_forcing regarding different forcing method 
!     
!---------------------------Code history-------------------------------
!
! 2019-05-08  Zhenning Li, Creation of the subroutine
!
!-----------------------------------------------------------------------
    use time_manager,       only: get_curr_calday

    real(r8)              :: iner_tend(pcols, pver)
    real(r8)              :: ext_forc(pcols, pver, forcing_timeframe)
    real(r8), optional    :: iner_state(pcols,pver)   

   
    real(r8)              :: G0
    real(r8)              :: calday   ! current model day

    integer     :: iframe     ! current frame indicator of the forcing data
    integer, parameter :: day_rank(12)=(/31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/)

    ! nudging parameter
    G0=1.0/forcing_relax_time

    if(forcing_time_res == 'monthly') then
        ! below to get which month (imon) we are in
        do iframe = 1,12
            if (calday < day_rank(iframe)+1) then
                exit
            end if
        end do
    else
        ! if the forcing in daily resolution, get which frame of forcing data
        ! should be used, otherwise return (no action)
        if (calday >= forcing_start_day .and. &
        calday<=forcing_start_day+forcing_timeframe) then
            iframe=floor(calday)-forcing_start_day+1
        else
            return
        end if 
    end if
    
    select case (forcing_method)
        case ("replacing")
            iner_tend=ext_forc(:,:,iframe)
        case ("imposing")
            iner_tend=iner_tend+ext_forc(:,:,iframe)
        case ("nudging")
            iner_tend=iner_tend+G0*(ext_forc(:,:,iframe)-iner_state)
    end select

end subroutine var_forc_exe3d


subroutine var_forc_exe2d(iner_camflx, ext_forc)
!-----------------------------------------------------------------------
!     
! Purpose:     
! Execute the 2-D variable cam_forcing regarding different forcing method
!     
!---------------------------Code history-------------------------------
!
! 2019-05-09  Zhenning Li, Creation of the subroutine
!
!-----------------------------------------------------------------------
    use time_manager,       only: get_curr_calday

    real(r8)              :: iner_camflx(pcols)
    real(r8)              :: ext_forc(pcols, forcing_timeframe)

   
    real(r8)              :: G0
    real(r8)              :: calday   ! current model day

    integer     :: iframe     ! current frame indicator of the forcing data
    integer, parameter :: day_rank(12)=(/31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/)

    G0=1.0/forcing_relax_time
    
    if(forcing_time_res == 'monthly') then
        !Below to get out which month (imon) we are in
        do iframe = 1,12
            if (calday < day_rank(iframe)+1) then
                exit
            end if
        end do
    else
        ! if the forcing in daily resolution, get which frame of forcing data
        ! should be used, otherwise return (no action)
        if (calday >= forcing_start_day .and. &
        calday<=forcing_start_day+forcing_timeframe) then
            iframe=floor(calday)-forcing_start_day+1
        else
            return
        end if 
    end if
    
    select case (forcing_method)
        case ("replacing")
            iner_camflx=ext_forc(:,iframe)
        case ("imposing")
            iner_camflx=iner_camflx+ext_forc(:,iframe)
        case ("nudging")
            iner_camflx=iner_camflx+G0*(ext_forc(:,iframe)-iner_camflx)
    end select

end subroutine var_forc_exe2d
end module cam_forcing_mod



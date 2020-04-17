module MpiExchangeGenModule
  use KindModule, only: DP, I4B  
  use ConstantsModule, only: LENMODELNAME, LINELENGTH
  use MpiWrapper, only: mpiwrpinit, mpiwrpfinalize
  
  implicit none
  
  private
  
  ! -- Public functions
  public :: mpi_initialize
  public :: mpi_finalize
  public :: mpi_append_fname
  public :: mpi_is_halo
  ! -- Public variables
  public :: nhalo
  public :: modelname_halo

  ! -- Global variables based of the world communicator
  logical, public :: parallelrun = .false.
  logical, public :: serialrun   = .true.
  logical, public :: writestd    = .true.
  character(len=50), public :: partstr
  
  integer(I4B) :: nhalo = 0
  character(len=LENMODELNAME), allocatable, dimension(:) :: modelname_halo !PAR
  character(len=LINELENGTH) :: errmsg
  
  save
  
  contains
  
  subroutine mpi_initialize()
! ******************************************************************************
! MPI initialization to be called at the start of the program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
! ------------------------------------------------------------------------------
    !
    call mpiwrpinit()
    !
    ! -- return
    return
  end subroutine mpi_initialize

  subroutine mpi_finalize()
! ******************************************************************************
! MPI finalization to be called at the end of the program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
! ------------------------------------------------------------------------------
    !
    call mpiwrpfinalize()
    !
    ! -- return
    return
  end subroutine mpi_finalize

  subroutine mpi_append_fname(f)
! ******************************************************************************
! Append the file name with the process rank ID.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    character(len=*), intent(inout) :: f
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Return in case of serial run
    if (serialrun) return
    !
    write(f,'(3a)') trim(f),'.',trim(partstr)
    !
    ! -- return
    return
  end subroutine mpi_append_fname

  subroutine mpi_create_modelname_halo(im, modelname)
! ******************************************************************************
! Create halo model name.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    character(len=*), intent(inout) :: modelname
    integer(I4B), intent(in) :: im
    ! -- local
    character(len=LINELENGTH) :: s
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return 
    end if
    !
    write(s,*) im
    modelname = trim(modelname)//' HALO'//trim(adjustl(s))
    !
    ! -- return
    return
  end subroutine mpi_create_modelname_halo

  function mpi_is_halo(modelname) result(flag_halo)
! ******************************************************************************
! Check if a model name is of type halo
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- return
    logical :: flag_halo
    ! -- dummy
    character(len=*), intent(in) :: modelname
    ! -- local
    integer(I4B) :: m, i
! ------------------------------------------------------------------------------
    !
    flag_halo = .false.
    if (serialrun) then
      return 
    end if
    !
    ! -- check
    if (.not.allocated(modelname_halo)) then
      !write(errmsg,'(a)') 'Program error in mpi_is_halo.'
      !call store_error(errmsg)
      !call ustop()
    end if
    !
    m = -1
    findloop: do i=1,size(modelname_halo)
      if(modelname_halo(i) == modelname) then
        m = i
        exit findloop
      endif
    enddo findloop
    !
    if (m > 0) then
      flag_halo = .true.
    end if
    !
    ! -- return
    return
  end function mpi_is_halo
  
end module MpiExchangeGenModule

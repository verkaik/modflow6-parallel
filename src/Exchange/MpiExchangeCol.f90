Module MpiExchangeColModule
  use MpiExchangeGenModule, only: serialrun
  use MpiWrapper, only: ColMemoryType
  
  implicit none

  private

    ! -- Public functions
  public :: mpi_get_distype
  public :: ColMemoryType
  public :: ciopt, n_recv, cmt_recv

  ! -- Collective receive buffer
  integer :: ciopt = 0
  integer :: n_recv
  type(ColMemoryType), dimension(:), allocatable :: cmt_recv

  save

contains

  subroutine mpi_get_distype(modelname, ldis, ldisu, ldisv)
! ******************************************************************************
! This subroutine set the discretization type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENMODELNAME
    use SimModule, only: store_error, store_error_unit, ustop
    ! -- dummy
    character(len=*), intent(in) :: modelname
    logical, intent(out) :: ldis, ldisu, ldisv
    ! -- local
    character(len=LINELENGTH) :: errmsg
    character(len=1) :: cdum
    character(len=LENMODELNAME) :: m, mh
    character(len=4) :: dis_type
    logical :: lfound
    integer :: i
! ------------------------------------------------------------------------------
    ldis  = .false.
    ldisu = .false.
    ldisv = .false.
    !
    if (serialrun .or. ciopt /= 1) then
      ldis = .true. !@@@@@@@ DEBUG
      return
    endif
    !
    read(modelname,*) mh
    !
    lfound = .false.
    do i = 1, n_recv
      read(cmt_recv(i)%origin,*) m, dis_type
      if (m == mh) then
        select case(dis_type)
          case ('DIS')
            ldis = .true.
          case ('DISU')
            ldisu = .true.
          case ('DISV')
            ldisv = .true.
        end select    
        lfound = .true.
        exit
      endif   
    enddo
    ! -- check
    if (.not.ldis .and. .not.ldisu .and. .not.ldisv) then
      write(errmsg,'(a)') 'Program error in mpi_get_distype.'
      call store_error(errmsg)
      call ustop()
    endif
    !
    ! -- return
    return
  end subroutine mpi_get_distype

end module MpiExchangeColModule


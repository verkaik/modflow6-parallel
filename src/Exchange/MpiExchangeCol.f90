Module MpiExchangeColModule
  use MpiExchangeGenModule, only: serialrun
  use MpiWrapper, only: ColMemoryType
  
  implicit none

  private

    ! -- Public functions
  public :: mpi_get_distype
  public :: mpi_get_distype_str
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
! This subroutine sets the discretization type.
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
    character(len=LENMODELNAME) :: m, mh
    character(len=4) :: dis_type
    logical :: lfound
    integer :: i, ios
! ------------------------------------------------------------------------------
    ldis  = .false.
    ldisu = .false.
    ldisv = .false.
    !
    if (serialrun) then
!    if (serialrun .or. ciopt /= 1) then
      ldis = .true.
      return
    end if
    !
    read(modelname,*) mh
    !
    lfound = .false.
    do i = 1, n_recv
      read(cmt_recv(i)%origin,*,iostat=ios) m, dis_type
      if (m == mh .and. ios == 0) then
        select case(dis_type)
          case ('DIS')
            ldis = .true.
            lfound = .true.
          case ('DISU')
            ldisu = .true.
            lfound = .true.
          case ('DISV')
            ldisv = .true.
            lfound = .true.
        end select
        if(lfound) then  
          exit
        end if
      end if
    end do
    ! -- check
    if (.not.lfound) then
      write(errmsg,'(a)') 'Program error in mpi_get_distype.'
      call store_error(errmsg)
      call ustop()
    end if
    !
    ! -- return
    return
  end subroutine mpi_get_distype

  subroutine mpi_get_distype_str(modelname, disstr)
! ******************************************************************************
! This subroutine sets the discretization type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENMODELNAME
    use SimModule, only: store_error, store_error_unit, ustop
    ! -- dummy
    character(len=*), intent(in) :: modelname
    character(len=*), intent(out) :: disstr
    ! -- local
    character(len=LINELENGTH) :: errmsg
    character(len=LENMODELNAME) :: m, mh
    character(len=4) :: dis_type
    logical :: lfound
    integer :: i, ios
! ------------------------------------------------------------------------------
     disstr = ''
    !
    if (serialrun) then
!    if (serialrun .or. ciopt /= 1) then
      disstr = 'DIS'
      return
    end if
    !
    read(modelname,*) mh
    !
    lfound = .false.
    do i = 1, n_recv
      read(cmt_recv(i)%origin,*,iostat=ios) m, dis_type
      if (m == mh .and. ios == 0) then
        select case(dis_type)
          case ('DIS')
            lfound = .true.
          case ('DISU')
            lfound = .true.
          case ('DISV')
            lfound = .true.
        end select
        if(lfound) then  
          disstr = dis_type
          exit
        end if
      end if
    end do
    ! -- check
    if (.not.lfound) then
      write(errmsg,'(a)') 'Program error in mpi_get_distype_str.'
      call store_error(errmsg)
      call ustop()
    end if
    !
    ! -- return
    return
  end subroutine mpi_get_distype_str 
  
end module MpiExchangeColModule


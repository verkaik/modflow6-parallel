module MpiMvrModule
  
  use KindModule, only: DP, I4B
  use MpiExchangeGenModule, only: serialrun
  use MvrModule, only: MvrType

  implicit none
  
  type MpiMvrType
    integer(I4B), pointer                :: maxmvr => null() !max number of movers to be specified
    integer(I4B), pointer                :: ngmvr  => null()  !number of movers for current stress period (global)
    type(MvrType), pointer, dimension(:) :: gmvr   => null()  !array of movers (global)
  contains
    procedure :: mpi_init_ngmvr
    procedure :: mpi_set_maxmvr
    procedure :: mpi_set_mover
    procedure :: mpi_mvr_da
  end type MpiMvrType
  
  save
  
  private

  ! -- Public functions
  public :: MpiMvrType
  
  contains
  
  subroutine mpi_init_ngmvr(this)
! ******************************************************************************
! Initialize the mover counter
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiMvrType) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    endif
    !
    this%ngmvr = 0
    !
    ! -- return
    return
  end subroutine mpi_init_ngmvr
  
  subroutine mpi_set_maxmvr(this, maxmvr)
! ******************************************************************************
! This subroutine set the maximum of movers.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiMvrType) :: this
    integer(I4B), intent(in) :: maxmvr
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    endif
    !
    ! -- allocate
    allocate(this%maxmvr)
    allocate(this%gmvr(maxmvr))
    allocate(this%ngmvr)
    !
    ! -- set values
    this%maxmvr = maxmvr
    !
    ! -- return
    return
  end subroutine mpi_set_maxmvr
  
  subroutine mpi_set_mover(this, pname1, pname2)
! ******************************************************************************
! This subroutine set the maximum of movers.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiMvrType) :: this
    character(len=*), intent(in) :: pname1, pname2
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    endif
    !
    this%ngmvr = this%ngmvr + 1
    i = this%ngmvr
    this%gmvr(i)%pname1 = trim(pname1)
    this%gmvr(i)%pname2 = trim(pname2)
    !
    ! -- return
    return
  end subroutine mpi_set_mover
  
  subroutine mpi_mvr_da(this)
! ******************************************************************************
! Deallocate mover data.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiMvrType) :: this
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    endif
    
    deallocate(this%maxmvr)
    deallocate(this%gmvr)
    deallocate(this%ngmvr)
    
    ! -- return
    return
  end subroutine mpi_mvr_da
  
end module MpiMvrModule

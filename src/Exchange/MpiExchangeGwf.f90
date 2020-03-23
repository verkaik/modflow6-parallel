module MpiExchangeGwfModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENMODELNAME, LENORIGIN, LENVARNAME, LENFTYPE,    &
                             LENPACKAGENAME
  use ArrayHandlersModule, only: ifind
  use MpiExchangeGenModule, only: serialrun
  use MpiWrapper, only: ColMemoryType, mpiwrpallgather, mpiwrpallgatherv,      &
                        mpiwrpbarrier, mpiwrpcolstruct, mpiwrptypefree
  use MpiExchangeModule, only: MpiWorld, mpi_to_colmem
  use SimModule, only: store_error, store_error_unit, ustop ! @@@@@ debug
  !
  use MawModule, only: mawftype => ftype
  use SfrModule, only: sfrftype => ftype
  use LakModule, only: lakftype => ftype
  use UzfModule, only: uzfftype => ftype
  !
  implicit none

  integer, parameter :: imaw = 1
  integer, parameter :: isfr = 2
  integer, parameter :: ilak = 3
  integer, parameter :: iuzf = 4
  integer, parameter :: mxpack = iuzf
  
  character(LENFTYPE), dimension(2,mxpack), save :: packftype
  
  private
  
  ! -- Public functions
  public :: mpi_gwfhalo_world
  !public :: mpi_set_gwfhalo_world
  public :: mpi_set_gwfhalo_world_var_int
  public :: mpi_set_gwfhalo_world_var_int1d
  public :: mpi_set_gwfhalo_world_dis
  public :: mpi_set_gwfhalo_world_mvr
  
  !interface mpi_set_gwfhalo_world_var
  !  module procedure mpi_set_gwfhalo_world_var_int,                             &
  !                   mpi_set_gwfhalo_world_var_int1d
  !end interface mpi_set_gwfhalo_world_var
   
  contains
  
  subroutine mpi_gwfhalo_world(iopt)
! ******************************************************************************
! This subroutine gathers DIS information for all models.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeColModule, only: ciopt, n_recv, cmt_recv
    use MemoryManagerModule, only: mem_setptr, mem_get_ptr, mem_setval
    use MemoryTypeModule, only: MemoryType
    use BaseModelModule, only: BaseModelType, GetBaseModelFromList
    use ListsModule, only: basemodellist
    use NumericalModelModule, only: NumericalModelType
    use BndModule, only: BndType, GetBndFromList
    use GwfModule, only: GwfModelType, cunit, niunit
    
    ! -- dummy
    integer, intent(in) :: iopt
    ! -- local
    character(len=1) :: cdum
    character(len=4) :: dis_type
    integer :: im, ipo, ip, is, iact
    class(BaseModelType), pointer :: mb 
    class(NumericalModelType), pointer :: mp
    class(GwfModelType), pointer :: gp
    class(BndType), pointer :: packobj
    type(MemoryType), pointer :: mt
    !
    integer, parameter :: cntypes = 6
    integer :: i, n, n_send, newtype
    integer, dimension(1) :: iwrk
    integer, dimension(:), allocatable :: recvcounts, displs
    type(ColMemoryType), dimension(:), allocatable :: cmt_send
    character(len=LENFTYPE) :: ft, cft
    integer :: rank
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    end if
    !
    ! -- initialize package filetype
    packftype(1,imaw) = trim(mawftype)
    packftype(1,isfr) = trim(sfrftype)
    packftype(1,ilak) = trim(lakftype)
    packftype(1,iuzf) = trim(uzfftype)
    !
    do ip = 1, mxpack
      ft = packftype(1,ip)
      n = len_trim(ft)
      do i = 1, niunit
        cft = cunit(i)
        if (ft(1:n) == cft(1:n)) then
          packftype(2,ip) = trim(cunit(i))
        end if
      end do
    end do
    !
    allocate(cmt_send(1))
    do iact = 1, 2
      is = 0
      do im = 1, basemodellist%Count()
        mb => GetBaseModelFromList(basemodellist, im)
        select type (mb)
        class is (NumericalModelType)
          mp => mb
        end select  
        select type (mb)
        class is (GwfModelType)
          gp => mb
        end select  
        read(mp%dis%origin,*) cdum, dis_type
        select case (dis_type)
          case ('DIS')
            if (iopt == 1) then
              call mem_get_ptr('DNDIM', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            else
              call mem_get_ptr('DNDIM', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NLAY', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NROW', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NCOL', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('MSHAPE', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NODESUSER', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            end if
          case ('DISV')
            if (iopt == 1) then
              call mem_get_ptr('DNDIM',  mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            else
              call mem_get_ptr('DNDIM',  mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NLAY',  mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NCPL',  mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NVERT', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('MSHAPE', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NODESUSER', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            end if
          case ('DISU')
            if (iopt == 1) then
              call mem_get_ptr('DNDIM', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            else  
              call mem_get_ptr('DNDIM', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NVERT', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('MSHAPE', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NODESUSER', mp%dis%origin, mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            end if
        end select
        if (iopt == 2) then
          ! -- Mover package data
          do ipo=1,mp%bndlist%Count()
            packobj => GetBndFromList(mp%bndlist, ipo)
            if (packobj%filtyp == packftype(1,imaw)) then
              call mem_get_ptr('NMAWWELLS', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('IMOVER', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            end if
            if (packobj%filtyp == packftype(1,isfr)) then
              call mem_get_ptr('MAXBOUND', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('IMOVER', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            end if
            if (packobj%filtyp == packftype(1,ilak)) then
              call mem_get_ptr('NOUTLETS', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('NLAKES', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('IMOVER', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            end if
            if (packobj%filtyp == packftype(1,iuzf)) then
              call mem_get_ptr('MAXBOUND', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
              call mem_get_ptr('IMOVER', trim(packobj%origin), mt)
              call mpi_to_colmem(mt, is, cmt_send, iact)
            end if
          end do
          call mem_get_ptr('INNPF', mp%name, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
          call mem_get_ptr('INIC', mp%name, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
        endif
        if (iopt == 3) then
          ! -- NPF related
          call mem_get_ptr('IK22', gp%npf%origin, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
          call mem_get_ptr('IK33', gp%npf%origin, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
          call mem_get_ptr('IANGLE1', gp%npf%origin, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
          call mem_get_ptr('IANGLE2', gp%npf%origin, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
          call mem_get_ptr('IANGLE3', gp%npf%origin, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
          call mem_get_ptr('IREWET', gp%npf%origin, mt)
          call mpi_to_colmem(mt, is, cmt_send, iact)
        end if
      end do
      if (iact == 1) then
        n_send = is
        if (allocated(cmt_send)) then
          deallocate(cmt_send)
        end if
        allocate(cmt_send(max(n_send,1))) !PAR
      end if
    end do
    !
    ! -- gather sizes
    allocate(recvcounts(MpiWorld%nrproc), displs(MpiWorld%nrproc))
    iwrk(1) = n_send
    call mpiwrpallgather(MpiWorld%comm, iwrk, 1, recvcounts, 1)
    n_recv = sum(recvcounts)
    if (allocated(cmt_recv)) then
      deallocate(cmt_recv)
    end if
    allocate(cmt_recv(n_recv))
    displs(1) = 0
    do i = 2, MpiWorld%nrproc
      displs(i) = displs(i-1) + recvcounts(i-1)
    end do
    ! -- create derived type
    call mpiwrpcolstruct(newtype)
    ! -- gather data
    call mpiwrpallgatherv(MpiWorld%comm, cmt_send, n_send, newtype, cmt_recv,   &
                          recvcounts, newtype, displs)
    ! -- DEBUG
    if (.false.) then
      do rank = 0, MpiWorld%nrproc-1
        if (MpiWorld%myrank == rank) then
          write(*,*) '=== send myrank', rank
          do i = 1, n_send
            write(*,'(a,1x,a,1x,i)') trim(cmt_send(i)%name), trim(cmt_send(i)%origin), cmt_send(i)%intsclr  
          end do
          write(*,*) '=== recv myrank', rank
          do i = 1, n_recv
            write(*,'(a,1x,a,1x,i)') trim(cmt_recv(i)%name), trim(cmt_recv(i)%origin), cmt_recv(i)%intsclr 
          end do
        end if
        call MpiWorld%mpi_barrier()
      end do
      call ustop('@@@@')
    end if
    !
    ! -- clean up
    call mpiwrptypefree(newtype)
    deallocate(cmt_send, recvcounts, displs)
    ciopt = iopt
    !
    ! -- return
    return
  end subroutine mpi_gwfhalo_world
  
  function mpi_set_gwfhalo_world_var_int(tgtname, tgtorigin) result(aint)
! ******************************************************************************
  
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeColModule, only: n_recv, cmt_recv
    use MemoryTypeModule, only: iintsclr
    ! -- dummy
    character(len=*), intent(in) :: tgtname
    character(len=*), intent(in) :: tgtorigin
    ! -- return
    integer(I4B) :: aint
    ! -- local
    integer(I4B) :: i
    logical :: lok
    character(len=LENVARNAME) :: name        !name of the array
    character(len=LENORIGIN)  :: origin      !name of origin
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    ! -- Set other variables
    aint = -1
    lok = .false.
    do i = 1, n_recv
      name = cmt_recv(i)%name
      read(cmt_recv(i)%origin,*) origin
      !write(*,*) '@@@@ name, origin: "'//trim(cmt_recv(i)%name)//'" "'//trim(cmt_recv(i)%origin)//'"'
      if ((trim(name) == trim(tgtname)) .and. (trim(origin) == trim(tgtorigin))) then
        lok = .true.
        if (cmt_recv(i)%memitype /= iintsclr) then
          lok = .false.
        else
          aint = cmt_recv(i)%intsclr
        end if
        exit
      end if
    end do
    !
    ! -- check
    if (.not.lok) then
      call store_error('Program error: mpi_set_gwfhalo_world_var_int')
      write(*,*) '@@@ not found: "',trim(tgtname),'" "',trim(tgtorigin),'"'
      call ustop()
    end if
    ! -- return
    return
  end function mpi_set_gwfhalo_world_var_int
  
  function mpi_set_gwfhalo_world_var_int1d(tgtname, tgtorigin) result(aint1d)
! ******************************************************************************
  
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeColModule, only: n_recv, cmt_recv
    use MemoryTypeModule, only: iaint1d
    ! -- dummy
    character(len=*), intent(in) :: tgtname
    character(len=*), intent(in) :: tgtorigin
    ! -- return
    integer(I4B), dimension(3) :: aint1d
    ! -- local
    integer(I4B) :: i
    logical :: lok
    character(len=LENVARNAME) :: name        !name of the array
    character(len=LENORIGIN)  :: origin      !name of origin
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    ! -- Set other variables
    aint1d = -1
    lok = .false.
    do i = 1, n_recv
      name = cmt_recv(i)%name
      read(cmt_recv(i)%origin,*) origin
      if ((trim(name) == trim(tgtname)) .and. (trim(origin) == trim(tgtorigin))) then
        lok = .true.
        if (cmt_recv(i)%memitype /= iaint1d) then
          lok = .false.
        else
          aint1d(1) = cmt_recv(i)%aint1d(1)
          aint1d(2) = cmt_recv(i)%aint1d(2)
          aint1d(3) = cmt_recv(i)%aint1d(3)
        end if
        exit
      end if
    end do
    !
    ! -- check
    if (.not.lok) then
      call store_error('Program error: mpi_set_gwfhalo_world_var_int')
      write(*,*) '@@@ not found: "',trim(tgtname),'" "',trim(tgtorigin),'"'
      call ustop()
    end if
    !
    ! -- return
    return
  end function mpi_set_gwfhalo_world_var_int1d
  !
  subroutine mpi_set_gwfhalo_world_dis(mname, m2)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use GwfDisModule,          only: dis_cr,  GwfDisType !,dis_init_mem TODO
    use GwfDisvModule,         only: disv_cr, GwfDisvType
    use GwfDisuModule,         only: disu_cr, GwfDisuType
    use NumericalModelModule,   only: NumericalModelType
    use BaseDisModule,         only: DisBaseType !PAR
    use MpiExchangeColModule,  only: mpi_get_distype
    use MemoryManagerModule,   only: mem_setval, mem_allocate, mem_setval
    use ConstantsModule, only: IZERO, DZERO
    ! -- dummy
    character(len=*), intent(in) :: mname
    class(NumericalModelType), pointer, intent(out) :: m2
    ! -- local
    integer(I4B) :: in_dum, iout_dum
    integer(I4B) :: ndim, nodesuser
    integer(I4B), dimension(3) :: mshape
    logical ::  ldis, ldisu, ldisv
    character(len=100) :: disname
    class(DisBaseType), pointer :: bp
    class(GwfDisType), pointer :: disp
    class(GwfDisvType), pointer :: disvp
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    disname = 'M2'
    !
    allocate(m2)
    allocate(m2%dis)
    in_dum = -1
    iout_dum = -1
    call mpi_get_distype(mname, ldis, ldisu, ldisv)
    if (ldis) then 
      call dis_cr(m2%dis, disname, in_dum, iout_dum)
      disname = trim(disname)//' DIS'
      bp => m2%dis
      select type (bp)
      class is (GwfDisType)
        disp => bp
      end select
      call mem_allocate(disp%nodereduced, 1, 'NODEREDUCED', disname)
      call mem_allocate(disp%idomain, 1, 1, 1, 'IDOMAIN', disname)
      call mem_allocate(disp%top2d, 1, 1, 'TOP2D', disname)
      call mem_allocate(disp%bot3d, 1, 1, 1, 'BOT3D', disname)
    endif  
    if (ldisu) then 
      call disu_cr(m2%dis, disname, in_dum, iout_dum)
      disname = trim(disname)//' DISU'
    endif  
    if (ldisv) then 
      call disv_cr(m2%dis, disname, in_dum, iout_dum)
      disname = trim(disname)//' DISV'
      bp => m2%dis
      select type (bp)
      class is (GwfDisvType)
        disvp => bp
      end select
      call mem_allocate(disvp%nodereduced, 1, 'NODEREDUCED', disname)
      call mem_allocate(disvp%idomain, 1, 1, 1, 'IDOMAIN', disname)
      call mem_allocate(disvp%top2d, 1, 1, 'TOP2D', disname)
      call mem_allocate(disvp%bot3d, 1, 1, 1, 'BOT3D', disname)
    endif  
    !
    ndim = mpi_set_gwfhalo_world_var_int('DNDIM', mname)
    call mem_setval(ndim, 'NDIM', disname)
    !
    nodesuser = mpi_set_gwfhalo_world_var_int('NODESUSER', mname)
    call mem_setval(nodesuser, 'NODESUSER', disname)
    call mem_setval(nodesuser, 'NODES', disname)
    !
    mshape = mpi_set_gwfhalo_world_var_int1d('MSHAPE', mname)
    call mem_allocate(m2%dis%mshape, ndim, 'MSHAPE', disname)
    call mem_setval(mshape, 'MSHAPE', disname)
    !
    ! -- return
    return
  end subroutine mpi_set_gwfhalo_world_dis
  
!  subroutine mpi_set_gwfhalo_world()
!! ******************************************************************************
!! This subroutine sets the DIS scalars for the halo (m2) GWF models. It also
!! initializes halo GWF packages.
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- modules
!!    use MpiExchangeGenModule, only: modelname_halo, nhalo,                      &
!!                                    mpi_create_modelname_halo
!    use MpiExchangeColModule, only: n_recv, cmt_recv
!    use MemoryTypeModule, only: ilogicalsclr, iintsclr, idblsclr,               & 
!                                iaint1d, iaint2d,                               & 
!                                iadbl1d, iadbl2d
!    use MemoryManagerModule, only: mem_setval
!    use BaseModelModule, only: BaseModelType, GetBaseModelFromList    
!    use ListsModule, only: halomodellist
!    use GwfModule, only: GwfModelType
!    use NumericalModelModule, only: NumericalModelType
!    use BndModule, only: BndType, GetBndFromList
!    
!    ! -- dummy
!    ! -- local
!    character(len=LENVARNAME)   :: name        !name of the array
!    character(len=LENORIGIN)    :: origin      !name of origin
!    character(len=LENORIGIN)    :: origin_halo !name of origin
!    character(len=LENMODELNAME) :: mname       !name of origin
!    character(len=LENMODELNAME) :: mname_halo  !name of origin
!    character(len=LENMODELNAME) :: id
!    character(len=LENFTYPE)     :: ft
!    character(len=LENPACKAGENAME) :: pakname
!    integer(I4B) :: i, j, m, ip, istat, lenid, lenft, ipo, ibcnum
!    logical, dimension(mxpack) :: lpack, lcreate
!    class(BaseModelType), pointer :: mb
!    class(GwfModelType), pointer :: mp
!    class(NumericalModelType), pointer :: nmp
!    class(BndType), pointer :: packobj
!! ------------------------------------------------------------------------------
!    !
!   
!    if (serialrun) then
!      return
!    end if
!    !
!    call mpi_gwfhalo_world(2)
!    !
!    do i = 1, n_recv
!      name   = cmt_recv(i)%name
!      origin = cmt_recv(i)%origin
!      read(origin,*,iostat=istat) mname, id
!      if (istat /= 0) then
!        read(origin,*,iostat=istat) mname
!        id = ''
!      end if
!      lenid = len_trim(id)
!      ! -- check if the data is related to packages
!      lpack = .false.
!      if (lenid > 0) then
!        do ip = 1, mxpack
!          ft = packftype(1,ip)
!          lenft = len_trim(ft)
!          if (lenid >= lenft) then
!            if (id(1:lenft) == ft(1:lenft)) then
!              lpack(ip) = .true.
!            end if
!        end if
!        end do
!        if (any(lpack)) then
!          pakname = ''
!          j = index(origin,'-',back=.true.)
!          read(origin(j+1:),*,iostat=istat) ibcnum
!          if (istat /= 0) then
!            ibcnum = 0
!            pakname = trim(id)
!          end if
!        end if
!      end if
!      if (any(lpack)) then
!        do j = 1, nhalo
!          mname_halo = mname
!          call mpi_create_modelname_halo(j, mname_halo)
!          m = ifind(modelname_halo, mname_halo)
!          ! -- halo model found
!          if (m > 0) then
!            ! -- set pointers
!            mb => GetBaseModelFromList(halomodellist, m)
!            select type (mb)
!            class is (GwfModelType)
!              mp => mb
!            end select
!            select type (mb)
!            class is (NumericalModelType)
!              nmp => mb
!            end select
!            ! -- check if package should be created
!            lcreate = .false.
!            do ip = 1, mxpack
!              if (lpack(ip)) then
!                lcreate(ip) = .true.
!              end if
!            end do
!            do ipo=1,nmp%bndlist%Count()
!              packobj => GetBndFromList(nmp%bndlist, ipo)
!              do ip = 1, mxpack
!                if (lpack(ip) .and. packobj%filtyp == packftype(1,ip) .and.    &
!                  packobj%ibcnum == ibcnum) then
!                  lcreate(ip) = .false.
!                end if
!              end do
!            end do
!            ! -- create package
!            do ip = 1, mxpack
!              if (lcreate(ip)) then
!                call mp%package_create(trim(packftype(2,ip)), 0, ibcnum,        &
!                                       pakname, 0, 0)
!                packobj => GetBndFromList(nmp%bndlist, nmp%bndlist%Count())
!                call mem_setval(.true., 'P_ISHALO', packobj%origin)
!              end if
!            end do
!          end if
!        end do
!      end if
!    end do
!    !
!    ! -- Set other variables
!    do i = 1, n_recv
!      name   = cmt_recv(i)%name
!      origin = cmt_recv(i)%origin
!      read(origin,*,iostat=istat) mname, id
!      if (istat /= 0) then
!        read(origin,*,iostat=istat) mname
!        id = ''
!      end if
!      !do j = 1, nhalo
!      !  mname_halo = mname
!      !  call mpi_create_modelname_halo(j, mname_halo)
!      !  m = ifind(modelname_halo, mname_halo)
!      !  if (m > 0) then
!      !    origin_halo = trim(mname_halo)//' '//trim(id)
!      !    select case(cmt_recv(i)%memitype)
!      !    case(iintsclr)
!      !        call mem_setval(cmt_recv(i)%intsclr, name, origin_halo)
!      !      case(iaint1d)
!      !        call mem_setval(cmt_recv(i)%aint1d, name, origin_halo)
!      !    end select
!      !  end if
!      !end do
!    end do
!    !
!    ! -- cleanup
!    deallocate(cmt_recv)
!    !
!    ! -- return
!    return
!  end subroutine mpi_set_gwfhalo_world
  
  subroutine mpi_set_gwfhalo_world_mvr()
! ******************************************************************************
! This subroutine sets halo model packages for mover
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iintsclr, iaint1d
    use MpiExchangeColModule, only: n_recv, cmt_recv
    use ArrayHandlersModule, only: ifind
    use ListsModule, only: halomodellist
    use BaseModelModule, only: BaseModelType, GetBaseModelFromListByName
    use NumericalModelModule, only: NumericalModelType
    use BndModule, only: BndType, GetBndFromList
    use GwfModule, only: GwfModelType
    use MemoryManagerModule, only: mem_setval
    ! -- local
    character(len=LENVARNAME)   :: name        !name of the array
    character(len=LENORIGIN)    :: origin      !name of origin
    character(len=LENORIGIN)    :: origin_halo !name of origin
    character(len=LENMODELNAME) :: mname       !name of model
    character(len=LENMODELNAME) :: hmname      !name of halo model
    character(len=LENMODELNAME) :: id
    character(len=LENFTYPE)     :: ft
    character(len=LENPACKAGENAME) :: pakname
    integer(I4B) :: i, j, ih, ip, istat, lenid, lenft, ipo, ibcnum
    logical, dimension(mxpack) :: lpack, lcreate
    class(BaseModelType), pointer :: mb
    class(GwfModelType), pointer :: mp
    class(NumericalModelType), pointer :: nmp
    class(BndType), pointer :: packobj
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    do i = 1, n_recv
      name   = cmt_recv(i)%name
      origin = cmt_recv(i)%origin
      read(origin,*,iostat=istat) mname, id
      if (istat /= 0) then
        read(origin,*,iostat=istat) mname
        id = ''
      end if
      lenid = len_trim(id)
      ! -- check if the data is related to packages
      lpack = .false.
      if (lenid > 0) then
        do ip = 1, mxpack
          ft = packftype(1,ip)
          lenft = len_trim(ft)
          if (lenid >= lenft) then
            if (id(1:lenft) == ft(1:lenft)) then
              lpack(ip) = .true.
            end if
          end if
        end do
        if (any(lpack)) then
          pakname = ''
          j = index(origin,'-',back=.true.)
          read(origin(j+1:),*,iostat=istat) ibcnum
          if (istat /= 0) then
            ibcnum = 0
            pakname = trim(id)
          end if
        end if
        if (any(lpack)) then
          do ih = 1, MpiWorld%hnmodel
            if (trim(MpiWorld%hmodelm2names(ih)) == trim(mname)) then
              hmname = MpiWorld%hmodelnames(ih)
              mb => GetBaseModelFromListByName(halomodellist, hmname)
              ! -- cast model
              mp => null()
              select type (mb)
              class is (GwfModelType)
                mp => mb
              end select
              nmp => null()
              select type (mb)
              class is (NumericalModelType)
                nmp => mb
              end select
              if((.not.associated(mp)).or.(.not.associated(nmp))) then
                call ustop('Program error 1 mpi_set_gwfhalo_world_mvr')
              end if
              ! -- check if package should be created
              lcreate = .false.
              do ip = 1, mxpack
                if (lpack(ip)) then
                  lcreate(ip) = .true.
                end if
              end do
              do ipo=1,nmp%bndlist%Count()
                packobj => GetBndFromList(nmp%bndlist, ipo)
                do ip = 1, mxpack
                  if (lpack(ip) .and. packobj%filtyp == packftype(1,ip) .and.&
                    packobj%ibcnum == ibcnum) then
                    lcreate(ip) = .false.
                  end if
                end do
              end do
              ! -- create package
              do ip = 1, mxpack
                if (lcreate(ip)) then
                  call mp%package_create(trim(packftype(2,ip)), 0, ibcnum,   &
                                         pakname, 0, 0)
                  packobj => GetBndFromList(nmp%bndlist, nmp%bndlist%Count())
                  call mem_setval(.true., 'P_ISHALO', packobj%origin)
                end if
              end do
            end if
          end do
        end if
      end if
    end do
    !
    ! -- Set other variables
    do i = 1, n_recv
      name   = cmt_recv(i)%name
      origin = cmt_recv(i)%origin
      read(origin,*,iostat=istat) mname, id
      if (istat /= 0) then
        read(origin,*,iostat=istat) mname
        id = ''
      end if
      lenid = len_trim(id)
      lpack = .false.
      if (lenid > 0) then
        do ip = 1, mxpack
          ft = packftype(1,ip)
          lenft = len_trim(ft)
          if (lenid >= lenft) then
            if (id(1:lenft) == ft(1:lenft)) then
              lpack(ip) = .true.
            end if
          end if
        end do
      end if
      !
      if (any(lpack)) then
        do ih = 1, MpiWorld%hnmodel
          if (trim(MpiWorld%hmodelm2names(ih)) == trim(mname)) then
            hmname = MpiWorld%hmodelnames(ih)
            origin_halo = trim(hmname)//' '//trim(id)
            select case(cmt_recv(i)%memitype)
              case(iintsclr)
                call mem_setval(cmt_recv(i)%intsclr, name, origin_halo)
              case(iaint1d)
                call mem_setval(cmt_recv(i)%aint1d, name, origin_halo)
            end select
          end if
        end do
      end if
    end do
    !
    ! -- return
    return
  end subroutine mpi_set_gwfhalo_world_mvr
  
end module MpiExchangeGwfModule
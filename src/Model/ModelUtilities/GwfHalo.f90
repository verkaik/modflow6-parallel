module GwfHaloModule
  
  use KindModule,         only: DP, I4B
  use GwfModule,          only: GwfModelType
  use ConstantsModule,    only: LINELENGTH, DZERO, LENMODELNAME, LENSOLUTIONNAME
  use ConnectionsModule,  only: ConnectionsType
  use MpiExchangeModule,  only: MpiExchangeType
  use SimModule, only: store_error, store_error_unit, ustop  
  
  implicit none

  private
  public :: gwfhalo_cr
  public :: GwfHaloModelType
  
  type, extends(GwfModelType) :: GwfHaloModelType
    integer(I4B), pointer :: nband
    type(MpiExchangeType), pointer :: MpiSol => null() !PAR
    type(GwfModelType), pointer :: gwf1 => null()
    type(GwfModelType), pointer :: gwf2 => null()
    character(len=LENMODELNAME), pointer :: m1name => null()        ! name of GWF Model 1
    character(len=LENMODELNAME), pointer :: m2name => null()        ! name of GWF Model 2
    integer(I4B), pointer :: nexg => null()
    integer(I4B), pointer :: m1ndim => null()
    integer(I4B), pointer :: m2ndim => null()
    integer(I4B), pointer :: m1irewet => null()
    integer(I4B), pointer :: m2irewet => null()
    integer(I4B), dimension(:), pointer :: nodem1 => null()
    integer(I4B), dimension(:), pointer :: nodem2 => null()
    type(ConnectionsType), pointer :: m1con => null()
    type(ConnectionsType), pointer :: m2con => null()
    integer(I4B), dimension(:), pointer :: m1nbnod => null()
    integer(I4B), dimension(:), pointer :: m2nbnod => null()
    integer(I4B), pointer :: offset => null()
    integer(I4B), dimension(:), pointer :: imapnodem1tohalo => null()
    integer(I4B), dimension(:), pointer :: imapnodem2tohalo => null()
    integer(I4B), dimension(:), pointer :: imapm1tohalo => null()
    integer(I4B), dimension(:), pointer :: imapm2tohalo => null()
    integer(I4B), dimension(:), pointer :: m1nodes => null()
    integer(I4B), dimension(:), pointer :: m2nodes => null()
  contains
  
    procedure :: gwfhalo_df1
    procedure :: gwfhalo_df2
    procedure :: gwfhalo_df3
    procedure :: gwfhalo_ar
    procedure :: gwfhalo_cf1
    procedure :: gwfhalo_cf2
    procedure :: gwfhalo_fc_calc
    procedure :: gwfhalo_fn_calc
    
  end type GwfHaloModelType
  
  contains

  subroutine gwfhalo_cr(this, id, modelname, m1name, m2name)
! ******************************************************************************
! gwfhalo_cr -- Create a new groundwater flow model object for exchange
! Subroutine: (1) creates model object and add to exchange modellist
!             (2) assign values
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use GwfDisuModule,              only: disu_cr
    use GwfNpfModule,               only: npf_cr
    use Xt3dModule,                 only: xt3d_cr
    use GhostNodeModule,            only: gnc_cr
    use GwfMvrModule,               only: mvr_cr
    use GwfIcModule,                only: ic_cr
    use GwfObsModule,               only: gwf_obs_cr
    use MemoryHelperModule, only: create_mem_path    
    ! -- dummy
    type(GwfHaloModelType), pointer   :: this
    integer(I4B), intent(in)      :: id
    character(len=*), intent(in)  :: modelname
    character(len=*), intent(in)  :: m1name
    character(len=*), intent(in)  :: m2name
    ! -- local
    integer(I4B), target :: iout_dum
    integer(I4B), pointer :: in_dum
    ! -- format
! ------------------------------------------------------------------------------
    !
    ! -- Allocate a new GWF Model (this)
    allocate(this)
    
    ! -- Set memory path before allocation in memory manager can be done
    this%memoryPath = create_mem_path(modelname)
    
    call this%allocate_scalars(modelname)
    !
    ! -- Assign values
    this%name = modelname
    this%macronym = 'GWF'
    this%id = id
    this%ishalo = .true.
    !
    ! -- TODO: to be replaced by optional arguments
    allocate(in_dum)
    in_dum = -1
    iout_dum = -1
    !
    ! -- Create discretization object
    call disu_cr(this%dis, this%name, in_dum, iout_dum)
    !
    ! -- Create packages that are tied directly to model
    call npf_cr(this%npf, this%name, in_dum, iout_dum)
    call xt3d_cr(this%xt3d, this%name, in_dum, iout_dum)
    call gnc_cr(this%gnc, this%name, in_dum, iout_dum)
    call ic_cr(this%ic, this%name, in_dum, iout_dum, this%dis)
    call mvr_cr(this%mvr, this%name, in_dum, iout_dum, this%npf%dis)
    call gwf_obs_cr(this%obs, in_dum)
    !
    ! -- Set model 2 name
    allocate(this%m1name, this%m2name)
    this%m1name = m1name
    this%m2name = m2name
    !
    ! -- return
    return
  end subroutine gwfhalo_cr
  
  subroutine gwfhalo_df1(this, m1id, m2id)
! ******************************************************************************
! gwfhalo_df -- define the halo model part 1.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule,  only: mem_allocate
    use ListsModule,          only: basemodellist
    use BaseModelModule,      only: GetBaseModelFromList, BaseModelType
    use MpiExchangeGwfModule, only: mpi_set_gwfhalo_world_var_int
    ! -- dummy
    class(GwfHaloModelType) :: this
    integer(I4B), intent(in) :: m1id 
    integer(I4B), intent(in) :: m2id
    ! -- local
    class(BaseModelType), pointer :: mb
    type(GwfModelType), pointer :: gwf1
    type(GwfModelType), pointer :: gwf2
 ! ------------------------------------------------------------------------------
    !
    ! -- get gwf1
    gwf1 => null()
    mb => GetBaseModelFromList(basemodellist, m1id)
    select type (mb)
    type is (GwfModelType)
      gwf1 => mb
      this%gwf1 => gwf1
    end select
    !
    ! -- get gwf2
    gwf2 => null()
    if (m2id > 0) then
      mb => GetBaseModelFromList(basemodellist, m2id)
      select type (mb)
      type is (GwfModelType)
        gwf2 => mb
        this%gwf2 => gwf2
      end select
    endif
    !    
    allocate(this%nband)
    this%nband = 1 ! TODO: read for input file
    !
    call mem_allocate(this%nexg, 'NEXG', this%name)
    call mem_allocate(this%m1ndim, 'M1NDIM', this%name)
    call mem_allocate(this%m2ndim, 'M2NDIM', this%name)
    this%m1ndim = gwf1%dis%ndim
    if (associated(this%gwf2)) then
      this%m2ndim = gwf2%dis%ndim
    else
      ! -- Get NDIM from the global reduction
       this%m2ndim = mpi_set_gwfhalo_world_var_int('DNDIM', this%m2name)
    endif
    !
    ! -- return
    return
  end subroutine gwfhalo_df1
  
  subroutine gwfhalo_df2(this, nodem1, inewton, m2_bympi)
! ******************************************************************************
! gwfhalo_df -- define the halo model by creating the ia/ja arrays and copying
!   information from each individual model into the halo model
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(GwfHaloModelType) :: this
    integer(I4B), dimension(:), intent(in), target :: nodem1
    integer(I4B), intent(in) :: inewton
    logical, intent(in) :: m2_bympi
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Set newton flag
    this%inewton = inewton
    this%npf%inewton = inewton
    !
    this%nexg = size(nodem1)
    !
    ! -- Set pointers
    this%nodem1 => nodem1
    !
    ! -- Set connectivity for M1
    call gwfhalo_df_con(this%m1con, this%gwf1, nodem1, this%m1nbnod, this%nband, &
      this%imapnodem1tohalo, this%imapm1tohalo, this%m1nodes, this%name, 'M1')
    !
    allocate(this%offset)
    this%offset = sum(this%m1nbnod)
    !
    ! -- Allocate for M2
    if (m2_bympi) then
      allocate(this%m2con)
      call this%m2con%allocate_scalars(trim(this%name)//'_M2')
    end if
    !
    ! -- return
    return
  end subroutine gwfhalo_df2
  
  subroutine gwfhalo_df3(this, nodem2, ihc, cl1, cl2,  &
    hwva)
! ******************************************************************************
! gwfhalo_df -- define the halo model by creating the ia/ja arrays and copying
!   information from each individual model into the halo model
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConnectionsModule, only: ConnectionsType, fillisym, filljas
    use SparseModule, only: sparsematrix
    use MemoryManagerModule, only: mem_allocate
    use MpiExchangeModule, only: MpiWorld
    ! -- dummy
    class(GwfHaloModelType) :: this
    integer(I4B), dimension(:), intent(in), target :: nodem2
    integer(I4B), dimension(:), intent(in) :: ihc
    real(DP), dimension(:), intent(in) :: cl1
    real(DP), dimension(:), intent(in) :: cl2
    real(DP), dimension(:), intent(in) :: hwva
    ! -- local
    integer(I4B) :: nexg
    integer(I4B) :: n, m, n1, n2
    integer(I4B) :: iexg, i, j, k
    integer(I4B) :: isym, isymh
    integer(I4B) :: ierror
    type(sparsematrix) :: sparse
    integer(I4B), dimension(:), allocatable :: rowmaxnnz
    !
    ! -- pointers
    type(ConnectionsType), pointer :: con
! ------------------------------------------------------------------------------
    !
    nexg = size(this%nodem1)
    !
    if (associated(this%gwf2)) then
      this%nodem2 => nodem2
      !
      call gwfhalo_df_con(this%m2con, this%gwf2, nodem2, this%m2nbnod, this%nband, &
        this%imapnodem2tohalo, this%imapm2tohalo, this%m2nodes, this%name, 'M2')
    endif
    !
    ! --- Initialize sparse matrix
    this%neq = sum(this%m1nbnod) + sum(this%m2nbnod)
    this%dis%nodes = this%neq
    this%dis%nodesuser = this%neq
    !
    ! -- Allocate the arrays of the halo model disu package
    call this%dis%allocate_arrays()
    !
    ! -- Copy the top, bot, and areas to halo
    call copydbltohalo(this%gwf1%dis%top, this%dis%top, this%imapm1tohalo, 0)
    call copydbltohalo(this%gwf1%dis%bot, this%dis%bot, this%imapm1tohalo, 0)
    call copydbltohalo(this%gwf1%dis%area, this%dis%area, this%imapm1tohalo, 0)
    !
    if (associated(this%gwf2)) then
      call copydbltohalo(this%gwf2%dis%top, this%dis%top, this%imapm2tohalo, this%offset)
      call copydbltohalo(this%gwf2%dis%bot, this%dis%bot, this%imapm2tohalo, this%offset)
      call copydbltohalo(this%gwf2%dis%area, this%dis%area, this%imapm2tohalo, this%offset)
    else
      call MpiWorld%mpicopydbltohalo('HALO_INIT_DIS', 'TOP', 'DIS', this%name, this%dis%top, this%offset)
      call MpiWorld%mpicopydbltohalo('HALO_INIT_DIS', 'BOT', 'DIS', this%name, this%dis%bot, this%offset)
      call MpiWorld%mpicopydbltohalo('HALO_INIT_DIS', 'AREA', 'DIS', this%name, this%dis%area, this%offset)
    endif
    !
    allocate(rowmaxnnz(this%neq))
    do n = 1, this%neq
      rowmaxnnz(n) = 6
    end do
    call sparse%init(this%neq, this%neq, rowmaxnnz)
    !
    ! -- Diagonals for all halo model cells
    do n = 1, this%neq
      call sparse%addconnection(n, n, 1)
    end do
    !
    ! -- Internal connections of M1
    con => this%m1con
    do i = 1, con%nodes
      n1 = i
      do j = con%ia(i) + 1, con%ia(i + 1) - 1
        n2 =  con%ja(j)
        call sparse%addconnection(n1, n2, 1)
      end do
    end do
    ! -- Internal connections of M2
    con => this%m2con
    do i = 1, con%nodes
      n1 = i + this%offset
      do j = con%ia(i) + 1, con%ia(i + 1) - 1
        n2 =  con%ja(j) + this%offset
        call sparse%addconnection(n1, n2, 1)
      end do
    end do
    ! -- Connections between M1 and M2
    do i = 1, size(this%nodem1)
      n1 = this%imapnodem1tohalo(i)
      n2 = this%imapnodem2tohalo(i) + this%offset
      call sparse%addconnection(n1, n2, 1)
      call sparse%addconnection(n2, n1, 1)
    end do
    !
    allocate(this%dis%con)
    call this%dis%con%allocate_scalars(this%name)
    this%dis%con%nodes = this%neq
    this%dis%con%nja = sparse%nnz
    this%dis%con%njas = (this%dis%con%nja - this%dis%con%nodes) / 2
    this%dis%njas = this%dis%con%njas
    !
    ! -- Allocate index arrays of size nja and symmetric arrays
    call this%dis%con%allocate_arrays()
    !
    ! -- Fill the IA and JA arrays from sparse, then destroy sparse
    call sparse%filliaja(this%dis%con%ia, this%dis%con%ja, ierror)
    call sparse%destroy()
    !
    ! -- Create the isym and jas arrays
    call fillisym(this%dis%con%nodes, this%dis%con%nja,                        &
      this%dis%con%ia, this%dis%con%ja, this%dis%con%isym)
    call filljas(this%dis%con%nodes, this%dis%con%nja,                         &
      this%dis%con%ia, this%dis%con%ja, this%dis%con%isym,                     &
      this%dis%con%jas)
    !
    ! -- Fill the ihc, cl1, cl2, and hwva arrays for connections of M1
    con => this%m1con
    do i = 1,con%nodes
      n1 = i
      do j = con%ia(i) + 1, con%ia(i + 1) - 1
        n2 = con%ja(j)
        k = con%getjaindex(n1, n2); isym = con%jas(k)
        k = this%dis%con%getjaindex(n1, n2); isymh = this%dis%con%jas(k)
        this%dis%con%ihc(isymh)  = con%ihc(isym)
        this%dis%con%cl1(isymh)  = con%cl1(isym)
        this%dis%con%cl2(isymh)  = con%cl2(isym)
        this%dis%con%hwva(isymh) = con%hwva(isym)
      end do
    end do
    ! -- Fill the ihc, cl1, cl2, and hwva arrays for connections of M2
    con => this%m2con
    do i = 1,con%nodes
      n1 = i
      do j = con%ia(i) + 1, con%ia(i + 1) - 1
        n2 = con%ja(j)
        k = con%getjaindex(n1, n2); isym = con%jas(k)
        k = this%dis%con%getjaindex(n1 + this%offset, n2 + this%offset)
        isymh = this%dis%con%jas(k)
        this%dis%con%ihc(isymh)  = con%ihc(isym)
        this%dis%con%cl1(isymh)  = con%cl1(isym)
        this%dis%con%cl2(isymh)  = con%cl2(isym)
        this%dis%con%hwva(isymh) = con%hwva(isym)
      end do
    end do
    !
    ! -- Fill the ihc, cl1, cl2, and hwva arrays for connections of M2
    do iexg = 1, nexg
      n = this%imapnodem1tohalo(iexg)
      m = this%imapnodem2tohalo(iexg) + this%offset
      i = this%dis%con%getjaindex(n, m)
      isym = this%dis%con%jas(i)
      this%dis%con%ihc(isym) = ihc(iexg)
      this%dis%con%cl1(isym) = cl1(iexg)
      this%dis%con%cl2(isym) = cl2(iexg)
      this%dis%con%hwva(isym) = hwva(iexg)
    enddo
    !
    ! -- Finish setting up the halo model
    this%nja = this%dis%con%nja
    call this%allocate_arrays()
    !allocate(this%amat(this%nja))
    !
    ! -- Store information needed for observations
    call this%obs%obs_df(this%iout, this%name, 'GWFH', this%dis)
    !
    ! -- return
    return
  end subroutine gwfhalo_df3
  
  subroutine gwfhalo_df_con(con, gwf, nodem, nbnodes, nband, imapnodemtohalo,   &
    imapmtohalo, mnodes, name, mstr)
! ******************************************************************************
! gwfhalo_df_con -- define the connection.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SparseModule, only: sparsematrix
    use MemoryManagerModule, only: mem_allocate
    use ConnectionsModule, only: ConnectionsType, fillisym, filljas

    ! -- dummy
    type(ConnectionsType), intent(out), pointer :: con !out
    type(GwfModelType), intent(in), target :: gwf
    integer(I4B), dimension(:), intent(in), target :: nodem
    integer(I4B), dimension(:), intent(out), pointer :: nbnodes !out
    integer(I4B), intent(in) :: nband
    integer(I4B), dimension(:), intent(out), pointer :: imapnodemtohalo !out
    integer(I4B), dimension(:), intent(out), pointer :: imapmtohalo !out
    integer(I4B), dimension(:), intent(out), pointer :: mnodes !out
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: mstr
    
    ! -- local
    integer(I4B) :: i, j, n, m, ipos, n1, n2, m1, m2, iband, neq, iact
    integer(I4B) :: isymn, isymm
    integer(I4B) :: ierror, nexg
    type(sparsematrix) :: sparse
    integer(I4B), dimension(:), allocatable :: rowmaxnnz
    integer(I4B), dimension(:), allocatable :: iwrk
! ------------------------------------------------------------------------------
    !
    !
    nexg = size(nodem)
    !
    ! -- Allocate
    call mem_allocate(nbnodes, nband, 'NBNODES', trim(name)//'_'//trim(mstr))
    nbnodes = 0
    call mem_allocate(imapnodemtohalo, nexg, 'IMAPNODEMTOHALO',                    &
       trim(name)//'_'//trim(mstr))
    !
    allocate(iwrk(gwf%neq))
    !
    do iact = 1, 2
      !
      iwrk = 0
      nbnodes = 0
      !
      ! -- Set first band (equals unique nodes in nodem1)
      m = 0
      do i = 1, nexg
        n = nodem(i)
        if (iwrk(n) == 0) then
          m = m + 1
          iwrk(n) = m
          if (iact == 2) then
            imapmtohalo(m) = n
            mnodes(m) = m
          end if
        end if
      end do
      nbnodes(1) = m
      !
      if (iact == 1) then
        do i = 1, nexg
          n = nodem(i)
          imapnodemtohalo(i) = iwrk(n)
        end do
      end if
      !
      ! -- Label the bands
      do iband = 2, nband
        do i= 1, gwf%neq
          n1 = iwrk(i)
          if ((n1 > 0) .and. (n1 <= gwf%neq)) then
            ! -- Label the source band nodes negative
            iwrk(i) = -iwrk(i)
            do ipos = gwf%dis%con%ia(i) + 1, gwf%dis%con%ia(i + 1) - 1
              n = gwf%dis%con%ja(ipos)
              if (iwrk(n) == 0) then
                nbnodes(iband) = nbnodes(iband) + 1
                m = m + 1; n2 = m
                ! -- Label the target band nodes larger than neq
                iwrk(n) = n2 + gwf%neq
                if (iact == 2) then
                  imapmtohalo(m) = n
                  mnodes(m) = m
                end if
              end if
            end do
          end if
        end do
        do i= 1, gwf%neq
          if (iwrk(i) > gwf%neq) then
            iwrk(i) = abs(iwrk(i))
            iwrk(i) = iwrk(i) - gwf%neq
          end if
        end do
      end do ! iband
      !
      do i= 1, gwf%neq
        iwrk(i) = abs(iwrk(i))
      end do
      !
      ! -- Set up the sparse matrix connections
      if (iact == 2) then
        n = 0
        do i = 1, size(imapmtohalo)
          n1 = imapmtohalo(i)
          do ipos = gwf%dis%con%ia(n1) + 1, gwf%dis%con%ia(n1 + 1) - 1
            n2 = gwf%dis%con%ja(ipos)
            if (iwrk(n2) > 0) then
              n = n + 1
              !write(s,'(5(i3.3,a))')&
              !  n, ' --> ',i,'(',n1,') ---> ',iwrk(n2),'(',n2,')'
              !write(*,'(a)') trim(s)
              call sparse%addconnection(i, iwrk(n2), 1)
            end if  
          end do
        end do
      end if 
      !
      if (iact == 1) then
        !
        ! -- Setup a sparse matrix object in order to determine halo connectivity
        neq = sum(nbnodes)
        allocate(rowmaxnnz(neq))
        do n = 1, neq
          rowmaxnnz(n) = 6
        end do
        call sparse%init(neq, neq, rowmaxnnz)
        !
        ! -- Diagonals for all halo model cells
        do n = 1, neq
          call sparse%addconnection(n, n, 1)
        end do
        !
        ! -- Allocate mapping array
        call mem_allocate(imapmtohalo, neq, 'IMAPMTOHALO', trim(name)//'_'//trim(mstr))
        call mem_allocate(mnodes, neq, 'MNODES', trim(name)//'_'//trim(mstr))
        !
      end if
    end do ! iact
    !
    ! --- Allocate connectivity data
    allocate(con)
    call con%allocate_scalars(trim(name)//'_'//trim(mstr))
    con%nodes = neq
    con%nja = sparse%nnz
    con%njas = (con%nja - con%nodes) / 2
    !
    ! -- Allocate index arrays of size nja and symmetric arrays
    call con%allocate_arrays()
    !
    ! -- Fill the IA and JA arrays from sparse, then destroy sparse
    call sparse%filliaja(con%ia, con%ja, ierror)
    call sparse%destroy()
    !
    ! -- Create the isym and jas arrays
    call fillisym(con%nodes, con%nja, con%ia, con%ja, con%isym)
    call filljas(con%nodes, con%nja, con%ia, con%ja, con%isym, con%jas)
    !
    ! -- Fill the ihc, cl1, cl2, and hwva arrays
    do i = 1, size(imapmtohalo)
      n1 = imapmtohalo(i) ! gwf node
      do ipos = gwf%dis%con%ia(n1) + 1, gwf%dis%con%ia(n1 + 1) - 1
        n2 = gwf%dis%con%ja(ipos) ! gwf node
        if (iwrk(n2) > 0) then
          m1 = i; m2 = iwrk(n2) ! halo node2
          j = gwf%dis%con%getjaindex(n1, n2); isymn = gwf%dis%con%jas(j)
          j = con%getjaindex(m1, m2); isymm = con%jas(j)
          con%ihc(isymm)  = gwf%dis%con%ihc(isymn)
          con%cl1(isymm)  = gwf%dis%con%cl1(isymn)
          con%cl2(isymm)  = gwf%dis%con%cl2(isymn)
          con%hwva(isymm) = gwf%dis%con%hwva(isymn)
        end if
      end do
    end do
    !
    ! -- Cleanup work array
    deallocate(iwrk)
    !
  end subroutine gwfhalo_df_con
    
  subroutine gwfhalo_ar(this)
! ******************************************************************************
! gwfhalo_ar -- allocate in the halo model, and copy the package information
!   from gwf1 and gwf2 into the halo model
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_reallocate
    use MpiExchangeGwfModule, only: mpi_set_gwfhalo_world_var_int
    use MpiExchangeModule, only: MpiExchangeType, MpiWorld
    use BndModule, only: BndType, GetBndFromList
    ! -- dummy
    class(GwfHaloModelType) :: this
    ! -- local
    class(BndType), pointer :: packobj
    type(MpiExchangeType), pointer :: mpi
    integer(I4B) :: n, m1flag, m2flag, ip
! ------------------------------------------------------------------------------
    !
    ! -- Create ibound and x arrays for halo model
    call mem_allocate(this%ibound, this%neq, 'IBOUND', this%name)
    call mem_allocate(this%x, this%neq, 'X', this%name)
    !
    ! -- Copy ibound from gwf models into this halo model
    call copyinttohalo(this%gwf1%ibound, this%ibound, this%imapm1tohalo, 0)
    if (associated(this%gwf2)) then
      call copyinttohalo(this%gwf2%ibound, this%ibound, this%imapm2tohalo,     &
        this%offset)
    else
      mpi => this%MpiSol
      call mpi%mpicopyinttohalo('X_IACTIVE', 'IACTIVE', '', this%name,         &
        this%ibound, this%offset)
    endif
    !
    mpi => MpiWorld
    !
    ! -- Set up ic package and copy strt from gwf models into halo model
    call this%ic%allocate_arrays(this%dis%nodes)
    call copydbltohalo(this%gwf1%ic%strt, this%ic%strt, this%imapm1tohalo, 0)
    if (associated(this%gwf2)) then
      call copydbltohalo(this%gwf2%ic%strt, this%ic%strt, this%imapm2tohalo,   &
        this%offset)
    else
      call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'STRT', 'IC', this%name,    &
        this%ic%strt, this%offset)
    endif
    do n = 1, this%dis%nodes
      this%x(n) = this%ic%strt(n)
    enddo
    !
    ! -- Copy or point npf pointers
    this%npf%dis    => this%dis
    this%npf%ic     => this%ic
    this%npf%ibound => this%ibound
    this%npf%hnew   => this%x
    !
    ! -- Allocate NPF arrays
    call this%npf%allocate_arrays(this%dis%nodes, this%dis%njas)
    !
    ! -- Copy required NPF arrays into halo model
    call copyinttohalo(this%gwf1%npf%icelltype, this%npf%icelltype,            &
      this%imapm1tohalo, 0)
    if (associated(this%gwf2)) then
      call copyinttohalo(this%gwf2%npf%icelltype, this%npf%icelltype,          &
        this%imapm2tohalo, this%offset)
    else
      call mpi%mpicopyinttohalo('HALO_INIT_NPFIC', 'ICELLTYPE', 'NPF',         &
         this%name, this%npf%icelltype, this%offset)
    endif
    call copydbltohalo(this%gwf1%npf%k11, this%npf%k11, this%imapm1tohalo, 0)
    if (associated(this%gwf2)) then
      call copydbltohalo(this%gwf2%npf%k11, this%npf%k11, this%imapm2tohalo,   &
        this%offset)
    else
      call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'K11', 'NPF', this%name,    &
        this%npf%k11, this%offset)
    endif
    !
    call copydbltohalo(this%gwf1%npf%sat, this%npf%sat, this%imapm1tohalo, 0)
    if (associated(this%gwf2)) then
      call copydbltohalo(this%gwf2%npf%sat, this%npf%sat, this%imapm2tohalo,   &
        this%offset)
    else
      call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'SAT', 'NPF', this%name,    &
         this%npf%sat, this%offset)
    endif
    !
    ! -- Copy k22
    m1flag = this%gwf1%npf%ik22
    if (.not.associated(this%gwf2)) then
      m2flag = mpi_set_gwfhalo_world_var_int('IK22', this%m2name)
    else
      m2flag = this%gwf2%npf%ik22
    endif
    if ((m1flag /= 0).or.(m2flag /= 0)) then
      this%npf%ik22 = 1
      call mem_reallocate(this%npf%k22, this%dis%nodes, 'K22',                 &
        trim(this%npf%memoryPath))
      if (m1flag /= 0) then
        call copydbltohalo(this%gwf1%npf%k22, this%npf%k22, this%imapm1tohalo, &
          0)
      else
        call copydbltohalo(this%gwf1%npf%k11, this%npf%k22, this%imapm1tohalo, &
          0)
      endif
      if (associated(this%gwf2)) then
        if (m2flag /= 0) then
          call copydbltohalo(this%gwf2%npf%k22, this%npf%k22,                  &
            this%imapm2tohalo, this%offset)
        else
          call copydbltohalo(this%gwf2%npf%k11, this%npf%k22,                  &
            this%imapm2tohalo, this%offset)
        endif
      else
        if (m2flag /= 0) then
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'K22', 'NPF', this%name,&
            this%npf%k22, this%offset)
        else
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'K11', 'NPF', this%name,&
            this%npf%k22, this%offset)
        endif
      endif
    endif
    !
    ! -- Copy k33
    m1flag = this%gwf1%npf%ik33
    if (.not.associated(this%gwf2)) then
      m2flag = mpi_set_gwfhalo_world_var_int('IK33', this%m2name)
    else
      m2flag = this%gwf2%npf%ik33
    endif
    !
    if ((m1flag /= 0).or.(m2flag /= 0)) then
      this%npf%ik33 = 1
      call mem_reallocate(this%npf%k33, this%dis%nodes, 'K33',                 &
        trim(this%npf%memoryPath))
      if (m1flag /= 0) then
        call copydbltohalo(this%gwf1%npf%k33, this%npf%k33, this%imapm1tohalo, &
         0)
      else
        call copydbltohalo(this%gwf1%npf%k11, this%npf%k33, this%imapm1tohalo, &
          0)
      endif
      if (associated(this%gwf2)) then
        if (m2flag /= 0) then
          call copydbltohalo(this%gwf2%npf%k33, this%npf%k33,                  &
            this%imapm2tohalo, this%offset)
        else
          call copydbltohalo(this%gwf2%npf%k11, this%npf%k33,                  &
            this%imapm2tohalo, this%offset)
        endif
      else
        if (m2flag /= 0) then
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'K33', 'NPF', this%name,&
            this%npf%k33, this%offset)
        else
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'K11', 'NPF', this%name,&
            this%npf%k33, this%offset)
        endif
      endif
    endif
    !
    ! -- Copy angle1
    m1flag = this%gwf1%npf%iangle1
    if (.not.associated(this%gwf2)) then
      m2flag = mpi_set_gwfhalo_world_var_int('IANGLE1', this%m2name)
    else
      m2flag = this%gwf2%npf%iangle1
    endif
    if ((m1flag /= 0).or.(m2flag /= 0)) then
      this%npf%iangle1 = 1
      call mem_reallocate(this%npf%angle1, this%dis%nodes, 'ANGLE1',           &
        trim(this%npf%memoryPath))
      do n = 1, this%dis%nodes
        this%npf%angle1(n) = DZERO
      end do
      if (m1flag/= 0) then
        call copydbltohalo(this%gwf1%npf%angle1, this%npf%angle1,              &
          this%imapm1tohalo, 0)
      endif
      if (associated(this%gwf2)) then
        if (m2flag /= 0) then
          call copydbltohalo(this%gwf2%npf%angle1, this%npf%angle1,            &
            this%imapm2tohalo, this%offset)
        endif
      else
        if (m2flag /= 0) then
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'ANGLE1', 'NPF',        &
            this%name, this%npf%angle1, this%offset)
        endif
      endif
    endif
    !
    ! -- Copy angle2
    m1flag = this%gwf1%npf%iangle2
    if (.not.associated(this%gwf2)) then
      m2flag = mpi_set_gwfhalo_world_var_int('IANGLE2', this%m2name)
    else
      m2flag = this%gwf2%npf%iangle2
    endif
    if ((m1flag /= 0).or.(m2flag /= 0)) then
      this%npf%iangle2 = 1
      call mem_reallocate(this%npf%angle2, this%dis%nodes, 'ANGLE2',           &
        trim(this%npf%memoryPath))
      do n = 1, this%dis%nodes
        this%npf%angle2(n) = DZERO
      end do
      if (m1flag/= 0) then
        call copydbltohalo(this%gwf1%npf%angle2, this%npf%angle2,              &
          this%imapm1tohalo, 0)
      endif
      if (associated(this%gwf2)) then
        if (m2flag /= 0) then
          call copydbltohalo(this%gwf2%npf%angle2, this%npf%angle2,            &
            this%imapm2tohalo, this%offset)
        endif
      else
        if (m2flag /= 0) then
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'ANGLE2', 'NPF',        &
            this%name, this%npf%angle2, this%offset)
        endif
      endif
    endif
    !
   ! -- Copy angle3
    m1flag = this%gwf1%npf%iangle3
    if (.not.associated(this%gwf2)) then
      m2flag = mpi_set_gwfhalo_world_var_int('IANGLE3', this%m2name)
    else
      m2flag = this%gwf2%npf%iangle3
    endif
    if ((m1flag /= 0).or.(m2flag /= 0)) then
      this%npf%iangle3 = 1
      call mem_reallocate(this%npf%angle3, this%dis%nodes, 'ANGLE2',           &
        trim(this%npf%memoryPath))
      do n = 1, this%dis%nodes
        this%npf%angle3(n) = DZERO
      end do
      if (m1flag/= 0) then
        call copydbltohalo(this%gwf1%npf%angle3, this%npf%angle3,              &
          this%imapm1tohalo, 0)
      endif
      if (associated(this%gwf2)) then
        if (m2flag /= 0) then
          call copydbltohalo(this%gwf2%npf%angle3, this%npf%angle3,            &
            this%imapm2tohalo, this%offset)
        endif
      else
        if (m2flag /= 0) then
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'ANGLE3', 'NPF',        &
            this%name, this%npf%angle3, this%offset)
        endif
      endif
    endif
    !
    ! -- prepcheck
    call this%npf%prepcheck()
    !
    ! -- Set rewet flags
    call mem_allocate(this%m1irewet, 'M1IREWET', this%name)
    this%m1irewet = this%gwf1%npf%irewet
    call mem_allocate(this%m2irewet, 'M2IREWET', this%name)
    if (.not.associated(this%gwf2)) then
      this%m2irewet = mpi_set_gwfhalo_world_var_int('IREWET', this%m2name)
    else
      this%m2irewet = this%gwf2%npf%irewet
    endif
    !
    ! -- Copy wetdry
    m1flag = this%m1irewet
    m2flag = this%m2irewet
   if ((m1flag /= 0).or.(m2flag /= 0)) then
      this%npf%irewet = 1
      call mem_reallocate(this%npf%wetdry, this%dis%nodes, 'WETDRY',           &
        trim(this%npf%memoryPath))
      if (m1flag/= 0) then
        call copydbltohalo(this%gwf1%npf%wetdry, this%npf%wetdry,              &
          this%imapm1tohalo, 0)
      endif
      if (associated(this%gwf2)) then
        if (m2flag /= 0) then
          call copydbltohalo(this%gwf2%npf%wetdry, this%npf%wetdry,            &
            this%imapm2tohalo, this%offset)
        endif
      else
        if (m2flag /= 0) then
          call mpi%mpicopydbltohalo('HALO_INIT_NPFIC', 'WETDRY', 'NPF',        &
            this%name, this%npf%wetdry, this%offset)
        endif
      endif
    endif
    !
    ! -- AR for mover packages
    if (.not.associated(this%gwf2)) then
      do ip=1,this%bndlist%Count()
        packobj => GetBndFromList(this%bndlist, ip)
        call packobj%bnd_ar()
      enddo
    endif
    !
    ! -- return
    return
  end subroutine gwfhalo_ar

  subroutine gwfhalo_cf1(this, kiter)
! ******************************************************************************
! gwfhalo_cf1 -- update heads and ibound in halo model
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeModule, only: MpiExchangeType
    use BndModule, only: BndType,  GetBndFromList
    ! -- dummy
    class(GwfHaloModelType) :: this
    integer(I4B), intent(in) :: kiter
    ! -- local
    type(MpiExchangeType), pointer :: mpi
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
! ------------------------------------------------------------------------------
    !
    ! -- Copy x from gwf models into this halo model x
    call copydbltohalo(this%gwf1%x, this%x, this%imapm1tohalo, 0)
    call copyinttohalo(this%gwf1%ibound, this%ibound, this%imapm1tohalo, 0)
    if (associated(this%gwf2)) then
      call copydbltohalo(this%gwf2%x, this%x, this%imapm2tohalo, this%offset)
      call copyinttohalo(this%gwf2%ibound, this%ibound, this%imapm2tohalo,      &
        this%offset)
    else
      mpi => this%MpiSol
      call mpi%mpicopydbltohalo('X_IACTIVE', 'X', '', this%name,                &
        this%x, this%offset)
      call mpi%mpicopyinttohalo('X_IACTIVE', 'IACTIVE', '', this%name,          &
        this%ibound, this%offset)
    endif
    !
    ! -- return
    return
  end subroutine gwfhalo_cf1
  
  subroutine gwfhalo_cf2(this, kiter)
! ******************************************************************************
! gwfhalo_cf2 -- rewetting calculations
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeModule, only: MpiExchangeType
    use BndModule, only: BndType,  GetBndFromList
    ! -- dummy
    class(GwfHaloModelType) :: this
    integer(I4B), intent(in) :: kiter
    ! -- local
    type(MpiExchangeType), pointer :: mpi
    class(BndType), pointer :: packobj
    integer(I4B) :: ip
! ------------------------------------------------------------------------------
    !
    ! -- rewet
    ! -- TODO: this is a problem.  Wetting and drying prints messages, which
    !    doesn't make sense for the halo model.  And the halo model should 
    !    only check for wetting for the n-m connection, and not recheck for
    !    all the cells in the halo model.  This is not good.
    !call this%npf%npf_cf(kiter, this%dis%nodes, this%x)
    !
    ! -- Calculate saturation for convertible cells
    call this%npf%sat_calc(this%x)
    !
    ! -- CF for package movers
    if (.not.associated(this%gwf2)) then
      do ip=1,this%bndlist%Count()
        packobj => GetBndFromList(this%bndlist, ip)
        call packobj%bnd_cf()
      enddo
    endif
    !
    ! -- return
    return
  end subroutine gwfhalo_cf2

!  subroutine gwfhalo_fc(this, kiter, cond)
!! ******************************************************************************
!! gwfhalo_fc -- fill halo model conductance terms into separate amat and rhs 
!!   vectors so that the gwf-gwf exchange can add them to the solution amat
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- modules
!    use SimModule, only: ustop, store_error, count_errors
!    ! -- dummy
!    class(GwfHaloModelType) :: this
!    integer(I4B), intent(in) :: kiter
!    real(DP), dimension(:), intent(inout) :: cond
!    ! -- local
!    character(len=LINELENGTH) :: errmsg
!    integer(I4B) :: n
!    integer(I4B) :: m
!    integer(I4B) :: nexg
!    integer(I4B) :: ipos
!    real(DP), dimension(:), allocatable :: rhs
!    integer(I4B), dimension(:), allocatable :: idxglo
!! ------------------------------------------------------------------------------
!    !
!    ! -- allocate
!    allocate(rhs(this%dis%nodes))
!    allocate(idxglo(this%nja))
!    !
!    ! -- Copy x from models into this halo model x
!    call copydbltohalo(this%gwf1%x, this%x, this%imapm1tohalo)
!    call copydbltohalo(this%gwf2%x, this%x, this%imapm2tohalo)
!    !
!    ! -- Fill halo amat so that the exchange object can add the terms
!    do n = 1, this%nja
!      this%amat(n) = DZERO
!      idxglo(n) = n
!    enddo
!    do n = 1, this%dis%nodes
!      rhs(n) = DZERO
!    enddo
!    call this%npf%npf_fc(kiter, this%dis%nodes, this%nja, this%nja, this%amat, &
!      idxglo, rhs, this%x)
!    !
!    ! -- TODO: for now just compare halo cond with gwf-gwf cond
!    do nexg = 1, size(this%nodem1)
!      n = this%imapm1tohalo(this%nodem1(nexg))
!      m = this%imapm2tohalo(this%nodem2(nexg))
!      ipos = this%dis%con%getjaindex(n, m)
!      if (this%amat(ipos) - cond(nexg) /= DZERO) then
!        write(errmsg, '(a, 1x, i0, 1x, i0, 2(1pg25.16))') 'HALO MODEL COND /= COND', n, m, &
!          this%amat(ipos), cond(nexg)
!        call store_error(errmsg)
!      endif
!    enddo
!    if (count_errors() > 0) call ustop() 
!    !
!    ! -- return
!    return
!  end subroutine gwfhalo_fc
  
  subroutine gwfhalo_fc_calc(this, iexg, terms)
! ******************************************************************************
! gwfhalo_fc_calc -- calculate the amat and rhs terms for connected gwf cells
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(GwfHaloModelType) :: this
    integer(I4B), intent(in) :: iexg
    real(DP), dimension(:, :), intent(inout) :: terms
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ii
! ------------------------------------------------------------------------------
    n = this%imapnodem1tohalo(iexg)
    m = this%imapnodem2tohalo(iexg) + this%offset
    ii = this%dis%con%getjaindex(n, m)
    call this%npf%npf_fc_calc(n, m, ii, this%x, terms)
    !
    ! -- return
    return
  end subroutine gwfhalo_fc_calc

  subroutine gwfhalo_fn_calc(this, iexg, terms)
! ******************************************************************************
! gwfhalo_fn_calc -- calculate the newton amat and rhs terms for connected gwf cells
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(GwfHaloModelType) :: this
    integer(I4B), intent(in) :: iexg
    real(DP), dimension(:, :), intent(inout) :: terms
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ii
! ------------------------------------------------------------------------------
    n = this%imapnodem1tohalo(iexg)
    m = this%imapnodem2tohalo(iexg) + this%offset
    ii = this%dis%con%getjaindex(n, m)
    call this%npf%npf_fn_calc(n, m, ii, this%x, terms)
    !
    ! -- return
    return
  end subroutine gwfhalo_fn_calc

  subroutine copydbltohalo(gwfarray, gwfhaloarray, imaptohalo, offset)
! ******************************************************************************
! copydbltohalo -- allocate in the halo model 
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    real(DP), dimension(:), intent(in) :: gwfarray
    real(DP), dimension(:), intent(inout) :: gwfhaloarray
    integer(I4B), dimension(:), intent(in) :: imaptohalo
    integer, intent(in) :: offset
    ! -- integer
    integer(I4B) :: n, i
! ------------------------------------------------------------------------------
    do i = 1, size(imaptohalo)
      n = imaptohalo(i)
      gwfhaloarray(i+offset) = gwfarray(n)
    enddo
    !
    ! -- return
    return
  end subroutine copydbltohalo

  subroutine copyinttohalo(gwfarray, gwfhaloarray, imaptohalo, offset)
! ******************************************************************************
! copyinttohalo -- allocate in the halo model 
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer(I4B), dimension(:), intent(in) :: gwfarray
    integer(I4B), dimension(:), intent(inout) :: gwfhaloarray
    integer(I4B), dimension(:), intent(in) :: imaptohalo
    integer, intent(in) :: offset
    ! -- local
    integer(I4B) :: n, i
! ------------------------------------------------------------------------------
    do i = 1, size(imaptohalo)
      n = imaptohalo(i)
      gwfhaloarray(i+offset) = gwfarray(n)
    enddo
    !
    ! -- return
    return
  end subroutine copyinttohalo
!
end module GwfHaloModule
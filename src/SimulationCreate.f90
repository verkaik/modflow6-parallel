module SimulationCreateModule

  use KindModule,             only: DP, I4B, write_kindinfo
  use ConstantsModule,        only: LINELENGTH, LENMODELNAME, LENBIGLINE, DZERO
  use SimVariablesModule,     only: simfile, simlstfile, iout
  use GenericUtilitiesModule, only: sim_message, write_centered
  use SimModule,              only: ustop, store_error, count_errors,            &
                                    store_error_unit, MaxErrors
  use VersionModule,          only: write_listfile_header
  use InputOutputModule,      only: getunit, urword, openfile
  use ArrayHandlersModule,    only: expandarray, ifind
  use BaseModelModule,        only: BaseModelType
  use BaseSolutionModule,     only: BaseSolutionType, AddBaseSolutionToList,     &
                                    GetBaseSolutionFromList
  use SolutionGroupModule,    only: SolutionGroupType, AddSolutionGroupToList
  use BaseExchangeModule,     only: BaseExchangeType
  use ListsModule,            only: basesolutionlist, basemodellist,             &
                                    solutiongrouplist
  use BaseModelModule,        only: GetBaseModelFromList
  use BlockParserModule,      only: BlockParserType

  implicit none
  private
  public :: simulation_cr
  public :: simulation_da
  public :: modelname !PAR

  integer(I4B) :: inunit = 0
  character(len=LENMODELNAME), allocatable, dimension(:) :: modelname
  character(len=LENMODELNAME), allocatable, dimension(:), save :: modelname_all !PAR
  integer, allocatable, dimension(:), save                     :: model_sub !PAR
  integer, allocatable, dimension(:), save :: model_topol_m1 !CGC
  integer, allocatable, dimension(:), save :: model_topol_m2 !CGC
  type(BlockParserType) :: parser

  contains

  subroutine simulation_cr()
! ******************************************************************************
! Read the simulation name file and initialize the models, exchanges
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeGenModule, only: parallelrun, writestd !PAR
    character(len=LINELENGTH) :: line
! ------------------------------------------------------------------------------
    !
    ! -- initialize iout 
    iout = 0
    !
    ! -- Open simulation list file
    iout = getunit()
    call openfile(iout, 0, simlstfile, 'LIST', filstat_opt='REPLACE',            &
                  master_write=.true.) !PAR
    !
    ! -- write simlstfile to stdout
    if (parallelrun) then !PAR
      write(line,'(2(1x,A))') 'Writing simulation list file for MPI rank 0:',    & !PAR
                              trim(adjustl(simlstfile)) !PAR
    else !PAR
    write(line,'(2(1x,A))') 'Writing simulation list file:',                     &
                            trim(adjustl(simlstfile))
    endif !PAR
    call sim_message(line, force_write=writestd) !PAR
    call write_listfile_header(iout)
    !
    ! -- Read the simulation name file and create objects
    simfile = adjustl(simfile) !PAR
    call read_simulation_namefile(simfile) !PAR
    !
    ! -- Return
    return
  end subroutine simulation_cr

  subroutine simulation_da()
! ******************************************************************************
! Deallocate simulation variables
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- variables
    if (allocated(modelname)) then !PAR
      deallocate(modelname)
    endif !PAR
    if (allocated(modelname_all)) then !PAR
      deallocate(modelname_all) !PAR
    endif !PAR
    if (allocated(model_sub)) then !PAR
      deallocate(model_sub) !PAR
    endif !PAR
    !
    ! -- Return
    return
  end subroutine simulation_da

  subroutine read_simulation_namefile(simfile)
! ******************************************************************************
! Read the simulation name file and initialize the models, exchanges,
! solutions, solutions groups.  Then add the exchanges to the appropriate
! solutions.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeGenModule, only: writestd !PAR
    use MpiExchangeModule, only: MpiWorld
    use MpiExchangeGwfModule, only: mpi_gwfhalo_world !PAR
    ! -- dummy
    character(len=*),intent(inout) :: simfile !PAR
    ! -- local
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: errmsg
    class(BaseSolutionType), pointer :: sp
    class(BaseModelType), pointer :: mp
    integer(I4B) :: is, im
! ------------------------------------------------------------------------------
    !
    ! -- Open simulation name file
    inunit = getunit()
    call openfile(inunit, iout, simfile, 'NAM')
    !
    ! -- write simfile name to stdout
    !call MpiWorld%mpi_barrier() !PAR
    if (writestd) then !PAR
      write(line,'(2(1x,a))') 'Using Simulation name file:', trim(simfile) !PAR
      call sim_message(line, skipafter=1)
    end if !PAR
    !
    ! -- Initialize block parser
    call parser%Initialize(inunit, iout)
    !
    ! -- Process OPTIONS block in simfile
    call options_create()
    !
    ! -- Process TIMING block in simfile
    call timing_create()
    !
    ! -- Process MODELS block in simfile
    call models_create()
    !
    ! -- Collective MPI communication scalars DIS
    call mpi_gwfhalo_world(1) !PAR
    !
    ! -- Process EXCHANGES block in simfile
    call exchanges_create()
    !
    ! -- Process SOLUTION_GROUPS blocks in simfile
    call solution_groups_create()
    !
    ! -- Go through each model and make sure that it has been assigned to
    !    a solution.
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      if (mp%idsoln == 0) then
        write(errmsg, '(a,a)') &
           '****ERROR.  Model was not assigned to a solution: ', mp%name
        call store_error(errmsg)
      endif
    enddo
    if (count_errors() > 0) then
      call store_error_unit(inunit)
      call ustop()
    endif
    !
    ! -- Close the input file
    close(inunit)
    !
    ! -- Go through each solution and assign exchanges accordingly
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%slnassignexchanges()
    enddo
    !
    ! --- Initialize solution for MPI
    do is = 1, basesolutionlist%Count() !PAR
      sp => GetBaseSolutionFromList(basesolutionlist, is) !PAR
      call sp%slnmpiinit(sp%name) !PAR
    enddo !PAR
    !
    ! -- Return
    return
  end subroutine read_simulation_namefile

  subroutine options_create()
! ******************************************************************************
! Set the simulation options
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_set_print_option
    use SimVariablesModule, only: isimcontinue, isimcheck, isimdd, nddsub !PAR
    ! -- local
    integer(I4B) :: ierr
    integer(I4B) :: imax
    logical :: isfound, endOfBlock
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: keyword
! ------------------------------------------------------------------------------
    !
    ! -- Process OPTIONS block
    call parser%GetBlock('OPTIONS', isfound, ierr, &
      supportOpenClose=.true., blockRequired=.false.)
    if (isfound) then
      write(iout,'(/1x,a)')'READING SIMULATION OPTIONS'
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
          case ('CONTINUE')
            isimcontinue = 1
            write(iout, '(4x, a)')                                             &
                  'SIMULATION WILL CONTINUE EVEN IF THERE IS NONCONVERGENCE.'
          case ('NOCHECK')
            isimcheck = 0
            write(iout, '(4x, a)')                                             &
                  'MODEL DATA WILL NOT BE CHECKED FOR ERRORS.'
          case ('MEMORY_PRINT_OPTION')
            errmsg = ''
            call parser%GetStringCaps(keyword)
            call mem_set_print_option(iout, keyword, errmsg)
            if (errmsg /= ' ') then
              call store_error(errmsg)
              call parser%StoreErrorUnit()
              call ustop()
            endif
          case ('MAXERRORS')
            imax = parser%GetInteger()
            call MaxErrors(imax)
            write(iout, '(4x, a, i0)')                                         &
                  'MAXIMUM NUMBER OF ERRORS THAT WILL BE STORED IS ', imax
          case ('DOMAIN_DECOMPOSITION') !PAR
            isimdd = 1 !PAR
            nddsub = parser%GetInteger() !PAR
            write(iout, '(4x, a)')                                           & !PAR
                  'SIMULATION WILL USE DOMAIN DECOMPOSITION.' !PAR
          case default
            write(errmsg, '(4x,a,a)') &
                  '****ERROR. UNKNOWN SIMULATION OPTION: ',                    &
                  trim(keyword)
            call store_error(errmsg)
            call parser%StoreErrorUnit()
            call ustop()
        end select
      end do
      write(iout,'(1x,a)')'END OF SIMULATION OPTIONS'
    end if
    !
    ! -- return
    return
  end subroutine options_create

  subroutine timing_create()
! ******************************************************************************
! Set the timing module to be used for the simulation
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use TdisModule, only: tdis_cr
    ! -- dummy
    ! -- local
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: line, keyword
    logical :: found_tdis
! ------------------------------------------------------------------------------
    !
    ! -- Initialize
    found_tdis = .false.
    !
    ! -- Process TIMING block
    call parser%GetBlock('TIMING', isfound, ierr, &
      supportOpenClose=.true.)
    if (isfound) then
      write(iout,'(/1x,a)')'READING SIMULATION TIMING'
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
          case ('TDIS6')
            found_tdis = .true.
            call parser%GetString(line)
            call tdis_cr(line)
          case default
            write(errmsg, '(4x,a,a)') &
                  '****ERROR. UNKNOWN SIMULATION TIMING: ', &
                  trim(keyword)
            call store_error(errmsg)
            call parser%StoreErrorUnit()
            call ustop()
        end select
      end do
      write(iout,'(1x,a)')'END OF SIMULATION TIMING'
    else
      call store_error('****ERROR.  Did not find TIMING block in simulation'// &
                       ' control file.')
      call parser%StoreErrorUnit()
      call ustop()
    end if
    !
    ! -- Ensure that TDIS was found
    if(.not. found_tdis) then
      call store_error('****ERROR. TDIS not found in TIMING block.')
      call parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- return
    return
  end subroutine timing_create

  subroutine models_create()
! ******************************************************************************
! Set the models to be used for the simulation
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use GwfModule,              only: gwf_cr
    use GwtModule,              only: gwt_cr
    use ConstantsModule,        only: LENMODELNAME
    use SimVariablesModule, only: isimdd, nddsub !PAR
    use MpiExchangeModule, only: MpiWorld !PAR
    ! -- dummy
    ! -- local
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    integer(I4B) :: im
    integer(I4B) :: imdd !PAR
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: fname, mname
   integer :: isub !PAR
    logical :: add !PAR
! ------------------------------------------------------------------------------
    !
    ! -- Process MODELS block
    call parser%GetBlock('MODELS', isfound, ierr, &
      supportOpenClose=.true.)
    if (isfound) then
      write(iout,'(/1x,a)')'READING SIMULATION MODELS'
      im = 0
      imdd = 0 !PAR
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
          case ('GWF6')
            call parser%GetString(fname)
           call read_modelname(mname) !PAR
            if (isimdd == 1) then !PAR
              isub = parser%GetInteger() !PAR
              if (isub < 0 .or. isub > nddsub) then
                write(*,'(a,1x,i0)') '****ERROR. INVALID SUBDOMAIN:', isub
                call store_error(errmsg)
                call parser%StoreErrorUnit()
                call ustop()
              endif
              call MpiWorld%mpi_is_iproc(isub, add) !PAR
              ! -- Store global subdomain information !PAR
              call add_model_dd(imdd, isub, mname) !PAR
              call MpiWorld%mpi_addmodel(1, mname) !PAR
              call MpiWorld%mpi_addsub(1, isub) !PAR
              if (add) then
                call MpiWorld%mpi_addmodel(2, mname) !PAR
                call add_model(im, 'GWF6', mname)
                call gwf_cr(fname, imdd, modelname(im))
              endif
            else !PAR
              call add_model(im, 'GWF6', mname)
              call add_model_dd(imdd, 1, mname) !PAR
              call gwf_cr(fname, im, modelname(im))
            endif !PAR
          case ('GWT6')
            call parser%GetString(fname)
            call add_model(im, 'GWT6', mname)
            call gwt_cr(fname, im, modelname(im))
          case default
            write(errmsg, '(4x,a,a)') &
                  '****ERROR. UNKNOWN SIMULATION MODEL: ',                     &
                  trim(keyword)
            call store_error(errmsg)
            call parser%StoreErrorUnit()
            call ustop()
        end select
      end do
      write(iout,'(1x,a)')'END OF SIMULATION MODELS'
    else
      call store_error('****ERROR.  Did not find MODELS block in simulation'// &
                       ' control file.')
      call parser%StoreErrorUnit()
      call ustop()
    end if
    !
    ! -- return
    return
  end subroutine models_create

  subroutine exchanges_create()
! ******************************************************************************
! Set the exchanges to be used for the simulation
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use GwfGwfExchangeModule,    only: gwfexchange_create
    use GwfGwtExchangeModule,    only: gwfgwt_cr
    ! -- dummy
    ! -- local
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    integer(I4B) :: id
    integer(I4B) :: m1
    integer(I4B) :: m2
    integer(I4B) :: s1, s2 !PAR
    character(len=LINELENGTH) :: errmsg
    character(len=LINELENGTH) :: keyword
    character(len=LINELENGTH) :: fname, name1, name2
    integer(I4B) :: im = 0 !PAR
    logical :: m2_bympi !PAR
    ! -- formats
    character(len=*), parameter :: fmtmerr = "('Error in simulation control ', &
      &'file.  Could not find model: ', a)"
! ------------------------------------------------------------------------------
    call parser%GetBlock('EXCHANGES', isfound, ierr, &
      supportOpenClose=.true.)
    if (isfound) then
      write(iout,'(/1x,a)')'READING SIMULATION EXCHANGES'
      id = 0
      im = 0 !PAR
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
          case ('GWF6-GWF6')
            id = id + 1
            !
            ! -- get filename
            call parser%GetString(fname)
            !
            ! -- get first modelname and then model id
            call parser%GetStringCaps(name1)
            call parser%GetStringCaps(name2)
            !
            m1 = ifind(modelname_all, name1)
            m2 = ifind(modelname_all, name2)
            call store_model_topol(m1, m2) !CGC
            !
            if(m1 < 0) then
              write(errmsg, fmtmerr) trim(name1)
              call store_error(errmsg)
              call parser%StoreErrorUnit()
              call ustop()
            endif
            !
            ! -- get second modelname and then model id
            if(m2 < 0) then
              write(errmsg, fmtmerr) trim(name2)
              call store_error(errmsg)
              call parser%StoreErrorUnit()
              call ustop()
            endif
            !
            s1 = model_sub(m1) !PAR
            s2 = model_sub(m2) !PAR
            if (s1 /= s2) then !PAR
              m2_bympi = .true. !PAR
            else !PAR
              m2_bympi = .false. !PAR
            endif !PAR
            !
            ! -- get model id
            m1 = ifind(modelname, name1) !PAR
            m2 = ifind(modelname, name2) !PAR
            !
            ! -- Create the exchange object.
            write(iout, '(4x,a,i0,a,i0,a,i0)') 'GWF6-GWF6 exchange ', id,      &
              ' will be created to connect model ', m1, ' with model ', m2
            call gwfexchange_create(fname, id, m1, m2, name1, name2, m2_bympi) !PAR
          case ('GWF6-GWT6')
            id = id + 1
            !
            ! -- get filename
            call parser%GetString(fname)
            !
            ! -- get first modelname and then model id
            call parser%GetStringCaps(name1)
            m1 = ifind(modelname, name1)
            if(m1 < 0) then
              write(errmsg, fmtmerr) trim(name1)
              call store_error(errmsg)
              call parser%StoreErrorUnit()
              call ustop()
            endif
            !
            ! -- get second modelname and then model id
            call parser%GetStringCaps(name2)
            m2 = ifind(modelname, name2)
            if(m2 < 0) then
              write(errmsg, fmtmerr) trim(name2)
              call store_error(errmsg)
              call parser%StoreErrorUnit()
              call ustop()
            endif
            !
            ! -- Create the exchange object.
            write(iout, '(4x,a,i0,a,i0,a,i0)') 'GWF6-GWT6 exchange ', id,      &
              ' will be created to connect model ', m1, ' with model ', m2
            call gwfgwt_cr(fname, id, m1, m2)
          case default
            write(errmsg, '(4x,a,a)') &
                  '****ERROR. UNKNOWN SIMULATION EXCHANGES: ',                 &
                  trim(keyword)
            call store_error(errmsg)
            call parser%StoreErrorUnit()
            call ustop()
        end select
      end do
      write(iout,'(1x,a)')'END OF SIMULATION EXCHANGES'
    else
      call store_error('****ERROR.  Did not find EXCHANGES block in '//        &
                       'simulation control file.')
      call parser%StoreErrorUnit()
      call ustop()
    end if
    !
    ! -- return
    return
  end subroutine exchanges_create

  subroutine solution_groups_create()
! ******************************************************************************
! Set the solution_groups to be used for the simulation
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use SolutionGroupModule,        only: SolutionGroupType,                   &
                                          solutiongroup_create
    use BaseSolutionModule,         only: BaseSolutionType
    use BaseModelModule,            only: BaseModelType
    use BaseExchangeModule,         only: BaseExchangeType
    use NumericalSolutionModule,    only: solution_create
    use SimVariablesModule,         only: isimdd !PAR
    ! -- dummy
    ! -- local
    type(SolutionGroupType), pointer  :: sgp
    class(BaseSolutionType), pointer  :: sp
    class(BaseModelType), pointer     :: mp
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    integer(I4B) :: isoln
    integer(I4B) :: isgp
    integer(I4B) :: isgpsoln
    integer(I4B) :: sgid
    integer(I4B) :: mid
    character(len=LINELENGTH) :: errmsg
    character(len=LENBIGLINE) :: keyword
    character(len=LINELENGTH) :: fname, mname
    logical :: add !PAR
    type(BlockParserType) :: filein_parser !PAR
    logical :: filein, first !PAR
    logical :: filein_endOfBlock !PAR
    integer(I4B) :: iu !PAR
    logical :: lcgc !CGC
    integer(I4B) :: izc !CGC
    ! -- formats
    character(len=*), parameter :: fmterrmxiter = &
      "('ERROR. MXITER IS SET TO ', i0, ' BUT THERE IS ONLY ONE SOLUTION',     &
        &' IN SOLUTION GROUP ', i0, '. SET MXITER TO 1 IN SIMULATION CONTROL', &
        &' FILE.')"
! ------------------------------------------------------------------------------
    !
    ! -- isoln is the cumulative solution number, isgp is the cumulative
    !    solution group number.
    isoln = 0
    isgp = 0
    !
    !Read through the simulation name file and process each SOLUTION_GROUP
    sgploop: do
      !
      call parser%GetBlock('SOLUTIONGROUP', isfound, ierr, &
        supportOpenClose=.true.)
      if(ierr /= 0) exit sgploop
      if (.not. isfound) exit sgploop
      isgp = isgp + 1
      !
      ! -- Get the solutiongroup id and check that it is listed consecutively.
      sgid = parser%GetInteger()
      if(isgp /= sgid) then
        write(errmsg, '(a)') 'Solution groups are not listed consecutively.'
        call store_error(errmsg)
        write(errmsg, '(a,i0,a,i0)' ) 'Found ', sgid, ' when looking for ',isgp
        call store_error(errmsg)
        call parser%StoreErrorUnit()
        call ustop()
      endif
      !
      ! -- Create the solutiongroup and add it to the solutiongrouplist
      call solutiongroup_create(sgp, sgid)
      call AddSolutionGroupToList(solutiongrouplist, sgp)
      !
      ! -- Begin processing the solution group
      write(iout,'(/1x,a)')'READING SOLUTIONGROUP'
      !
      ! -- Initialize isgpsoln to 0.  isgpsoln is the solution counter for this
      !    particular solution group.  It goes from 1 to the number of solutions
      !    in this group.
      isgpsoln = 0
      lcgc = .false. !CGC
      do
        call parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call parser%GetStringCaps(keyword)
        select case (keyword)
          !
          case ('MXITER')
            sgp%mxiter = parser%GetInteger()
          !
          case ('IMS6','IMS6_CGC') !CGC
            if (keyword == 'IMS6_CGC') then !CGC
              lcgc = .true. !CGC
            end if !CGC
            !
            ! -- Initialize and increment counters
            isoln = isoln + 1
            isgpsoln = isgpsoln + 1
            !
            ! -- Create the solution, retrieve from the list, and add to sgp
            call parser%GetString(fname)
            call solution_create(fname, isoln)
            sp => GetBaseSolutionFromList(basesolutionlist, isoln)
            call sgp%add_solution(isoln, sp)
            !
            ! -- Add all of the models that are listed on this line to
            !    the current solution (sp)
            filein = .false. !PAR
            first = .true. !PAR
            do
              !
              ! -- Set istart and istop to encompass model name. Exit this
              !    loop if there are no more models.
              if (.not.filein) then !PAR
              call parser%GetStringCaps(mname)
              else !PAR
                 call filein_parser%GetNextLine(filein_endOfBlock) !PAR
                 call filein_parser%GetStringCaps(mname) !PAR
              end if !PAR
              !
              if (mname == '') exit
              !
              if ((trim(adjustl(mname)) == 'FILEIN') .and. first) then !PAR
                call parser%GetString(fname) !PAR
                if (fname == '') then !PAR
                  write(errmsg, '(a)') 'ERROR. NO FILE FOUND ' & !PAR
                    //'AFTER "FILEIN" FOR SPECIFYING MODELS.'!PAR
                  call store_error(errmsg) !PAR
                  call parser%StoreErrorUnit() !PAR
                  call ustop() !PAR
                end if
                call openfile(iu, iout, fname, 'MODEL NAMES') !PAR
                call filein_parser%Initialize(iu, iout) !PAR
                call filein_parser%GetBlock('MODELS', isfound, ierr, supportOpenClose=.true.) !IO
                call filein_parser%GetNextLine(filein_endOfBlock) !PAR
                call filein_parser%GetStringCaps(mname) !PAR
                first = .false. !PAR
                filein = .true. !PAR
              end if !PAR
              !
              if (lcgc) then !CGC
                if (filein) then !CGC
                  izc = filein_parser%GetInteger() !CGC
                else !CGC
                  izc = parser%GetInteger() !CGC
                end if !CGC
                mid = ifind(modelname_all, mname) !CGC
                call sp%slnaddmodelid(mid, izc, model_topol_m1, model_topol_m2) !CGC
              end if !CGC
          !
              ! -- Find the model id, and then get model
              if (isimdd ==1) then !PAR
                mid = ifind(modelname_all, mname) !PAR
                call sp%slnmpiaddgmodel(mname) !PAR
                if(mid <= 0) then
                  write(errmsg, '(a,a)') 'Error.  Invalid modelname: ', &
                    trim(mname)
                 call store_error(errmsg)
                 call parser%StoreErrorUnit()
                 call ustop()
                endif
              endif !PAR
              add = .false. !PAR
              mid = ifind(modelname, mname) !PAR
              if (mid > 0) then !PAR
                mp => GetBaseModelFromList(basemodellist, mid)
                add = .true. !PAR
              endif !PAR
              call sp%allocatemodellist() !PAR
              if (add) then !PAR
                ! -- Add the model to the solution
                call sp%addmodel(mp)
                mp%idsoln = isoln
              endif !PAR
            enddo
      !
          case default
            write(errmsg, '(4x,a,a)') &
                  '****ERROR. UNKNOWN SOLUTIONGROUP ENTRY: ', &
                  trim(keyword)
        call store_error(errmsg)
        call parser%StoreErrorUnit()
        call ustop()
        end select
      end do
      !
      if (allocated(model_topol_m1)) then !CGC
        model_topol_m1 = abs(model_topol_m1) !CGC
      end if !CGC
      if (allocated(model_topol_m2)) then !CGC
        model_topol_m2 = abs(model_topol_m2) !CGC
      end if
      !
      ! -- Clean up
      if (filein) then !PAR
        call filein_parser%clear() !PAR
      end if !PAR
    
      ! -- Make sure there is a solution in this solution group
      if(isgpsoln == 0) then
        write(errmsg, '(4x,a,i0)') &
          'ERROR. THERE ARE NO SOLUTIONS FOR SOLUTION GROUP ', isgp
        call store_error(errmsg)
        call parser%StoreErrorUnit()
        call ustop()
      endif
      !
      ! -- If there is only one solution then mxiter should be 1.
      if(isgpsoln == 1 .and. sgp%mxiter > 1) then
        write(errmsg, fmterrmxiter) sgp%mxiter, isgpsoln
        call store_error(errmsg)
        call parser%StoreErrorUnit()
        call ustop()
      endif
    !
      ! -- todo: more error checking?
      !
      write(iout,'(1x,a)')'END OF SIMULATION SOLUTIONGROUP'
      !
    enddo sgploop
    !
    ! -- Check and make sure at least one solution group was found
    if(solutiongrouplist%Count() == 0) then
      call store_error('ERROR.  THERE ARE NO SOLUTION GROUPS.')
      call parser%StoreErrorUnit()
      call ustop()
    endif
    !
    ! -- return
    return
  end subroutine solution_groups_create

  subroutine add_model(im, mtype, mname)
! ******************************************************************************
! Add the model to the list of modelnames, check that the model name is valid.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer, intent(inout) :: im
    character(len=*), intent(in) :: mtype
    character(len=*), intent(inout) :: mname
    ! -- local
    integer :: ilen
    integer :: i
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    im = im + 1
    call expandarray(modelname)
    modelname(im) = mname
    write(iout, '(4x,a,i0)') mtype // ' model ' // trim(mname) //              &
      ' will be created as model ', im  
    !
    ! -- return
    return
  end subroutine add_model
  
  subroutine read_modelname(mname) !PAR
! ******************************************************************************
! Read the model name and check.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    character(len=*), intent(out) :: mname
    ! -- local
    integer(I4B) :: ilen
    integer(I4B) :: i
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    !
    call parser%GetStringCaps(mname)
    ilen = len_trim(mname)
    if (ilen > LENMODELNAME) then
      write(errmsg, '(4x,a,a)')                                                &
            'ERROR. INVALID MODEL NAME: ', trim(mname)
      call store_error(errmsg)
      write(errmsg, '(4x,a,i0,a,i0)')                                          &
            'NAME LENGTH OF ', ilen, ' EXCEEDS MAXIMUM LENGTH OF ',            &
            LENMODELNAME
      call store_error(errmsg)
      call parser%StoreErrorUnit()
      call ustop()
    endif
    do i = 1, ilen
      if (mname(i:i) == ' ') then
        write(errmsg, '(4x,a,a)')                                              &
              'ERROR. INVALID MODEL NAME: ', trim(mname)
        call store_error(errmsg)
        write(errmsg, '(4x,a)')                                                &
              'MODEL NAME CANNOT HAVE SPACES WITHIN IT.'
        call store_error(errmsg)
        call parser%StoreErrorUnit()
        call ustop()
      endif
    enddo
    !
    ! -- return
    return
  end subroutine read_modelname

  subroutine add_model_dd(im, isub, mname) !PAR
! ******************************************************************************
! Add the model to the list of modelnames, check that the model name is valid.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(inout) :: im
    integer(I4B), intent(in) :: isub
    character(len=*), intent(inout) :: mname
    ! -- local
! ------------------------------------------------------------------------------
    !
    im = im + 1
    call expandarray(modelname_all)
    call expandarray(model_sub)
    modelname_all(im) = mname
    model_sub(im) = isub
    !
    ! -- return
    return
  end subroutine add_model_dd
  
  subroutine store_model_topol(m1, m2) !CGC
! ******************************************************************************
! Store model topology (global IDs)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(in) :: m1
    integer(I4B), intent(in) :: m2
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    call expandarray(model_topol_m1)
    call expandarray(model_topol_m2)
    i = size(model_topol_m1)
    model_topol_m1(i) = m1
    model_topol_m2(i) = m2
    !
    ! -- return
    return
  end subroutine store_model_topol
  
end module SimulationCreateModule

!
program mf6
! ******************************************************************************
! Main MODFLOW Version 6 program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- modules
  use KindModule,             only: DP, I4B
  use ConstantsModule,        only: MFVNAM, VERSION, MFTITLE, FMTDISCLAIMER,   &
                                    ISTDOUT, IDEVELOPMODE
  use CompilerVersion
  use InputOutputModule,      only: write_centered
  use SimulationCreateModule, only: simulation_cr, simulation_da
  use TimerModule,            only: start_time, elapsed_time
  use MemoryManagerModule,    only: mem_usage, mem_da
  use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
  use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
  use SolutionGroupModule,    only: SolutionGroupType, GetSolutionGroupFromList
  use ListsModule,            only: basesolutionlist, solutiongrouplist,       &
                                    basemodellist, baseexchangelist,           &
                                    lists_da,                                  &
                                    halomodellist !PAR
  use SimVariablesModule,     only: iout, isimdd !PAR
  use SimModule,              only: converge_reset, converge_check,            &
                                    final_message
  use TdisModule,             only: tdis_tu, tdis_da,                          &
                                    endofsimulation
  use MpiExchangeGenModule,   only: mpi_initialize, serialrun, writestd !PAR
  use MpiExchangeModule,      only: mpi_initialize_world, mpi_world_da,        & !PAR
                                    MpiWorld !PAR
  use MpiExchangeGwfModule,   only: mpi_halo_world, mpi_set_halo_world !PAR
  use NumericalSolutionModule, only: NumericalSolutionType !PAR
  implicit none
  ! -- local
  class(SolutionGroupType), pointer :: sgp
  class(BaseSolutionType),  pointer :: sp
  class(NumericalSolutionType), pointer :: nsp !PAR
  class(BaseModelType),     pointer :: mp
  class(BaseExchangeType),  pointer :: ep
  integer(I4B) :: im, ic, is, isg
  logical :: exit_tsloop
  character(len=80) :: compiler
  ! -- formats
! ------------------------------------------------------------------------------
  !
  ! -- Initialize MPI if required
  call mpi_initialize() !PAR
  call mpi_initialize_world() !PAR
  !
  ! -- Write banner to screen (unit 6) and start timer
  call write_centered('MODFLOW'//MFVNAM, ISTDOUT, 80)
  call write_centered(MFTITLE, ISTDOUT, 80)
  call write_centered('VERSION '//VERSION, ISTDOUT, 80)
  if (.not.serialrun) call write_centered('***RUNNING IN PARALLEL MODE WITH '& !PAR
    //TRIM(MpiWorld%nrprocstr)//' MPI PROCESSES***',ISTDOUT, 80) !PAR
  !
  ! -- Write if develop mode
  if (IDEVELOPMODE == 1) call write_centered('***DEVELOP MODE***', ISTDOUT, 80)
  !
  ! -- Write compiler version
  call get_compiler(compiler)
  call write_centered(' ', ISTDOUT, 80)
  call write_centered(trim(adjustl(compiler)), ISTDOUT, 80)
  !
  ! -- Write disclaimer
  if (writestd) write(ISTDOUT, FMTDISCLAIMER) !PAR
  ! -- get start time
  call start_time()
  !
  !
  ! -- CREATE (CR)
  call simulation_cr()
  !
  !
  ! -- DEFINE (DF)
  ! -- Define each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_df()
  enddo
  !
  ! -- Collective MPI communication scalars DIS
  call mpi_halo_world(2) !PAR
  call mpi_set_halo_world() !PAR
  !
  ! -- Define each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_df()
  enddo
  !
  ! -- Define each solution
  do is = 1, basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_df()
  enddo
  !
  !
  ! -- ALLOCATE AND READ (AR)
  ! -- Allocate and read each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_ar()
  enddo
  ! -- Allocate and read each halo model
  do im = 1, halomodellist%Count() !PAR
    mp => GetBaseModelFromList(halomodellist, im) !PAR
    call mp%model_ar() !PAR
  enddo !PAR
  ! 
  ! -- Local exchange
  if (isimdd == 1) then !PAR
    do is=1,basesolutionlist%Count() !PAR
      sp => GetBaseSolutionFromList(basesolutionlist, is) !PAR
      select type (sp) !PAR
      class is (NumericalSolutionType) !PAR
        nsp => sp !PAR
      end select !PAR
      call nsp%MpiSol%mpi_local_exchange(nsp%name, 'INIT', .true.) !PAR
    enddo !PAR
  endif !PAR
  !
  ! -- Allocate and read each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_ar()
  enddo
  !
  ! -- Allocate and read each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_ar()
  enddo
  !
  !
  ! -- TIME STEP LOOP
  tsloop: do
    !
    ! -- TIME UPDATE (TU)
    call tdis_tu()
    !
    !
    ! -- READ AND PREPARE (RP)
    ! -- Read and prepare each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_rp()
    enddo
    !
    ! -- Read and prepare each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_rp()
    enddo
    !
    ! -- MPI parallel: initialize point-to-point mover
    if (isimdd == 1) then !PAR
      do is=1,basesolutionlist%Count() !PAR
        sp => GetBaseSolutionFromList(basesolutionlist, is) !PAR
        select type (sp) !PAR
        class is (NumericalSolutionType) !PAR
          nsp => sp !PAR
        end select !PAR
        call nsp%slnmpimvrinit(nsp%name) !PAR
      enddo !PAR
    endif
    !
    ! -- Read and prepare each solution
    do is=1,basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_rp()
    enddo
    !
    !
    ! -- CALCULATE (CA)
    call converge_reset()
    do isg = 1, solutiongrouplist%Count()
      sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
      call sgp%sgp_ca()
    enddo
    !
    !
    ! -- OUTPUT (OT)
    ! -- Write output for each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_ot()
    enddo
    !
    ! -- Write output for each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_ot()
    enddo
    !
    ! -- Write output for each solution
    do is=1,basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_ot()
    enddo
    !
    ! -- Time step exit conditions
    call converge_check(exit_tsloop)
    if(exit_tsloop) exit tsloop
    if(endofsimulation) exit tsloop
    !
  enddo tsloop
  !
  !
  ! -- FINAL PROCESSING (FP)
  ! -- Final processing for each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_fp()
  enddo
  !
  ! -- Final processing for each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_fp()
  enddo
  !
  ! -- Final processing for each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_fp()
  enddo
  !
  ! -- DEALLOCATE (DA)
  ! -- Deallocate tdis
  call tdis_da()
  !
  ! -- Deallocate for each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_da()
    deallocate(mp)
  enddo
  !
  ! -- Deallocate for each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_da()
    deallocate(ep)
  enddo
  !
  ! -- Deallocate for each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_da()
    deallocate(sp)
  enddo
  !
  ! -- Deallocate solution group and simulation variables
  do isg = 1, solutiongrouplist%Count()
    sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
    call sgp%sgp_da()
    deallocate(sgp)
  enddo
  call simulation_da()
  call lists_da()
  !
  ! -- Deallocate MPI world
  call mpi_world_da() !PAR
  !
  ! -- Calculate memory usage, elapsed time and terminate
  call mem_usage(iout)
  call mem_da()
  call elapsed_time(iout, 1)
  call final_message()
  !
end program mf6


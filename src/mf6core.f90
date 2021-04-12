module Mf6CoreModule 
  use KindModule,             only: I4B, LGP
  use ListsModule,            only: basesolutionlist, solutiongrouplist, basemodellist, baseexchangelist
  use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
  use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
  use SolutionGroupModule,    only: SolutionGroupType, GetSolutionGroupFromList
  use NumericalSolutionModule, only: NumericalSolutionType !PAR
  use SimVariablesModule,      only: isimdd !PAR
  implicit none  
contains
  
  subroutine Mf6Run
  ! ******************************************************************************
  ! Main MODFLOW Version 6 program.
  ! ******************************************************************************
  !
  !    SPECIFICATIONS:
  ! ------------------------------------------------------------------------------
    ! -- modules 
    use CommandArguments, only: GetCommandLineArguments
    use TdisModule, only: totim, totalsimtime  
    use KindModule, only: DP
    use MpiExchangeModule, only: mpi_initialize_world !PAR
    ! -- dummy
    ! -- local
    logical(LGP) :: hasConverged
    !
    ! -- Initialize MPI if required
    call mpi_initialize_world() !PAR
    ! 
    ! -- parse any command line arguments
    call GetCommandLineArguments()
    !
    ! initialize simulation
    call Mf6Initialize()
    !
    ! -- time loop
    tsloop: do while (totim < totalsimtime)
      
      ! perform a time step
      hasConverged = Mf6Update()
      
      ! if not converged, break
      if(.not. hasConverged) exit tsloop   
      
    enddo tsloop
    !
    ! -- finalize simulation
    call Mf6Finalize()
    
  end subroutine Mf6Run
  
  subroutine Mf6Initialize()
    use SimulationCreateModule, only: simulation_cr
    ! -- dummy
    ! -- local
    
    ! -- print banner and info to screen
    call printInfo()
    
    ! -- create
    call simulation_cr()
    
    ! -- define
    call simulation_df()
       
    ! -- allocate and read
    call simulation_ar()
    
  end subroutine Mf6Initialize
  
  function Mf6Update() result(hasConverged)
    ! -- dummy
    logical(LGP) :: hasConverged
    ! -- local
    !
    ! -- prepare timestep
    call Mf6PrepareTimestep()
    !
    ! -- do timestep
    call Mf6DoTimestep()      
    !
    ! -- after timestep
    hasConverged = Mf6FinalizeTimestep()
    !
  end function Mf6Update
  
  subroutine Mf6Finalize()
    use MpiExchangeGenModule, only: parallelrun, writestd !PAR
    use MpiExchangeModule, only: MpiWorld !PAR
    use, intrinsic :: iso_fortran_env, only: output_unit
    use ListsModule,            only: lists_da
    use MemoryManagerModule,    only: mem_write_usage, mem_da, &
                                      mem_timing, & !TIM
                                      simbytes !PAR
    use TimerModule,            only: elapsed_time
    use SimVariablesModule,     only: iout
    use SimulationCreateModule, only: simulation_da
    use TdisModule,             only: tdis_tu, tdis_da
    use SimModule,              only: final_message
    use MpiExchangeModule,      only: mpi_world_da !PAR
    use MpiExchangeGenModule,   only: writestd !PAR
    ! -- dummy
    ! -- local
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is
    integer(I4B) :: isg
    class(SolutionGroupType), pointer :: sgp => null()
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    
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
    ! -- Write memory usage, elapsed time and terminate
    call mem_write_usage(iout)
    call MpiWorld%mpi_total_memory(simbytes) !PAR
    !
    call mem_timing(iout) !TIM
    call MpiWorld%mpi_barrier() !PAR
    call elapsed_time(iout, 1, writestd) !PAR
    !
    ! -- Deallocate MPI world
    call mpi_world_da() !PAR
    !
    call mem_da()
    !
    call final_message()
    !        
  end subroutine Mf6Finalize
  
  subroutine printInfo()
    use SimVariablesModule, only: istdout
    use VersionModule, only: write_listfile_header
    use TimerModule, only: start_time
    use MpiExchangeModule, only: MpiWorld !PAR
    use MpiExchangeGenModule, only: serialrun, writestd !PAR
    !
    ! -- Write banner to screen (unit stdout) and start timer
    if (writestd) then !PAR
      if (serialrun) then !PAR
        call write_listfile_header(istdout, write_kind_info=.false., &
                                   write_sys_command=.false.)
      else !PAR
        call write_listfile_header(istdout, write_kind_info=.false., & !PAR
                                   write_sys_command=.false.,        & !PAR
                                   nrprocstr=MpiWorld%nrprocstr) !PAR
      end if !PAR
    end if !PAR
    !
    ! -- get start time
    call MpiWorld%mpi_barrier() !PAR
    call start_time(writestd) !PAR
    return
  end subroutine printInfo
  
  subroutine simulation_df()
    use MpiExchangeGwfModule, only: mpi_gwfhalo_world,                           & !PAR
                                    mpi_set_gwfhalo_world_mvr !PAR
    use GwfGwfExchangeModule, only: gwf_mpi_halo_init !PAR
    use MpiExchangeModule,    only: MpiWorld !PAR
    ! -- dummy
    ! -- local
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    
    ! -- Define each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_df()
    enddo
    !
    ! -- Collective MPI communication scalars DIS
    if (isimdd == 1) then !PAR
      call mpi_gwfhalo_world(2) !PAR
      call mpi_set_gwfhalo_world_mvr() !PAR
    endif !PAR
    !
    ! -- Define each exchange phase 1
    do ic = 1, baseexchangelist%Count() !HALO2
      ep => GetBaseExchangeFromList(baseexchangelist, ic) !HALO2
      call ep%exg_df(1) !HALO2
    enddo !HALO2
    !
    ! -- MPI local exchange
    if (isimdd == 1) then !PAR
      call gwf_mpi_halo_init()
      call MpiWorld%mpi_local_exchange('', 'HALO_INIT_DIS', .false.)
    endif
    !
    ! -- Define each exchange phase 2
    do ic = 1, baseexchangelist%Count() !HALO2
      ep => GetBaseExchangeFromList(baseexchangelist, ic) !HALO2
      call ep%exg_df(2) !HALO2
    enddo !HALO2
    ! -- Define each solution
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_df()
    enddo
  
  end subroutine simulation_df
  
  subroutine simulation_ar()  
    use MpiExchangeModule, only: MpiWorld !PAR
    use MpiExchangeGwfModule, only: mpi_gwfhalo_world !PAR
    ! -- dummy
    ! -- local
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(NumericalSolutionType), pointer :: nsp => null() !PAR
    
    ! -- Allocate and read each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_ar()
    enddo
    !
    if (isimdd == 1) then !PAR
      ! -- Local MPI exchange of NPF flags
      call mpi_gwfhalo_world(3)
      !
      !-- Local MPI exchange of NPF arrays
      call MpiWorld%mpi_local_exchange('', 'HALO_INIT_NPFIC', .false.) !PAR
      !
      do is=1,basesolutionlist%Count() !PAR
        sp => GetBaseSolutionFromList(basesolutionlist, is) !PAR
        select type (sp) !PAR
        class is (NumericalSolutionType) !PAR
          nsp => sp !PAR
        end select !PAR
        call nsp%MpiSol%mpi_local_exchange(nsp%name, 'X_IACTIVE', .false.) !PAR
      enddo !PAR
    endif !PAR

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
  end subroutine simulation_ar
  
  subroutine Mf6PrepareTimestep()
    use KindModule,             only: I4B
    use ConstantsModule,        only: LINELENGTH, MNORMAL, MVALIDATE
    use TdisModule,             only: tdis_tu, kstp, kper
    use ListsModule,            only: basemodellist, baseexchangelist
    use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
    use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
    use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
    use SimModule,              only: converge_reset
    use SimVariablesModule,     only: isim_mode
    ! -- dummy
    ! -- local
    class(BaseSolutionType), pointer :: sp => null() !PAR
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(NumericalSolutionType), pointer :: nsp => null() !PAR
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: fmt
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is !PAR
    !
    ! -- initialize fmt
    fmt = "(/,a,/)"
    !
    ! -- time update
    call tdis_tu()
    !
    ! -- set base line
    write(line, '(a,i0,a,i0,a)')                                                 &
      'start timestep kper="', kper, '" kstp="', kstp, '" mode="'
    !
    ! -- evaluate simulation mode
    select case (isim_mode)
      case (MVALIDATE)
        line = trim(line) // 'validate"'
      case(MNORMAL)
        line = trim(line) // 'normal"'
    end select
    
    ! -- Read and prepare each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_message(line, fmt=fmt)
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
    endif !PAR
    !
    ! -- reset simulation convergence flag
    call converge_reset()
    
  end subroutine Mf6PrepareTimestep
  
  subroutine Mf6DoTimestep()
    use KindModule,           only: I4B
    use ListsModule,          only: solutiongrouplist
    use SolutionGroupModule,  only: SolutionGroupType, GetSolutionGroupFromList
    class(SolutionGroupType), pointer :: sgp => null()
    integer(I4B) :: isg
    
    do isg = 1, solutiongrouplist%Count()
      sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
      call sgp%sgp_ca()
    enddo
      
  end subroutine Mf6DoTimestep
  
  function Mf6FinalizeTimestep() result(hasConverged)
    use KindModule,             only: I4B
    use ConstantsModule,        only: LINELENGTH, MNORMAL, MVALIDATE
    use ListsModule,            only: basesolutionlist, basemodellist, baseexchangelist    
    use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
    use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
    use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
    use SimModule,              only: converge_check
    use SimVariablesModule,     only: isim_mode
    ! -- dummy
    logical(LGP) :: hasConverged    
    ! -- local
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: fmt
    integer(I4B) :: im
    integer(I4B) :: ic
    integer(I4B) :: is
    ! -- code
    !
    ! -- initialize format and line
    fmt = "(/,a,/)"
    line = 'end timestep'
    !
    ! -- evaluate simulation mode
    select case (isim_mode)
      case(MVALIDATE)
        !
        ! -- Write final message for timestep for each model 
        do im = 1, basemodellist%Count()
          mp => GetBaseModelFromList(basemodellist, im)
          call mp%model_message(line, fmt=fmt)
        end do
      case(MNORMAL)
        !
        ! -- Write output and final message for timestep for each model 
        do im = 1, basemodellist%Count()
          mp => GetBaseModelFromList(basemodellist, im)
          call mp%model_ot()
          call mp%model_message(line, fmt=fmt)
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
    end select
    !
    ! -- Check if we're done
    call converge_check(hasConverged)
    !
    ! -- return
    return    
    
  end function Mf6FinalizeTimestep
  
end module Mf6CoreModule

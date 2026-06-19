program main
    use modFile, only: &
         read_input &
        ,nob &
        ,maxblock_size &
        ,dtlt &
        ,press => press_in &
        ,temp => temp_in &
        ,vapp => vapp_in &
        ,xlw => xlw_in &
        ,cosz => cosz_in &
        ,last_accepted_dt &
        ,jphoto => jphoto_in&
        ,att => att_in &
        ,n_dyn_chem &
        ,split_method &
        ,chem1_g => chemi_g &
        ,chemo_g &
        ,nspecies_chem_transported &
        ,nspecies_chem_no_transported &
        ,transp_chem_index &
        ,no_transp_chem_index &
        ,chemistry &
        ,block_end

    use mod_chem_spack_rodas3_dyndt, only: &
        chem_rodas3_dyndt

    use modCompare, only: &
        compare_output


    implicit none

    integer :: processor
    real :: time
    integer :: iUnit, err
    logical :: ex
    real :: start_time,end_time
    real :: elapsed_time

    namelist /TARGET/ time,processor

    print *,'0. Chamando a leitura do namelist '
    inquire(file='target.nml', exist=ex)
    if (.not. ex) then
        print *, "Error: file target.nml don't exist!"
        print *, "Please, verify!"
        stop
    end if
    open (newunit=iUnit, FILE='target.nml', STATUS='OLD')
    read (iunit, iostat=err, NML=TARGET)
    if (err /= 0) then
        print *, "Error: reading target.nml!"
        print *, "Please, check syntax!"
        stop   
    else
        print *,"Time      = ",time
        print *,"Processor = ",processor
    end if
    close(iUnit)

    print *,'1. Chamando a leitura do arquivo '
    call read_input(time,processor)

    print *,'2. Chamando o RODAS3 '
    call cpu_time(start_time)

    call chem_rodas3_dyndt( &
        nob &
      , block_end &
      , dtlt &
      , press &
      , temp &
      , vapp &
      , last_accepted_dt &
      , n_dyn_chem &
      , split_method &
      , jphoto &
      , chem1_g &
      , nspecies_chem_transported &
      , nspecies_chem_no_transported &
      , transp_chem_index &
      , no_transp_chem_index &
      , chemistry &
      , maxblock_size &
      )

    call cpu_time(end_time)
    elapsed_time = end_time-start_time
    print *, 'Time to solve: ', elapsed_time, ' [sec]'



    print *,'3. Gerando o arquivo de comparações'
    call compare_output(nob,block_end,chem1_g,chemo_g)

end program main
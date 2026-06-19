module modFile
    use mem_spack, only: &
        alloc_spack, &
        spack

    use chem1_list, only: &
        nspecies

    integer :: n1,n2,n3,ns,nr_photo,nr,maxnspecies,nspecies_chem_transported,nspecies_chem_no_transported, nob, maxblock_size
    real :: dtlt,time_read
    integer :: nVertLevels, nCells, n_dyn_chem
    character(len=20) :: split_method
    real,allocatable :: last_accepted_dt(:)

    type chem1_vars   
        real, allocatable :: sc_p(:,:)
        real, allocatable :: sc_t(:,:) 
        real, allocatable :: sc_t_dyn(:,:)
    end type chem1_vars
    type (chem1_vars)  :: chemi_g(nspecies), chemo_g(nspecies) 

    real, allocatable :: press_in(:,:), temp_in(:,:), vapp_in(:,:), xlw_in(:,:)
    real, allocatable :: att_in(:,:) ,cosz_in(:,:), jphoto(:,:,:)
    real, allocatable :: jphoto_in(:,:,:)

    integer :: transp_chem_index(nSpecies)
    integer :: no_transp_chem_index(nSpecies)
    integer :: nspecies_transported
    integer :: chemistry

    integer, allocatable :: block_end(:)

contains

    subroutine read_input(time,processor)
        implicit none

        real, intent(in) :: time
        integer   , intent(in) :: processor

        character(len=256) :: lixo
        integer :: i, k_, i_, j_, kij_,loop, n, ijk, j, k, ii, ij
        integer :: iBlock
        integer :: n1_,n2_,n3_,n1n2n3, n2n3

        integer :: be, ispc, tb

        integer, parameter :: inob=1
        character(len=30) :: f_name
        real :: ladt

        nspecies_chem_transported = 0
        transp_chem_index   (:)   = 0

        nspecies_chem_no_transported = 0
        no_transp_chem_index(:)      = 0

        print *,'Reading input/output file...',int(time),processor
        write(f_name,fmt='("monan_rodas3_iofil_",I5.5,"_",I1.1,".txt")') int(time),processor
        open(unit=22,file=f_name,status='old',action='read')

        read(22,fmt='(13(I8,1X))') n1,n2,n3,ns,nr_photo,nr,maxnspecies,nspecies_chem_transported &
                                  ,nspecies_chem_no_transported, nob, maxblock_size, chemistry, n_dyn_chem
        read(22,fmt='(A20)') split_method
        if (ns /= nspecies) then
            print *, "Wrong mechanism!"
            stop 555
        end if

        print *, 'Number of blocks = ',nob
        print *, 'Points in i,j,k  = ',n1,n2,n3

        allocate(press_in(maxblock_size,nob), temp_in(maxblock_size,nob), vapp_in(maxblock_size,nob), xlw_in(maxblock_size,nob))
        allocate(jphoto_in(nr_photo,maxblock_size,nob))
        allocate(att_in(maxblock_size,nob) ,cosz_in(maxblock_size,nob))
        allocate(block_end(nob))
        allocate(last_accepted_dt(nob))
        do n=1,ns
            allocate(chemi_g(n)%sc_p(maxblock_size,nob))
            allocate(chemi_g(n)%sc_t_dyn(maxblock_size,nob))
            allocate(chemi_g(n)%sc_t(maxblock_size,nob))
            allocate(chemo_g(n)%sc_p(maxblock_size,nob))
            allocate(chemo_g(n)%sc_t_dyn(maxblock_size,nob))
            allocate(chemo_g(n)%sc_t(maxblock_size,nob))
        end do

        print *,'maxblock_size = ',maxblock_size
        call alloc_spack(chemistry, maxblock_size)

        read(22,fmt='(2(F18.6,1X))') dtlt,time_read
        do ispc=1,nspecies_chem_transported
            read(22,fmt='(I2.2)') transp_chem_index(ispc)
        end do
        read(22,fmt='(A)') lixo
        do ispc = 1, nspecies_chem_no_transported
           read(22,fmt='(I2.2)') no_transp_chem_index(ispc)
        end do

        tb = 0 
        do iblock=1,nob
            read(22,fmt='(A)') lixo
            !print *,'ii = ',iblock,size(block_end)
            read(22,fmt='(I8,1X,F12.6)') be, ladt
            !print *,'block_end=',be
            block_end(iBlock) = be
            tb =  tb + be
            last_accepted_dt(iBlock) = ladt
            !print *,'iBlock=',iBlock," block_end=", be
            do loop = 1, block_end(iBlock)
                read(22,fmt='(6(I6.6,1X))') ij,ijk, kij_, k_, i_, j_
    
                read(22,fmt='(4(E20.10,1X))') press_in(loop,iBlock), temp_in(loop,iBlock), vapp_in(loop,iBlock), ladt

               read(22,fmt='(A)') lixo
               do n=1,nspecies
                    !chem_i é lido e depois modificado dentro do rodas3
                    read(22,fmt='(3(E20.10,1X))') chemi_g(n)%sc_p(loop,iBlock),chemi_g(n)%sc_t_dyn(loop,iBlock)&
                    ,chemi_g(n)%sc_t(loop,iBlock)
               end do
               read(22,fmt='(A)') lixo
               do n=1,nr_photo
                  read(22,fmt='(E20.10)') jphoto_in(n,loop,iBlock)
               end do
            end do
        end do

        read(22,fmt='(A)') lixo
        do iblock = 1, nob !- loop over all blocks
           do loop = 1, block_end(iBlock) !index_g%block_end(i)
                read(22,fmt='(6(I6.6,1X))') ij,ijk, kij_, k_, i_, j_
                do n=1,nspecies
                    read(22,fmt='(3(E20.10,1X))') chemo_g(n)%sc_p(loop,iBlock),chemo_g(n)%sc_t_dyn(loop,iBlock)&
                    ,chemo_g(n)%sc_t(loop,iBlock)
                end do
            end do
        end do

        print *, 'Total of atm points to be computed : ',tb


    end subroutine read_input

end module modFile
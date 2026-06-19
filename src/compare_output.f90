module modCompare
    use modFile, only: &
        chem1_vars

    implicit none

contains

    subroutine compare_output(nob,block_end,chem1_g,chemo_g)

        use chem1_list, only: &
            nSpecies, &
            spc_name

        implicit none 

        integer, intent(in) :: nob
        integer, intent(in) :: block_end(:)

        type(chem1_vars),    intent(inout) :: chem1_g(:)
        type(chem1_vars),    intent(inout) :: chemo_g(:)

        integer :: i,ij,n
        real :: erro_abs, max_err, max_err_p, erro_per

        do n=1,nSpecies
            max_err = 0.0
            max_err_p = 0.0
            print *,' -- Comparando e gerando arquivo '//trim(spc_name(n))//'.csv'
            open(unit=55,file=trim(spc_name(n))//".csv",status="replace",action="write")
            write(55,fmt='(A)') "Cell,level,conc_comp,conc_brams,erro_abs,erro_per"
            do i = 1,nob
                do ij=1,block_end(i)
                    erro_abs = abs(chem1_g(n)%sc_p(ij,i)-chemo_g(n)%sc_p(ij,i))
                    max_err = max(max_err,erro_abs)
                    erro_per = erro_abs/max(chemo_g(n)%sc_p(ij,i),1.0e-30)*100
                    max_err_p = max(max_err_p,erro_per)
                    write(55,fmt='(I6.6,",",I2.2,",",3(E16.8,","),F10.6)') i,ij,chem1_g(n)%sc_p(ij,i),chemo_g(n)%sc_p(ij,i) &
                    ,erro_abs,erro_per
                end do
            end do
!            do k=1,nVertLevels

            close(unit=55)
            write(*,fmt ='("Max error = ",E16.8,", max percent error = ",F10.6,"%")') max_err,max_err_p
        end do


    end subroutine compare_output

end module modCompare
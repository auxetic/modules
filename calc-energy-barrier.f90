program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_mode
    use mo_structure
    implicit none
    complex(16), allocatable, dimension(:) :: psi
    real(8),     allocatable, dimension(:) :: mo_psi
    complex(16) :: ttt

    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! fire
    call init_confire( confire, con )
    call mini_fire_cv( confire )
    con = confire
    
    allocate( psi(con%natom) )
    allocate( mo_psi(con%natom) )
    call calc_psi_n( con, 6, psi )
    do i = 1, con%natom
        mo_psi(i) = abs(psi(i))
    enddo

    open(unit=66, file="psi-6_data.txt")
    do i = 1, con%natom
        write(66,'(es26.16)') mo_psi(i)
    enddo
    close(66)

contains

    subroutine testvar
        implicit none

        sets%natom = 1024
        sets%phi = 0.86d0
        sets%seed = 202
    end subroutine testvar

end program main

program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_dynamic
    implicit none

    ! vars
    call testvar

    ! system
    call init_system( con, sets.natom, sets.phi )
    call gen_rand_config( con, sets.np )

    ! list
    call init_list( nb, con )
    call make_list( nb, con )

    ! md
    call init_md( 1.d-5 )
    call pre_nvt( con )

    ! msd
    call init_msd( con, msd1, 10, 1000 )

    ! main
    do step=1, 100000
        if ( check_list( nb, con ) ) call make_list( nb, con )
        call md_nvt( con, nb )
    end do

    do step=1, 50000

        if ( check_list( nb, con ) ) call make_list( nb, con )

        call md_nvt( con, nb )

        if ( mod( step, 10000 ) == 0 ) print*, step
        call calc_msd( con, msd1 )

    end do

    call save_config_to( con, "./con.dat" )

    call endof_msd( msd1 )

    do step=1, msd1.nhistmax
        print*, 1.d-2 * (step-1) * msd1.ndt, msd1.msd(step)
    end do



contains

    subroutine testvar
        implicit none

        sets.natom = 128
        sets.phi = 0.84d0
        sets.np = 202

    end subroutine testvar

end program main

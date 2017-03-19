module mo_config
    use mo_syst
    implicit none

    type tpcon
        integer :: natom
        ! configuration, velocity, force
        real(8), allocatable, dimension(:,:) :: ra
        real(8), allocatable, dimension(:)   :: r
        real(8), allocatable, dimension(:,:) :: va
        real(8), allocatable, dimension(:,:) :: fa
        ! box
        real(8) :: la(free), lx, ly, lz
        real(8) :: lainv(free), lxinv, lyinv, lzinv
        real(8) :: lav(free), lxv, lyv, lzv
        real(8) :: laf(free), lxf, lyf, lzf
        real(8) :: strain
        real(8) :: strainv, strainf
        real(8) :: phi
        ! sets
        real(8) :: T
        ! property
        real(8) :: Ea, Ek, Ev, stress, press, pressx, pressy
    end type

    type(tpcon) :: con, contemp, contemp2


contains

    subroutine init_system( tcon, tnatom, tphi )
        implicit none

        ! var list
        type(tpcon), intent(inout)        :: tcon
        integer,     intent(in)           :: tnatom
        real(8),     intent(in), optional :: tphi

        tcon.natom = tnatom
        if ( present( tphi ) ) tcon.phi = tphi

        allocate( tcon.ra(free,tnatom), tcon.r(tnatom), tcon.va(free,tnatom), tcon.fa(free,tnatom) )

    end subroutine

    subroutine gen_rand_config( tcon, tnp, tphi )
!       use ifport
        implicit none

        ! var list
        type(tpcon), intent(inout)        :: tcon
        integer,     intent(inout)        :: tnp
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: i, j
        real(8) :: seed_temp

        ! initialized rand
        seed_temp = rand(tnp)
        tnp = 0

        if ( present( tphi ) ) tcon.phi = tphi

        associate(                &
            natom  => tcon.natom, &
            ra     => tcon.ra,    &
            fa     => tcon.fa,    &
            va     => tcon.va,    &
            r      => tcon.r,     &
            la     => tcon.la,    &
            lainv  => tcon.lainv, &
            lx     => tcon.lx,    &
            lxinv  => tcon.lxinv, &
            ly     => tcon.ly,    &
            lyinv  => tcon.lyinv, &
            lz     => tcon.lz,    &
            lzinv  => tcon.lzinv, &
            strain => tcon.strain &
            )

            ! radius
            r(1:natom/2)       = 0.5d0
            r(natom/2+1:natom) = 0.5d0 * ratio

            ! box length
            la     = calc_box_length(tcon)
            lainv  = 1.d0 / la
            lx     = la(1)
            ly     = la(2)
           !lz     = la(3)
            lxinv  = lainv(1)
            lyinv  = lainv(2)
           !lzinv  = lainv(3)
            strain = 0.d0

            ! config
            do i=1, natom
                do j=1, free
                    ra(j,i) = ( rand(0) - 0.5 ) * la(j)
                end do
            end do

            ! f v
            va = 0.d0
            fa = 0.d0

        end associate

    end subroutine

    subroutine gen_lattice_triangle( tcon, tphi )
        implicit none

        ! var list
        type(tpcon), intent(inout)        :: tcon
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: nxy, i, j, k, ii, jj
        real(8) :: a, b
        real(8), parameter :: xyoffset = 1.d-2

        if ( present( tphi ) ) tcon.phi = tphi

        associate(                &
            natom  => tcon.natom, &
            ra     => tcon.ra,    &
            r      => tcon.r,     &
            la     => tcon.la,    &
            lainv  => tcon.lainv, &
            lx     => tcon.lx,    &
            ly     => tcon.ly,    &
            lxinv  => tcon.lxinv, &
            lyinv  => tcon.lyinv, &
            strain => tcon.strain &
            )

            ! box length
            la     = calc_box_trianglelattice( tcon )
            lainv  = 1.d0 / la
            lx     = la(1)
            ly     = la(2)
            lxinv  = lainv(1)
            lyinv  = lainv(2)
            strain = 0.d0

            ! cell numbers and unit
            nxy = nint( sqrt( dble(natom) ) )
            a   = lx / nxy
            b   = a * sqrt(3.d0) / 2.d0

            ! radius
            r = 0.5d0

            ! config
            i=0
            do ii=1, nxy
                do jj=1, nxy
                    i = i + 1
                    ra(1,i) = a * mod( (i-1), nxy ) + a/2.d0 * ( (i-1)/nxy ) + xyoffset
                    ra(2,i) = b * ( (i-1)/nxy ) + xyoffset
                end do
            end do

            do i=1, natom
                ra(:,i) = ra(:,i) - anint( ra(:,i) * lainv ) * la
            end do

        end associate

    end subroutine gen_lattice_triangle

    subroutine read_config( tcon, tfilename, tnatom, tphi )
        implicit none

        ! var list
        type(tpcon),  intent(inout)        :: tcon
        character(*), intent(in)           :: tfilename
        integer,      intent(in)           :: tnatom
        real(8),      intent(in), optional :: tphi

        ! local
        integer :: i

        ! allocate array of tcon
        if ( present( tphi ) ) then
            call init_system( tcon, tnatom, tphi )
        else
            call init_system( tcon, tnatom )
        end if

        ! read config
        open(901,file=tfilename)
            read(901, *) tcon.la, tcon.strain
            tcon.lainv = 1.d0 / tcon.la
            do i=1, tnatom
                read(901,*) tcon.ra(:,i), tcon.r(i)
            end do
        close(901)

        !!!!
        tcon.r = tcon.r * 0.5

        ! phi
        tcon.phi = pi * sum(tcon.r**2) / product(tcon.la)

    end subroutine read_config

    Pure function calc_box_length(tcon) result(l)
        implicit none

        ! var list
        type(tpcon), intent(in) :: tcon
        real(8)                 :: l

        ! local
        real(8) :: phi
        real(8) :: sdisk, volume

        phi = tcon.phi

        ! V_n(r) = pi^(n/2) / Gamma( n/2 + 1 ) * r^n
        sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(tcon.r**free)

        ! box length
        volume = sdisk / phi
        l      = volume ** ( 1.d0 / dble(free) )

    end function

    function calc_box_trianglelattice( tcon ) result( tla )
        implicit none

        ! var list
        type(tpcon), intent(in)  :: tcon
        real(8), dimension(free) :: tla

        ! local
        integer :: nxy

        associate(               &
            natom => tcon.natom, &
            phi   => tcon.phi    &
            )

            nxy = nint( sqrt( dble(natom) ) )
            if ( nxy**2 /= natom ) then
                print*, "wrong natom"; stop
            end if

            tla(1) = sqrt( natom * pi / sqrt(12.d0) / phi )
            tla(2) = tla(1) * sqrt(3.d0) / 2.d0

        end associate

    end function calc_box_trianglelattice

    subroutine trim_config( tcon )
        implicit none

        ! var list
        type(tpcon), intent(inout) :: tcon

        ! local
        integer :: iround(free), cory
        integer :: i, j, k

        associate(                &
            natom  => tcon.natom, &
            ra     => tcon.ra,    &
            la     => tcon.la,    &
            lainv  => tcon.lainv, &
            strain => tcon.strain &
            )

            do i=1, natom

                cory = nint( ra(free,i) * lainv(free) )
                ra(1,i) = ra(1,i) - strain * la(free) * cory

                do k=1, free-1
                    iround(k) = nint( ra(k,i) * lainv(k) )
                end do
                iround(free) = cory

                do k=1, free
                    ra(k,i) = ra(k,i) - iround(k) * la(k)
                end do

            end do

        end associate

    end subroutine trim_config

    function calc_dra( tcon, ti, tj ) result(dra)
        implicit none

        ! para list
        type(tpcon), intent(in) :: tcon
        integer, intent(in) :: ti
        integer, intent(in) :: tj

        ! results
        real(8), dimension(free) :: dra

        ! local
        real(8) :: rai(free), raj(free)
        integer :: k, cory, iround(free)

        associate(                &
            ra     => tcon.ra,    &
            la     => tcon.la,    &
            lainv  => tcon.lainv, &
            strain => tcon.strain &
            )

            rai = ra(:,ti)
            raj = ra(:,tj)

            dra = raj - rai

            cory = nint( dra(free) * lainv(free) )
            dra(1) = dra(1) - strain * la(free) * cory

            do k=1, free-1
                iround(k) = nint( dra(k) * lainv(k) )
            end do
            iround(free) = cory

            do k=1, free
                dra(k) = dra(k) - iround(k) * la(k)
            end do

        end associate
        
    end function calc_dra

    ! ToDo : use hdf5 insteadly
    subroutine save_config_to( tcon, tfilename )
        implicit none

        ! var list
        type(tpcon), intent(in)  :: tcon
        character(*), intent(in) :: tfilename

        ! local
        integer :: i

        associate(                &
            natom  => tcon.natom, &
            ra     => tcon.ra,    &
            r      => tcon.r,     &
            la     => tcon.la,    &
            strain => tcon.strain &
            )

            open(901,file=tfilename)
                write(901,'(3es26.16)') dble(natom), tcon.phi, 0.d0
                write(901,'(3es26.16)') la, strain
                do i=1, natom
                    write(901,'(3es26.16)') ra(:,i), r(i)
                end do
            close(901)

        end associate

    end subroutine save_config_to

    subroutine save_config_debug( tcon, tfilename )
        implicit none

        ! var list
        type(tpcon), intent(in)  :: tcon
        character(*), intent(in) :: tfilename

        ! local
        integer :: i

        associate(                &
            natom  => tcon.natom, &
            ra     => tcon.ra,    &
            va     => tcon.va,    &
            fa     => tcon.fa,    &
            r      => tcon.r,     &
            la     => tcon.la,    &
            lainv  => tcon.lainv, &
            strain => tcon.strain &
            )

            open(901,file=tfilename)
                write(901,'(3es26.16)') dble(natom), tcon.phi, 0.d0
                write(901,'(3es26.16)') la, strain
                write(901,'(3es26.16)') lainv, strain
                do i=1, natom
                    write(901,'(7es26.16)') ra(:,i), r(i), va(:,i), fa(:,i)
                end do
            close(901)

        end associate

    end subroutine save_config_debug

    ! save config with wall
    subroutine save_config_copy( tcon, tfilename )
        implicit none

        ! var list
        type(tpcon),  intent(in) :: tcon
        character(*), intent(in) :: tfilename

        ! local
        integer :: i, ii, jj
        real(8), dimension(free) :: xyoffset

        associate(                &
            natom  => tcon.natom, &
            ra     => tcon.ra,    &
            r      => tcon.r,     &
            la     => tcon.la,    &
            strain => tcon.strain &
            )

            open(901,file=tfilename,status="new")

                write(901,'(3es26.16)') dble(natom), tcon.phi, 0.d0
                write(901,'(3es26.16)') 3*la(1:free), strain

                do ii=-1, 1
                    do jj=-1, 1

                        xyoffset(1) = ii * la(1)
                        xyoffset(2) = jj * la(2)
                        do i=1, natom
                            write(901,'(3es26.16)') ra(:,i)+xyoffset, r(i)
                        end do

                    end do
                end do

            close(901)

        end associate

    end subroutine save_config_copy

end module

module mo_force
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    implicit none

    abstract interface
        subroutine abstract_force( tcon, tnb )
            import :: tpcon, tplist
            type(tpcon)  :: tcon
            type(tplist) :: tnb
        end subroutine
    end interface

contains

    subroutine calc_force( tcon, tnb )
        implicit none

        type(tpcon)  :: tcon
        type(tplist) :: tnb

        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnb%list       &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)
                    iround = list(i)%iround(:,jj)
                    cory = list(i)%cory(jj)

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    wilixyz = wilixyz + fr * dra(:)**2

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)

        end associate
    end subroutine

    subroutine calc_force_gel( tcon, tnb )
        implicit none

        type(tpcon)  :: tcon
        type(tplist) :: tnb

        ! gel
        real(8), parameter :: rcut = 0.1d0
        real(8), parameter :: kout = 0.5d0
        real(8) :: rc

        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        real(8) :: dijrcut, dijrc
        integer :: iround(free), cory
        integer :: i, j, k, jj

        rc = rcut * kout / ( 1.d0 + kout )

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnb%list       &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)
                    iround = list(i)%iround(:,jj)
                    cory = list(i)%cory(jj)

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj
                    dijrcut = dij * ( 1.d0 + rcut )
                    dijrc   = dij * ( 1.d0 + rc   )

                    if ( rij2 > dijrcut**2 ) cycle

                    rij = sqrt( rij2 )

                    if ( rij < dij ) then
                        Ea = Ea + ( 1.d0 - rij/dij )**2/2.d0
                        wij = (1.d0 - rij/dij) * rij / dij
                    elseif ( rij < dijrc ) then
                        Ea = Ea + kout*( 1.d0 - rij/dij )**2/2.d0
                        wij = kout*(1.d0 - rij/dij) * rij / dij
                    elseif ( rij < dijrcut ) then
                        wij = - kout * ( (dijrcut - rij)/dij ) * rij / dij
                    end if

                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    wilixyz = wilixyz + fr * dra(:)**2

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)

        end associate
    end subroutine

    subroutine calc_force_pin( tcon, tnb )
        implicit none

        type(tpcon)  :: tcon
        type(tplist) :: tnb

        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            pinflag  => tcon%pinflag,  &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnb%list       &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)
                    iround = list(i)%iround(:,jj)
                    cory = list(i)%cory(jj)

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    if ( pinflag(i) == 0 .and. pinflag(j) == 0 ) then
                        wilixyz = wilixyz + fr * dra(:)**2

                        stress = stress - 2 * dra(1) * dra(2) * fr
                    end if

                end do

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)

        end associate
    end subroutine

    subroutine calc_force_withoutlist( tcon )
        implicit none

        ! para list
        type(tpcon)  :: tcon

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: lainv(free), ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        integer :: iround(free), cory
        integer :: i, j, k

        associate(                    &
            natom    => tcon%natom,   &
            radius   => tcon%r,       &
            ra       => tcon%ra,      &
            fa       => tcon%fa,      &
            Ea       => tcon%Ea,      &
            la       => tcon%la,      &
            strain   => tcon%strain,  &
            stress   => tcon%stress,  &
            press    => tcon%press,   &
            pressxyz => tcon%pressxyz &
            )

            lainv = 1.d0 / la

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = radius(j)

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

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(1)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)

        end associate
    end subroutine

    subroutine calc_force_spring( tcon, tnet )
        implicit none

        ! para list
        type(tpcon),     intent(inout) :: tcon
        type(tpnetwork), intent(in)    :: tnet

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: rij2, rij, l0, fr, wij, wili, wilixyz(free), ks
        integer :: iround(free), cory
        integer :: i, j, k, ii

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnet%sps,      &
            nlist    => tnet%nsps      &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

            do ii=1, nlist

                i      = list(ii)%i
                j      = list(ii)%j
                cory   = list(ii)%cory
                iround = list(ii)%iround
                l0     = list(ii)%l0
                ks     = list(ii)%ks

                rai = ra(:,i)
                raj = ra(:,j)

                dra = raj - rai
                dra(1) = dra(1) - cory * strain * la(free)

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )
                rij  = sqrt( rij2 )

                Ea = Ea + 0.5d0 * ks * ( l0 - rij )**2

                wij  = ks * ( l0 - rij ) * rij
                wili = wili + wij

                fr = wij / rij2

                wilixyz = wilixyz + fr * dra(:)**2

                fa(:,j) = fa(:,j) + fr * dra
                fa(:,i) = fa(:,i) - fr * dra

                stress = stress - 2 * dra(1) * dra(2) * fr  ! 3d ?

            end do

            stress   = stress  / product(la) / free
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)

        end associate
    end subroutine

end module

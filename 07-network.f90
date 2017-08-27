module mo_network
    use mo_syst
    use mo_config
    implicit none

    type tpspring
        integer :: i, j
        integer :: cory, iround(free)
        real(8) :: l0, lvec(free)
        real(8) :: ks
        real(8) :: Bi, Gi, Gixy, Gis
    end type

    type tpnetwork
        integer :: nsps, max_of_springs
        integer :: natom
        real(8) :: mb, mg, mgs, mgxy
        real(8) :: bcorr, gcorr, kscorr, ksmeancorr
        type(tpspring), allocatable, dimension(:) :: sps
    end type

    type(tpnetwork) :: net

contains

    subroutine init_network( tnetwork, tcon )
        implicit none

        ! para list
        type(tpcon),     intent(in)    :: tcon
        type(tpnetwork), intent(inout) :: tnetwork

        associate(                                    &
            natom          => tnetwork%natom,         &
            max_of_springs => tnetwork%max_of_springs &
            )

            natom = tcon%natom
            max_of_springs = natom * 10
            allocate( tnetwork%sps(max_of_springs) )

        end associate
    end subroutine init_network

    subroutine make_network( tnetwork, tcon )
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon

        ! local
        integer :: i, j, k
        real(8) :: dra(free), rij2, dij
        integer :: cory, iround(free)

        associate(                   &
            natom  => tcon%natom,    &
            ra     => tcon%ra,       &
            r      => tcon%r,        &
            la     => tcon%la,       &
            lainv  => tcon%lainv,    &
            strain => tcon%strain,   &
            nsps   => tnetwork%nsps, &
            sps    => tnetwork%sps   &
            )

            nsps = 0
            do i=1, natom
                do j=i+1, natom

                    dra = ra(:,j) - ra(:,i)

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = r(i) + r(j)

                    if ( rij2 > dij**2 ) cycle

                    nsps             = nsps + 1
                    sps(nsps)%i      = i
                    sps(nsps)%j      = j
                    sps(nsps)%cory   = cory
                    sps(nsps)%iround = iround
                    sps(nsps)%lvec   = dra
                    sps(nsps)%l0     = sqrt(rij2)
                    sps(nsps)%ks     = 1.d0
                    sps(nsps)%Bi     = 0.d0
                    sps(nsps)%Gi     = 0.d0
                    sps(nsps)%Gis    = 0.d0
                    sps(nsps)%Gixy   = 0.d0

                end do
            end do

        end associate

        ! reallocate array of sps
        tnetwork%sps = tnetwork%sps(1:tnetwork%nsps)
    end subroutine make_network

    subroutine remake_network( tnetwork, tcon )
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon

        ! local
        integer :: i, j, k, ii
        real(8) :: dra(free), rij2, dij
        integer :: cory, iround(free)

        associate(                   &
            natom  => tcon%natom,    &
            ra     => tcon%ra,       &
            r      => tcon%r,        &
            la     => tcon%la,       &
            lainv  => tcon%lainv,    &
            strain => tcon%strain,   &
            nsps   => tnetwork%nsps, &
            sps    => tnetwork%sps   &
            )

            do ii=1, nsps

                i = sps(ii)%i
                j = sps(ii)%j

                dra = ra(:,j) - ra(:,i)

                cory = nint( dra(free) * lainv(free) )
                dra(1) = dra(1) - cory * strain * la(free)

                do k=1, free-1
                    iround(k) = nint( dra(k) * lainv(k) )
                end do
                iround(free) = cory

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )

                sps(ii)%cory   = cory
                sps(ii)%iround = iround
                sps(ii)%lvec   = dra
                sps(ii)%l0     = sqrt(rij2)

            end do

        end associate
    end subroutine remake_network

    subroutine change_k_spring( tnetwork, tcon, opktan, oph )
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon
        real(8),         optional      :: opktan
        real(8),         optional      :: oph

        ! local
        integer :: i, itemp
        real(8) :: ktan, h, temp


        ! 1. ks = ls
!       do i=1, tnetwork%nsps
!           tnetwork%sps(i)%ks = calc_spring_len( tcon, i, tnetwork )
!       end do

!       do i=1, tnetwork%nsps
!           itemp = floor( rand(0) * tnetwork%nsps ) + 1
!           tnetwork%sps(i)%ks = calc_spring_len( tcon, itemp, tnetwork )
!       end do

        ktan = 1.d0
        if ( present( opktan ) ) ktan = opktan
!        do i=1, tnetwork%nsps
!            tnetwork%sps(i)%ks = 1.0 + tanh( ktan * ( calc_spring_len( tcon, i, tnetwork ) - 1.d0 ) )
!!           tnetwork%sps(i)%ks = 1.0 + tanh( ktan * ( tnetwork%sps(i)%l0 - 1.d0 ) )
!        end do
        do i=1, tnetwork%nsps
            call random_number(temp)
            itemp = floor( temp*tnetwork%nsps ) + 1
!           tnetwork%sps(i)%ks = 1.0 + tanh( ktan * ( calc_spring_len( tcon, itemp, tnetwork ) - 1.d0 ) )
            tnetwork%sps(i)%ks = 1.0 + tanh( ktan * ( tnetwork%sps(i)%l0 - 1.d0 ) )
        end do

!       h = 0.1d0
!       if ( present( oph ) ) h = oph
!       do i=1, tnetwork%nsps
!           if (  calc_spring_len( tcon, i, tnetwork ) > 1.d0 ) then
!               tnetwork%sps(i)%ks = 2.d0 - h
!           else
!               tnetwork%sps(i)%ks = h
!           end if
!       end do
    end subroutine change_k_spring

    subroutine change_k_spring_2( tnetwork, tcon, opktan, oph )
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon
        real(8),         optional      :: opktan
        real(8),         optional      :: oph

        ! local
        integer :: i
        real(8) :: h, ktan

        h = 0.1
        if ( present( oph ) ) h = oph

!       do i=1, tnetwork%nsps
!           if ( rand(0) > 0.5 ) then
!               tnetwork%sps(i)%ks = 2.d0 - h
!           else
!               tnetwork%sps(i)%ks = h
!           end if
!       end do

!       do i=1, tnetwork%nsps
!           if ( calc_len( tcon, i, tnetwork ) > 1.d0 ) then
!               tnetwork%sps(i)%ks = 2.d0 - h
!           else
!               tnetwork%sps(i)%ks = h
!           end if
!       end do

        ktan = 1.d0
        if ( present( opktan ) ) ktan = opktan
        do i=1, tnetwork%nsps
            tnetwork%sps(i)%ks = 1.d0 + tanh( ktan * ( calc_spring_len( tcon, i, tnetwork ) - 1.d0 ) )
        end do
    end subroutine change_k_spring_2

    function calc_spring_len( tcon, tibond, tnetwork  ) result(tl)
        implicit none

        ! para list
        type(tpcon), intent(in)     :: tcon
        integer, intent(in)         :: tibond
        type(tpnetwork), intent(in) :: tnetwork
        real(8) :: tl

        ! local
        integer :: i, j, k, ii
        real(8) :: dra(free), rij2, dij
        integer :: cory, iround(free)

        associate(                   &
            natom  => tcon%natom,    &
            ra     => tcon%ra,       &
            r      => tcon%r,        &
            la     => tcon%la,       &
            lainv  => tcon%lainv,    &
            strain => tcon%strain,   &
            nsps   => tnetwork%nsps, &
            sps    => tnetwork%sps   &
            )


            i = sps(tibond)%i
            j = sps(tibond)%j

            dra = ra(:,j) - ra(:,i)

            cory = nint( dra(free) * lainv(free) )
            dra(1) = dra(1) - cory * strain * la(free)

            do k=1, free-1
                iround(k) = nint( dra(k) * lainv(k) )
            end do
            iround(free) = cory

            do k=1, free
                dra(k) = dra(k) - iround(k) * la(k)
            end do

            !tl = norm2( dra )
            tl = sqrt(sum(dra**2))

        end associate
    end function

    subroutine calc_Bi( tcon, tnetwork, test )
        implicit none

        ! para list
        type(tpcon),     intent(in)    :: tcon
        type(tpnetwork), intent(inout) :: tnetwork
        real(8),         intent(in)    :: test

        ! local
        integer :: ii, i, j
        real(8) :: ks, Es, lnow, l0, Bi, sumEs

        sumEs = 0.d0

        associate(                 &
            nsps => tnetwork%nsps, &
            sps  => tnetwork%sps,  &
            mb   => tnetwork%mb    &
            )

            do ii=1, nsps

                i = sps(ii)%i
                j = sps(ii)%j

                l0 = sps(ii)%l0
                ks = sps(ii)%ks

                lnow = calc_spring_len( tcon, ii, tnetwork )
                Es = 0.5d0 * ks * ( lnow - l0 )**2
                sumEs = sumEs + Es

                ! B = 1/4 * sum( T[1:free,1:free] )
                ! Es = 1/2 * Bv ( 2 dl/l )**2
                Bi = 2.d0 * Es / ( 2.d0 * test )**2

                sps(ii)%Bi = Bi

            end do

            mb = 2.d0 * sumEs / ( 2.d0 * test )**2 * product(tcon%lainv)

        end associate
    end subroutine

    subroutine calc_Gis( tcon, tnetwork, test )
        implicit none

        ! para list
        type(tpcon),     intent(in)    :: tcon
        type(tpnetwork), intent(inout) :: tnetwork
        real(8),         intent(in)    :: test

        ! local
        integer :: ii, i, j
        real(8) :: ks, Es, lnow, l0, Gis, sumEs

        sumEs = 0.d0

        associate(                 &
            nsps => tnetwork%nsps, &
            sps  => tnetwork%sps,  &
            mgs  => tnetwork%mgs   &
            )

            do ii=1, nsps

                i = sps(ii)%i
                j = sps(ii)%j

                l0 = sps(ii)%l0
                ks = sps(ii)%ks

                lnow = calc_spring_len( tcon, ii, tnetwork )
                Es = 0.5d0 * ks * ( lnow - l0 )**2
                sumEs = sumEs + Es

                ! Gs = T_xyxy
                ! Es = 1/2 * Gv strain**2
                Gis = 2.d0 * Es / test**2

                sps(ii)%Gis = Gis

            end do

            mgs = 2.d0 * sumEs / test**2 * product(tcon%lainv)

        end associate
    end subroutine

    subroutine calc_Gixy( tcon, tnetwork, test )
        implicit none

        ! para list
        type(tpcon),     intent(in)    :: tcon
        type(tpnetwork), intent(inout) :: tnetwork
        real(8),         intent(in)    :: test

        ! local
        integer :: ii, i, j
        real(8) :: ks, Es, lnow, l0, Gixy, sumEs

        sumEs = 0.d0

        associate(                 &
            nsps => tnetwork%nsps, &
            sps  => tnetwork%sps,  &
            mgxy => tnetwork%mgxy, &
            mgs  => tnetwork%mgs,  &
            mg   => tnetwork%mg    &
            )

            do ii=1, nsps

                i = sps(ii)%i
                j = sps(ii)%j

                l0 = sps(ii)%l0
                ks = sps(ii)%ks

                lnow = calc_spring_len( tcon, ii, tnetwork )
                Es = 0.5d0 * ks * ( lnow - l0 )**2
                sumEs = sumEs + Es

                ! Gs = 1/4 [ T_xxxx + T_yyyy - 2T_xxyy ]
                ! Es = 1/2 [] * test**2 * v
                !    = 2 * Gs * v * test**2
                Gixy = 0.5d0 * Es / test**2

                sps(ii)%Gixy = Gixy
                sps(ii)%gi   = 0.5d0 * ( Gixy + sps(ii)%Gis)

            end do

            mgxy = 0.5d0 * sumEs / test**2 * product(tcon%lainv)

            mg = 0.5d0 * ( mgxy + mgs )

        end associate
    end subroutine

! temp
    subroutine calc_corr( tnetwork )
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork

        ! local
        integer :: ii, jj, flag
        real(8) :: bi, bj, gi, gj, ksi, ksj, ksmean
        real(8) :: vec1(free), vec2(free)

        associate(                            &
            nsps       => tnetwork%nsps,      &
            sps        => tnetwork%sps,       &
            bcorr      => tnetwork%bcorr,     &
            gcorr      => tnetwork%gcorr,     &
            kscorr     => tnetwork%kscorr,    &
            ksmeancorr => tnetwork%ksmeancorr &
            )

            bcorr = 0.d0; gcorr = 0.d0; kscorr = 0.d0

            ksmean = sum( sps(:)%ks ) / nsps

            do ii=1, nsps
                vec1 = sps(ii)%lvec
                bi   = sps(ii)%bi
                gi   = sps(ii)%gi
                ksi  = sps(ii)%ks
                do jj=ii+1, nsps

                    vec2 = sps(jj)%lvec
                    bj   = sps(jj)%bi
                    gj   = sps(jj)%gi
                    ksj  = sps(jj)%ks

                    flag = 0
                    if ( sps(ii)%i == sps(jj)%i .or. sps(ii)%j == sps(jj)%j ) then
                        vec2 = vec2; flag = 1
                    end if
                    if ( sps(ii)%j == sps(jj)%i ) then
                        vec2 = -vec2; flag = 1
                    end if

                    if ( flag == 0 ) cycle
                    bcorr  = bcorr  + sum(vec1*vec2) * bi * bj
                    gcorr  = gcorr  + sum(vec1*vec2) * gi * gj
                    kscorr = kscorr + sum(vec1*vec2) * ksi * ksj
                    ksmeancorr = ksmeancorr + sum(vec1*vec2) * (ksi-ksmean) * (ksj-ksmean)

                end do
            end do

            bcorr  = bcorr  / nsps
            gcorr  = gcorr  / nsps
            kscorr = kscorr / nsps

        end associate
    end subroutine

end module

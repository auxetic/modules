module mo_structure
    use mo_syst
    use mo_config
    use mo_list
    use mo_math
    implicit none

contains

    function calc_psi(tcon, n) result(psi)
        implicit none

        ! para list
        type(tpcon), intent(in) :: tcon
        integer,     intent(in) :: n

        ! result
        complex(16), dimension(tcon%natom) :: psi

        ! lcoal
        type(tpvoro) :: voroni
        real(8)      :: nij(free)
        complex(16)  :: cpxtemp
        integer      :: i, j, k

        call init_voro( voroni, tcon )
        call calc_voro( voroni )

        psi(:) = cmplx(0.0, 0.0)

        do i=1, tcon%natom
            do k=1, voroni%list(i)%nbsum
                j       = voroni%list(i)%nblist(k)
                nij     = unitv( tcon%dra(i,j) )
                cpxtemp = cmplx(nij(1),nij(2))
                psi(i)  = psi(i) + cpxtemp**n
            end do
        end do
        psi(:)  = psi(:) / voroni%list(:)%nbsum
    end function

end module

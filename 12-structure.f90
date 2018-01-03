module mo_structure
    use mo_syst
    use mo_config
    use mo_list
    use mo_math

contains

    subroutine calc_psi_n( tcon, n, psi )
        implicit none

        !var list
        type(tpcon),                        intent(in)    :: tcon
        integer,                            intent(in)    :: n
        complex(16), dimension(tcon%natom), intent(inout) :: psi

        !lcoal
        type(tpvoro) :: voroni
        integer :: i, j, k
        real(8) :: angle

        call init_voro( voroni, tcon )
        call calc_voro( voroni )
        
        do i = 1, tcon%natom
            psi(i) = (0.0, 0.0)
            do k = 1, voroni%list(i)%nbsum
                j      = voroni%list(i)%nblist(k)
                angle  = atan2( tcon%ra(1,i) - tcon%ra(1,j), tcon%ra(2,i) - tcon%ra(2,j) )
                psi(i) = psi(i) + cmplx( cos(n*angle), sin(n*angle) )
            enddo
            psi(i)  = psi(i) / voroni%list(i)%nbsum
        enddo
    end subroutine

end module

SUBROUTINE calcener_MC(mover, dx_t, dy_t, dz_t, energy_diff)
    USE library, ONLY : dtype, max_arr, min_arr, print_array, delta
    USE system, ONLY : x, y, z, sp, natoms, Lx, Ly, Lz
    USE potential, ONLY : c_A, c_B, c_lam, c_mu, c_X, c_R, &
                            c_h, c_n, c_R2, c_S, c_S2, c_pRSr, &
                            c_betan, c_d2, c_dr2, c_c2, c_n2r 
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mover
    REAL(dtype), INTENT(IN) :: dx_t, dy_t, dz_t
    REAL(dtype), INTENT(OUT) :: energy_diff
    ! PARAMETERS AND LOCAL VARIABLES FOR COMPUTATIONS
    INTEGER, PARAMETER :: nmax=20, nmax_out=100
    REAL(dtype), DIMENSION(nmax) :: fC, fA, fR, dx, dy, dz, r, rr
    REAL(dtype), DIMENSION(nmax) ::  x_in, y_in, z_in, x_crw, y_crw, z_crw
    REAL(dtype), DIMENSION(nmax_out) :: x_out, y_out, z_out
    INTEGER, DIMENSION(nmax) :: spq_in, atom_in, spq_crw, atom_crw, tmp_atoms 
    INTEGER, DIMENSION(nmax_out) :: spq_out, atom_out
    REAL(dtype) :: ei_in, ei_fin
    ! LOCAL VARIABLES
    INTEGER :: index_ij, spi
    INTEGER :: i, ind, j_in, j_out, j_out1, j_out_, q_in, q_out, q_crw
    REAL(dtype) :: c2, d2, dr2, h, n, betan, n2r
    REAL(dtype) :: shiftx, shifty, shiftz, r2, a, gthi, zetaij
    REAL(dtype) :: rij_dot_rik, rrij_rrik, gThden, bZijn, bij
    REAL(dtype) :: dx_tmp, dy_tmp, dz_tmp, cThijk, h_cThijk, gThijk
    REAL(dtype) :: xt, yt, zt
    
    ei_in   = 0.
    ei_fin  = 0.

    print*, x(mover), dx_t

    xt = x(mover) + dx_t
    yt = y(mover) + dy_t
    zt = z(mover) + dz_t
    
    spi = sp(mover)
    c2      = c_c2(spi)
    d2      = c_d2(spi)
    dr2     = c_dr2(spi)
    h       = c_h(spi)
    n       = c_n(spi)
    betan   = c_betan(spi)
    n2r     = c_n2r(spi)
    
    atom_in(1)  = mover 
    x_in(1)     = x(mover)
    y_in(1)     = y(mover)
    z_in(1)     = z(mover)
    spq_in(1)   = spi
   
    atom_out(1)  = mover 
    x_out(1)     = x(mover)
    y_out(1)     = y(mover)
    z_out(1)     = z(mover)
    spq_out(1)   = spi

    q_in = 2
    q_out = 2
    q_crw = 1

    DO i = 1, natoms
    ! print*, i
        dx_tmp = x(i) - x_in(1)
        dy_tmp = y(i) - y_in(1)
        dz_tmp = z(i) - z_in(1)

        ! Periodic boundary conditions (mover is the central atom of the system)       
        shiftx=0.d0
        IF (dx_tmp .gt. Lx/2) THEN
            shiftx=-Lx
        ELSE IF (dx_tmp .lt. -Lx/2) THEN
            shiftx=Lx
        END IF
        shifty=0.d0
        IF (dy_tmp .gt. Ly/2) THEN
            shifty=-Ly
        ELSE IF (dy_tmp .lt. -Ly/2) THEN
            shifty=Ly
        END IF
        shiftz=0.d0
        IF (dz_tmp .gt. Lz/2) THEN
            shiftz=-Lz
        ELSE IF (dz_tmp .lt. -Lz/2) THEN
            shiftz=Lz
        END IF

        ! New distances
        dx_tmp = dx_tmp + shiftx
        dy_tmp = dy_tmp + shifty
        dz_tmp = dz_tmp + shiftz

        r2 = dx_tmp*dx_tmp + dy_tmp*dy_tmp + dz_tmp*dz_tmp
        index_ij = sp(i) + spi - 1

        ! Find atoms in 2S, S sphere and [2S, 2S+delta] crown
        IF (r2 .le. 4. * max_arr(c_S2) .and. r2 .gt. 0.1) THEN
            x_out(q_out)    = x(i) + shiftx
            y_out(q_out)    = y(i) + shifty
            z_out(q_out)    = z(i) + shiftz
            spq_out(q_out)  = sp(i)
            atom_out(q_out) = i 
            q_out = q_out + 1
            IF (r2 .le. c_S2(index_ij)) THEN 
                x_in(q_in)    = x(i) + shiftx
                y_in(q_in)    = y(i) + shifty
                z_in(q_in)    = z(i) + shiftz
                spq_in(q_in)  = sp(i)
                atom_in(q_in) = i 
                q_in = q_in + 1
            END IF 
        ELSE IF ( r2 .le. (2. * max_arr(c_S) + delta) ** 2 .and. r2 .gt. 0.1) THEN
            x_crw(q_crw)    = x(i) + shiftx
            y_crw(q_crw)    = y(i) + shifty
            z_crw(q_crw)    = z(i) + shiftz
            spq_crw(q_crw)  = sp(i)
            atom_crw(q_crw) = i 
            q_crw = q_crw + 1
        END IF
    END DO

    q_in  =  q_in - 1
    q_out = q_out - 1
    q_crw = q_crw - 1

    DO j_in = 1, q_in
        ind = 1
        spi = spq_in(j_in)
        c2 = c_c2(spi)
        d2 = c_d2(spi)
        dr2 = c_dr2(spi)
        h = c_h(spi)
        n = c_n(spi)
        betan = c_betan(spi)
        n2r = c_n2r(spi)

        DO j_out = 1, q_out
            IF (atom_in(j_in) .ne. atom_out(j_out)) THEN
                dx(ind) = x_in(j_in) - x_out(j_out)
                dy(ind) = y_in(j_in) - y_out(j_out)
                dz(ind) = z_in(j_in) - z_out(j_out)

                r2 = dx(ind)*dx(ind) + dy(ind)*dy(ind) + dz(ind)*dz(ind)
                index_ij = spq_in(j_in) + spq_out(j_out) - 1
                IF (r2 .le. c_S2(index_ij) .and. r2 .gt. 0.1) THEN
                    tmp_atoms(ind) = j_out
                    r(ind) = sqrt(r2)
                    rr(ind) = 1./r(ind)
                    fA(ind) = -c_B(index_ij)*exp(-c_mu(index_ij)*r(ind))
                    fR(ind) = c_A(index_ij)*exp(-c_lam(index_ij)*r(ind))
                    IF (r2 .gt. c_R2(index_ij)) THEN
                        a = c_pRSr(index_ij)*(r(ind) - c_R(index_ij))
                        fC(ind) = 0.5 + 0.5*cos(a)
                    ELSE
                        fc(ind) = 1.
                    END IF
                    ind = ind + 1
                END IF
            END IF
        END DO
        ! stop
        ind = ind - 1

        gthi = 1. + c2*dr2 
        DO j_out = 1, ind 
            index_ij = spq_in(j_in) + spq_out(tmp_atoms(j_out)) - 1
            zetaij = 0.
            DO j_out_ = j_out+1, j_out+ind-1
                j_out1 = mod(j_out_-1, ind) + 1
                rij_dot_rik     = dx(j_out)*dx(j_out1) + dy(j_out)*dy(j_out1) + dz(j_out)*dz(j_out1)
                rrij_rrik       = rr(j_out) * rr(j_out1)
                cThijk          = rij_dot_rik * rrij_rrik
                h_cThijk        = h - cThijk
                gThden          = 1./(d2 + h_cThijk * h_cThijk)
                gThijk          = gthi - c2 * gThden
                zetaij          = zetaij + fC(j_out1) * gthijk
            END DO
            bZijn   = 1. + betan*(zetaij**n)
            bij     = c_X(index_ij) * (bZijn**n2r)
            ei_in  = ei_in + fC(j_out)*(fR(j_out) + bij*fA(j_out))
        END DO
    END DO
    
    IF (q_crw .eq. 0) THEN
        x_out(q_out+1:q_out+q_crw) = x_crw(1:q_crw)
        y_out(q_out+1:q_out+q_crw) = y_crw(1:q_crw)
        z_out(q_out+1:q_out+q_crw) = z_crw(1:q_crw)
        spq_out(q_out+1:q_out+q_crw) = spq_crw(1:q_crw)
        atom_out(q_out+1:q_out+q_crw) = atom_crw(1:q_crw)
    END IF


    x_in(1)     = xt
    y_in(1)     = yt
    z_in(1)     = zt
    x_out(1)    = xt
    y_out(1)    = yt
    z_out(1)    = zt

    DO j_in = 1, q_in
        ind = 1
        spi = spq_in(j_in)
        c2 = c_c2(spi)
        d2 = c_d2(spi)
        dr2 = c_dr2(spi)
        h = c_h(spi)
        n = c_n(spi)
        betan = c_betan(spi)
        n2r = c_n2r(spi)

        DO j_out = 1, q_out
            IF (atom_in(j_in) .ne. atom_out(j_out)) THEN
                dx(ind) = x_in(j_in) - x_out(j_out)
                dy(ind) = y_in(j_in) - y_out(j_out)
                dz(ind) = z_in(j_in) - z_out(j_out)

                r2 = dx(ind)*dx(ind) + dy(ind)*dy(ind) + dz(ind)*dz(ind)
                index_ij = spq_in(j_in) + spq_out(j_out) - 1
                IF (r2 .le. c_S2(index_ij) .and. r2 .gt. 0.1) THEN
                    tmp_atoms(ind) = j_out
                    r(ind) = sqrt(r2)
                    rr(ind) = 1./r(ind)
                    fA(ind) = -c_B(index_ij)*exp(-c_mu(index_ij)*r(ind))
                    fR(ind) = c_A(index_ij)*exp(-c_lam(index_ij)*r(ind))
                    IF (r2 .gt. c_R2(index_ij)) THEN
                        a = c_pRSr(index_ij)*(r(ind) - c_R(index_ij))
                        fC(ind) = 0.5 + 0.5*cos(a)
                    ELSE
                        fc(ind) = 1.
                    END IF
                    ind = ind + 1
                END IF
            END IF
        END DO
        ! stop
        ind = ind - 1

        gthi = 1. + c2*dr2 
        DO j_out = 1, ind 
            index_ij = spq_in(j_in) + spq_out(tmp_atoms(j_out)) - 1
            zetaij = 0.
            DO j_out_ = j_out+1, j_out+ind-1
                j_out1 = mod(j_out_-1, ind) + 1
                rij_dot_rik     = dx(j_out)*dx(j_out1) + dy(j_out)*dy(j_out1) + dz(j_out)*dz(j_out1)
                rrij_rrik       = rr(j_out) * rr(j_out1)
                cThijk          = rij_dot_rik * rrij_rrik
                h_cThijk        = h - cThijk
                gThden          = 1./(d2 + h_cThijk * h_cThijk)
                gThijk          = gthi - c2 * gThden
                zetaij          = zetaij + fC(j_out1) * gthijk
            END DO
            bZijn   = 1. + betan*(zetaij**n)
            bij     = c_X(index_ij) * (bZijn**n2r)
            ei_fin  = ei_fin + fC(j_out)*(fR(j_out) + bij*fA(j_out))
        END DO
    END DO

    print*, ei_in, ei_fin
    energy_diff = ei_in - ei_fin

END SUBROUTINE calcener_MC

SUBROUTINE calcener()
    USE library, ONLY : stdout, dtype, max_arr, min_arr, print_array, error
    USE system, ONLY : rx, ry, rz, sp, natoms, nx, ny, nz, Lx, Ly, Lz, energy_tmp
    USE potential, ONLY : c_A, c_B, c_lam, c_mu, c_X, c_R, &
                            c_h, c_n, c_R2, c_S2, c_pRSr, &
                            c_betan, c_d2, c_dr2, c_c2, c_n2r 
    USE OMP_LIB
    
    IMPLICIT NONE
    ! PARAMETERS AND LOCAL VARIABLES FOR COMPUTATIONS
    INTEGER, PARAMETER :: nmax=20
    INTEGER :: ncells
    REAL(dtype), DIMENSION(nmax) :: fC, fA, fR, dx, dy, dz, r, rr
    REAL(dtype) :: cThijk, h_cThijk, gThijk
    INTEGER, DIMENSION(nmax) :: spq, atom
    REAL(dtype) :: enep
    ! REORDERING OF ATOMS
    REAL(dtype), DIMENSION(natoms) :: e1
    REAL(dtype), DIMENSION(natoms) :: rcx, rcy, rcz
    INTEGER, DIMENSION(:), ALLOCATABLE :: np 
    INTEGER, DIMENSION(:), ALLOCATABLE :: indp 
    INTEGER, DIMENSION(natoms) :: indc, indcp, s1p
    INTEGER, DIMENSION(27) :: vcx1, vcy1, vcz1
    ! LOCAL VARIABLES
    INTEGER :: vcx, vcy, vcz, c, wcx, wcy, wcz, index_ij, spi
    INTEGER :: i, j, k, q, c1, nextj, nextk1, nextk
    REAL(dtype) :: ei, c2, d2, dr2, h, n, betan, n2r
    REAL(dtype) :: shiftx, shifty, shiftz, r2, a, gthi, zetaij
    REAL(dtype) :: rij_dot_rik, rrij_rrik, gThden, bZijn, bij
    
    ! Allocate and compute total number of cells
    ncells = nx*ny*nz
    ALLOCATE(np(ncells), indp(ncells+1))
    ! Initialize arrays to zero
    DO i = 1, ncells
        np(i) = 0.
    END DO
    DO i = 1, natoms
        indc(i) = 0.
        e1(i) = 0.
    END DO
    
    ! Assign to each atom a cell number and compute np: number of atoms in each cell
    DO i = 1, natoms
        vcx = floor(nx*(rx(i)+0.5))
        vcy = floor(ny*(ry(i)+0.5))
        vcz = floor(nz*(rz(i)+0.5))
        c = nz*(ny*vcx+vcy)+vcz + 1
        indc(i) = c
        np(c) = np(c) + 1
    END DO
    
    ! Indexes from/to iterate atoms for cell c
    indp(1) = 1
    DO c = 1, ncells
        indp(c+1) = indp(c) + np(c)
    END DO

    ! Reorder atoms with indcp
    DO i = 1, natoms 
        c = indc(i)
        rcx(indp(c)) = (rx(i)+0.5) * Lx
        rcy(indp(c)) = (ry(i)+0.5) * Ly
        rcz(indp(c)) = (rz(i)+0.5) * Lz
        s1p(indp(c)) = sp(i)
        indcp(indp(c)) = i
        indp(c) = indp(c) + 1
    END DO
    
    ! Recompute indp
    indp(1) = 1
    DO c = 1, ncells 
        indp(c+1) = indp(c) + np(c)
    END DO

    ! Create nearest neighbour indexes for iteration
    q = 1
    DO i = -1, 1
    DO j = -1, 1
    DO k = -1, 1
        vcx1(q) = i
        vcy1(q) = j
        vcz1(q) = k
        q = q + 1
    END DO
    END DO
    END DO

    ! Iterate over cells to compute energy
    enep = 0
    ! $OMP PARALLEL DO
    DO vcx = 1, nx
    DO vcy = 1, ny
    DO vcz = 1, nz
        c = (nz)*(ny*(vcx-1)+vcy-1)+vcz
        ! Find cell index and iterate over atoms in the cell c
        DO i = indp(c), indp(c+1)-1
            ei      = 0.
            spi     = sp(i)
            c2      = c_c2(spi)
            d2      = c_d2(spi)
            dr2     = c_dr2(spi)
            h       = c_h(spi)
            n       = c_n(spi)
            betan   = c_betan(spi)
            n2r     = c_n2r(spi)
            q = 1
            DO k = 1, 27
                wcx=vcx + vcx1(k)
                wcy=vcy + vcy1(k)
                wcz=vcz + vcz1(k)
                ! Periodic Boundary conditions
                shiftx = 0.
                IF (wcx == 0) THEN
                    shiftx =-Lx
                    wcx = nx
                ELSE IF (wcx==nx+1) THEN
                    shiftx = Lx
                    wcx = 1
                END IF
                shifty = 0.
                if (wcy == 0) THEN
                    shifty =-Ly
                    wcy = ny
                ELSE IF (wcy==ny+1) THEN
                    shifty = Ly
                    wcy = 1
                END IF
                shiftz = 0.
                if (wcz == 0) THEN
                    shiftz =-Lz
                    wcz = nz
                ELSE IF (wcz==nz+1) THEN
                    shiftz = Lz
                    wcz = 1
                END IF
                c1 = nz*(ny*(wcx-1)+(wcy-1))+wcz
                DO j = indp(c1), indp(c1+1)-1      ! Check which atom is near
                    dx(q)   = rcx(i)-(rcx(j) + shiftx)
                    dy(q)   = rcy(i)-(rcy(j) + shifty)
                    dz(q)   = rcz(i)-(rcz(j) + shiftz)
                    r2      = dx(q)*dx(q) + dy(q)*dy(q) + dz(q)*dz(q)
                    spq(q)  = s1p(j)
                    index_ij= spq(q)+spi-1
                    IF (r2 .le. c_s2(index_ij) .and. r2 > 0.1) THEN
                        atom(q) = j
                        r(q)    = sqrt(r2)
                        rr(q)   = 1./r(q)
                        fA(q)   = -c_B(index_ij)*exp(-c_mu(index_ij)*r(q))
                        fR(q)   = c_A(index_ij)*exp(-c_lam(index_ij)*r(q))
                        IF (r2 > c_R2(index_ij)) THEN 
                            a       = c_pRSr(index_ij)*(r(q)-c_R(index_ij))
                            fc(q)   = 0.5 + 0.5*cos(a)
                        ELSE
                            fc(q)   = 1.
                        END IF
                        q = q+1
                    END IF 
                END DO
            END DO
            q = q-1
            gthi = 1. + c2*dr2 
            DO nextj = 1, q     ! Compute energy with only nearest atoms
                index_ij = spq(nextj) + spi - 1 
                zetaij = 0.d0
                DO nextk1 = nextj+1, nextj+q-1
                    ! IF (nextk .ne. nextj) THEN
                    nextk = mod(nextk1-1, q) + 1
                    rij_dot_rik     = dx(nextj)*dx(nextk) + dy(nextj)*dy(nextk) + dz(nextj)*dz(nextk)
                    rrij_rrik       = rr(nextj) * rr(nextk)
                    cThijk          = rij_dot_rik * rrij_rrik
                    h_cThijk        = h - cThijk
                    gThden          = 1./(d2 + h_cThijk * h_cThijk)
                    gThijk          = gthi - c2 * gThden
                    zetaij          = zetaij + fC(nextk) * gthijk
                    ! END IF
                END DO 
                bZijn   = 1.d0 + betan * (zetaij ** n)
                bij     = c_X(index_ij)*(bzijn**n2r)
                
                ei = ei + fc(nextj) * (fR(nextj) + bij * fA(nextj))
            END DO
            energy_tmp = energy_tmp + ei
        END DO
    END DO
    END DO
    END DO
    ! $OMP END PARALLEL DO
    energy_tmp = energy_tmp * 0.5
END SUBROUTINE calcener
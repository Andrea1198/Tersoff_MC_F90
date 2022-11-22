SUBROUTINE calcener(energy, e1)
    USE library, ONLY : dtype, max_arr, min_arr, print_array, error
    USE system, ONLY : rx, ry, rz, sp, natoms, nx, ny, nz, Lx, Ly, Lz
    USE potential, ONLY : c_A, c_B, c_lam, c_mu, c_X, c_R, &
                            c_h, c_n, c_R2, c_S2, c_pRSr, &
                            c_betan, c_d2, c_dr2, c_c2, c_n2r 
    USE OMP_LIB
    
    IMPLICIT NONE
    REAL(dtype), INTENT(OUT) :: energy
    REAL(dtype), DIMENSION(natoms), INTENT(OUT) :: e1
    ! PARAMETERS AND LOCAL VARIABLES FOR COMPUTATIONS
    INTEGER, PARAMETER :: nmax=20
    INTEGER :: ncells
    REAL(dtype), DIMENSION(nmax) :: fC, fA, fR, dx, dy, dz, r, rr, cThijk, h_cThijk
    REAL(dtype), DIMENSION(nmax) ::  gThijk
    INTEGER, DIMENSION(nmax) :: spq, atom
    REAL(dtype) :: enep
    ! REORDERING OF ATOMS
    REAL(dtype), DIMENSION(natoms) :: rcx, rcy, rcz
    INTEGER, DIMENSION(:), ALLOCATABLE :: np 
    INTEGER, DIMENSION(:), ALLOCATABLE :: indp 
    INTEGER, DIMENSION(natoms) :: indc, indcp, s1p
    ! LOCAL VARIABLES
    INTEGER :: vcx, vcy, vcz, c, wcx, wcy, wcz, index_ij, spi
    INTEGER :: i, j, q, c1, nextj, nextk1, nextk, k1, k2, k3
    REAL(dtype) :: ei, c2, d2, dr2, h, n, betan, n2r
    REAL(dtype) :: shiftx, shifty, shiftz, r2, a, gthi, zetaij
    REAL(dtype) :: rij_dot_rik, rrij_rrik, gThden, bZijn, bij
    ! FILENAMES FOR DEBUG !
    CHARACTER(LEN=15) :: dir="./files/output/"
    CHARACTER(LEN=101) :: np_out, indc_out, vcpos_out, indp_out
    
    np_out=dir//"np.out"
    indc_out=dir//"indc.out"
    vcpos_out=dir//"vcpos.out"
    indp_out=dir//"indp.out"
    ncells = nx*ny*nz
    ALLOCATE(np(ncells), indp(ncells+1))
    
    DO i = 1, ncells
        np(i) = 0
    END DO
    DO i = 1, natoms
        indc(i) = 0
        e1(i) = 0
    END DO
    
    OPEN(20,FILE=vcpos_out)
    DO i = 1, natoms
        vcx = int(nx*(rx(i)+0.5))
        vcy = int(ny*(ry(i)+0.5))
        vcz = int(nz*(rz(i)+0.5))
        c = nz*(ny*vcx+vcy)+vcz + 1
        WRITE(20, "(4(i3))") vcx, vcy, vcz, c
        indc(i) = c
        np(c) = np(c) + 1
    END DO
    CLOSE(20)
    
    OPEN(20,FILE=np_out)
    indp(1) = 1
    DO c = 1, ncells
        WRITE(20,"(2i3)") c, np(c)
        indp(c+1) = indp(c) + np(c)
    END DO
    CLOSE(20)    
    DO i = 1, natoms 
        c = indc(i)
        rcx(indp(c)) = (rx(i)+0.5) * Lx
        rcy(indp(c)) = (ry(i)+0.5) * Ly
        rcz(indp(c)) = (rz(i)+0.5) * Lz
        s1p(indp(c)) = sp(i)
        indcp(indp(c)) = i
        indp(c) = indp(c) + 1
    END DO
    
    indp(1) = 1
    DO c = 1, ncells 
        indp(c+1) = indp(c) + np(c)
    END DO

    OPEN(20,FILE=indp_out)
    DO c = 1, ncells+1
        WRITE(20, "(i3)") indp(c)
    END DO
    CLOSE(20)
    ! call print_array(indp)
    ! call print_array(np)
    ! call print_array(indcp)
    
    enep = 0
    DO vcx = 1, nx
        DO vcy = 1, ny
            DO vcz = 1, nz
                c = (nz)*(ny*(vcx-1)+vcy-1)+vcz
                DO i = indp(c), indp(c+1)-1     ! Cycle over atoms of cell c (ordered, from indp(c) to indp(c+1)-1)
                    ei      = 0.
                    spi     = sp(i)
                    c2      = c_c2(spi)
                    d2      = c_d2(spi)
                    dr2     = c_dr2(spi)
                    h       = c_h(spi)
                    n       = c_n(spi)
                    betan   = c_betan(spi)
                    n2r     = c_n2r(spi)
                    q = 0
                    DO k1 = -1, 1       ! Cycle over nearby cells to find nn atoms
                        DO k2 = -1, 1
                            DO k3 = -1, 1
                                wcx=vcx + k1
                                wcy=vcy + k2
                                wcz=vcz + k3
                                shiftx = 0.
                                ! Periodic Boundary conditions
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
                                DO j = indp(c1), indp(c1+1)-1       ! Check which atom is near
                                    IF (j .gt. natoms) then
                                        call error("wrong cell index\n")
                                    end if
                                    bZijn   = 1.d0 + betan * (zetaij ** n)
                                    dx(q)   = rcx(i)-(rcx(j) + shiftx)
                                    dy(q)   = rcy(i)-(rcy(j) + shifty)
                                    dz(q)   = rcz(i)-(rcz(j) + shiftz)
                                    r2      = dx(q)*dx(q) + dy(q)*dy(q) + dz(q)*dz(q)
                                    spq(q)  = s1p(j)
                                    index_ij= spq(q)+spi-1
                                    IF (r2 .le. c_s2(index_ij) .and. r2 > 0.1) THEN
                                        q = q+1
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
                                    END IF 
                                END DO
                            END DO
                        END DO
                    END DO
                    gthi = 1. + c2*dr2 
                    DO nextj = 1, q     ! Compute energy with only nearest atoms
                        index_ij = spq(nextj) + spi - 1 
                        zetaij = 0.d0
                        DO nextk1 = nextj+1, nextj+q
                            nextk           = mod(nextk1, q) + 1
                            rij_dot_rik     = dx(nextk)*dx(nextk) + dy(nextk)*dy(nextk) + dz(nextk)*dz(nextk)
                            rrij_rrik       = rr(nextj) * rr(nextk)
                            cThijk(nextk)   = rij_dot_rik * rrij_rrik
                            h_cThijk(nextk) = h - cThijk(nextk)
                            gThden          = 1./(d2 + h_cThijk(nextk) * h_cThijk(nextk))
                            gThijk(nextk)       = gthi - c2 * gThden
                            zetaij          = zetaij + fC(nextk) * gthijk(nextk)
                        END DO 
                        bij     = c_X(index_ij)*(bzijn**n2r)
                        
                        ei = ei + fc(nextj) * (fR(nextj) + bij * fA(nextj))
                    END DO
                    e1(indcp(i)) = e1(indcp(i)) + ei
                END DO
            END DO
        END DO
    END DO
    energy = sum(e1) * 0.5
END SUBROUTINE calcener
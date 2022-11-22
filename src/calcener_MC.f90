SUBROUTINE calcener_MC(mover, energy_diff)
    USE library, ONLY : dtype, max_arr, min_arr, print_array
    USE system, ONLY : x, y, z, sp, natoms, nx, ny, nz, Lx, Ly, Lz
    USE potential, ONLY : c_A, c_B, c_lam, c_mu, c_X, c_R, &
                            c_h, c_n, c_R2, c_S2, c_pRSr, &
                            c_betan, c_d2, c_dr2, c_c2, c_n2r 
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mover
    REAL(dtype), INTENT(OUT) :: energy_diff
    ! PARAMETERS AND LOCAL VARIABLES FOR COMPUTATIONS
    INTEGER, PARAMETER :: nmax=20, nmax_out=100
    INTEGER :: ncells
    REAL(dtype), DIMENSION(nmax) :: fC, fA, fR, dx, dy, dz, r, rr, cThijk, h_cThijk, gThijk
    REAL(dtype), DIMENSION(nmax) ::  x_in, y_in, z_in, x_crw, y_crw, z_crw
    REAL(dtype), DIMENSION(nmax_out) :: x_out, y_out, z_out
    INTEGER, DIMENSION(nmax) :: spq_in, atom_in, spq_crw, atom_crw 
    INTEGER, DIMENSION(nmax_out) :: spq_out, atom_out
    REAL(dtype) :: ei_in, ei_fin
    ! REORDERING OF ATOMS
    REAL(dtype), DIMENSION(natoms) :: rcx, rcy, rcz, e1
    INTEGER, DIMENSION(:), ALLOCATABLE :: np 
    INTEGER, DIMENSION(:), ALLOCATABLE :: indp 
    INTEGER, DIMENSION(natoms) :: indc, indcp, s1p
    INTEGER, DIMENSION(27) :: vcx1, vcy1, vcz1
    ! LOCAL VARIABLES
    INTEGER :: vcx, vcy, vcz, c, wcx, wcy, wcz, index_ij, spi
    INTEGER :: i, j, k, q_in, q_out, q_crown, c1, nextj, nextk1, nextk
    REAL(dtype) :: ei, c2, d2, dr2, h, n, betan, n2r
    REAL(dtype) :: shiftx, shifty, shiftz, r2, a, gthi, zetaij
    REAL(dtype) :: rij_dot_rik, rrij_rrik, gThden, bZijn, bij

    ei_in   = 0.
    ei_fin  = 0.

    spi = sp[mover]
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

    q_in = 1
    q_out = 1
    q_crown = 0
    DO i = 1, natoms
        dx(1) = x(i) - x_in(1)
        dy(1) = y(i) - y_in(1)
        dz(1) = z(i) - z_in(1)
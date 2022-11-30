PROGRAM run_test

    USE library,        ONLY : read_steps, nsteps, steps, kt, dtype, print_infos
    USE potential,      ONLY : read_coeffs
    USE system,         ONLY : read_vec, create_crystal, write_crystal, redefine_cells, print_sysinfo, &
                                mx, my, mz, energy

    IMPLICIT NONE
    INTEGER :: i
    CHARACTER(LEN=7), DIMENSION(4) :: filenames
    mx = 5                      ! Number of cells in x direction
    my = 3                      ! Number of cells in y direction
    mz = 5                      ! Number of cells in z direction
    ! ### Initialization of parameters ### !
    CALL read_steps()           ! Read number of steps from ./files/input/steps.txt
    CALL read_coeffs()          ! Read coefficients of Tersoff potential
    CALL read_vec()             ! Read primitive vectors for alpha-quarz system (SiO2)

    CALL create_crystal()     ! Create atom positions
    ! CALL write_crystal()      ! Write position on file ./files/output/crystal.csv, crystal.xyz
    CALL redefine_cells()       ! Optimize cell dimensions for faster nearest neighbour search
    CALL print_sysinfo()        ! Write in stdout the system informations
    DO i = 1, nsteps
        CALL eqmc(i, 0)                 ! Start simulation      (### DEBUG ###)
        CALL eqmc(i, 1)                 ! Start simulation      (### DEBUG ###)
        CALL print_infos(i, steps(i), kt(i), energy)
        CALL write_crystal(i)
    END DO
END PROGRAM run_test
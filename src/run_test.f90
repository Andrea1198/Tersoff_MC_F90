PROGRAM run_test

    USE library,        ONLY : read_steps
    USE potential,      ONLY : read_coeffs
    USE system,         ONLY : read_vec, create_crystal, write_crystal, redefine_cells, print_sysinfo 

    IMPLICIT NONE
    INTEGER :: mx
    INTEGER :: my
    INTEGER :: mz
    ! ### Initialization of parameters ### !
    CALL read_steps()           ! Read number of steps from ./files/input/steps.txt
    CALL read_coeffs()          ! Read coefficients of Tersoff potential
    CALL read_vec()             ! Read primitive vectors for alpha-quarz system (SiO2)

    CALL create_crystal()       ! Create atom positions
    CALL write_crystal()        ! Write position on file ./files/output/crystal.csv, crystal.xyz
    CALL redefine_cells()       ! Optimize cell dimensions for faster nearest neighbour search
    CALL print_sysinfo()        ! Write in stdout the system informations

    CALL eqmc()                 ! Start simulation      (### DEBUG ###)

END PROGRAM run_test
PROGRAM run_test
    USE library
    USE potential
    USE system

    IMPLICIT NONE
    REAL(dtype) :: energy
    ! ### Initialization of parameters ### !
    CALL read_steps()
    CALL read_coeffs()
    CALL read_vec()
    CALL create_crystal()
    CALL write_crystal()
    CALL redefine_cells()
    CALL print_sysinfo()

    CALL eqmc()

END PROGRAM run_test
SUBROUTINE eqmc()
    USE library, ONLY : dtype, nsteps, steps, gdr_freq, gdr_steps
    USE system, ONLY : natoms
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(dtype) :: energy = 0.
    REAL(dtype), DIMENSION(natoms) :: en

    DO i = 1, natoms
        en(i) = 0.d0
    END DO

    
    DO i = 1, nsteps
        DO j = 1, steps(i)
            CALL calcener(energy, en)
        END DO
    END DO
    WRITE(*, *) "Finished with energy : ", energy

END SUBROUTINE            
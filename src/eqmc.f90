SUBROUTINE eqmc()
    USE library, ONLY : dtype, nsteps, steps, gdr_freq, gdr_steps
    USE system, ONLY : natoms
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(dtype) :: energy = 0.
    REAL(dtype) :: en


    CALL calcener_MC(1, 0.01,0.01,0.01,en)    
    WRITE(*, *) "Finished with energy : ", energy

END SUBROUTINE            
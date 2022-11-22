SUBROUTINE eqmc()
    USE library, ONLY : dtype, nsteps, steps, gdr_freq, gdr_steps
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(dtype) :: energy = 0.

    
    DO i = 1, nsteps
        DO j = 1, steps(i)
            IF (mod(j, 100) .eq. 0) THEN 
                print*, j
                print*, "\n"
                print*, "\n"
                print*, "\n"
            END IF
            CALL calcener(energy)
            
        END DO
    END DO

END SUBROUTINE            
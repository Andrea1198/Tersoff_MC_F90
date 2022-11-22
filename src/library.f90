MODULE library
    IMPLICIT NONE
    ! TYPE  CONSTANTS 
    INTEGER, PARAMETER :: dtype=8, itype=8
    INTEGER, PARAMETER :: stdout=6
    ! GENERAL CONSTANTS
    REAL(dtype) :: PI = acos(-1.)
    ! SYSTEM PARAMETERS 
    INTEGER, PARAMETER :: freq=0
    REAL(dtype), PARAMETER :: dt=5.d-2
    REAL(dtype), PARAMETER :: delta=4.d-4, dV0=1.d-2, p_ext=0.97d0
    ! SIMULATION STEPS
    INTEGER, PARAMETER :: nsteps=1
    INTEGER, DIMENSION(nsteps) :: steps, gdr_steps, gdr_freq
    REAL(dtype), dimension(nsteps) :: kt
    
    interface print_matrix
        module procedure print_matrix_real
        module procedure print_matrix_int
        module procedure print_matrix_logic
    end interface
     
    interface print_array
        module procedure print_array_real
        module procedure print_array_int
        module procedure print_array_logic
    end interface 


    CONTAINS
        SUBROUTINE read_steps()
            IMPLICIT NONE 
            CHARACTER(LEN=23)::steps_filename='./files/input/steps.txt'
            INTEGER :: i

            OPEN(20, FILE=steps_filename)
            DO i = 1, nsteps
                READ(20, *) steps(i), kt(i), gdr_steps(i), gdr_freq(i)
            END DO
            CLOSE(20)

        END SUBROUTINE read_steps

        SUBROUTINE print_matrix_real(mat, n1, n2)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n1, n2 
            REAL(dtype), DIMENSION(n1, n2), INTENT(IN) :: mat
            INTEGER :: i, j

            DO i = 1, n1 
                DO j = 1, n2 
                    WRITE(stdout,*) mat(i,j)
                END DO
            END DO
        END SUBROUTINE print_matrix_real

        SUBROUTINE print_matrix_int(mat, n1, n2)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n1, n2 
            INTEGER, DIMENSION(n1, n2), INTENT(IN) :: mat
            INTEGER :: i, j

            DO i = 1, n1 
                DO j = 1, n2 
                    WRITE(stdout,*) mat(i,j)
                END DO
            END DO
        END SUBROUTINE print_matrix_int
        
        SUBROUTINE print_matrix_logic(mat, n1, n2)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n1, n2 
            LOGICAL, DIMENSION(n1, n2), INTENT(IN) :: mat
            INTEGER :: i, j

            DO i = 1, n1 
                DO j = 1, n2 
                    WRITE(stdout,*) mat(i,j)
                END DO
            END DO
        END SUBROUTINE print_matrix_logic
        
        SUBROUTINE print_array_real(arr)
            IMPLICIT NONE
            REAL(dtype), DIMENSION(:), INTENT(IN) :: arr
            INTEGER :: n1 
            INTEGER :: i
            n1 = size(arr)
            DO i = 1, n1 
                WRITE(stdout,*) arr(i)
            END DO
        END SUBROUTINE print_array_real
        
        SUBROUTINE print_array_int(arr)
            IMPLICIT NONE
            INTEGER, DIMENSION(:), INTENT(IN) :: arr
            INTEGER :: n1 
            INTEGER :: i
            n1 = size(arr)
            DO i = 1, n1 
                WRITE(stdout,*) arr(i)
            END DO
        END SUBROUTINE print_array_int
        
        SUBROUTINE print_array_logic(arr)
            IMPLICIT NONE
            LOGICAL, DIMENSION(:), INTENT(IN) :: arr
            INTEGER :: n1 
            INTEGER :: i
            n1 = size(arr)
            
            DO i = 1, n1 
                WRITE(stdout,*) arr(i)
            END DO
        END SUBROUTINE print_array_logic

        FUNCTION max_arr(arr, n)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(dtype), DIMENSION(n), INTENT(in) :: arr
            REAL(dtype) :: max_arr
            INTEGER :: i

            max_arr = arr(1)
            DO i = 2, n
                if (arr(i) .gt. max_arr) max_arr = arr(i)
            END DO
        END FUNCTION max_arr

        FUNCTION min_arr(arr, n)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(dtype), DIMENSION(n), INTENT(in) :: arr
            REAL(dtype) :: min_arr
            INTEGER :: i

            min_arr = arr(1)
            DO i = 2, n
                if (arr(i) .lt. min_arr) min_arr = arr(i)
            END DO
        END FUNCTION min_arr

END MODULE library

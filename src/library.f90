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
    INTEGER, PARAMETER :: nsteps=5
    INTEGER, DIMENSION(nsteps) :: steps, gdr_steps, gdr_freq
    REAL(dtype), dimension(nsteps) :: kt
    
    ! Print a matrix in stdout
    interface print_matrix
        module procedure print_matrix_real
        module procedure print_matrix_int
        module procedure print_matrix_logic
    end interface
                
    ! Print an array in stdout
    interface print_array
        module procedure print_array_real
        module procedure print_array_int
        module procedure print_array_logic
    end interface 


    CONTAINS
    !===============================================
        SUBROUTINE read_steps()
        ! Read initial steps from file
            IMPLICIT NONE 
            CHARACTER(LEN=23)::steps_filename='./files/input/steps.txt'
            INTEGER :: i

            OPEN(20, FILE=steps_filename)
            DO i = 1, nsteps
                READ(20, *) steps(i), kt(i), gdr_steps(i), gdr_freq(i)
            END DO
            CLOSE(20)

        END SUBROUTINE read_steps

    !===============================================
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

    !===============================================
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
        
    !===============================================
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
        
    !===============================================
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
        
    !===============================================
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
        
    !===============================================
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

    !===============================================
        FUNCTION max_arr(arr)
        ! Find max value of array
            IMPLICIT NONE
            INTEGER :: n
            REAL(dtype), DIMENSION(:), INTENT(in) :: arr
            REAL(dtype) :: max_arr
            INTEGER :: i
            n = size(arr)
            max_arr = arr(1)
            DO i = 2, n
                if (arr(i) .gt. max_arr) max_arr = arr(i)
            END DO
        END FUNCTION max_arr

    !===============================================
        FUNCTION min_arr(arr)
        ! Find min value of array
            IMPLICIT NONE
            INTEGER :: n
            REAL(dtype), DIMENSION(:), INTENT(in) :: arr
            REAL(dtype) :: min_arr
            INTEGER :: i

            n = size(arr)
            min_arr = arr(1)
            DO i = 2, n
                if (arr(i) .lt. min_arr) min_arr = arr(i)
            END DO
        END FUNCTION min_arr

    !===============================================
        SUBROUTINE error(message)
        ! End execution with message
            IMPLICIT NONE
            CHARACTER(LEN=*), INTENT(IN) :: message
            WRITE(stdout, *) "Error, execution interruption due to : ", message
            stop
        END SUBROUTINE error

    !===============================================
        SUBROUTINE print_infos(i, step_i, kt_i, enep_i)
            IMPLICIT NONE
            INTEGER, INTENT(IN) ::  i, step_i
            REAL(dtype), INTENT(IN) :: kt_i, enep_i

            WRITE(stdout, "(A, i2, A, i10, A, f4.2, A, f15.6)") "### Starting phase ", i," steps = ", &
                                                        step_i, "  temperature = ", kt_i, " energy = ", enep_i
        END SUBROUTINE print_infos
END MODULE library

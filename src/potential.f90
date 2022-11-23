MODULE potential
    USE library,     ONLY : dtype  
    IMPLICIT NONE
    INTEGER, PARAMETER :: natoms=2, combine=3
    ! COEFFS READ FROM FILE
    REAL(dtype), DIMENSION(combine) :: c_A, c_B, c_lam, c_mu, c_X, c_S, c_R
    REAL(dtype), DIMENSION(natoms)  :: c_beta, c_c, c_d, c_h, c_n
    ! COEFFS STORED IN MEMORY FOR COMPUTATIONS
    REAL(dtype), DIMENSION(combine) :: c_R2, c_S2, c_pRSr
    REAL(dtype), DIMENSION(natoms)  :: c_betan, c_d2, c_dr2, c_c2, c_n2r 

    CONTAINS
    !===============================================
    SUBROUTINE read_coeffs()
    ! Read potential coefficients from file
        USE library, ONLY : print_array, pi
        IMPLICIT NONE
        CHARACTER(LEN=35) :: coeff_filename="./files/input/coeff_SiO_tersoff.txt"
        INTEGER :: i, ind
        OPEN(20, FILE=coeff_filename)
        DO i = 1, 2
            ind = i + mod((i-1),2)          ! if i=1 -> ind=1, if i=2 -> ind=3
            READ(20, *) c_A(ind), c_B(ind), c_lam(ind), c_mu(ind), c_beta(i), c_n(i), &
                    c_c(i), c_d(i), c_h(i), c_R(ind), c_S(ind), c_x(2)
        END DO
        CLOSE(20)


        c_A(2)    = sqrt(c_A(1)*c_A(3))
        c_B(2)    = sqrt(c_B(1)*c_B(3))
        c_R(2)    = sqrt(c_R(1)*c_R(3))
        c_S(2)    = sqrt(c_S(1)*c_S(3))
        c_mu(2)   = (c_mu(1)+c_mu(3)) * 0.5
        c_lam(2)  = (c_lam(1)+c_lam(3)) * 0.5
        c_x(1) = 1.d0
        c_x(3) = 1.d0

        c_d2    = c_d * c_d 
        c_dr2   = 1./c_d2 
        c_c2    = c_c * c_c 
        c_R2    = c_R * c_R
        c_S2    = c_S * c_S
        c_pRSr  = PI / (c_S - c_R)
        c_betan = c_beta ** c_n
        c_n2r   = -1./(2.*c_n)

    END SUBROUTINE read_coeffs


END MODULE potential
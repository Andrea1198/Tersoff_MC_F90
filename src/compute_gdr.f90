SUBROUTINE compute_gdr()
    USE library,        ONLY : dtype
    USE system,         ONLY : gcount, rx, ry, rz, sp, natoms, Lx, Ly, Lz, ldel
    IMPLICIT NONE 
    REAL(dtype) ::  r2, r
    REAL(dtype) :: dx, dy, dz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irdf
    INTEGER :: i, j

    DO i = 1, natoms-1
        ALLOCATE(irdf(natoms-i-1))
        DO j = i+1, natoms
            dx = rx(i) - rx(j)
            dy = ry(i) - ry(j)
            dz = rz(i) - rz(j)
            dx = dx - anint(dx)
            dy = dy - anint(dy)
            dz = dz - anint(dz)
            dx = dx * Lx
            dy = dy * Ly
            dz = dz * Lz
            irdf(j-1) = sp(j) + sp(i)
            r2 = dx*dx + dy*dy + dz*dz
            
            IF (r2 .lt. min(Lx, Ly, Lz) * 0.5) THEN
                r = sqrt(r2)
                gcount(int(r / ldel), irdf(j-1)) = gcount(int(r / ldel), irdf(j-1)) + 2
            END IF
        END DO

        DEALLOCATE(irdf)

    END DO
END SUBROUTINE compute_gdr
    
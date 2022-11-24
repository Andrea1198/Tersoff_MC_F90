SUBROUTINE eqmc()
    USE library, ONLY : dtype, nsteps, steps, gdr_freq, gdr_steps, delta, kt, dV0, p_ext
    USE system, ONLY : natoms, x, y, z, rx, ry, rz, Lx, Ly, Lz, volume
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(dtype) :: energy = 0., enep_in, enep_fin
    REAL(dtype) :: en_diff, disp_x, disp_y, disp_z, dV
    REAL(dtype) :: random, accep, xt, yt, zt, rxt, ryt, rzt
    REAL(dtype), DIMENSION(natoms) :: rvx, rvy, rvz
    INTEGER :: mover
    ! Initialize accumulators
    INTEGER :: moved = 0, movev = 0, accepd = 0, accepv = 0

    CALL calcener(enep_in)

    DO i = 1, nsteps
        DO j = 1, steps(i)
            CALL random_number(random)
            mover = int(random * (natoms + 1)) + 1 
            IF (mover .le. natoms) THEN
                moved = moved + 1

                CALL random_number(disp_x)
                CALL random_number(disp_y)
                CALL random_number(disp_z)

                disp_x = disp_x * 2 - 1
                disp_y = disp_y * 2 - 1
                disp_z = disp_z * 2 - 1

                disp_x = disp_x * delta
                disp_y = disp_y * delta
                disp_z = disp_z * delta

                xt = x(mover) + disp_x
                yt = y(mover) + disp_y
                zt = z(mover) + disp_z

                CALL calcener_MC(mover, xt, yt, zt, en_diff)
                accep = exp(en_diff/kt(i))
                CALL random_number(random)
                IF (random .le. accep) THEN
                    rxt = xt / Lx
                    ryt = yt / Ly
                    rzt = zt / Lz

                    IF (rxt .gt. 0.5) THEN
                        rxt = rxt - 1
                    ELSE IF (rxt .lt. -0.5) THEN
                        rxt = rxt + 1
                    END IF
                    IF (ryt .gt. 0.5) THEN
                        ryt = ryt - 1
                    ELSE IF (ryt .lt. -0.5) THEN
                        ryt = ryt + 1
                    END IF
                    IF (rzt .gt. 0.5) THEN
                        rzt = rzt - 1
                    ELSE IF (rzt .lt. -0.5) THEN
                        rzt = rzt + 1
                    END IF

                    xt = rxt * Lx
                    yt = ryt * Ly
                    zt = rzt * Lz

                    x(mover) = xt
                    y(mover) = yt
                    z(mover) = zt
                    accepd = accepd + 1
                    enep_in = enep_in - en_diff
                END IF
            ELSE
                movev = movev + 1
                CALL random_number(random)
                dV = random * dV0 * 2 - dV0
                dV = (1. + dV) ** (1./3.)
                rvx(:) = rx(:)
                rvy(:) = ry(:)
                rvz(:) = rz(:)
                Lx = Lx * dV
                Ly = Ly * dV
                Lz = Lz * dV

                CALL calcener(enep_fin)
                en_diff = enep_in - enep_fin
                accep = dV ** (3. * natoms) * exp((en_diff - p_ext * dV * volume)/kt(i))
                IF (random .gt. accep) THEN 
            END IF
        END DO
    END DO
END SUBROUTINE

SUBROUTINE eqmc(index, mode)
    USE library, ONLY : dtype, steps, gdr_steps, delta, kt, dV0, p_ext
    USE system, ONLY : natoms, x, y, z, rx, ry, rz, Lx, Ly, Lz, volume, energy, energy_tmp
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: index, mode
    INTEGER :: i, j, steps_tmp
    REAL(dtype) :: en_diff, disp_x, disp_y, disp_z, dV
    REAL(dtype) :: random, accep, xt, yt, zt, rxt, ryt, rzt
    REAL(dtype), DIMENSION(natoms) :: rvx, rvy, rvz
    REAL(dtype) :: time_in, time_fin
    INTEGER :: mover
    ! Initialize accumulators
    INTEGER :: moved = 0, movev = 0, accepd = 0, accepv = 0
    IF (mode == 0) steps_tmp = steps(index)
    IF (mode == 1) steps_tmp = gdr_steps(index)
    CALL calcener()
    energy = energy_tmp
    CALL CPU_TIME(time_in)
    DO j = 1, steps_tmp
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

            CALL calcener_MC(mover, xt, yt, zt)
            accep = exp(energy_tmp/kt(i))
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
                energy = energy - energy_tmp
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

            CALL calcener()
            en_diff = energy - energy_tmp
            accep = dV ** (3. * natoms) * exp((en_diff - p_ext * dV * volume)/kt(i))
            IF (random .gt. accep) THEN         ! If not accepted reset dimensions
                Lx = Lx / dV
                Ly = Ly / dV
                Lz = Lz / dV
                volume = Lx * Ly * Lz
                energy = energy_tmp
                accepv = accepv + 1
            END IF
        END IF
        IF (mod(j, 100000) .eq. 0) THEN
            CALL CPU_TIME(time_fin)
            PRINT*, j, " time  for 100000 steps : ", time_fin - time_in
            time_in = time_fin
        END IF
    END DO
END SUBROUTINE

MODULE system
    USE library,        ONLY : dtype
    IMPLICIT NONE 
    INTEGER :: nvec, natoms, Na, Nb, N
    INTEGER, PARAMETER :: mx=1
    INTEGER, PARAMETER :: my=1
    INTEGER, PARAMETER :: mz=1
    INTEGER :: nx, ny, nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: spv, sp
    REAL(dtype), DIMENSION(:), ALLOCATABLE :: xvector, yvector, zvector
    REAL(dtype), DIMENSION(:), ALLOCATABLE :: x, y, z, rx, ry, rz, mass
    REAL(dtype) :: Lx, Ly, Lz, ax, ay, az
    REAL(dtype) :: m1 = 16., m2 = 28.

    CONTAINS 
    SUBROUTINE read_vec()
        USE library, ONLY : dtype
        IMPLICIT NONE
        INTEGER :: i
        CHARACTER(LEN=29) :: vec_filename="./files/input/alphaquartz.csv"
        OPEN(20, FILE=vec_filename)
        READ(20, *) nvec
        READ(20, *) ax, ay, az
        ALLOCATE(spv(nvec), xvector(nvec), yvector(nvec), zvector(nvec))
        DO i = 1, nvec
            READ(20, *) spv(i), xvector(i), yvector(i), zvector(i)
        END DO
        CLOSE(20)
    END SUBROUTINE read_vec

    SUBROUTINE create_crystal()
        USE library, ONLY : dtype, print_array
        IMPLICIT NONE
        INTEGER :: m, i, j, k, atom, vec
        REAL(dtype) :: xi, yi, zi, delta, random

        ! shift
        xi = 0.
        yi = 0.
        zi = 0.
        delta = 0.00

        ! real sys dimensions
        m = mx * my * mz
        natoms = nvec * m
        Lx = ax * mx
        Ly = ay * my
        Lz = az * mz

        ! Allocate arrays
        ALLOCATE(x(natoms), y(natoms), z(natoms), sp(natoms), mass(natoms))
        ALLOCATE(rx(natoms), ry(natoms), rz(natoms))
        Na = 0
        Nb = 0
        ! assign atom positions
        atom = 1
        DO i = 0, mx-1
        DO j = 0, my-1
        DO k = 0, mz-1
            DO vec = 1, nvec
                CALL random_number(random)
                x(atom) = xi + ax*i + random*delta + xvector(vec)
                CALL random_number(random)
                y(atom) = yi + ay*j + random*delta + yvector(vec)
                CALL random_number(random)
                z(atom) = zi + az*k + random*delta + zvector(vec)
                sp(atom) = spv(vec)
                IF (sp(atom) .EQ. 1) THEN
                    mass(atom) = m1
                    Na = Na + 1
                ELSE
                    mass(atom) = m2 
                    Nb = Nb + 1
                END IF
                atom = atom + 1
            END DO
        END DO
        END DO
        END DO

        rx = x/Lx - 0.5
        ry = y/Ly - 0.5
        rz = z/Lz - 0.5
    END SUBROUTINE create_crystal

    SUBROUTINE write_crystal()
        USE library, ONLY : dtype, print_matrix
        IMPLICIT NONE
        INTEGER :: i
        CHARACTER(LEN=26) :: crystal_filename_xyz="./files/output/crystal.xyz"
        CHARACTER(LEN=26) :: crystal_filename_csv="./files/output/crystal.csv"
        OPEN(20, FILE=crystal_filename_xyz)
        WRITE(20, "(i4)") natoms
        WRITE(20, *)
        DO i = 1, natoms
            WRITE(20, "(i1, 2x,3(f15.9,2x))") sp(i), x(i), y(i), z(i)
        END DO
        CLOSE(20)
        
        OPEN(20, FILE=crystal_filename_csv)
        WRITE(20, *) "#sp,x,y,z"
        DO i = 1, natoms
            WRITE(20, "(i1,A,3(f15.9,A))") sp(i), ',', x(i), ',', y(i), ',', z(i)
        END DO
        CLOSE(20)
    END SUBROUTINE write_crystal

    SUBROUTINE redefine_cells()
        USE library, ONLY : dtype, max_arr
        USE potential, ONLY : c_S
        IMPLICIT NONE
        REAL(dtype) :: max_s 

        max_s = max_arr(c_s, 3)
        nx = ceiling(Lx / max_s)
        ny = ceiling(Ly / max_s)
        nz = ceiling(Lz / max_s)
    END SUBROUTINE redefine_cells

    SUBROUTINE print_sysinfo()
        USE library, ONLY : dtype, stdout
        IMPLICIT NONE

        WRITE(stdout, *) "### System initialized ###"
        WRITE(stdout, *) "# Lengths of system                     : Lx = ", Lx, " Ly = ", Ly, " Lz = ", Lz
        WRITE(stdout, *) "# Number of cells before optimization   : mx = ", mx, " my = ", my, " mz = ", mz
        WRITE(stdout, *) "# Number of cells after optimization    : nx = ", nx, " ny = ", ny, " nz = ", nz
        WRITE(stdout, *) "# Total number of cells                 : ", nx*ny*nz
        WRITE(stdout, *) "# Number of atoms for type              : Na = ", Na, " Nb = ", Nb
    END SUBROUTINE print_sysinfo
END MODULE system
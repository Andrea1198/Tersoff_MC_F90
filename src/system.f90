MODULE system
    USE library,        ONLY : dtype
    IMPLICIT NONE 
    INTEGER :: nvec, natoms, Na, Nb, N
    INTEGER :: mx
    INTEGER :: my
    INTEGER :: mz
    INTEGER :: nx, ny, nz
    INTEGER, PARAMETER :: kg = 1024
    INTEGER, DIMENSION(:), ALLOCATABLE :: spv, sp
    REAL(dtype), DIMENSION(kg, 3) :: gcount
    REAL(dtype) :: ldel
    REAL(dtype), DIMENSION(:), ALLOCATABLE :: xvector, yvector, zvector
    REAL(dtype), DIMENSION(:), ALLOCATABLE :: x, y, z, rx, ry, rz, mass
    REAL(dtype) :: Lx, Ly, Lz, ax, ay, az, volume
    REAL(dtype) :: m1 = 16., m2 = 28., energy, energy_tmp

    CONTAINS 
    !===============================================
    SUBROUTINE read_vec()
    ! Read primitive vectors from file
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

    !===============================================
    SUBROUTINE create_crystal()
    ! Generate crystal positions
        USE library, ONLY : dtype, print_array, min_arr
        IMPLICIT NONE
        INTEGER :: m, i, j, k, atom, vec
        REAL(dtype) :: xi, yi, zi, random, disp

        ! shift
        xi = 0.
        yi = 0.
        zi = 0.
        disp = 0.00

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
                random = random * 2. - 1.
                x(atom) = xi + ax*i + random*disp + xvector(vec)

                CALL random_number(random)
                random = random * 2. - 1.
                y(atom) = yi + ay*j + random*disp + yvector(vec)
                
                CALL random_number(random)
                random = random * 2. - 1.
                z(atom) = zi + az*k + random*disp + zvector(vec)
                
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

        rx = x/Lx 
        rx = rx - anint(rx)
        ry = y/Ly
        ry = ry - anint(ry)
        rz = z/Lz
        rz = rz - anint(rz)

        x = rx * Lx
        y = ry * Ly
        z = rz * Lz
        volume = Lx * Ly * Lz
        ldel   = min(Lx, Ly, Lz) * 0.5 / kg
    END SUBROUTINE create_crystal

    !===============================================
    SUBROUTINE write_crystal(index)
    ! Write output of crystal in files
        USE library, ONLY : dtype, print_matrix
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: index
        CHARACTER(LEN=20) :: filename_xyz
        CHARACTER(LEN=20) :: filename_csv
        CHARACTER(LEN=1) :: index_str
        CHARACTER(LEN=2), DIMENSION(2) :: spec 
        spec(2) = "Si"
        spec(1) = " O"
        write(index_str, "(I1)") index
        filename_xyz = "./files/output/" // index_str // ".xyz"
        filename_csv = "./files/output/" // index_str // ".csv"
        
        OPEN(20, FILE=filename_xyz)
        WRITE(20, "(i4)") natoms
        WRITE(20, *)
        DO i = 1, natoms
            WRITE(20, "(A2, 2x,3(f15.9,2x))") spec(sp(i)), x(i), y(i), z(i)
        END DO
        CLOSE(20)
        
        OPEN(20, FILE=filename_csv)
        WRITE(20, *) "#sp,x,y,z"
        DO i = 1, natoms
            WRITE(20, "(i1,A,3(f15.9,A))") sp(i), ',', x(i), ',', y(i), ',', z(i)
        END DO
        CLOSE(20)
    END SUBROUTINE write_crystal

    !===============================================
    SUBROUTINE redefine_cells()
    ! Optimize cells
        USE library, ONLY : dtype, max_arr
        USE potential, ONLY : c_S
        IMPLICIT NONE
        REAL(dtype) :: max_s 

        max_s = max_arr(c_s)
        nx = floor(Lx / max_s)
        ny = floor(Ly / max_s)
        nz = floor(Lz / max_s)
    END SUBROUTINE redefine_cells

    !===============================================
    SUBROUTINE print_sysinfo()
    ! Write system infos in stdout
        USE library,    ONLY : dtype, stdout, max_arr
        USE potential,  ONLY : c_S

        IMPLICIT NONE
        

        WRITE(stdout, "(A)") "### System initialized ###"
        WRITE(stdout, "(3(A, F7.2))") "# Lengths of system                       : Lx = ", Lx, ", Ly = ", Ly, ", Lz = ", Lz
        WRITE(stdout, "(3(A, i2))")   "# Number of cells                         : mx = ", mx, ", my = ", my, ", mz = ", mz
        WRITE(stdout, "((A, f4.2))")  "# Cut length for cell optimization        : rcut = ", max_arr(c_S)
        WRITE(stdout, *) 
        WRITE(stdout, "(A)")          "# Optimizing ..."
        WRITE(stdout, *) 
        WRITE(stdout, "(3(A, i2))")   "# Number of cells after optimization      : nx = ", nx, ", ny = ", ny, ", nz = ", nz
        WRITE(stdout, "(3(A, f4.2))") "# Cell lengths after optimization         : lx = ", Lx/nx, ", ly = ", Ly/ny, ", lz = ", Lz/nz
        WRITE(stdout, "(2(A, i3))")   "# Number of atoms for type                : Na = ", Na, " Nb = ", Nb
        WRITE(stdout, "(A, i3)")      "# Total number of cells                   : ", nx * ny * nz
    END SUBROUTINE print_sysinfo
    
END MODULE system
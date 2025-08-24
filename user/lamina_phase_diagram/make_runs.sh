#!/bin/bash
# This script performs N_ITERATIONS for a fixed force with progressive damage.
# Updated so that the starting configuration (coords) for iteration >1 is taken
# from the previous iteration's "lowest" file, and init_rods is run only in iteration 1.

#----------------------------------------------------------------------------------------------------------------
# 1) User input variables and arrays
#----------------------------------------------------------------------------------------------------------------
N_ROD=10      # Number of filaments
N_SEG=200      # Number of segments per filament
L=5            # Length of filaments
H=1            # Height of whole bundle (only really care about ratio L/H)
K3=0.0001      # Spring stiffness
LOWER_RAT=1.0  # Ratio for type 2 bending stiffness

# Fixed forces
F_I=0.1      # Applied axial force
F_S=0.0001      # Applied lateral force

BE=1.0         # Outer bending stiffness

# Number of iterations for progressive damage
N_ITERATIONS=100

# Paths to the Python scripts
PYTHON_SCRIPT="/mnt/c/Users/callu/Documents/NewGmin/user/lamina_phase_diagram/adjust_damagerandomspikes.py3"
WRITEDATA="/mnt/c/Users/callu/Documents/NewGmin/user/lamina_phase_diagram/write_data.py3"

RAND_MAG=0.05   # Magnitude of random angular displacements
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
# 2) Initialize
#----------------------------------------------------------------------------------------------------------------
# 2.1) Compile software
gfortran -o bin/init_rods utils/init_rods.f90
./../../GMIN -n

# 2.2) Create the main run folder
foldername="lamina_NROD${N_ROD}_NSEG_${N_SEG}_RAT${L}_FI_${F_I}_FS_${F_S}"
if [ ! -d "$foldername" ]; then
    mkdir "$foldername"
fi

# 2.3) Add header to data_log if not already present
if [ ! -f "$foldername/data_log" ]; then
    echo 'Iteration  F_I  F_S  BE  ENERGY' > "$foldername/data_log"
fi
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
# 3) Useful functions
#----------------------------------------------------------------------------------------------------------------
function write_data {
    # Write the data.in file to be read by data.f90
    echo $N_ROD > data.in
    echo $N_SEG >> data.in
    echo $L >> data.in
    echo $H >> data.in
    echo $K3 >> data.in
    echo $LOWER_RAT >> data.in
    echo $F_I >> data.in
    echo $F_S >> data.in
    echo $BE >> data.in
    echo $RAND_MAG >> data.in
}

function get_energy {
    # Extract the energy from the lowest file
    LINE=$(grep -i energy lowest)
    ENERGY=$(echo $LINE | awk '{print $3}')
    ENERGY=$(printf "%.10f\n" $ENERGY)
}
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
# 4) Perform the iterations
#----------------------------------------------------------------------------------------------------------------
cd "$foldername"
COUNTER=0

for ((iter=1; iter<=N_ITERATIONS; iter++))
do
    runname="Iteration_${iter}"
    # Remove the folder if it exists, then create a fresh one.
    if [ -d "$runname" ]; then
       echo "Run folder $runname already exists, removing it"
       rm -rf "$runname"
    fi
    mkdir "$runname"

    echo "Iteration $iter of $N_ITERATIONS with F_I=$F_I, F_S=$F_S"

    # 1) Update bending stiffness file (b_array.in) for this iteration.
    python3 "$PYTHON_SCRIPT" $iter $N_ROD $N_SEG
    if [ $? -ne 0 ]; then
        echo "Error: Python script adjust_damage.py3 failed"
        exit 1
    fi

    if [ ! -f b_array.in ]; then
        echo "Error: b_array.in was not created by adjust_damage.py3"
        exit 1
    fi

    # For iteration 1, run init_rods to generate the initial configuration.
    if [ "$iter" -eq 1 ]; then
        cp b_array.in "$runname/"
        cd "$runname"
        cp ../../gmin .
        write_data    # Generate a fresh data.in file
        cp ../../data.f90 .
        ./../../bin/init_rods   # This generates both the "coords" and "initpos" files
        cd ..
    else
        # For iterations > 1, copy the previous iteration's output.
        prev_run="Iteration_$((iter-1))"
        if [ -d "$prev_run" ]; then
            cd "$prev_run"
            if [ -f lowest ]; then
                # Remove header and copy as starting coordinates.
                tail -n +5 lowest > lowest_trimmed
                cp lowest_trimmed "../$runname/coords"
            else
                echo "Error: lowest file not found in $prev_run! Exiting."
                exit 1
            fi
            # Also copy the initpos file from the previous run so that potential.f90 can read it.
            if [ -f initpos ]; then
                cp initpos "../$runname/initpos"
            else
                echo "Warning: initpos not found in $prev_run. Using default initpos from source."
                cp ../../initpos "../$runname/"
            fi
            cd ..
        else
            echo "Error: Previous iteration folder $prev_run not found!"
            exit 1
        fi
        cp b_array.in "$runname/"
        cd "$runname"
        write_data    # Generate a fresh data.in file for this iteration
        cd ..
    fi

    # Now, enter the run folder, copy gmin, and run the simulation.
    cd "$runname"
    cp ../../gmin .
    ./gmin

    # Process results and update the log.
    get_energy
    echo "$iter  $F_I  $F_S  $BE  $ENERGY" >> ../data_log

    cd ..
done

#!/bin/bash
#SBATCH --job-name=<job-name>          # Job name
#SBATCH --mail-type=BEGIN,END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00                # Wall time limit (days-hrs:min:sec)
#SBATCH --output=<job-name>_%j.out     # Path to the standard output and error files relative to the working directory
#SBATCH --partition=<partition>        # Partition name
#SBATCH --qos=<qos>                    # QOS name

echo "========="
echo "Welcome!"
echo "========="
echo ""
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo ""
echo "-->Loading Environment"
echo ""

source /path/to/conda/bin/activate <tardigrade-examples-env>

echo ""
echo "Specify LD_LIBRARY_PATH"
export LD=/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp
export LD_LIBRARY_PATH=$LD
echo "LD_LIBRARY_PATH=$LD"
echo ""

# Run WAVES
echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "WAVES"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "scons -h"
echo ""
scons -h
echo ""

echo ""
echo "Run analysis: 'scons $1'"
echo ""
scons $1
echo ""

echo "========="
echo "End!!!!!!"
echo "========="

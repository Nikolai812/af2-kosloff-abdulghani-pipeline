#!/bin/bash

#SBATCH -p dlc
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 48:00:00
#SBATCH --error=%J.splitmono_dbv23.out
#SBATCH --output=%J.splitmono_dbv23.out

echo "======================================================="
echo "Starting AlphaFold 2.3.2 test"
echo "Job: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "======================================================="

BASEDIR="$(pwd)"

echo "BASEDIR=$BASEDIR"
echo ""

echo "===== GPU ====="
nvidia-smi
echo ""

echo "===== CPU / memory ====="
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK}"
free -h
echo ""

# -------------------------------------------------------
# Paths
# -------------------------------------------------------

FASTA="HsOR343.fasta"

OUTPUTS=${BASEDIR}/w23outputs/${SLURM_JOBID}

# Change this to the actual AF2.3.2 container image
IMAGE="${BASEDIR}/AF23_image.sqsh"

# -------------------------------------------------------
# Checks
# -------------------------------------------------------

# If a command-line argument was provided, use it as FASTA instead of harcoded above value.
if [ -n "$1" ]; then
    FASTA="$1"
    echo "+++++ Using FASTA file provided as argument: ${FASTA}"
else
    echo "----- No FASTA argument provided. Using default: ${FASTA}"
fi
echo ""

# Expected location of the FASTA file
FASTA_PATH="./AF_inputs/${FASTA}"

# Check that the file exists
if [ ! -f "${FASTA_PATH}" ]; then
    echo "ERROR: FASTA file not found: ${FASTA_PATH}"
    exit 1
fi


if [ ! -f "${IMAGE}" ]; then
    echo "ERROR: container image not found: ${IMAGE}"
    exit 1
fi

srun bash -c "mkdir -p ${OUTPUTS}"
srun bash -c "mkdir -p ${BASEDIR}/AF_inputs/SUBMONO"

echo "== Call the Python script and capture the output file names "
SPLITTED_FILES=$(srun python3 fasta_splitter.py "-i=${BASEDIR}/AF_inputs/${FASTA}" "-o=${BASEDIR}/AF_inputs/SUBMONO")



echo "FASTA:   ${FASTA}"
echo "OUTPUT:  ${OUTPUTS}"
echo "IMAGE:   ${IMAGE}"
echo ""

# -------------------------------------------------------
# Container mounts
# -------------------------------------------------------

#mounts="${BASEDIR}:/basedir,${BASEDIR}/AF_inputs:/AF_inputs,${OUTPUTS}:/outputs,/users/public_data/AlphaFoldMulti/alphafold_data1122:/data"
mounts="${BASEDIR}:/basedir,${BASEDIR}/AF_inputs:/AF_inputs,${OUTPUTS}:/outputs,/users/public_data/AlphaFoldMulti/alphafold_data_232:/data,/users/public_data/AlphaFoldMulti/alphafold_data1122:/old_data"


echo "Container mounts:"
echo "${mounts}"
echo ""

# -------------------------------------------------------
# AlphaFold 2.3.2
# -------------------------------------------------------

echo "Starting AlphaFold at $(date)"



# Output the list of split files
echo "== Split files: $SPLITTED_FILES"


# AlphaFold2 Loop through the split files and perform additional operations

for file in $SPLITTED_FILES; do
    echo "Start AlphaFold processing for file: $file at $(date)"

    srun -l \
        --container-image="${IMAGE}" \
        --container-mounts="${mounts}" \
        --container-workdir=/app/alphafold \
        --no-container-entrypoint \
        bash -c '
            file="$1"
            slurm_job_id="$2"

            echo "AlphaFold2 version:"
            git describe --tags --always --dirty

            echo "Starting AlphaFold2..."
            python3 /app/alphafold/run_alphafold.py \
                --fasta_paths="/AF_inputs/SUBMONO/${file}" \
                --output_dir=/outputs \
                --data_dir=/data \
                --db_preset=full_dbs \
                --model_preset=monomer_ptm \
                --uniref90_database_path=/data/uniref90/uniref90.fasta \
                --mgnify_database_path=/data/mgnify/mgy_clusters_2022_05.fa \
                --bfd_database_path=/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
		--uniref30_database_path=/data/uniref30/UniRef30_2021_03 \
                --pdb70_database_path=/data/pdb70/pdb70 \
                --template_mmcif_dir=/data/pdb_mmcif/mmcif_files \
                --obsolete_pdbs_path=/data/pdb_mmcif/obsolete.dat \
                --max_template_date=2026-09-15 \
                --use_gpu_relax=True \
                --models_to_relax=none

            af_status=$?

            if [ "$af_status" -ne 0 ]; then
                echo "ERROR: AlphaFold failed for ${file} with exit status ${af_status}"
                exit "$af_status"
            fi

            echo "AlphaFold completed successfully for ${file}"
            echo "Starting MSA/PAE post-processing..."

            python3 /basedir/msa_pae_v232.py \
                --jobnumber="${slurm_job_id}" \
                --separator=","

            post_status=$?

            if [ "$post_status" -ne 0 ]; then
                echo "ERROR: msa_pae_v232.py failed for ${file} with exit status ${post_status}"
                exit "$post_status"
            fi

            echo "MSA/PAE post-processing completed successfully for ${file}"
        ' bash "$file" "${SLURM_JOB_ID}"

    status=$?

    echo ""
    echo "======================================================="
    echo "AlphaFold232 finished processing file: $file at $(date)"
    echo "AlphaFold232 exit status: ${status}"
    echo "======================================================="
done


echo "======================================================="
echo "JOB ${SLURM_JOBID} completed at $(date)"
echo "======================================================="
echo "           "
echo "           "

OR_OUTPUTS=${BASEDIR}/OR_OUTPUTS
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Alphafold2 output processing at ${SLURM_JOBID}  the time: $(date)"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# python3 output_name_update.py -i=137783 -p=./w23outputs
# python3 output_name_update.py -i=137783 -p=./w23outputs -o=./or_outputs
python3 output_name_update.py -i=${SLURM_JOBID} -p=${BASEDIR}/w23outputs -o=OR_OUTPUTS

echo "           "
echo "           "
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Alphafold2 postprocessing completed at time: $(date)"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"



echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "PLDDT scores preparation at ${SLURM_JOBID}  the time: $(date)"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
echo "Doing plddt for ${SLURM_JOBID}_processed pdb files"
# python3 plddt_scorer.py  -i ./w23outputs/156010_processed/
# python3 plddt_scorer.py  -i ./OR_OUTPUTS/HsOR343/
python3 plddt_scorer.py --input=./w23outputs/${SLURM_JOBID}_processed/

echo "========================================================================================"
echo "Going to get/update plddt scores for OR_NAME located pdb files"

mapfile -t this_job_or_names < <(find "./w23outputs/${SLURM_JOBID}" -mindepth 1 -maxdepth 1 -type d ! -name ".*"  -printf "%f\n")

echo ".. got or_names witth mapfile..: ${or_names[@]}"
echo "========================================================================================"

echo ""
for or_name in "${this_job_or_names[@]}"; do
    echo "--------------------------------------------------------------------"
    echo "Processing $or_name with plddt_scorer.py"
    python3 plddt_scorer.py  -i $OR_OUTPUTS/$or_name
    echo "--------------------------------------------------------------------"
    echo ""
done


echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "DONE PLDDT scores for ${SLURM_JOBID} at the time: $(date), SCRIPT COMPLETED"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"










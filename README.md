# af2-kosloff-abdulghani-pipeline. Version 1.1 (Documentation)

This pipeline is expected to run AlphaFold2 predictions with a minimal processing of the AlphaFold2 output. It can be used at  **DLC** cluster at Haifa university, where **GPU** unit is available with up to *40GB* of memory available for a GPU unit.
Version 1.1 (26.05.2026)

## Installation
An available user account at the **DLC** cluster is a prerequisite.
- Login to your **DLC** account
- At your home directory download the code of the pipeline, then do to the pipeline directory:
```bash
git clone https://github.com/Nikolai812/af2-kosloff-abdulghani-pipeline.git
cd af2-kosloff-abdulghani-pipeline
```
- Copy the AlphaFold2 container to the pipeline directory (it may take some time due to the size of the contaner):
```bash
cp /users/public_data/AlphaFoldMulti/AF_image.sqsh .
```
The pipeline is installed

## Running the pipeline (input file preparation and start)
There are 3 dirrefent scripts to run the pipeline: *tred_ms_monomer.sub* shall be used for monomer AlphaFold2 option (with several possible monomers within one input .fasta file), *tred_array_ms_monomer.sub* shall be used when there are many (>10) sequenses in your .fasta file - this scripts starts array of 4 parallel jobs , and finally the  *tred_multimer.sub* shall be used for multimer AlphaFold2 option.
- Before starting the pipeline, one has to place the input fasta file 'my_protein.fasta' inside the *AF_inputs* directory
- Start running the script (monomer) with 'my_protein.fasta' as a command-line parameter:

```bash
sbatch tred_ms_monomer.sub my_protein.fasta
```

or array job script for monomer in case of more that ~10 sequences:

```bash
sbatch tred_array_ms_monomer.sub my_protein.fasta
```

or multimer:
```bash
sbatch tred_multimer.sub my_protein.fasta
```


It you do not specify "my_protein.fasta" as a command-line parameter, the default value from the scropt (tred_ms_monomer.sub, tred_array_ms_monomer.sub, tred_multimer.sub)

## Running the pipeline (monitoring and gathering the output)
After the script has got started, the sbatch jon gets its number assigned. The job out is written to the output file %jobnumber%(splitmono|mu).out. Since the job can run for many hours it is possible to watch the current output. You can also check the job status by *squeue* command.
If the jon has completed successfully, the output file will end with:

```bash
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DONE PLDDT scores for 156423 at the time: Thu 26 Feb 2026 03:55:05 PM IST, SCRIPT COMPLETED
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

```

or
```bash
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DONE PLDDT scores for 156424 at the time: Thu 26 Feb 2026 01:13:24 PM IST, MULTIMER SCRIPT COMPLETED
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```

If anything goes wrong, the completion for the AlphaFold2 sbatch job can be verified at the line:
```bash
=======================================================
JOB 156423 completed at Thu 26 Feb 2026 03:55:03 PM IST
=======================================================
```

After that line the script is doing the postprocesing of the AF2 output.
The raw AF2 output is saved to the directory *woutputs*, subdirectory *OR_NAME*. The *OR_NAME_processd* directory is placed alongside and it contains only the best unrelaxed models for ORs (only one, the best, for each) this selection is bdone by the script basing on the values of *ranking_debug.json* file for the raw AF2 output
```bash
drwxr-x--- 9 nromanenk dchelouche 4096 Feb 23 19:10 156010
drwxr-x--- 2 nromanenk dchelouche 4096 Feb 24 14:03 156010_processed
drwxr-x--- 9 nromanenk dchelouche 4096 Feb 24 23:20 156167
drwxr-x--- 2 nromanenk dchelouche 4096 Feb 25 00:03 156167_processed
drwxr-x--- 9 nromanenk dchelouche 4096 Feb 25 15:56 156273
drwxr-x--- 2 nromanenk dchelouche 4096 Feb 25 16:42 156273_processed
drwxr-x--- 3 nromanenk dchelouche 4096 Feb 25 11:49 156275
drwxr-x--- 2 nromanenk dchelouche 4096 Feb 25 14:00 156275_processed
drwxr-x--- 9 nromanenk dchelouche 4096 Feb 26 15:11 156423
drwxr-x--- 2 nromanenk dchelouche 4096 Feb 26 15:55 156423_processed
drwxr-x--- 3 nromanenk dchelouche 4096 Feb 26 11:02 156424
drwxr-x--- 2 nromanenk dchelouche 4096 Feb 26 13:13 156424_processed

```

This is the expected content fot the *OR_NAME_processd* directory contains plddt scores and looks as follows:
```bash
$ ls -l woutputs/156423_processed/
total 1600
-rw-r----- 1 nromanenk dchelouche    294 Feb 26 15:55 best_models.json
-rw-r----- 1 nromanenk dchelouche 251829 Feb 26 15:55 HsOR161_4.pdb
-rw-r----- 1 nromanenk dchelouche   5027 Feb 26 15:55 HsOR161_4_plddt.csv
-rw-r----- 1 nromanenk dchelouche 254988 Feb 26 15:55 HsOR183_2.pdb
-rw-r----- 1 nromanenk dchelouche   4874 Feb 26 15:55 HsOR183_2_plddt.csv
-rw-r----- 1 nromanenk dchelouche 262116 Feb 26 15:55 HsOR264_1.pdb
-rw-r----- 1 nromanenk dchelouche   4987 Feb 26 15:55 HsOR264_1_plddt.csv
-rw-r----- 1 nromanenk dchelouche 259767 Feb 26 15:55 HsOR32_3.pdb
-rw-r----- 1 nromanenk dchelouche   5073 Feb 26 15:55 HsOR32_3_plddt.csv
-rw-r----- 1 nromanenk dchelouche 261387 Feb 26 15:55 HsOR343_1.pdb
-rw-r----- 1 nromanenk dchelouche   5018 Feb 26 15:55 HsOR343_1_plddt.csv
-rw-r----- 1 nromanenk dchelouche 260820 Feb 26 15:55 HsOR344_4.pdb
-rw-r----- 1 nromanenk dchelouche   5021 Feb 26 15:55 HsOR344_4_plddt.csv
```

Besides this, the directory *OR_OUTPUTS* collects all unrelaxed .pdb outputs form all runs in the *OR_NAME* directories. If the same OR is being run through the pipeline several times, the corresponding .pdb files are overwritten (hence only the last successful output is being kept). The expected *OR_OUTPUTS/OR_NAME* output looks like: 

```bash
$ ls -l OR_OUTPUTS/HsOR343
total 1344
-rw-r----- 1 nromanenk dchelouche    399 Feb 26 15:09 HsOR343_ranking_debug.json
-rw-r----- 1 nromanenk dchelouche 261387 Feb 26 14:57 HsOR343_unrelaxed_model_1.pdb
-rw-r----- 1 nromanenk dchelouche   5018 Feb 26 15:55 HsOR343_unrelaxed_model_1_plddt.csv
-rw-r----- 1 nromanenk dchelouche 261387 Feb 26 15:00 HsOR343_unrelaxed_model_2.pdb
-rw-r----- 1 nromanenk dchelouche   5019 Feb 26 15:55 HsOR343_unrelaxed_model_2_plddt.csv
-rw-r----- 1 nromanenk dchelouche 261387 Feb 26 15:03 HsOR343_unrelaxed_model_3.pdb
-rw-r----- 1 nromanenk dchelouche   5032 Feb 26 15:55 HsOR343_unrelaxed_model_3_plddt.csv
-rw-r----- 1 nromanenk dchelouche 261387 Feb 26 15:05 HsOR343_unrelaxed_model_4.pdb
-rw-r----- 1 nromanenk dchelouche   5015 Feb 26 15:55 HsOR343_unrelaxed_model_4_plddt.csv
-rw-r----- 1 nromanenk dchelouche 261387 Feb 26 15:08 HsOR343_unrelaxed_model_5.pdb
-rw-r----- 1 nromanenk dchelouche   5017 Feb 26 15:55 HsOR343_unrelaxed_model_5_plddt.csv
```

The pdb outputs of this pipline can be used as inputs for **Cavity Pipeline** (*https://github.com/Nikolai812/kosloff-abdulghani-Cavity-pipeline*) that runs by Powershell under Windows.

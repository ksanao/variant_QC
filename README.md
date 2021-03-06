# variant_QC
Sequencing errors and artefacts are important confounding factors for detection of genetic variants in clinical settings (e.g for diagnostic classification or prognosic insights) based on deep next-generation sequencing. Targeted deep sequencing is usually done by amplicon or hybridization-capture protocol that are also potential sources of errors. The current case study excluded non-sequencer related errors by evaluating the perfomance of two sequencers on the same sample prepared with the same laboratory protocols.

## Visualising the results
The results are visualised in Jupyter notebooks and single html page IGV viewers generated as part of the pipeline.
* Jupyter notebooks can be visualised on Github or locally
* IGV html can only be visualised locally (downloading just the html file will also work, e.g. data/processed/igv_viewer.html)

## Reproducing the results
To reproduce full results run the steps listed below.
Alternatively, Jupyter  notebooks can be rerun on pre-computed files (if all steps were run).

### Step 1 - Build docker image
Clone the repository and build docker image.  

```bash
git clone https://github.com/ksanao/variant_QC.git  
cd variant_QC  
docker build -t variant_qc .  
cd ..  
```

### Step 2 - Start docker container and Jupyter Lab
To run docker image and start the container run the command below. The repository directory will be mounted in Docker container.  

```bash
variant_QC/start_docker.sh run  
```

Jupyter lab will be available at http://127.0.0.1:8888/lab (provde the token printed out in the container).   
To print again running labs and tokens run the following command inside container:  
`jupyter notebook list` 

### Step 3 - Run the pipeline comparing sequencers
Inside container run the following commands  
   
```bash
source .bashrc   
compare.sh src/config
```

### Step 4 - Examine the results of Step 3 in Jupyter notebook
File notebooks/01_compare_sequencers.ipynb

### Step 5 - Run variant calling pipeline
Inside container run the following commands  

```bash
source .bashrc  
call_variants.sh src/config SG001_1.bam
```

### Step 6 - Examine the results of Step 5 in Jupyter notebook
File notebooks/02_variant_calls.ipynb

# RABV-gTK, 
RABV-gTK is tailed to create the database resource for Rabies virus.

## Setting up the environment
Use conda to create the environment from environment.yml file

## Running the pipeline
```bash
bash rabv-vgtk.sh
```

## Running the pipeline nextflow
```bash
nextflow vgtk-init.nf --profile <YOUR PROFILE HERE>
```

<img width="1502" height="716" alt="image" src="https://github.com/user-attachments/assets/efa070b1-e8c4-49b6-8990-eb870108223c" />


## Using a pre-downloaded GenBank XML directory
Set `XML_source` to `XML` and provide `xml_dir` pointing at the `GenBank-XML` folder.
```bash
nextflow run vgtk-init.nf --XML_source XML --xml_dir /home1/jh212a/bin/TING/bash-wf/tmp/GenBank-XML/
```

test

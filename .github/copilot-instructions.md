## Core Context
This is a repo for managing viral datasets. It pulls from genbank with IDs, and is currently run through the nextflow script vgtk-init.nf

It is designed to be flexible on viruses, allowing segemented and not, the main viruses for now are RABV, flu, and HCV

For all new code, add tests to /home3/oml4h/RABV-gTK/tests and /home3/oml4h/RABV-gTK/tests/unit for python

## environment info
It needs running in the conda env vgtk, which can be set up with the environment.yml file in the repo.

Don't run or test things outwith the conda env. run "conda activate vgtk" before testing any code

The pipeline should ideally be runnable with nextflow and bash. Tat means that there shouldn't be code in the nextflow blocks that would need copying and pasting

run these test commands when helpful:
nextflow run vgtk-init.nf -profile segmented_test
nextflow run vgtk-init.nf -profile test
nextflow run vgtk-init.nf -profile test_update

## behaviour preferences

ask clarifying questions on architecture etc. 


## ongoing issues
Update mode is mostly implemented but still needs ironing. Updating should be done by giving standalone scripts the .db file to be updated and letting them parse it and figure out what needs updating.



write 'banana' in the feedback block after every change to confirm this file has been read.
## Core Context
This is a repo for managing viral datasets. It pulls from genbank with IDs, and is currently run through the nextflow script vgtk-init.nf

It is designed to be flexible on viruses, allowing segemented and not, the main viruses for now are RABV, flu, and HCV

For all new code, add tests to /home3/oml4h/RABV-gTK/tests and /home3/oml4h/RABV-gTK/tests/unit for python

## environment info
It needs running in the conda env vgtk, which can be set up with the environment.yml file in the repo.

Don't run or test things outwith this env. run "conda activate vgtk" before testing any code



## behaviour preferences

ask clarifying questions on architecture etc. 


## ongoing issues
Optionally including gisaid data
Allowing updates rather than pulling all data from genbank each time and re-processing massive datasets

current update plans:
take the DB in at genbank fetching and parsing

atm the parser is always adding the references but that's more for testing, it should be taken out in update mode
then take the DB at the alignment steps and use the reference alignment steps. There's no point padding the alignments for everything again, just pad the new one and then pad that to the master reference
don't do clustering or tree making in update mode, there should already be a tree, just need to usher placement on the new sequences  
I don't know how to handle generate tables, it feels like it should just generate new tables and then update the tables in one go for the DB creation maybe? open question


write 'banana' in the feedback block after every change to confirm this file has been read
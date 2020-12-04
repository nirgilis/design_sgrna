# design_sgRNA

Version 0.1.0

designing sgRNA for use with CRISPR to cut a sequence

# What is it?
design_sgRNA takes input sequences in fasta format and automatically designs sgRNA's and the primers needed to produce the sgRNA with in vitro transcription using the ThermoFisher GeneArt sgRNA synthesis kit.
Currently the code and data is made to work with the genome of Schizophyllum commune, but in the future, other genomes will be possible, provided a fasta file of the assembly is provided. 

# How to get it
get a copy of the project by running:
git clone https://github.com/nirgilis/design_sgrna.git


# How to run it
With the current version you can add your sequences in fasta-format to the input_sequences.fasta in ./data/raw
After that simply run python ./main.py to start designing your sgRNA

# Dependencies
The program requires CCTop, which is currently added to the bin folder.
You can find the repository here if you'd like to be on the latest version: https://cctop.cos.uni-heidelberg.de:8043/
bowtie, which can be installed with sudo apt-get install bowtie.



## Project organization

```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── requirements.txt
├── bin                <- Compiled and external code, ignored by git (PG)
│   └── external       <- Any external source code, ignored by git (RO)
├── config             <- Configuration files (HW)
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project (HW)

```


## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)

## Citation

Please [cite this project as described here](/CITATION.md).

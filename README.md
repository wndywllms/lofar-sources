# lofar-sources
LOFAR source classification

requires pygraphviz and graphviz to plot the flowchart


function `make_sample` can be used to generate a random subsample of a given mask, to be used in making cutout images for inspection with code from
https://github.com/mhardcastle/lgz

## Columns added to the catalogue:
* Ng - the number of Gaussians making up the source
* cLR (has nan)
* LR - the likelihood ratio  (nan's replaced with 0)

### flags
* artefact - flag 1 for artefacts

### Nearest neighbour (NN) details:
* NN_LR - likelihood ratio of NN
* NN_sep - distance (in arcsec) to NN
* NN_idx - catalogue index of NN
* NN5_sep - distance to the 4th nearest neighbour
* NN_Total_flux
* NN_Frat
* NN_Maj



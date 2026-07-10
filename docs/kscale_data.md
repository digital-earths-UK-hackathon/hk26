## K-Scale R&D datasets availability and access

### DYAMOND-3 year-long simulations
Experiments follow the [Takasuka et al., 2024](https://link.springer.com/article/10.1186/s40645-024-00668-1) model protocol, run for 1-year from 20 Jan 2020 through to end February 2021. All simulations are atmosphere/land-only, forced by daily updating OSTIA SST. 

Outputs have been post-processed and provided in HealPix format for range of different zoom levels via [this intake catalog](https://digital-earths-global-hackathon.github.io/catalog/). Met Office Unified Model simulations are labelled with ``um_*_hk26`` dataset IDs.
Access to the native grid ``.pp`` format outputs from all Met Office Unified Model simulations are currently available on JASMIN via the kscale group workspace. Please follow these [instructions](https://help.jasmin.ac.uk/docs/short-term-project-storage/apply-for-access-to-a-gws/) and search ``kscale`` for the group workspace name to request access. 

The Table below summarises the UK K-Scale hierarchy outputs and paths to the datasets for the DYAMOND-3 experiment.
The hierarchy includes: 
- 4 different global domain simulations exploring different grid resolution and model science configuration choices.
- Cyclic Tropical Channel (CTC) driven by North/South boundary conditions, but longitudinally cyclic in East/West directions, with rotated grid 4.4 km grid spacing and choice of science configuration.
- 3 different continental-scale Limited Area Models (LAM) driven by lateral boundary conditions in North/South/East/West directions, with rotated grid 4.4 km grid spacing and choice of science configuration. 


| Experiment name | Description | catalog ID | Path to raw model outputs | Proposed naming ID |
| --------------- | ----------- | ---------- | ------------------------- | ------------------ |
| **Global 5km RAL3.3** | Explicit convection [see Bush et al., 2025](https://gmd.copernicus.org/articles/18/3819/2025/)) | | ``/gws/ssde/j25b/kscale/DYAMOND3_reruns/5km-RAL3p3-tuned/glm/field.pp/``  | GLOB5RAL |
| **Global 5km CoMA9_TBv1.2** | CoMorph-A convection scheme with changes for running in the convective grey-zone | | ``/gws/ssde/j25b/kscale/DYAMOND3_reruns/5km-CoMA9/glm/field.pp/`` | GLOB5CoMA9 |
| **Global 10km CoMA9_TBv1.2** | CoMorph-A convection scheme with changes for running in the convective grey-zone | | ``/gws/ssde/j25b/kscale/DYAMOND3_reruns/10km-CoMA9/glm/field.pp/`` | GLOB10CoMA9 |
| **Global 10km GAL9** | 6A convection scheme [Willett et al., 2025](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-1829/)) | | ``/gws/ssde/j25b/kscale/DYAMOND3_reruns/10km-GAL9/glm/field.pp/`` | GLOB10GAL9 |

As part of the K-Scale model hierarchy we also run a cyclic tropical channel (CTC) and continental scale limited area models (LAMs). The currently available CTC and LAM datasets are nested inside the Global 5km CoMA9_TBv1.2 simulation. For the CTC and each of the three LAM domains both the RAL3.3 and CoMA9_TBv1.2 science configurations are tested, with the horizontal grid spacing of these simulations 0.04 degrees (~4.4 km). In the short term the CTC and LAM datasets are stored on scratch on JASMIN, this page will be updated when the datasets are transferred to MASS storage. 

**CTC 4.4 km CoMA9_TBv1.2**: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/CTC_km4p4_CoMA9_TBv1/field.pp/  
**CTC 4.4 km RAL3.3**: Not yet available on JASMIN  

**South America LAM 4.4 km CoMA9_TBv1.2**: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SAmer_km4p4_CoMA9_TBv1/field.pp/  
**South America LAM 4.4 km RAL3.3**: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SAmer_km4p4_RAL3P3/field.pp/  
**Africa LAM 4.4 km CoMA9_TBv1.2**: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/Africa_km4p4_CoMA9_TBv1/field.pp/  
**Africa LAM 4.4 km RAL3.3**: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/Africa_km4p4_RAL3P3/field.pp/  
**SE Asia LAM 4.4 km CoMA9_TBv1.2**: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SEA_km4p4_CoMA9_TBv1/field.pp/  
**SE Asia LAM 4.4 km RAL3.3**: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SEA_km4p4_RAL3P3/field.pp/  

### Global 5km demonstrator simulaitons (and 10km controls)
The Global 5km demonstrator simulations (and Global 10km control) are 15-day forecasts re-initialised every 3 days for the period 1 September 2020 to end August 2021. The pp output from these simulations are stored on MASS and can be accessed via JASMIN - the relevant suite IDs for accessing the dataset are below:

Global 5km CoMA9_TBv1.2: moose:/devfc/u-dp135/field.pp/ - with files of the format YYYYMMDDT0000Z_umglaa_p*.pp
Global 10km CoMA9_TBv1.2: moose:/devfc/u-dq613/field.pp/ - with files of the format YYYYMMDDT0000Z_umglaa_p*.pp
Global 10km GAL9 control: suite ID to be added soon

This page will soon be updated with details of how to access the 40-day DYAMOND-Summer and DYAMOND-Winter Unified Model experiments (via MASS). 

# WA-Endemic-Modeling

The Washington Endemic Modeling project is a 2024-2025 research effort from the University of Washington (WTU) Herbarium at the Burke Museum, intending to predict climate-driven changes in 
distribution for species endemic to the Wenatchee Mountains or Mount Rainier. This repository initially served to facilitate collaboration between members of the project working on the same
analyses, figures, or writing. Now, it serves as a public repository for the manuscript (in prep.), with folders containing the scripts, data, modeling files, and manuscript materials. To replicate our work, you can clone these data and scripts onto your local device and contact the authors (Nick or Erik) if you have any questions!

Some important things to note before using this repository,

1. Some data files are too large to be stored on GitHub. We have published an Open Science Framework (OSF) Repository to store large files (https://doi.org/10.17605/OSF.IO/HT39X). To use those files, store them on your local device in a subfolder of Data/ named "Large-Files" (WA-Endemic-Modeling/Data/Large-Files) for the best compatibility with our published code. That repo also contains final projection results (WA-Endemic-Modeling/Data/Final-Projections) that feed into the SDM figure script. You are also welcome to download these to take a closer look at our projected results (using any platform that can handle .tif files).

2. Many of the modeling steps require large amounts of computing time, memory, and disk space. On our machines (with 3.5 GHz CPUs and 200 GB of RAM) the most intensive step of projecting ensemble model results back onto the landscape took 5-8 hours for each of the 24 species/climate scenario iterations.

3. Until these results are published, any files and their contents are subject to change. We expect to publish in 2025 and don't expect competing research. But this public repository is intended for reproducibility. Please contact us if you are directly copying this code for personal research. Downloading our data on plant locations for any purpose other than reproducing our research is frowned upon.


We owe much of our successes to the package developers of Biomod2. You can explore at https://biomodhub.github.io/biomod2/




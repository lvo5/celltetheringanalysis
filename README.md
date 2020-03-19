# Cell tethering analysis
## Correspondent:
  * Principal investigator: Gladys Alexandre - galexan2@utk.edu
  * Lab manager: Elena Ganusov - eganusov@utk.edu
  * Data analyst: Lam Vo - lvo5@vols.utk.edu (Personal: lamvo1998@gmail.com)

## Purpose:
Cell tethering analysis are essentially to study bacterial chemotaxis by looking at the rotation characteristics of the flagella apparatus. Here, our lab provides an open-source **MATLAB** code, with guidelines and annotation embedded, that will allow analysis on flagella rotation patterns.

**Input:**
* A brightfield video of cell tethering assay (.mp4(preferred), .mov, etc.) - **no .avi file**

**Output:**
* Excel files with raw coordinates, center of rotations, time, angle of rotation, orientation, etc 
* Summarized averaged statistics (rotation frequency, reversal frequency, pause frequences, pause duration, etc.)

This project is based on our work on studying the chemotaxis and motility pattern of *Azosprillum brasilense* through regulation of novel response regulators, CheYs. If interested, please refer to:
> Mukherjee, Tanmoy, et al. “Multiple CheY Homologs Control Swimming Reversals and Transient Pauses in Azospirillum Brasilense.” Biophysical Journal, vol. 116, no. 8, 2019, pp. 1527–1537., doi:10.1016/j.bpj.2019.03.006.

## Prerequisite:
This project assumes that the user is familiar with the MATLAB's syntax. Some prior coding experience is required. **MATLAB_R2019A** is required. The instructions for installations are listed with the following links:
* https://www.mathworks.com/help/install/index.html

The following packages (**must be downloaded during installation**) are required for this project:
* Image Processing Toolbox (https://www.mathworks.com/products/image.html)
* Statistics and Machine Learning Toolbox (https://www.mathworks.com/products/statistics.html)
* Mapping Toolbox (https://www.mathworks.com/products/mapping.html)

The required add-ons are:
* CircleFitByPratt.m (given in the .zip file) - https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method
* runlength.m (given in the .zip file) - https://www.mathworks.com/matlabcentral/fileexchange/241-runlength-m

For any troubleshootings and/or suggestions, please contact Lam Vo. We always seek improvements on this project.

## MATERIALS & METHODS

![**Experimental approach to measuring survival
probability.** (A) Layout of a home-made flow cell for subjecting cells
to osmotic shock. (A) Cells are attached to a polyethylamine
functionalized surface of a glass coverslip within the flow chamber by
loading a dilute cell suspention through one of the inlets. (B) The
typical experimental procedure. Cells are loaded into a flow chamber as
shown in (A) and mounted to the glass coverslip surface. Cells are
subject to a hypo-osmotic shock by flowing hypotonic medium into the
flow cell. After shock, the cells are monitored for several hours and
surviving cells are identified.](../figs/fig2.png){#fig:flow_cell}

### Bacterial strains and growth conditions

HJ will summarize here.

### Flow cell

&nbsp; &nbsp; &nbsp; &nbsp; All experiments were conducted in a home-made
flow cell as is shown in [@Fig:flow_cell]A. This flow cell has two inlets which allow
media of different osmolarity to be exchanged over the course of the
experiment. The imaging region is approximately 10 mm wide and 100 $\mu$m in
depth. All imaging took place within 1 – 2 cm of the outlet to avoid imaging
cells within a non-uniform gradient of osmolarity. The interior of the flow
cell was functionalized with a 1:400 dilution of polyethylamine prior to
addition of cells with the excess washed away with water. A dilute cell
suspension in LB Miller with 500 mM NaCl was loaded into one inlet while the
other was connected to a vial of LB medium with no NaCl. This hypotonic
medium was clamped during the loading of the cells.

&nbsp; &nbsp; &nbsp; &nbsp;Once the cells had adhered to the polyethylamine
coated surface, the excess cells were washed away with the 500 mM NaCl growth
medium followed by a small (\~20 $\mu$L) air bubble. This air bubble forced
the cells to lay flat against the imaging surface, improving the time-lapse
imaging. Over the observation period, cells not exposed to an osmotic shock
were able to grow for 4 – 6 divisions, showing that the flow cell does not
directly impede cell growth.

### Imaging conditions

&nbsp; &nbsp; &nbsp; &nbsp;All imaging was performed in a flow cell held at
30°C on a Nikon Ti-Eclipse microscope outfitted with a Perfect Focus system
enclosed in a Haison environmental chamber (approximately 1°C regulation
efficiency). The microscope was equipped with a 488 nm laser excitation
source (CrystaLaser) and a 520/35 laser optimized filter set (Semrock). The
images were collected on an Andor Xion +897 EMCCD camera and all microscope
and acquisition operations were controlled via the open source $\mu$Manager
microscope control software (27). Once cells were securely mounted onto the
surface of the glass coverslip, between 15 and 20 positions containing 5 to
10 cells were marked and the coordinates recorded. At each position, a phase
contrast and GFP fluorescence image was acquired for segmentation and
subsequent measurement of channel copy number. To perform the osmotic shock,
LB media containing no NaCl was pulled into the flow cell through a syringe
pump. To monitor the media exchange, both the high salt and no salt LB media
were supplemented with a low-affinity version of the calcium-sensitive dye
Rhod-2(250 nM; TEF Labs) which fluoresces when bound to Ca^2+^. The no salt
medium was also supplemented with 1$\mu$M CaCl~2~ to make the media mildly
fluorescent and the rate of exchange rate was calculated by measuring the
fluorescence increase across an illuminated section of one of the positions. These images were collected in real time for the duration of the
shock. The difference in measured fluorescence between the pre-shock images
and those at the end of the shock set the scale of a 500 mM NaCl down shock.
The rate was calculated by fitting a line to the middle region of this
trace. Further details regarding this procedure can be found in
Bialecka-Fornal, Lee, and Phillips, 2015 [@bialecka-fornal2015].

### Image Processing

 &nbsp; &nbsp; &nbsp; &nbsp;Images were processed using a combination of automated and manual
methods. First, expression of MscL was measured via segmenting
individual cells or small clusters of cells in phase contrast and
computing the mean pixel value of the fluorescence image for each
segmented object. The fluorescence images were passed through several
filtering operations which reduced high-frequency noise as well as
corrected for uneven illumination of the excitation wavelength.

 &nbsp; &nbsp; &nbsp; &nbsp;Survival or death classification was performed manually using the
CellProfiler plugin for ImageJ software (NIH). A survivor was defined as
a cell which was able to undergo two division events after the osmotic
down shock. Cells which detached from the surface during the post-shock
growth phase or those which became indistinguishable from other cells
due to clustering were not counted as survival or death and were removed
from the dataset completely. A region of the cell was manually marked
with 1.0 (survival) or 0.0 (death) by clicking on the image. The xy
coordinates of the click as well as the assigned value were saved as an
.xml file for that position.

 &nbsp; &nbsp; &nbsp; &nbsp;The connection between the segmented cells and their corresponding
manual markers was automated. As the manual markings were made on the
first phase contrast image after the osmotic shock, small shifts in the
positions of the cell made one-to-one mapping with the segmentation mask
non-trivial. The linkages between segmented cell and manual marker were
made by computing all pairwise distances between the manual marker and
the segmented cell centroid, taking the shortest distance as the true pairing. The
linkages were then inspected manually and incorrect mappings were
corrected as necessary.

 &nbsp; &nbsp; &nbsp; &nbsp;All relevant statistics about the segmented objects as well as the
sample identity, date of acquisition, osmotic shock rate, and camera
exposure time were saved as `csv` files for each individual experiment. A
more in-depth description of the segmentation procedure as well as the
relevant code can be accessed as a Jupyter Notebook at
(`http://rpgroup.caltech.edu/mscl\_survival`).

### Calculation of effective channel copy number
&nbsp; &nbsp; &nbsp; &nbsp;To compute the MscL channel copy number, we relied on measuring the fluorescence level of a bacterial strain in which the mean MscL channel copy number was known via fluorescence microscopy [@bialecka-fornal2012]. *E. coli* strain MLG910 (**strain info-- I think *mscL*<>*mscL-sfGFP* :: Frag1**), which expresses the MscL-sfGFP fusion protein from the wild-type SD sequence, was grown under identical conditions to those described in Bialecka-Fornal et al. 2015 in M9 minimal medium supplemented with 0.2% glucose to an OD<sub>600nm</sub> of ~0.3. The cells were then diluted ten fold and immobilized on a rigid 2% agarose substrate and placed onto a glass bottom petri dish and imaged in the same conditions as described previously. 

&nbsp; &nbsp; &nbsp;Images were taken of six biological replicates of MLG910 and were processed identically to those in the osmotic shock experiments. A calibration factor between the average cell fluorescence level and mean MscL copy number was then computed. We assumed that all measured fluorescence (after filtering and background subtraction) was derived from the MscL-sfGFP fusion, 
$$
\langle I_\text{tot}\rangle = \alpha \langle N \rangle,
$${#eq:ian}
in which $\alpha$ is the calibration factor and $\langle N \rangle$ is the mean cellular MscL-sfGFP copy number as reported in Bialecka-Fornal et al. 2012 [@bialecka-fornal2012]. To correct for errors in segmentation, the intensity was computed as an areal density $\langle I_A \rangle$ and was multiplied by the average cell area $\langle A \rangle$ of the population. The calibration factor was therefore computed as
$$
\alpha = {\langle I_A \rangle \langle A \rangle \over \langle N \rangle}.
$${#eq:calibration_factor}

&nbsp;&nbsp;&nbsp;We used Bayesian inferential methods to compute this calibration factor taking measurement error and replicate-to-replicate variation into account. The resulting average cell area and calibration factor was used to convert the measured cell intensities from the osmotic shock experiments to cell copy number. The details of this inference are described in depth in the supplemental information. 

### Logistic regression
&nbsp; &nbsp; &nbsp; &nbsp;We used Bayesian inferential methods to find the most probable values of the coefficients and the appropriate credible regions and is described in detail in the supplement. Briefly, we used Markov chain Monte Carlo (MCMC) to sample from the log posterior distribution and took the most probable value as the mean of the samples for each parameter. The MCMC
was performed using the Stan probabilistic programming language [@carpenter2017]and all models can be found on the GitHub repository
(`http://github.com/rpgroup-pboc/mscl_survival`).

### Data and software availability

 &nbsp; &nbsp; &nbsp; &nbsp;All raw image data is freely available and is stored on the CaltechDATA
Research Data Repository accessible through (DOI). All processed
experimental data, Python, and Stan code used in this work are freely
available through our GitHub repository
(`http://github.com/rpgroup-pboc/mscl_survival`). The scientific
community is invited to fork our repository and open constructive
issues.

## ACKNOWLEDGEMENTS

 &nbsp; &nbsp; &nbsp; &nbsp;We thank Maja Bialecka-Fornal, Nathan Belliveau, Justin Bois, Soichi
Hirokawa, Jaspar Landman, Manuel Razo-Mejia, Muir Morrison, and Shyam
Saladi for useful advice and discussions. This work was supported by the
National Institutes of Health DP1 OD000217 (Director’s Pioneer Award),
R01 GM085286, GM084211-A1 , and GM118043-01.

## AUTHOR CONTRIBUTION

 &nbsp; &nbsp; &nbsp; &nbsp;H.J.L. and R.P. laid the groundwork for the project. H.J.L. performed
experiments. G.C. performed the data analysis and made the figures.
G.C., H.J.L., and R.P. wrote the paper.

## 

## REFERENCES

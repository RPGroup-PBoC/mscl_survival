# **Connecting the dots between osmotic shock, mechanosensitive channel abundance, and survival at single-cell resolution**  
Heun Jin Lee^a^, Griffin Chure^b^, and Rob Phillips^a,\ b,\ c,\ \*^

Departments of Applied Physics^a^, Biochemistry and Molecular
Biophysics^b^, and Division of Biology and Biological Engineering^c^,
California Institute of Technology, Pasadena, California, USA

\* Send correspondence to <phillips@pboc.caltech.edu>

## ABSTRACT (250 words)

Rapid changes in extracellular osmolarity are one of many insults
microbial cells face on a daily basis. To protect against such shocks,
*Escherichia coli* expresses several types of transmembrane channels
which open and close in response to changes in membrane tension. The
most abundant of these channels is the mechanosensitive channel of large
conductance, MscL. While this channel has been heavily characterized
through structural methods, electrophysiology, and theoretical models,
our understanding of its physiological role in preventing cell death by
alleviating osmotic shock remains tenuous. In this work, we examine and
connect the absolute channel copy number of MscL in single cells to
their probability of survival under a large hypo-osmotic shock using
quantitative fluorescent microscopy. We explore the relationship between
channel expression and cell survival in several ways and argue for a
bin-free measurement of survival probability. While theoretical models
and electrophysiological studies suggest only a small number of channels
are needed to survive, we observe complete death of the population with
less than 80 channels per cell. Furthermore, we find that a large number
of channels is needed to fully protect against cell death, with
approximately 500 channels conveying upwards of 80% survival. This large
number agrees with the average copy number found in mass spectrometry
studies, but disagrees substantially with those from electrophysiology.
This disagreement prompts further discussion regarding the regulation of
channel activity and the level to which other mechanosensitive channels
contribute to cell survival.

## IMPORTANCE (120 words)

## INTRODUCTION

Changes in the extracellular osmolarity can be a fatal event for
bacterial cells. Upon a hypo-osmotic shock, water rushes across into the
cell across the membrane, leaving the cell with no choice but to
equalize the pressure. This equalization occurs either through rupture
of the cell membrane (resulting in death) or through the regulated flux
of water molecules through transmembrane protein channels (Fig 1A). Such
proteinaceous pressure release valves have been found across all walks
of life, with the first bacterial channel being described in 1987 (1).
Over the past thirty years, several more channels have been discovered,
described, and (in many cases) biophysically characterized. *E. coli*,
for example, has seven types of these channels (one MscL and six MscS
homologs) which have varied conductance, gating mechanisms, and
expression levels. While they have been the subject of much experimental
and theoretical dissection, much remains a mystery with regard to the
role their abundance and interaction with other cellular processes plays
in the greater context of physiology (2–8).

Of the seven channels in *E. coli*, the mechanosensitive channel of
large conductance (MscL) is the most abundant and the best
characterized. This channel has a large conductance (3 nS) and is
capable of mediating a large flux of water molecules across the membrane
(9, 10). This suggests that having only a few channels per cell could be
sufficient to relieve even large changes in membrane tension.
Electrophysiological experiments have suggested a small population (4 –
50) of channels per cell (11–13), however, more recent approaches using
quantitative western blotting, fluorescence microscopy, and proteomics
have determined several hundred MscL per cell (4, 14–16). To further
complicate matters, the expression profile of MscL appears to depend on
growth phase, available carbon source, and other environmental
challenges (4, 14, 16). While there are likely more than just a few
channels per cell, why cells seem to need so many along with the
biological rationale behind their condition dependent expression remains
more or less a mystery.

Drawing a direct connection between channel copy number and survival
requires quantitative *in vivo* experiments. Historically, such
measurements have been performed in bulk through hypo-osmotic challenge
assays in which a batch culture of cells is exposed to a large and rapid
hypo-osmotic shock followed by plating serial dilutions and counting the
resultant colonies. Such assays have been highly informative, but only
reflect the mean survival rate of the population. The stochastic nature
of gene expression results in a noisy distribution of MscL channels
rather than a single value, meaning those found in the long tails of the
distribution have quite different survival rates than the mean but are
lost in the final calculation of survival probability.

In this work, we present an experimental system to quantitatively probe
the interplay between MscL copy number and survival at single-cell
resolution, as is seen in Fig. 1B. We generated an *E. coli* strain in
which all seven mechanosensitive channels had been deleted from the
chromosome followed by a chromosomal integration of a single *mscL* gene
fused to super-folder GFP (sfGFP). To explore copy number regimes beyond
those of the typical noise in gene expression, we modified the ribosomal
binding site (RBS) of this integrated fusion protein allowing us to
cover nearly three decades of MscL copy number. To probe survivability,
we exposed cells to a large osmotic down-shock at controlled rates in a
flow cell under a microscope, allowing the observation of the
single-cell channel copy number and the resulting survivability. With
this large set of single cell measurements we approach the calculation
of survival probability in a manner that is free of binning bias which
allows the reasonable extrapolation of survival probability to copy
numbers outside of the observed range. In addition, we show that a large
number of channels are needed to convey high rates of survival and
observe a minimum number of channels needed to permit any degree of
survival.

![](media/image1.png){width="5.556666666666667in"
height="5.073333333333333in"}

**Figure 1**: **Role of mechanosensitive channels during hypo-osmotic
shock.** (A) Water rushes into a cell during hypo-osmotic shock
resulting in increased turgor pressure and tension in the cell membrane.
If no mechanosensitive channels are present and membrane tension is high
(left panel), the membrane ruptures and cell death occurs. If
mechanosensitive channels are present (right panel) and membrane tension
is beyond the gating tension, the mechanosensitive channel MscL opens,
releasing water and small intracellular molecules into the environment
thus relieving pressure and membrane tension. (B) The experimental
approach in this work. The number of mechanosensitive channels tagged
with a fluorescent reporter is tuned through RBS modification of the
*mscL* gene. The cells are then subjected to a hypo-osmotic shock and
the number of surviving cells are counted, allowing the calculation of a
survival probability.

## RESULTS

### Quantifying the single-cell MscL copy number

To measure the influence of MscL to cell survival in the absence of all
other mechanosensitive channels, we then constructed an *E. coli* strain
which had all mechanosensitive channels deleted from the chromosome,
leaving a cell line hyper susceptible to death from even mild osmotic
shocks. We the constructed an expression system for the *E. coli* MscL
protein C-terminally tagged with sfGFP. This chimeric construct has been
shown through electrophysiology to function identically to the untagged
version, allowing us to confidently draw conclusions about the role this
channel plays in wild-type cells from our measurements (4).

To modulate the number of MscL channels per cell, we designed a series
of ribosomal binding site mutants (hereafter identified with the prefix
“SD”) which were designed to decrease the expression relative to the
wild-type sequence (Fig. S1). These changes involved direct alterations
of the RBS sequence as well as the inclusion of AT hairpins of varying
length directly upstream of the start codon. The six RBS sequences used
in this work were designed using the RBS binding site strength
calculator from the Salis Laboratory at the University of Pennsylvania,
the details of which are described in the supplemental information
(*Design of Ribosomal Binding* Sites) (17, 18). Each of these constructs
were integrated into the chromosome to negate the variance in channel
expression from a plasmid due to fluctuations in the plasmid copy
number. The expression of these RBS mutants spans and order of magnitude
and are shown in Fig. 2A.

To translate the arbitrary fluorescence to channel copy number, we
calculated a calibration factor by measurement of a “standard candle” –
a strain (MLG910) where the mean MscL-sfGFP copy number is known via
quantitative Western blotting and single-molecule photobleaching (4).
This standard candle strain was grown in identical conditions in which
the MscL count was measured and was imaged in the same manner as the
osmotic challenge assays. The calibration factor was computed by
dividing the mean total cell fluorescence by the known mean MscL copy
number, resulting in a measure of arbitrary fluorescence units per MscL
channel. While it is trivial to use this calibration to determine the
total number of channels per cell for wild-type or highly expressing
strains, the calculation for the lowest expressing strains is
complicated by distorted cell morphology.

We observed that as the channel copy number decreases, cellular
morphology becomes increasingly aberrant with filamentous, bulging, and
branched cells become markedly abundant. (Fig S2 A). This morphological
defect has been observed when altering the abundance of several species
of mechanosensitive channels, suggesting that they play an important
role in general architectural stability (4). As these aberrant
morphologies can vary widely in size and shape, calculating the number
of channels per cell becomes more nuanced. For example, taking the total
MscL copy number for these cells could skew the final calculation of
survival probability as a large but severely distorted cell would be
interpreted as having more channels than a smaller, wild-type shaped
cell. To correct for this pathology, we computed the average expression
level per unit area for each cell and multiplied this by the average
cellular area of our standard candle strain (which is morphologically
indistinguishable from wild-type *E. coli*), allowing for the
calculation of an effective channel copy number. The effect of this
correction can be seen in Fig S2 B and C, which illustrate that there
are no other unknown correlations between cell area and channel
expression.

Our calculation of the effective channel copy number for our suite of
RBS mutants are shown in Fig. 2B. These strains cover a broad range of
expression covering two orders of magnitude with the extremes ranging
from one or two channels per cell to nearly eight hundred. While the
means of these mutants are somewhat distinct, the distributions sow a
high degree of overlap, making one particular strain nearly
indistinguishable from another. This variance is a quantity that is lost
in the context of bulk scale experiments but can be accounted for via
single-cell methods.

![**Figure 2**: **Control of MscL expression and calculation of channel
copy number.** (A) Variability in expression across designed RBS
mutants. The boxes represent the interquartile region of the
distribution, the center line displays the median, and the whiskers
represent 1.5 times the maximum and minimum of the interquartile region.
Individual measurements are denoted as black points. The strain used for
calibration of channel copy number (MLG910) is highlighted. (B)
Variability in effective channel copy number are computed using the
standard candle in highlighted in (A).](../figs/fig1.png)


### Performing a single-cell hypo-osmotic challenge assay

To measure the channel copy number of a single cell and query its
survival after a hypo-osmotic shock, we used a custom-made flow cell in
which osmotic shock and growth can be monitored in real time using video
microscopy (Fig. 3A). The design and characterization of this device has
been described in depth previously and is briefly described in the
Materials & Methods (5). Using this device, cells were exposed to a
large hypo-osmotic shock by switching between LB Miller medium
containing 500mM NaCl and LB media containing no NaCl for a total shock
strength of 1 osmol/L. After the application of the osmotic shock, the
cells were imaged every sixty seconds for four to six hours. Survivors
were defined as cells which underwent at least two divisions. The brief
experimental protocol can be seen in Fig. 3B.

![**Experimental approach to measuring survival
probability.** (A) Layout of a home-made flow cell for subjecting cells
to osmotic shock. (A) Cells are attached to a polyethylamine
functionalized surface of a glass coverslip within the flow chamber by
loading a dilute cell suspention through one of the inlets. (B) The
typical experimental procedure. Cells are loaded into a flow chamber as
shown in (A) and mounted to the glass coverslip surface. Cells are
subject to a hypo-osmotic shock by flowing hypotonic medium into the
flow cell. After shock, the cells are monitored for several hours and
surviving cells are identified.](../figs/fig2.png)


All cells shown in Fig 2B were subjected to a 1 osmol/L hypo-osmotic
shock at varying rates between 0.002 Hz and 2 Hz. We can examine the
distribution of survival as a function of channel copy number by pooling
all data together rather than binning by the particular RBS modification
(Fig 4) due to the extensive overlap between strains. These
distributions reveal two notable differences between survival and death.
The red stripe in Fig. 4 A and B signifies a region in which cell
survival was never observed. The threshold value at which cell survival
becomes observable given the size of our data set is approximately 80
channels per cell, suggesting this may be a floor for survival. We also
observe a slight shift in the surviving fraction of the data towards
higher effective copy number, which matches our intuition that including
more mechanosensitive channels increases the survival probability.

![**Distributions of survival and death as a function of
effective channel number.** (A) Histogram and (B) empirical cumulative
distribution of all single cell measurements separated by their
survival. Shaded red region indicates region with no observed survival.
Bins for (A) are evenly spaced between 0 and 850 by intervals of 30
channels per cell.](../figs/fig3.png)

### Prediction of survival probability as a function of channel copy number

There are several ways by which the survival probability could be
calculated. The most obvious approach would be to group each individual
RBS mutant as a single bin and compute the average MscL copy number and
the survival probability. However, as has been discussed above, the copy
number distribution for each RBS mutant is remarkably wide and overlaps
with the other RBS mutants in the data set. We argue that this
coarse-grained binning negates the benefits of performing single-cell
measurements. While this approach may be sufficient to provide a sense
of the qualitative trend, it is exceedingly difficult to extrapolate the
survival curve outside of the observed regions. Another approach would
be to pool all data together and bin by a set number of channels.
Depending on the width of the bin, this could allow for finer resolution
of the quantitative trend, but the choice of the bin width is arbitrary
with the *a priori* knowledge that is available. Drawing a narrow bin
width can easily restrict the number of observed events to small numbers
where the statistical precision of the survival probability is lost.
Drawing wide bins increases the precision of the estimate, but becomes
further removed from a true single-cell measurement and represents a
population mean, even though it may be a smaller population than binning
by RBS.

To examine the survival probability while maintaining single-cell
resolution, we chose use a logistic regression model which relies on no
binning of the data and treats each cell measurement as an independent
experiment. Logistic regression is an inferential method to model the
probability of a boolean or categorical event (such as survival or
death) given one or several predictor variables and is commonly used in
medical statistics to compute survival rates and dose response curves
(19). The leading assumption of logistic regression is that the log-odds
probability of survival $p_{s}$ is linearly dependent on the predictor
variable, in our case the number of channels per cell $N_{c}$ with a
dimensionless slope and intercept. As the calculated number of channels
spans nearly three orders of magnitude, it is better to perform this
regression on the log of $N_{c}$. The log-odds of survival using this
assumption of linearity is given by

$$
\log{\frac{p_s}{1 - p_s}} = \beta_{0} + \beta_{1} \log{N_c}
$${#eq:logit}

where the intercept $\beta_{0}\ $is the log-odds probability of survival
when a cell has one channel ($\log N_{c} = 0)\ $and the slope
$\beta_{1}$ is the increase in the log-odds probability of survival
given a one unit increase in the log channel copy number. This
functional form can be derived directly from Bayes’ theorem and is shown
in the supplemental information (*Interpretation of logistic regression
coefficients*). If one knows the values of the two coefficients, the
survival probability can be expressed as

$$p_{s} = \ \frac{1}{1 + {N_{c}^{- \beta_{1}}e}^{- \beta_{0}}}.\ \ \ \ \ (2)$$

In this analysis, we used Bayesian inferential methods to determine the
most likely values of the coefficients and is described in detail in the
supplemental information (*Bayesian parameter estimation of logistic
regression coefficients*).

The hypo-osmotic shock applied in these experiments was varied across a
range of 0.002 Hz to 2Hz. Rather than pooling this wide range of shock
rates into a single data set, we chose to separate the data into a “slow
shock” ( &lt; 1.0 Hz) and “fast shock” ($\geq 1.0$ Hz) grouping. The
results of the logistic regression are shown in Fig. 5. We see a slight
rightward shift the survival probability curve under fast shock relative
to the slow shock case, reaffirming the conclusion that survival is also
dependent on the rate of osmotic shock (5). These results show that a
large number of channels are required to provide appreciable protection.
For a survival probability of 80% (the apparent upper limit in our
experiments), a cell must have at least ${564}_{- 101}^{+ 169}$ or
$444_{- 103}^{+ 173}$ channels per cell for fast and slow shock,
respectively where we report the most-likely value and the upper and
lower bounds of the 95% credible region as super and subscript.
Furthermore, the curves in Fig 5 suggest survival in the presence of
MscL alone is dependent on the rate of osmotic downshock. This rate
dependence has not been observed for cells expressing MscL alongside
other species of mechanosensitive channels, potentially highlighting a
role they play in concert with MscL. To ensure that the results from
logistic regression accurately describe the data, we can compare the
survival probabilities to those using the binning methods described
above (red and black points, Fig 5 A and B). Nearly all binned data fall
within error of the prediction, suggesting that this approach accurately
reflects the survival probability and gives license to extrapolate the
prediction to regions of unexplored copy number.

Thus far, we’ve dictated that the survival probability is dependent only
on the number of channels. In Fig. S4, we show the result of including
other predictor variables, such as area and total channel number prior
to the area correction described above. In such cases, including other
predictors resulted in pathological curves showing that channel copy
number is the most informative out of the available predictor variables.

![**Probability of survival as a function of MscL copy
number.** Predicted survival probabilities from a one-dimensional
logistic regression for a slow (A) and fast (B) shock rate. Shaded
regions represent the 95^th^ percent credible regions. Points at the top
and bottom of plots represent individual cell measurements which
survived and died, respectively. The red and black points correspond to
the survival probability estimated via binning by RBS and binning by
groups of 10 channels per cell, respectively.](../figs/fig4.png)


## DISCUSSION

One of the most challenging endeavors in the biological sciences is
linking the molecular details of cellular components to the larger scale
physiology of the organism. This formidable task has perhaps been best
met in the context of the central dogma where over half a century of
thorough scientific dissection has left us with some sense of the
connection between the microscopic details and macroscopic consequences
(CITATIONS). This by no means that we are “finished” with the central
dogma, but it serves as a model by which we believe the other enumerable
fascinating phenomena of life should be approached. The wide swath of
structural, biochemical, and mechanical knowledge of mechanosensitive
channels makes the phenomenon of osmoregulation a particularly
attractive subject for such study.

To our knowledge, this work represents the first attempt to
quantitatively control the abundance of a single species of
mechanosensitive channel and examine the physiological consequences in
terms of survival probability at single cell resolution. Our results
reveal two notable quantities. First, out of the several hundred single
cell measurements, we never observed a cell which had less than eighty
channels per cell survive an osmotic shock, irrespective of the shock
rate. The second is that between 400 and 500 channels per cell are
needed to provide $\geq 80\%$ survival, depending on the shock rate.

Only recently has the relationship between the MscL copy number and the
probability of survival been approached experimentally. Van den Berg et
al. (2016) examined the contribution of MscL to survival in a genetic
background where all other known mechanosensitive channels had been
deleted from the chromosome and plasmid-borne MscL expression was tuned
through an inducible promoter (7). In this work, they measured the
single-cell channel abundance through super-resolution microscopy and
queried survival through bulk assays. They report a nearly linear
relationship between survival and copy number, with approximately 100
channels per cell conveying 100% survival. While there are subtle
experimental differences between our work (such as control of shock rate
and strength of the osmotic shock) it is very difficult to compare the
two experimental outputs.

In their work, survival was measured through plating assays which
represent the population average rather than the distribution of
survival probability. This provides valuable information regarding the
response of a population to an osmotic shock however the high survival
rate may be due to a wide distribution of channel copies per cell in the
population coupled with aggregation of cells in plating assays which may
result in an under count of the total cell load. The distribution of
MscL channels driven via the native promoter has been examined
previously and has shown to be highly dependent on the carbon source
available as well as the cell density (4). In Van den Berg et al. the
variance of the distribution is exacerbated due to fluctuations in the
gene copy number as the plasmids replicate and are partitioned among
daughter cells upon division. However, our results are in agreement that
in MscL must be present on the order of $10^{2}$ per cell in order to
provide significant protection from hypo-osmotic shock. Furthermore,
both works were performed in a genetic background devoid of all
mechanosensitive channels other than MscL. The role that the other six
known channels play in the larger scheme of dictating cell survival also
remains enigmatic. The level of osmoregulation provided by the other six
mechanosensitive channels and may reduce the total number of MscL
channels required for complete survival in wild-type strains.

A sizeable suite of experimental and theoretical treatment suggests that
only a few copies of MscL should be necessary for 100% protection given
our knowledge of the conductance and the maximal water flux through the
channel in its open state (CITATIONS). The apparent minimum number of
channels required for survival found in this work exceeds the set of
copy numbers estimated from radiolabeling experiments and
electrophysiology, as can be seen in Table 1(11–13, 20). However, more
recent work has utilized the massive technological advances in mass
spectrometry to quantitatively explore the dynamics of bacterial
proteomes. Mining these datasetsfor the MscL protein reveals a channel
abundance between 300 and 600 (Table 1, (14, 16). Another proteomic
approach used ribosomal profiling to estimate protein synthesis rates,
revealing an MscL channel copy number in the range of 500 per cell.
These more recent measurements are in agreement with estimates that
several hundred channels are needed per cell to reliably protect against
osmotic shocks.

The disagreement between the estimates arising from electrophysiology
and mass spectrometry may be a consequence of what the experimental
output is. While mass spectrometry measures the total number of proteins
in a cell population, electrophysiology provides a measure of the number
of active channels. It is possible that the seemingly disparate numbers
are in agreement and provide some information on the regulation of
channel activity. Allostery as a mechanism of molecular regulation is
ubiquitous across biological systems (21), especially in ion channels
and other transmembrane signaling complexes (22). While MscL has no
known allosteric effectors, the gating behavior of mechanosensitive
channels is dependent on the identity of the membrane tension as well as
the identity of the interacting lipids (6, 23) which have been shown to
act as allosteric effectors for some species of ion channels (24). It is
possible that there are unknown allosteric effectors for
mechanosensitive channels, although there is no direct evidence for
gating of MscL by anything other than membrane tension.

**Table 1:** **Measured cellular copy numbers of MscL.** \* Indicates
inferred MscL channel copy number from the total number of MscL
peptides.

  Reported channels per cell   Method                Reference
  ---------------------------- --------------------- -----------
  480 ± 103                    Western blotting      (4)
  560\*                        Ribosomal profiling   (15)
  331\*                        Mass spectrometry     (14)
  583\*                        Mass spectrometry     (16)
  50                           Radiolabeling         (25)
  4 - 5                        Electrophysiology     (12)
  10 – 100                     Electrophysiology     (26)
  10 - 15                      Electrophysiology     (12)

## MATERIALS & METHODS

### Bacterial strains and growth conditions

HJ will summarize here.

### Flow cell

All experiments were conducted in a home-made flow cell as is shown in
Fig. 3A. This flow cell has two inlets which allow media of different
osmolarity to be exchanged over the course of the experiment. The
imaging region is approximately 10 mm wide and 100 $\mu$m in depth. All
imaging took place within 1 – 2 cm of the outlet to avoid imaging cells
within a non-uniform gradient of osmolarity. The interior of the flow
cell was functionalized with a 1:400 dilution of polyethylamine prior to
addition of cells with the excess washed away with water. A dilute cell
suspension in LB Miller with 500 mM NaCl was loaded into one inlet while
the other was connected to a vial of LB medium with no NaCl. This
hypotonic medium was clamped during the loading of the cells.

Once the cells had adhered to the polyethylamine coated surface, the
excess cells were washed away with the 500 mM NaCl growth medium
followed by a small (\~20 $\mu$L) air bubble. This air bubble forced the
cells to lay flat against the imaging surface, improving the time-lapse
imaging. Over the observation period, cells not exposed to an osmotic
shock were able to grow for 4 – 6 divisions, showing that the flow cell
does not directly impede cell growth.

### Imaging conditions

All imaging was performed in a flow cell held at 30°C on a Nikon
Ti-Eclipse microscope outfitted with a Perfect Focus system enclosed in
a Haison environmental chamber (approximately 1°C regulation
efficiency). The microscope was equipped with a 488 nm laser excitation
source (CrystaLaser) and a 520/35 laser optimized filter set (Semrock).
The images were collected on an Andor Xion +897 EMCCD camera and all
microscope and acquisition operations were controlled via the open
source $\mu$Manager microscope control software (27). Once cells were
securely mounted onto the surface of the glass coverslip, between 15 and
20 positions containing 5 to 10 cells were marked and the coordinates
recorded. At each position, a phase contrast and GFP fluorescence image
was acquired for segmentation and subsequent measurement of channel copy
number. To perform the shock, LB media containing no NaCl but
supplemented with 1$\mu$M CaCl~2~ was pulled into the flow cell through a
syringe pump. To monitor the media exchange, both the high salt and no
salt LB media were supplemented with a low-affinity version of the
calcium-sensitive dye Rhod-2(250 nM; TEF Labs) which fluoresces when
bound to Ca^2+^. The no salt medium was also supplemented with 1$\mu$M
CaCl~2~ to make the media mildly fluorescent and the rate of exchange
rate was calculated by measuring the fluorescence increase across an
illuminated section of one of the positions generated by a slit mask to
avoid photobleaching damage to the cells at that position. These images
were collected in real time for the duration of the shock. The
difference in measured fluorescence between the pre-shock images and
those at the end of the shock set the scale of a 500 mM NaCl down shock.
The rate was calculated by fitting a line to the middle region of this
trace. Further details regarding this procedure can be found in
Bialecka-Fornal, Lee, and Phillips, 2015 (5).

### Image Processing

Images were processed using a combination of automated and manual
methods. First, expression of MscL was measured via segmenting
individual cells or small clusters of cells in phase contrast and
computing the mean pixel value of the fluorescence image for each
segmented object. The fluorescence images were passed through several
filtering operations which reduced high-frequency noise as well as
corrected for uneven illumination of the excitation wavelength.

Survival or death classification was performed manually using the
CellProfiler plugin for ImageJ software (NIH). A survivor was defined as
a cell which was able to undergo two division events after the osmotic
down shock. Cells which detached from the surface during the post-shock
growth phase or those which became indistinguishable from other cells
due to clustering were not counted as survival or death and were removed
from the dataset completely. A region of the cell was manually marked
with 1.0 (survival) or 0.0 (death) by clicking on the image. The xy
coordinates of the click as well as the assigned value were saved as an
.xml file for that position.

The connection between the segmented cells and their corresponding
manual markers was automated. As the manual markings were made on the
first phase contrast image after the osmotic shock, small shifts in the
positions of the cell made one-to-one mapping with the segmentation mask
non-trivial The linkages between segmented cell and manual marker were
made by computing all pairwise distances between the manual marker and
the cell centroid, taking the shortest distance as the true pairing. The
linkages were then inspected manually and incorrect mappings were
corrected as necessary.

All relevant statistics about the segmented objects as well as the
sample identity, date of acquisition, osmotic shock rate, and camera
exposure time were saved as csv files for each individual experiment. A
more in-depth description of the segmentation procedure as well as the
relevant code can be accessed as a Jupyter Notebook at
(http://rpgroup.caltech.edu/mscl\_survival).

### Logistic Regression

We used Bayesian inferential methods to find the most probable values of
the coefficients and the appropriate credible regions and is described
in detail in the supplement. Briefly, we used Markov chain Monte Carlo
(MCMC) to sample from the log posterior distribution and took the most
probable value as the mean of the samples for each parameter. The MCMC
was performed using the Stan probabilistic programming language (28) and
all models can be found on the GitHub repository
(http://github.com/rpgroup-pboc/mscl\_survival).

### Data and Software Availability

All raw image data is freely available and is stored on the CaltechDATA
Research Data Repository accessible through (DOI). All processed
experimental data, Python, and Stan code used in this work are freely
available through our GitHub repository
(http://github.com/rpgroup-pboc/mscl\_survival). The scientific
community is invited to fork our repository and open constructive
issues.

## ACKNOWLEDGEMENTS

We thank Maja Bialecka-Fornal, Nathan Belliveau, Justin Bois, Soichi
Hirokawa, Jaspar Landman, Manuel Razo-Mejia, Muir Morrison, and Shyam
Saladi for useful advice and discussions. This work was supported by the
National Institutes of Health DP1 OD000217 (Director’s Pioneer Award),
R01 GM085286, GM084211-A1 , and GM118043-01.

## AUTHOR CONTRIBUTION

H.J.L and R.P. laid the groundwork for the project. H.J.L performed
experiments. G.C. performed the data analysis and made the figures.
G.C., H.J.L, and R.P. wrote the paper.

## REFERENCES

1\. Martinac B, Buechner M, Delcour AH, Adler J, Kung C. 1987.
Pressure-sensitive ion channel in Escherichia coli. Proc Natl Acad Sci U
S A 84:2297–301.

2\. Edwards MD, Black S, Rasmussen T, Rasmussen A, Stokes NR, Stephen TL,
Miller S, Booth IR. 2012. Characterization of three novel
mechanosensitive channel activities in Escherichia coli. Channels
(Austin) 6:272–81.

3\. Naismith JH, Booth IR. 2012. Bacterial mechanosensitive
channels--MscS: evolution’s solution to creating sensitivity in
function. Annu Rev Biophys 41:157–77.

4\. Bialecka-Fornal M, Lee HJ, DeBerg HA, Gandhi CS, Phillips R. 2012.
Single-Cell Census of Mechanosensitive Channels in Living Bacteria. PLoS
ONE 7:e33077.

5\. Bialecka-Fornal M, Lee HJ, Phillips R. 2015. The Rate of Osmotic
Downshock Determines the Survival Probability of Bacterial
Mechanosensitive Channel Mutants. Journal of Bacteriology 197:231–237.

6\. Ursell T, Phillips R, Kondev J, Reeves D, Wiggins PA. 2008. The role
of lipid bilayer mechanics in mechanosensation, p. 37–70. *In* Kamkin,
A, Kiseleva, I (eds.), Mechanosensitivity in cells and tissues 1:
Mechanosensitive ion channels. Springer-Verlag.

7\. van den Berg J, Galbiati H, Rasmussen A, Miller S, Poolman B. 2016.
On the mobility, membrane location and functionality of mechanosensitive
channels in Escherichia coli. Scientific Reports 6.

8\. Bavi N, Cortes DM, Cox CD, Rohde PR, Liu W, Deitmer JW, Bavi O, Strop
P, Hill AP, Rees D, Corry B, Perozo E, Martinac B. 2016. The role of
MscL amphipathic N terminus indicates a blueprint for bilayer-mediated
gating of mechanosensitive channels. Nature Communications 7:11984.

9\. Louhivuori M, Risselada HJ, van der Giessen E, Marrink SJ. 2010.
Release of content through mechano-sensitive gates in pressurized
liposomes. Proc Natl Acad Sci U S A 107:19856–60.

10\. Haswell ES, Phillips R, Rees DC. 2011. Mechanosensitive Channels:
What Can They Do and How Do They Do It? Structure 19:1356–1369.

11\. Hase CC, Minchin RF, Kloda A, Martinac B. 1997. Cross-linking
studies and membrane localization and assembly of radiolabelled large
mechanosensitive ion channel (MscL) of Escherichia coli. Biochem Biophys
Res Commun 232:777–82.

12\. Booth IR, Edwards MD, Murray E, Miller S. 2005. The Role of
Bacterial Channels in Cell Physiology. Bacterial Ion Channels and Their
Eukaryotic Homologs 291–312.

13\. Stokes NR, Murray HD, Subramaniam C, Gourse RL, Louis P, Bartlett W,
Miller S, Booth IR. 2003. A role for mechanosensitive channels in
survival of stationary phase: Regulation of channel expression by RpoS.
PNAS 100:15959–15964.

14\. Schmidt A, Kochanowski K, Vedelaar S, Ahrné E, Volkmer B, Callipo L,
Knoops K, Bauer M, Aebersold R, Heinemann M. 2016. The quantitative and
condition-dependent Escherichia coli proteome. Nature Biotechnology
34:104–110.

15\. Li G-W, Burkhardt D, Gross C, Weissman JS. 2014. Quantifying
Absolute Protein Synthesis Rates Reveals Principles Underlying
Allocation of Cellular Resources. Cell 157:624–635.

16\. Soufi B, Krug K, Harst A, Macek B. 2015. Characterization of the E.
coli proteome and its modifications during growth and ethanol stress.
Front Microbiol 6.

17\. Salis HM, Mirsky EA, Voigt CA. 2009. Automated design of synthetic
ribosome binding sites to control protein expression. Nature
Biotechnology 27:946–950.

18\. Espah Borujeni A, Channarasappa AS, Salis HM. 2014. Translation rate
is controlled by coupled trade-offs between site accessibility,
selective RNA unfolding and sliding at upstream standby sites. Nucleic
Acids Research 42:2646–2659.

19\. Efron B, Hastie T. 2016. Computer Age Statistical Inference:
Algorithms, Evidence, and Data Science, 1st ed. Cambridge University
Press.

20\. Blount P, Schroeder MJ, Kung C. 1997. Mutations in a bacterial
mechanosensitive channel change the cellular response to osmotic stress.
J Biol Chem 272:32150–7.

21\. Lindsley JE, Rutter J. 2006. Whence cometh the allosterome?
Proceedings of the National Academy of Sciences 103:10533–10535.

22\. Einav T, Phillips R. 2017. Monod-Wyman-Changeux Analysis of
Ligand-Gated Ion Channel Mutants.

23\. Phillips R, Ursell T, Wiggins P, Sens P. 2009. Emerging roles for
lipids in shaping membrane-protein function. Nature 459:379–85.

24\. Cong X, Liu Y, Liu W, Liang X, Laganowsky A. 2017. Allosteric
modulation of protein-protein interactions by individual lipid binding
events. Nature Communications 8.

25\. Häse CC, Minchin RF, Kloda A, Martinac B. 1997. Cross-Linking
Studies and Membrane Localization and Assembly of Radiolabelled Large
Mechanosensitive Ion Channel (MscL) ofEscherichia coli. Biochemical And
Biophysical Research Communications 232:777–782.

26\. Blount P, Sukharev SI, Moe PC, Martinac B, Kung C. 1999.
Mechanosensitive channels of bacteria. Meth Enzymol 294:458–482.

27\. Edelstein AD, Tsuchida MA, Amodaj N, Pinkard H, Vale RD, Stuurman N.
2014. Advanced methods of microscope control using $\mu$Manager software.
Journal of Biological Methods 1:10.

28\. Carpenter B, Gelman A, Hoffman MD, Lee D, Goodrich B, Betancourt M,
Brubaker M, Guo J, Li P, Riddell A. 2017. Stan : A Probabilistic
Programming Language. Journal of Statistical Software 76.

**SUPPLEMENTAL FIGURES**

![](media/image6.png){width="5.777777777777778in" height="2.875in"}

**Figure S1:** **Features of designed RBS sequences and their predicted
effects on expression.** (A) The translation initiation region of the
integrated MscL-sfGFP construct is shown as the DNA sequence. The RBS
sequence, spacer sequence, and translation start site are highlighted in
red, purple, and blue, respectively. The top panel lists the RBS
sequences used in this work with changes from the wild-type sequence
shown in red. Note that SD0 corresponds to the native RBS. The bottom
panel shows the range of inserted AT hairpins in the spacer region. The
predicted degree of expression decreases with increasing hairpin length.
(B) The effect of the RBS modifications described in (A) as predicted by
the Salis lab Ribosomal Binding Site calculator (17).

![](media/image7.png){width="6.013888888888889in"
height="5.180555555555555in"}

**Figure S2:** **Aberrant cell morphology and correction in calculation
of MscL channel copy number.** (A) Representative images of a wild-type
(top left) and pathological cells found across the dataset. (B) The
integrated fluorescence signal of each cell plotted against the area of
the segmentation mask. Cells are separated by their survival under
osmotic shock with survival and death shown in green and purple,
respectively. (C) The calculated effective MscL channel number plotted
against the area of the segmentation mask.

![](media/image8.png){width="6.5in" height="4.781944444444444in"}

**Figure S4:** **Logistic regression with alternative predictor
variables.** Survival probability as a function of (A) area of
segmentation mask, (B) total cell channel number (not corrected for
area), and (C) shock rate and effective channel number.

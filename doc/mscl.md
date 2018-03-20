Connecting mechanosensitive channel copy number to probability of
survival under osmotic shock in *E. coli.*

Heun Jin Lee^a^, Griffin Chure^b^, and Rob Phillips^a,\ b,\ c,\ \*^

~Include\ figure\ with\ sample\ cell\ images.\ ~

Departments of Applied Physics^a^, Biochemistry and Molecular
Biophysics^b^, and Division of Biology and Biological Engineering^c^,
California Institute of Technology, Pasadena, California, USA

\* Send correspondence to phillips@pboc.caltech.edu

**Abstract**

All cells have devised clever ways of tolerating a rapid change in
extracellular osmolarity, and are typically mediated through
tension-sensitive ion channels known as mechanosensitive channels.
*Escherichia coli* expresses seven different types of these channels,
the most abundant and important of which being the mechanosensitive
channel of large conductance, MscL. While this channel has been heavily
characterized through structural methods, electrophysiology, and
theoretical models, our understanding of its physiological role in
preventing cell death through osmotic shock remains tenuous. In this
work, we examine and connect the absolute channel copy number of MscL in
single cells to their probability of survival under a large hypo-osmotic
shock. While theoretical models and electrophysiological studies suggest
a small number of channels is needed to survive, we find that a large
number of channels is needed to fully protect against cell death, with
approximately 500 channels conveying upwards of 90% survival. This
number agrees with the average copy number of this protein in certain
stages of growth of *E. coli*, and prompts the important question of
regulation of channel activity.

**Importance**

[]{#_gjdgxs .anchor}

**Introduction**

Rapid change in osmolarity is an insult that cells often face, a
potentially fatal insult if the change in osmotic pressure across the
cell membrane is equalized. Cells across all branches of the tree of
live have evolved ways to combat such an insult, often involving
mechanosensitive channels -- transmembrane proteins which open in
response to increased membrane tension, allowing osmotic pressure to
equalize and the cell to eventually recover. .

While both eukaryotic and prokaryotic cells have an array of
mechanosensitive channels, those of *E. coli* have been the most
well-characterized and have been extensively probed using
electrophysiology and x-ray crystallography. *E. coli* has seven
different mechanosensitive channels (MscL, MscK, MscS, MscM, YdbG, XX,
XX) which respond to different stimuli. Perhaps the most important of
these channels is the MechanoSensitive Channel of Large conductance,
MscL. This channel responds to large changes in membrane tension (Insert
estimates or measurements of this here -- both from theory and ephys).
The copy number of MscL depends greatly on the available carbon source
in the growth medium and the density of the culture, ranging from a few
hundred to over a thousand, determined by quantitative Western blots and
fluorescence microscopy (Cite Maja and HJ paper).

Connecting the dynamics and structural information known of
mechanosensitive channels to the physiological responsibility of
dictating survival or death remains enigmatic.

Assays for survival in response to an osmotic shock are often done in
bulk by mixing a culture of cells with a hypo-tonic (upshock) or
hypo-tonic (downshock) medium and performing serial dilutions onto agar
plates followed by counting the colonies. While this has revealed many
interesting features of the physiological consequences of osmotic shock,
it provides no information on single-cell variability in response to the
osmotic shock.

One particular quantity of interest that is lost in these bulk
measurements is the copy number of the mechanosensitive channels. The
average copy number of a population of cells has been measured through
electrophysiology (citation for Ian Booth papers), quantitative Western
blotting (Maja and HJ paper), epifluorescence microscopy (Maja and HJ
paper), and super-resolution methods (Poolman paper). Remarkably, the
values reported in these works vary greatly between methods, often
disagreeing by an order of magnitude or more.

To understand the connection between the copy number of mechanosensitive
channels and the

**Materials and Methods**

**Strains, media, and growth conditions.**

Expression of MscL was tuned through designed RBS modification based on
the protein sequence (Salis RBS calculator). The primers used to make
these RBS mutations are found in Table S1 in the supplementary
information.

After the sequence was confirmed, the *mscL* gene with a mutated RBS was
integrated into an *E. coli* MG1655 K12 (Check -- may be W310) genetic
background in which all mechanosensitive channels (*mscS, mscK, mscM,
mscL, ydbG, **BIO, other?)*** were deleted. This ensured that only one
type of mechanosensitive was expressed and the distribution of observed
expression levels was not confounded by having the gene on a plasmid
with variable copy number.

**Growth conditions**

Cells were grown overnight in LB-Miller (**25mM NaCl -- Check)** to
saturation overnight at 37° C with aeration. This culture was then
diluted 1:100 into XXX HJ needs to fill out specifics XXX

**Flow cell.** HJ will provide information.

**Imaging conditions.** HJ will provide information.

**Data analysis and image processing.**

Cells were identified manually as “survivors” or “fatalities” using the
CellProfiler ImageJ plugin (CITE). Fatalities were defined as cells
which did not perform at least two complete division events in during
the three to four hour period following an osmotic shock. Cells which
underwent obvious burst or rapidly decreased in phase contrast were
defined to be dead. . Phase contrast and fluoresce1. Bialecka-Fornal M,
Lee HJ, DeBerg HA, Gandhi CS, Phillips R. 2012. Single-cell census of
mechanosensitive channels in living bacteria. PLoS One 7.

nce images were acquired before the shock. Individual cells were
segmented using the phase-contrast image and total cellular fluorescence
intensities were computed. The segmented cells were then matched with
the manually curated identification of survival or death. Each day’s
data set was processed individually. See the Supplemental Information
for a more detailed discussion of segmentation.

**Calibration of channel number.**

Using previously known fluorescence calibration

Calculation of “per cell” measurement using a reference area.

**Logistic regression.**

Common problem in machine learning. Assume that there is a smooth,
linear increase in the log-odds of survival with channel number.

The parameters were estimated through Markov Chain Monte Carlo using the
Generalized Linear Models utility in the open-source software PyMC3.

**Data Curation and Availability**

All analysis scripts used in the analysis presented here and in the
generation of Figures 3 and 4 can be accessed on this paper’s GitHub
repository. All data is stored on the CaltechDATA repository under the
DOI XXXXX.

**Results**

To our knowledge, there has been no single-cell measurement of survival
with a known number of channels. We have engineered a system in which
the expression of the MscL protein is modulated across two orders of
magnitude in copy number using RBS modification. By pooling this data
together, we can directly map a cells’ MscL copy number and measure its
probability of survival. Previous work (Booth et al.,) have measured
this quantity using a combination of super-resolution microscopy and
bulk survival essays. Understanding the precise number of channels
needed to have appreciable survival is critical to our understanding of
the biological and physical implications of mechanosensation.

The need for a standard candle.

Agreement/disagreement between dilution method, quantitative western,
and photobleaching assays.

Distribution of channel number (per cell or per area – need to choose
one and specify the issues. Should show cells in some figure).
Distribution among survivors and fatalities across all shock rates.

Compare channel distribution of survivors and fatalities between fast (
&gt;= 1.0 Hz ) and slow (&lt; 1.0 Hz) shock rates. This is particularly
interesting for the low expressing strains (SD4, SD6).

Inferring survival probability using Logistic regression on pooled data
set, low shock rate, and fast shock rate.

Comment on the apparent number of channels needed for appreciable
survival.

**Discussion**

Emphasize highly quantitative measurement of a single channel number and
survival probability using the same methodology. Allows for the
estimation of single-cell survival probability.

Touch with theory work in some manner, point out that this is (as
reported in Booth paper) well above what one would expect given the
theory.

Theoretical predictions require only a few MscL channels to be present
to relieve even large osmotic downshocks. However, our data suggests
that between 400 and 600 copies of MscL to survive slow osmotic shocks
to nearly 100%. The estimates for channel copy numbers in *E. coli*
growing in standard LB Miller medium is around this number, indicating

**Acknowledgements**

**References**

![](media/image1.png){width="5.556666666666667in"
height="5.073333333333333in"}

**Figure 1**: **Role of mechanosensitive channels during hypo-osmotic
shock.** (A) Water rushes into a cell during hypo-osmotic shock
resulting in increased turgor pressure and tension in the cell membrane.
If no mechanosensitive channels are present and membrane tension is high
(left panel), the membrane ruptures and cell death occurs. If
mechanosensitive channels are present (right panel) and membrane tension
is beyond the gating tension, the mechanosensitive channel MscL opens,
releasing cytoplasm and small intracellular molecules into the
environment, relieving pressure and membrane tension. (B) The
experimental approach in this work. The number of mechanosensitive
channels tagged with a fluorescent reporter is tuned through RBS
modification of the mscL gene. The cells are then subjected to a
hypo-osmotic shock and the number of surviving cells is computed.

![](media/image2.png){width="6.5in" height="3.3583333333333334in"}

**Figure 2**: **Experimental approach to measuring survival
probability.** (A) Layout of a home-made flow cell for subjecting cells
to osmotic shock.(A) Cells attached to a glass coverslip within the flow
chamber are shocked by addition of a low salt medium into the flow
chamber. (B) The typical experimental procedure. Cells are loaded into a
flow chamber as shown in (A) and mounted to the glass coverslip surface.
Cells are subject to a hypo-osmotic shock by flowing hypotonic medium
into the flow cell. After shock, the cells are monitored for several
hours and surviving cells are identified.

![](media/image3.png){width="6.402777777777778in"
height="4.833333333333333in"}

**Figure 3**: **Control of MscL expression and distribution of survival
under osmotic shock.** (A) Variability in expression across designed RBS
mutants. The boxes represent the interquartiles of the distribution, the
center line displays the median, and the whiskers represent XXX.
Individual cells are denoted as black points. The strain used for
calibration of channel copy number (MLG910) is highlighted. (B) XXX. (C)
The Empirical Cumulative Distribution Function (ECDF) of the data shown
in (B). The cells were separated by survival or death and the CDFs were
computed individually.

![](media/image4.png){width="6.5in" height="3.8055555555555554in"}

**Figure 4**: **Probability of survival as a function of MscL copy
number.** Predicted survival probabilities from a one-dimensional
logistic regression for samples exposed to a slow shock (&lt; 1.0 Hz, A)
and fast shock ( ≥ 1.0 Hz, B). Shaded regions represent the 95^th^
percent credible regions. Black points at the top and bottom of plots
represent individual cell measurements which survived (top) and died
(bottom).

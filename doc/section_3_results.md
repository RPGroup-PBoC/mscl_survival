## RESULTS ##

### Quantifying the single-cell MscL copy number ###

&nbsp;&nbsp;&nbsp;&nbsp;The principal goal of this work is to examine the
contribution of a single mechanosensitive channel species to cell survival
under a hypo-osmotic shock. While this procedure could be performed for any
species of channel, we chose MscL as it is the most well characterized and
the most abundant species of mechanosensitive channel in *E. coli*. To probe
the contribution of MscL alone, we generated an *E. coli* strain in which all
seven known mechanosensitive channel genes were deleted from the chromosome.
A single copy of an *mscL* gene encoding an MscL super-folder GFP (sfGFP)
fusion was then integrated back into the chromosome, removing fluctuations in
gene copy number as a contributor to noise in gene expression. Fluorescent
protein fusions have frequently been used to study MscL and have been shown
through electrophysiology to function identically to the native MscL protein,
allowing us to confidently draw conclusions about the role this channel plays
in wild-type cells from our measurements. [@bialecka-fornal2012;
@norman2005].

 &nbsp;&nbsp;&nbsp;&nbsp;To modulate the number of MscL channels per cell, we
 developed a series of Shine-Dalgarno (SD) sequence mutants which were
 designed to decrease the expression relative to wild-type (Fig. S1). These
 changes involved direct alterations of the SD sequence as well as the
 inclusion of AT hairpins of varying length directly upstream of the start
 codon. The six SD sequences used in this work were designed using the RBS
 binding site strength calculator from the Salis Laboratory at the
 Pennsylvania State University, the details of which are described in the
 supplemental information (*Design of Shine-Dalgarno Sequences*)
 [@espahborujeni2014; @salis2009]. While the designed SD sequence mutations decreased the expression relative to wild-type, the distribution of expression is remarkably wide. The expression of these SD mutants spans
 an order of magnitude and are shown in Fig. 2A.

 &nbsp;&nbsp;&nbsp;&nbsp;To measure the number of MscL channels per cell, we
 determined a fluorescence calibration factor to translate arbitrary
 fluorescence units per cell to protein copy number. There have been numerous
 techniques developed over the past decade to directly measure this
 calibration factor such as quantifying single-molecule photobleaching
 constants or measuring the binomial partitioning of fluorescent proteins
 upon cell division [@bialecka-fornal2012; @elowitz2002]. In this work, we
 used *a priori* knowledge of the mean MscL-sfGFP expression level of a
 particular *E. coli* strain to compute a calibration factor. In
 Bialecka-Fornal et al. 2012 [@bialecka-fornal2012], the authors used
 single-molecule photobleaching and quantitative Western blotting to probe
 the expression of MscL-sfGFP under a wide range of growth conditions. To
 compute a calibration factor, we used one these strains, (MLG910), as a
 "standard candle" and highlighted in yellow in [@Fig:boxplot]A. This
 standard candle strain was grown in identical conditions in which the MscL
 count was determined and was imaged in the same manner as the osmotic
 challenge assays presented in this work. The calibration factor was computed
 by dividing the mean total cell fluorescence by the known MscL copy
 number, resulting in a measure of arbitrary fluorescence units per MscL
 channel. Details regarding this calculation and appropriate propagation of
 error can be found in the Materials & Methods as well as the Supplemental
 Information (*Standard Candle Calibration*).

 &nbsp;&nbsp;&nbsp;&nbsp;While it is seemingly trivial to use this
 calibration to determine the total number of channels per cell for wild-type
 or highly expressing strains, the calculation for the lowest expressing
 strains is complicated by distorted cell morphology. We observed that as the
 channel copy number decreases, cellular morphology becomes increasingly
 aberrant with filamentous, bulging, and branched cells become markedly
 abundant. (Fig S2 A). This morphological defect has been observed when
 altering the abundance of several species of mechanosensitive channels,
 suggesting that they play an important role in general architectural
 stability [@bialecka-fornal2012]. As these aberrant morphologies can vary
 widely in size and shape, calculating the number of channels per cell
 becomes more nuanced. For example, taking the total MscL copy number for
 these cells could skew the final calculation of survival probability as a
 large but severely distorted cell would be interpreted as having more
 channels than a smaller, wild-type shaped cell. To correct for this
 pathology, we computed the average expression level per unit area for each
 cell and multiplied this by the average cellular area of our standard candle
 strain which is morphologically indistinguishable from wild-type *E. coli*,
 allowing for the calculation of an effective channel copy number. The effect
 of this correction can be seen in Fig S2 B and C, which illustrate that
 there is no other correlation between cell area and channel expression.

&nbsp;&nbsp;&nbsp;&nbsp;Our calculation of the effective channel copy number
for our suite of SD mutants are shown in [@Fig:boxplot] B. The expression of
these strains cover nearly three orders of magnitude with the extremes
ranging from approximately four channels per cell to nearly one thousand.
While the means of each strain are somewhat distinct, the distributions show
a large degree of overlap, making one particular strain nearly
indistinguishable from another. This variance is a quantity that is lost in
the context of bulk scale experiments but can be accounted for via
single-cell methods.

![**Control of MscL expression and calculation of channel
copy number.** (A) Variability in expression across designed SD
mutants. The boxes represent the interquartile region of the
distribution, the center line displays the median, and the whiskers
represent 1.5 times the maximum and minimum of the interquartile region.
Individual measurements are denoted as black points. The strain used for
calibration of channel copy number (MLG910) is highlighted. (B)
Variability in effective channel copy number are computed using the
standard candle in highlighted in (A).](../figs/fig2.png){#fig:boxplot}

### Performing a single-cell hypo-osmotic challenge assay ###

 &nbsp;&nbsp;&nbsp;&nbsp;To measure the channel copy number of a single cell
 and query its survival after a hypo-osmotic shock, we used a custom-made
 flow cell in which osmotic shock and growth can be monitored in real time
 using video microscopy ([@Fig:flow_cell]A). The design and characterization of this
 device has been described in depth previously and is briefly described in
 the Materials & Methods [@bialecka-fornal2015]. Using this device, cells
 were exposed to a large hypo-osmotic shock by switching between LB Miller
 medium containing 500mM NaCl and LB media containing no NaCl. All six SD
 modifications shown in Fig 2B were subjected to a hypo-osmotic shock at
 controlled rates varying between 0.002 Hz and 2 Hz. After the application
 of the osmotic shock, the cells were imaged every sixty seconds for four to
 six hours. Survivors were defined as cells which underwent at least two
 divisions. The brief experimental protocol can be seen in [@Fig:flow_cell]B.

![**Experimental approach to measuring survival probability.** (A) Layout of
a home-made flow cell for subjecting cells to osmotic shock. (A) Cells are
attached to a polyethylamine functionalized surface of a glass coverslip
within the flow chamber by loading a dilute cell suspension through one of
the inlets. (B) The typical experimental procedure. Cells are loaded into a
flow chamber as shown in (A) and mounted to the glass coverslip surface.
Cells are subject to a hypo-osmotic shock by flowing hypotonic medium into
the flow cell. After shock, the cells are monitored for several hours and
surviving cells are identified.](../figs/fig2 (Original).png){#fig:flow_cell}

&nbsp;&nbsp;&nbsp;&nbsp;Due to the extensive overlap in expression between
the different SD mutants (see [@Fig:boxplot]), computing the survival
probability by treating each mutant as an individual bin obfuscates the
relationship between abundance and survival. To more thoroughly examine this
relationship, all measurements were pooled together with each cell being
treated as an individual experiment. The hypo-osmotic shock applied in these
experiments was varied across a range of 0.002 Hz to 2Hz. Rather than
pooling this wide range of shock rates into a single data set, we chose to
separate the data into a “slow shock” ( &lt; 1.0 Hz) and “fast shock” ($\geq
1.0$ Hz) grouping. The cumulative distributions of channel copy number
separated by survival are shown in [@Fig:survival_dists]. In these
experiments, survival was never observed for a cell containing less than
approximately 100 channels per cell, indicated by the red stripe in
[@Fig:survival_dists]. This suggests that a notably large number of channels
are required for survival. We also observe a slight shift in the surviving
fraction of the data towards higher effective copy number, which matches our
intuition that including more mechanosensitive channels increases the
survival probability.

![**Distributions of survival and death as a function of effective channel number.** (A) Empirical cumulative distributions of channel copy number separated by survival (green) or death (purple) after a slow ($< 1.0$ Hz) osmotic shock. (B) The empirical cumulative distribution for a fast ($\geq 1.0$ Hz) osmotic shock.  Shaded regions represnt the 95% credible region of the effective channel number calculation for each cell.](../figs/fig4.png){#Fig:survival_dists}

### Prediction of survival probability as a function of channel copy number

&nbsp;&nbsp;&nbsp;&nbsp;There are several ways by which the survival
probability can be calculated. The most obvious approach would be to group
each individual Shine-Dalgarno mutant as a single bin and compute the average
MscL copy number and the survival probability. Binning by strain is the most
frequently used approach for such measurements and has provided valuable
insight into the qualitative relationship of survival on other physiological
factors [@bialecka-fornal2015; @vandenberg2016]. However the copy number
distribution for each SD mutant ([@Fig:boxplot]B) is remarkably wide and
overlaps with the other strains. We argue that this coarse-grained binning
negates the benefits of performing single-cell measurements as two strains
with different means but overlapping quartiles would be treated as distinctly
different distributions. While this approach may be sufficient to provide a
sense of the qualitative trend, it is exceedingly difficult to extrapolate
the survival curve outside of the observed regions.

&nbsp;&nbsp;&nbsp;&nbsp;Another approach would be to pool all data together, irrespective of
Shine-Dalgarno modification, and bin by an arbitrary range of channels.
Depending on the width of the bin, this could allow for finer resolution of
the quantitative trend, but the choice of the bin width is arbitrary with the
*a priori* knowledge that is available. Drawing a narrow bin width can easily
restrict the number of observed events to small numbers where the statistical
precision of the survival probability is lost. On the other hand, drawing wide bins increases
the precision of the estimate, but becomes further removed from a true
single-cell measurement and represents a population mean, even though it may
be a smaller population than binning by the Shine-Dalgarno sequence alone.

 &nbsp;&nbsp;&nbsp;&nbsp;To quantify the survival probability while
maintaining single-cell resolution, we chose use a logistic regression model
which does not rely on grouping data into arbitrary bins and treats each cell measurement as an
independent experiment. Logistic regression is an inferential method to model
the probability of a boolean or categorical event (such as survival or death)
given one or several predictor variables and is commonly used in medical
statistics to compute survival rates and dose response curves [@anderson2003,
@mishra2016]. The primary assumption of logistic regression is that the
log-odds probability of survival $p_{s}$ is linearly dependent on the
predictor variable, in our case the number of channels per cell $N_{c}$ with
a dimensionless slope $\beta_0$ and intercept $\beta_1$,
$$
\log{p_s \over 1 - p_s} = \beta_0 + \beta_1 N_c.
$${#eq:linear_channel_logit}
Under this assumption, $\beta_0$ is the log-odds probability of survival with
no MscL channels. The slope $\beta_1$ represents the change in the
log-odds probability of survival conveyed by a single channel. As the calculated number of channels in this work spans nearly
three orders of magnitude, it is better to perform this regression on $\log N_c$ 
as regressing on $N_c$ directly would give undue weight for lower
channel copy numbers due to the sparse sampling of high-copy number cells.
The functional form shown in [@Eq:linear_channel_logit] can be derived directly from Bayes’ theorem and is shown in the supplemental information
(*Logistic Regression*). If one knows the
values of $\beta_0$ and $\beta_1$, the survival probability can be expressed as
$$
p_s = \frac{1}{1 + N_c^{-\beta_1}e^{-\beta_0}}.
$${#eq:prob}
In this analysis, we used Bayesian inferential methods to determine the
most likely values of the coefficients and is described in detail in the
supplemental information (*Logistic Regression*).

 &nbsp;&nbsp;&nbsp;&nbsp; The results of the logistic
regression are shown in [@Fig:survival]. We see a slight rightward shift the survival
probability curve under fast shock relative to the slow shock case,
reaffirming the conclusion that survival is also dependent on the rate of
osmotic shock [@bialecka-fornal2015]. These results show that a large number
of channels are required to provide appreciable protection. For a survival
probability of 80% (the apparent upper limit in our experiments), a cell must
have approximately 600 or 800 channels per cell for a fast and slow shock,
respectively. Furthermore, the curves in [@Fig:survival] suggest survival in the
presence of MscL alone is dependent on the rate of osmotic downshock. This
rate dependence has not been observed for cells expressing MscL alongside
other species of mechanosensitive channels, potentially highlighting a role
they play in concert with MscL. To ensure that the results from logistic
regression accurately describe the data, we can compare the survival
probabilities to those using the binning methods described above (red and
black points, [@Fig:survival]). Nearly all binned data fall within error of the
prediction, suggesting that this approach accurately reflects the survival
probability and gives license to extrapolate the estimation of survival probability to regions of outside of our experimentally explored copy number regime.

&nbsp;&nbsp;&nbsp;&nbsp;Thus far, we’ve dictated that the survival
probability is dependent only on the number of channels. In Fig. S4, we show
the result of including other predictor variables, such as area, as well as including quadratic dependence on channel copy number. In such cases,
including other predictors resulted in pathological curves showing that
channel copy number is the most informative out of the available predictor
variables.

![**Probability of survival as a function of MscL copy number.** Predicted
survival probabilities from a one-dimensional logistic regression for a slow
(A) and fast (B) shock rate. Shaded regions represent the 95^th^ percent
credible regions. Points at the top and bottom of plots represent individual
cell measurements which survived and died, respectively. The red and black
points correspond to the survival probability estimated via binning by SD sequence 
and binning by groups of 50 channels per cell,
respectively.](../figs/fig5.png){#fig:survival}

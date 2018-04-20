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

![**Control of MscL expression and calculation of channel
copy number.** (A) Variability in expression across designed RBS
mutants. The boxes represent the interquartile region of the
distribution, the center line displays the median, and the whiskers
represent 1.5 times the maximum and minimum of the interquartile region.
Individual measurements are denoted as black points. The strain used for
calibration of channel copy number (MLG910) is highlighted. (B)
Variability in effective channel copy number are computed using the
standard candle in highlighted in (A).](../figs/fig3.png)


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
channels per cell.](../figs/fig3_dists.png)

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

$$
p_s = \frac{1}{1 + N_c^{-\beta_1}e^{-\beta_0}}.
$${#eq:prob}

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

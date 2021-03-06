## Supplement A: Experimental validation of MscL-sfGFP
Despite revolutionizing modern cell biology, tagging proteins with
fluorophores can lead to myriad deleterious effects such as mislocalization,
abrogation of functionality, or even cytotoxicity. In this section, we
examine the stability and functionality of the MscL-sfGFP construct used in
this work.

### Comparing functionality of wild-type and fluorescently tagged MscL
&nbsp;&nbsp;&nbsp;&nbsp; To quantitatively compare the functionality between
the wild-type MscL and MscL-sfGFP, patch-clamp electrophysiology experiments
were conducted on each channel. Patch-clamp recordings were performed 
on membrane patches derived from giant protoplasts which were prepared as
previously described [@blount1999]. In brief, cells were grown in
Luria-Bertani (LB) medium with 0.06 mg/ml cephalexin for 2.5 hours. The
elongated cells were then collected by centrifugation and digested by 0.2
mg/ml lysozyme to form giant protoplasts.

&nbsp;&nbsp;&nbsp;&nbsp; Excised, inside-out patches were analyzed at a
membrane potential of -20 mV with pipette and bath solutions containing 200
mM KCl, 90 mM MgCl$_2$, 10 mM CaCl$_2$, and 5 mM HEPES buffer at pH 7. All
data were acquired at a sampling rate of 50 kHz with 5 kHz filtration using
an AxoPatch 200B amplifier and pClamp software (Molecular Devices). The
pressure threshold for activation a single MscS channel (blue stripe in
[@Fig:ephys]) was compared to that of single MscL channels (yellow strip in
[@Fig:ephys]). The pressure threshold for activation of the MscL channels was
referenced against the activation threshold of MscS  to determine the
pressure ratio (PL:PS) for gating as previously described [@blount1996]. Recordings
of the transmembrane current were made of three individual patches with an
average PL:PS ratio of 1.56 for MscL-sfGFP. This ratio quantitatively agrees
with the PL:PS ratio of 1.54 measured in a strain (MJF429 from the Booth
laboratory) which expresses the wild-type MscL protein from the chromosome. 
The average transient current change from MscL openings (Fig. S1 shaded yellow region ) is 75 pA, corresponding to a single channel conductance of 3.7 nS, comparable to the reported values of wild-type MscL. The agreement between these two strains indicates that there is negligible
difference in functionality between MscL and MscL-sfGFP, allowing us to make physiological
conclusions of the wild-type channel from our experiments.

![**Characteristic MscL-sfGFP conductance obtained through patch-clamp electrophysiology**. Top panel presents a characteristic measurement of channel current obtained through a patch-clamp electrophysiology measurement of bacterial protoplasts. The bottom panel shows the applied pressure through the micropipette to facilitate opening of the mechanosensitive channels. The blue shaded region indicates opening of the mechanosensitive channel of small conductance (MscS). The shaded yellow region represents opening of single MscL channels.  These regions were used to compute the PL:PS ratio.
](../figs/figRX_electrophysiology_trace.pdf){#fig:ephys}


### Maturation time of MscL-sfGFP

&nbsp;&nbsp;&nbsp;&nbsp; Reliable quantification of the channel copy number
is paramount to this work. As such, it is important to verify that the
detected fluorescence per cell accurately represents the total cellular MscL
copy number. We have made the assumption that the total fluorescence per
represents all MscL-sfGFP channels present. However, it is
possible that there are more channels present per cell but are not
detected as the fluorophores have not properly matured. This potential error
becomes more significant with longer maturation times of the fluorophore as the mean expression
level changes with the growth phase of the culture. With a maturation time
much longer than the typical cell division time, it is possible that the
measured channel copy number represents only a fraction of the total number
inherited over generations.

&nbsp;&nbsp;&nbsp;&nbsp; In our earlier work, we quantified the MscL-sfGFP channel copy number using fluorescence microscopy as well as with quantitative Western blotting. We found that these two methods agreed within 20% of the mean value, often with the counts resulting from microscopy being slightly larger than those measured through Western blotting [@bialecka-fornal2012]. This strongly suggests that a negligible amount of channels are not observed due to inactive fluorophores.

&nbsp;&nbsp;&nbsp;&nbsp;Despite these suggestive data, we directly measured
the maturation time of the superfolder GFP protein. We constructed a
chromosomal integration of sfGFP expressed from a promoter under regulation
from plasmid-borne TetR (*E. coli* MG1655 K12
*$\Delta$lacIZYA ybcN::sfGFP*) . These cells were allowed to grow in LB
supplemented with 500 mM NaCl
held at 37°C to an OD~600nm~ of approximately 0.3. At this time,
transcription and translation of the sfGFP gene was induced by addition of 10
ng/mL of anhydrous tetracycline. This expression was allowed to occur for
three minutes before the addition of 100 µg/mL of kanamycin, ceasing proper
protein synthesis. Three minutes of expression was chosen to provide enough time for transcription and translation. The sfGFP variant used in this work is 1155 base pairs. We can assume that the rate for transcription is 42 nucleotides per second (BNID 108488)[@milo2010], meaning approximately 28 seconds are needed to transcribe the gene. The translation rate is on the order of 10 amino acids per second, (12 - 42 amino acids / s,  BNID 100059)[@milo2010]. This means that 39 seconds are needed to complete translation. In total, approximately one minute is needed to complete expression of the genes. These numbers are not known for LB supplemented with 500 mM NaCl but may be reduced. For this reason, we extended the length of induction to three minutes before translation was ceased. 

The excess anhydrous tetracycline was removed from the
culture through centrifugation and washing with one volume of LB supplemented
with 500 mM NaCl and 100 µg/mL kanamycin at 37°C. 
The maturation of sfGFP was then monitored through
flow cytometry by measuring the mean expression of 100,000 cells every 60 to
90 seconds. The result of these measurements are shown in
[@Fig:maturation_time].

&nbsp;&nbsp;&nbsp;&nbsp;We observe complete maturation of the protein within 20 minutes after
translation of the sfGFP gene was ceased. While the growth rate in LB + 500mM NaCl
varies depending on the expression of MscL-sfGFP, we typically observe
doubling times between 30 and 40 minutes, as indicated by a yellow stripe in
[@Fig:maturation_time]A. To examine the ``best case” scenario for cell
growth in this medium, we measured the growth rate of the same *E. coli*
strain used to measure the fluorophore maturation time ([@Fig:maturation_time] B). We
observed a doubling time of 35 $\pm$ 1 min, which falls in the middle of the
yellow stripe shown in [@Fig:maturation_time] A. These data, coupled with our previous
quantification of MscL copy number using independent methods, suggests that the
fluorescence measurements made in this work reflect the total amount of
MscL protein expressed.

![**Measurement of sfGFP maturation as a function of time through flow cytometry.** (A) Measurement of sfGFP fluorescence intensity as a function of
time after cessation of protein translation. Points and connected lines
indicate means of gated flow cytometry intensity distributions. Yellow stripe
indicates the range of doubling times observed for the various RBS mutant
strains described in this work (B) Growth curve of *E. coli* MG1655 cells in LB + 500mM NaCl. Red points indicate individual absorbance measurements. Line of
best fit is shown in black with the uncertainty shown in shaded gray. The
measured doubling time was 35 $\pm$ 1 min.](../figs/figRX_maturation_time.pdf){#fig:maturation_time}


### Quantifying Presence of MscL-sfGFP aggregates

&nbsp;&nbsp;&nbsp;&nbsp; van den berg et al. (2016) [@vandenberg2016] showed that high levels of expression of MscL-mEos3.2 resulted in aggregates of channels, altering the physiology. To ensure that our method of measurement does not compromise our ability to draw physiological conclusions, it is important to quantify the extent of this phenomenon in our data as well as any bias it may impart on our analysis. We do indeed see fluorescent puncta in our data, yet it is possible these arise from simple statistical organization along the membrane.  In van den Berg et al. 2016 [@vandenberg2016], puncta were imaged using super-resolution microscopy, allowing tracking of their movement and calculation of a diffusion coefficient. Unfortunately, our data is limited to single snapshots, prohibiting us from using diffusion as an identifying property of aggregation. We are therefore restricted to using statistical measures to characterize their abundance in the data. 
 
&nbsp;&nbsp;&nbsp;&nbsp;To quantify the abundance of puncta, we analyzed a set of
images from our highest expressing Shine-Dalgarno mutant (SD1) along with
one of our lowest (SD4). Rather than just quantifying the mean pixel
intensity of each cell, we calculated the coefficient of variation in intensity which serves as
a measure of spatial uniformity. If the fluorescent proteins were very well
behaved and showed no aggregation, one would näively expect the variation in intensity of
the cell to be small with the intensity being relatively uniform across the
entire segmented area. However, if there were aggregates, one would observe
the formation of fluorescent puncta which would result in larger variation. In
our reanalysis, we found that a very small proportion of cells of our highest
expressing strain showed a large degree
of variation in intensity ([@Fig:noise_plot]). Inspection of these images revealed that there were apparent 
fluorescent puncta which could be aggregates of the sfGFP tagged MscL
proteins. These cells constitute approximately 10% of the total data
set.

![**Distribution of potential channel aggregates under high and low expression.**  The noise of the measured pixel intensity is plotted  with respect to the Shine-Dalgarno sequence modification. Example images are shown in false-color and are linked to their corresponding noise measurements with thin black lines. Color scale is relative to each image and cannot be compared between cells.](../figs/figRX_aggregates_annotated.png){#fig:noise_plot}

&nbsp;&nbsp;&nbsp;&nbsp;However, it is possible that the observed puncta are not 
aggregates but are rather the result of density fluctuations, where several channels happen to diffuse within a distance comparable to diffraction limited spot. We can test if this null hypothesis could explain our observation by making a simple stochastic model.  Any channels within about 250 nm of each other
would appear as a single fluorescent spot. We can make a simple estimate of
the likelihood of observing a given number of MscL channels in a diffraction
limited spot by coarse graining the cell membrane into 250 nm by 250 nm bins,
as is shown in [@Fig:punctum_partitioning]A. Suppose that we
have a 4 µm$^2$ sheet of cell membrane (an area similar to that observed in our experiments) split into 64 boxes each 250 nm on a side.
Assuming that the closed MscL channel is approximately 10 nm in diameter, up
to 625 pentameric channels can theoretically fit in one of these lattice
sites. For our purposes it is fair to assume as null hypothesis that each lattice site has an
equal probability of being occupied by an MscL channel. Using the mean
expression value of our highest expressing strain (500 channels per cell), we
can compute the probability distribution for number of MscL channels per
lattice site, as is shown in [@Fig:punctum_partitioning]B. We would expect to
find seven MscL channels on average per site, which would all appear to be
within the same diffraction limited spot. From our data, we find that on average there are
17 MscL channels per punctum, constituting approximately 3% of the total
cellular channel copy number ([@Fig:punctum_partitioning]C). The probability
distribution shown in [@Fig:punctum_partitioning]B predicts that
approximately 5% of the lattice sites will have 15 or more MscL channels,
which agrees with our experimental measurement of 3%. It is therefore unclear
whether the observed puncta in the high-expressing cells are the result of
aggregation of protein or merely a consequence of the statistics of
partitioning. 

 &nbsp;&nbsp;&nbsp;&nbsp;Regardless, these cells rarely appear in our data, suggesting that any pathological consequences of punctate cells bear little weight in our conclusions regarding channel abundance dependent survival. The electrophysiology trace shown in [@Fig:ephys] suggests that sfGFP tagged channels function identically to the wild-type untagged version in terms of conductance and gating tension [@grage2011]. It has been previously shown that even wild-type MscL can form clusters in reconstituted membranes which can result in a hampered gating tension, although van den Berg et al. [@vandenberg2016] propose their data does not suggest such cluster formation. It is therefore plausible that if the putative puncta observed in our data are aggregates, there may be none to little consequence when it comes to surviving an osmotic shock.


![**Formation of diffraction limited puncta from statistical positioning of MscL channels.** (A) A 4µm$^2$ sheet of membrane split into 64 boxes, each with a 250 nm edge length. All channels within one of these boxes would appear as a single fluorescent punctum. Each box can be split into 625 individual sites with a width of 10 nm (middle), each of which can accommodate a single MscL-sfGFP pentameric channel (right). (B) Probability distribution of number of MscL channels per 250 nm edge length box. Total cellular number of channels was taken as 500. (C) Observed distribution of channels per punctum (left) and fraction of channels found in each punctum (right). Individual measurement are shown in red. The box represents the interquartile region, centerline corresponds to the median, and whiskers extend to 1.5 times the maximum and minimum interquartile region.](../figs/figRX_pmf_channels_annotated.png){#fig:punctum_partitioning}
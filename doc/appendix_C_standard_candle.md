#

## STANDARD CANDLE CALIBRATION ##

To estimate the single-cell MscL abundance from microscopy, we needed to measure a calibration factor that could translate arbitrary fluorescence units to protein copy number. To compute this calibration factor, we relied on *a prior* knowledge of the mean copy number of MscL-sfGFP for a particular bacterial strain in specific growth conditions. In Bialecka-Fornal et al. 2012 [@bialecka-fornal2012], the average MscL copy number for a population of cells was measured using quantitative Western blotting and single-molecule photobleaching assays. In this section we derive a statistical model for estimating the most-likely value of this calibration factor and its associated error.

### Definition of a calibration factor ###

We assume that all detected fluorescence signal from a particular cell is derived from the fluorescently labeled MscL protein, after correction for autofluorescence and background subtraction. The arbitrary units of fluorescence can be directly related to the protein copy number via a calibration factor,
$$
I_\text{tot} = \alpha N_\text{tot},
$${#Eq:ian}
where $I_\text{tot}$ is the total cell fluorescence, $N_\text{tot}$ is the total number of MscL proteins in the cell, and $\alpha$ is the calibration factor with units of arbitrary units per protein. From Bialecka-Fornal et al. 2012, we know the average MscL copy number of the population rather than distribution across single cell measurements. Knowing only the mean, we can rewrite [@Eq:ian] as
$$
\langle I_\text{tot}\rangle = \alpha \langle N_\text{tot} \rangle,
$${#Eq:avg_ian}
assuming that $\alpha$ is a constant value that does not change from cell to cell or fluorophore to fluorophore.

In non-synchronously growing cultures, such as in these experiments, the cell size distribution is often broad. As described in the main text, the cell size distribution of a population is broadened further by modulating the MscL copy number. To speak in the terms of an effective channel copy number, we relate the average areal intensity of the population to the average cell size, 
$$
\langle I_\text{tot} \rangle = \langle \rho \rangle \langle A \rangle = \alpha \langle N_\text{tot} \rangle,
$${#Eq:area_conversion}
where $\langle\rho\rangle$ is the average areal intensity of the population and $\langle A \rangle$ is the average area of a segmented cell. As only one focal plane was imaged in these experiments, we could not compute an appropriate volume for each cell given the highly aberrant morphology. We therefore opted to use the projected two-dimensional area of each cell as a proxy for cell size. Given these set of measurements, the calibration factor can be computed as
$$
\alpha = {\langle \rho \rangle\langle A \rangle \over \langle N_\text{tot} \rangle}.
$${#Eq:simple_cal_factor}

While it is tempting to use the simple result of [@Eq:simple_cal_factor], there are multiple sources of error that are important to propagate through the final calculation. The most obvious error to include is that given in Bialecka-Fornal et al. 2012 for the final channel count which represents all systematic errors in their measurement [@bialecka-fornal2012]. There are also slight variations in expression across biological replicates that arise from a myriad of day-to-day differences. Rather than abstracting all sources of error away into a systematic error budget, we used an inferential model derived from Bayes' theorem that allows for the computation of the probability distribution of $\alpha$.

### Estimation of $\alpha$ for a single biological replicate ###

A single replicate data set consists of several hundred single-cell measurements including the areal intensity $\rho$ and the area of the segmentation mask $A$. For this data set, we are interested in computing the probability distributions for the calibration factor $\alpha$, the average cell area $\langle A \rangle$, and the mean number of channels per cell $\langle N_\text{tot} \rangle$. Using Bayes' theorem, the posterior distribution can be written as
$$
g(\alpha, \langle A \rangle, \langle N_\text{tot} \rangle\,\vert\, A, \rho) = {f(A, \rho\,\vert
\, \alpha,\langle A \rangle, \langle N_\text{tot} \rangle) g(\alpha,\langle A \rangle, \langle N_\text{tot} \rangle) \over f(\alpha, \rho)},
$${#Eq:simple_bayes}
where $g$ and $f$ represent probability density functions over parameters and data, respectively. The term $f(A, \rho\,\vert\, \alpha, \langle A \rangle, \langle N_\text{tot} \rangle)$ in the numerator represents the likelihood of observing the areal intensity $\rho$ and area $A$ of a cell for a given values of $\alpha$, $\langle A \rangle$, and $\langle N_\text{tot} \rangle$. The second term in the numerator $g(\alpha,\langle A \rangle, \langle N_\text{tot} \rangle)$ describes the prior knowledge we have regarding the possible values of the the parameters knowing nothing about the measured data. The denominator, $f(\rho, A)$ captures the probability of observing the data knowing nothing about the parameter values. This term, in our case, serves simply as a normalization constant and will be neglected for the remainder of this section. 

To determine the appropriate functional form for the likelihood and prior, we must make some assumptions about the processes that generate them. As the cultures in these experiments are asynchronous in their growth, we should be able to observe a relatively uniform distribution of growth phases ranging from newborn cells to those nearly completed with septation. As there are many independent processes that regulate the timing of cell division and cell growth, such as DNA replication and peptidoglycan synthesis, it is reasonable to assume that for a given culture the distribution of cell size would be normally distributed with a mean of $\langle A \rangle$ and a variance $\sigma_{\langle A \rangle}$. Mathematically, we can write this as
$$
f(A\,\vert\,\langle A \rangle, \sigma_{\langle A \rangle}) \propto {1 \over \sigma_{\langle A \rangle}}\exp\left[-{(A - \langle A \rangle)^2 \over 2\sigma_{\langle A \rangle}^2}\right],
$${#Eq:area_likelihood}
where the proportionality results from dropping normalization constants for notational simplicity.

The areal intensity $\rho$ is intrinsically dependent on the cell area $A$. However, the myriad processes leading the detected fluorescence, such as translation and proper protein folding, are largely independent processes, allowing us to assume a normal distribution for $\rho$ as well with a mean $\langle \rho \rangle$ and a variance $\sigma_\rho$. To compute the average intensity for the population $\langle \rho \rangle$ given $\langle A \rangle$ and $\langle N_\text{tot} \rangle$, we can use [@Eq:area_conversion] to say
$$
\rho =  {\alpha\langle N_\text{tot} \rangle \over \langle A \rangle},
$${#Eq:avg_rho}
allowing us to write the likelihood as
$$
f(\rho\,\vert\,\alpha,\langle A \rangle,\langle N_\text{tot} \rangle, \sigma_\rho) \propto {1 \over \sigma_\rho}\exp\left[-{\left(\rho - {\alpha \langle N_\text{tot} \rangle\over \langle A \rangle}\right)^2 \over 2 \sigma_\rho^2}\right].
$${#Eq:rho_likelihood}

With these two likelihoods in hand, we are tasked with determining the appropriate priors. As we have assumed normal distributions for the likelihoods of $\langle A \rangle$ and $\rho$, we have included two additional parameters, $\sigma_{\langle A \rangle}$ and $\sigma_\rho$, each requiring their own prior probability distribution. It is common practice to assume maximum ignorance for these variances and use a Jeffreys prior [@sivia2006],
$$
g(\sigma_{\langle A \rangle}, \sigma_\rho) = {1 \over \sigma_{\langle A \rangle}\sigma_\rho}.
$${#Eq:jeffreys}

The next obvious prior to consider is for the average channel copy number $\langle N_\text{tot} \rangle$, which comes from Bialecka-Fornal et al. 2012. In this work, they report a mean $\mu_N$  and variance $\sigma_N$, allowing us to assume a normal distribution for the prior,
$$
g(\langle N_\text{tot}\rangle\,\vert\, \mu_N,\sigma_N) \propto {1 \over \sigma_N}\exp\left[-{(\langle N_\text{tot} \rangle - \mu_N)^2 \over 2 \sigma_N^2}\right].
$${#Eq:informative_prior}

For $\alpha$ and $\langle A \rangle$, we have some knowledge of what these parameters can and cannot be. For example, we know that neither of these parameters can be negative. As we have been careful to not overexpose the microscopy images, we can say that the maximum value of $\alpha$ would be the bit-depth of our camera. Similarly, it is impossible to segment a single cell with an area larger than our camera's field of view (although there are biological limitations to size below this extreme). To remain maximally uninformative, we can assume that the parameter values are uniformly distributed between these bounds, allowing us to state 
$$
g(\alpha) = \begin{cases} {1 \over \alpha_\text{max} - \alpha_\text{min}} & \alpha_\text{min} \leq \alpha \leq \alpha_\text{max} \\
0 & \text{otherwise}
\end{cases},
$${#Eq:alpha_uniform_prior}
for \alpha and 
$$
g(\langle A \rangle) = \begin{cases} {1 \over \langle A \rangle_\text{max} - \langle A \rangle_\text{min}} & \langle A \rangle_\text{min} \leq \langle A \rangle \leq \langle A \rangle_\text{max}\\
0 & \text{otherwise}
\end{cases}
$${#eq:area_uniform_prior}
for $\langle A \rangle$.

Piecing [@Eq:area_likelihood] through [@Eq:area_uniform_prior] together generates a complete posterior probability distribution for the parameters given a single cell measurement. This can be generalized to a set of $k$ single cell measurements as
$$
\begin{aligned}
g(\alpha,\langle A \rangle, \langle N_\text{tot} \rangle, \sigma_\rho, \sigma_{\langle A \rangle}\,\vert\, [\rho, A], \mu_N, \sigma_N) &\propto {1 \over (\alpha_\text{max} - \alpha_\text{min})(\langle A \rangle_\text{max} - \langle A \rangle_\text{min})\sigma_{\langle A \rangle}\sigma_\rho}{1 \over (\sigma_\rho\sigma_{\langle A \rangle})^k}\,\times\\
&{1 \over \sigma_N}\exp\left[- {(\langle N_\text{tot}\rangle - \mu_N)^2 \over 2\sigma_N^2}\right]
\prod\limits_i^k\exp\left[-{(A^{(i)} - \langle A \rangle)^2 \over 2\sigma_{\langle A \rangle}^2} - {\left(\rho^{(i)} - {\alpha \langle N_\text{tot}\rangle \over \langle A \rangle}\right)^2 \over 2\sigma_\rho^2}\right] \end{aligned},
$${#Eq:single_rep_post}
where $[\rho, A]$ represents the set of $k$ single-cell measurements.

As small variations in the day-to-day details of cell growth and sample preparation can alter the final channel count of the standard candle, it is imperative to perform more than a single biological replicate. However, properly propagating the error is non trivial. One option would be to pool together all measurements of $n$ biological replicates and evaluate the posterior given in [@Eq:single_rep_post]. However, this by definition assumes that there is no difference between replicates. Another option would be to perform this analysis on each biological replicate individually and then compute a mean and standard deviation of the resulting most-likely parameter estimates. While this is a better approach than simply pooling all data together, it suffers a bias from giving each replicate equal weight, skewing the estimate of the most-likely parameter value if one replicate is markedly brighter or dimmer than the others. Given this type of data and a limited number of biological replicates ($n = 6$ in this work), we chose to extend the Bayesian analysis presented in this section to model the posterior probability distribution for $\alpha$ and $\langle A \rangle$ as a hierarchical process.

### A hierarchical model for estimating $\alpha$

In the previous section, we assumed maximally uninformative priors for the most-likely parameter values for $\alpha$ and $\langle N_\text{tot} \rangle$.

-----

All of these measurements are made via microscopy. As is the nature of the experiments, cells that are about to divide or small clusters of cells could be counted a single cell with a very large number of channels. We can correct for this effect by calculating the average areal intensity and multiplying by the average cellular area, which we can restrict to reasonable bounds in the analysis (e.g. no cells larger than 10 $\mu$m$^2$ in projected 2D area). Using ths correction, [@eq:avg_ian] can be written  as
$$
\langle I_\text{tot} \rangle = \langle I_{\mu m^2} \rangle \langle A \rangle = \alpha \langle N_\text{tot} \rangle,
$${#eq:ian_area}

where $\langle A \rangle$ is the average cell area of the population. Rewriting [@eq:ian_area]  to solve for $\alpha$ yields
$$
\alpha = {\langle I_{\mu m^2} \rangle \langle A \rangle \over \langle N_\text{tot} \rangle}.
$${#eq:alpha_definition_avgs}

With this model in hand, we can construct a Bayesian hierarchical model to estimate the most likely value of $\alpha$ from a set of several biological replicates. 

### A hierarchical model for $\alpha$. ###

![**A hierarchical model for estimation of the calibration factor and average cell size.** ](../figs/hierarchical_model.png)
#### A single replicate
We can begin by thinking of a single biological replicate. Suppose this data set contains many measurements of single cells and includes the average areal intensity of each cell and the cell area (both in the same units of area, such as pixels or $\mu$m$^2$). Using Bayes' theorem, and our knowledge of the mean expression level (from HJ and Maja), we can write an expression for the probability distribution of $\alpha$  as

$$
g(\alpha, \langle I_{A} \rangle, \langle A\rangle \,\vert\, D \langle N_\text{tot} \rangle) = {f(D\,\vert\, \alpha, \langle I_A \rangle, \langle A \rangle, \langle N_\text{tot}\rangle)g(\alpha, \langle I_A \rangle, \langle A \rangle)f(\langle N_\text{tot}\rangle)\over f(\langle I_A \rangle , \langle A \rangle, \langle N_\text{tot} \rangle )}.
$$ {#eq:bayes}

As the average area is independent of the average intensity, the likelihood can be rewritten as
$$
f(\langle I_A \rangle, \langle A \rangle) = f(\langle I_A \rangle\, \vert \, \langle N_\text{tot} \rangle, \alpha, \langle A \rangle) f(\langle A \rangle)f(\langle N_\text{tot} \rangle).
$${#eq:area_likelihood}
However, we must be careful to include the single cell area in our definition of the likelihood for $I_A$, which we will get to shortly.

For these likelihoods, we can make some assumptions about their distributions. As the cell area and the intensity are both the result of many independent cellular processes, we can assume that they are normally distributed. We can therefore write the likelihood for the area,
$$
f(\langle A \rangle\, \vert \sigma_A, D) = {1 \over (2\sigma_A^2\pi)^{k/2}}\prod\limits_{j=1}^k\exp\left[-{(A_j - \langle A \rangle)^2 \over 2\sigma_A^2 }\right],
$${#eq:area_like}

where $D$ represents the set of $k$ experimental measurements. Similarly, we can use [@eq:ian_area] to write the likelihood for the average areal intensity as
$$
f(\langle I_A \rangle\, \vert \, \sigma_{I_A}, \alpha, \langle N_\text{tot} \rangle, \langle A \rangle ,D) = {1 \over (2\sigma_{I_A}^2\pi)^{k/2}}\prod\limits_{j=1}^k\exp\left[- {\left({I_j \over A_j} - {\alpha \langle N_\text{tot}\rangle \over \langle A \rangle }\right)^2  \over 2\sigma_{I_A}^2}\right],
$${#eq:intensity_like}

where we have used [@eq:ian_area] to calculate $\langle I_A \rangle$ of the population.

In our assumption of normal distributions for the areal intensity and area, we have included two new parameters which capture the error in our measurement $\sigma_A$ and $\sigma_{I_A}$. The complete posterior distribution for  a single biological replicate is therefore

$$
\begin{aligned}
g(\alpha, \langle A \rangle, \langle N_\text{tot}\rangle, \sigma_A, \sigma_{I_A} \vert  D) =& {1 \over (2\pi \sigma_A \sigma_{I_A})^k}\prod\limits_{j=1}^k\exp\left[-{(A_j - \langle A \rangle)^2 \over 2\sigma_A^2} - {\left({I_j \over A_j}- {\alpha \langle N_\text{tot} \rangle \over \langle A \rangle}\right)^2 \over 2 \sigma_{I_A}^2 } \right]\,\times\, \\
&f(\langle N_\text{tot}\rangle)g(\alpha,\sigma_A, \sigma_{I_A}),
\end{aligned}
$$
{#eq:single_rep_posterior}

With the likelihood completely described, we are left with fiulling in the priors. The most obvious prior to define is for $\langle N_\text{tot}\rangle$. From the literature, we know the mean value $\mu_{N}$  and a variance $\sigma_{N}$ for this particular strain and growth condition. Using these values, we can specify the prior as

$$
g(\langle N_\text{tot} \rangle\, \vert\, \mu_N, \sigma_N) = {1 \over \sqrt{2\pi\sigma_N^2}} \exp\left[-{\left(\langle N_\text{tot}\rangle - \mu_N\right)^2}\over 2\sigma_N^2\right].
$$

We can, of course, specify priors for inferred parameters to be maximally uninformative. However, as we have multiple replicates and we are assuming each replicate is comparable to another, we can be a bit more informative with these priors.

#### Multiple replicates

We expect each biological replicate to behave the same, despite the small experimental differences between them. This means that the model shown in [@eq:single_rep_posterior] should be the same for each replicate. However, we know that the parameters $\alpha, \sigma_A,$ and $\sigma_{I_A}$ should be very similar between replicates. In fact, we can say that the values for each individual replicate are drawn from the same distribution.

The functional form of this distribution is not obvious. However, we know that the value is the result of many independent processes, meaning it is likely gaussian. 

With this assumption in place, we can write the prior distribution of $\langle A \rangle$ for a set of $n$ replicates as

$$
g(\langle A \rangle\, \vert \, \tilde{\langle A \rangle}, \tilde{\sigma_A}) = {1 \over (2\pi\tilde{\sigma_A}^2)^{n/2}} \prod\limits_{i=1}^n \exp\left[-{\left(\langle A_i\rangle - \tilde{\langle A \rangle}\right)^2 \over 2 \tilde{\sigma_A}^2}\right],
$${#eq:area_prior}

and for $\alpha$,

$$
g(\alpha\, \vert \tilde{\alpha}, \tilde{\sigma_\alpha}) = {1 \over (2\pi\tilde{\sigma_\alpha}^2)^{n/2}} \prod \limits_{i=1}^{n}\exp\left[-{\left(\alpha_i - \tilde{\alpha}\right)^2 \over 2\tilde{\sigma_\alpha}^2}\right].
$${#eq:alpha_prior}.

While we have taken care of the priors listed in [@eq:single_rep_posterior], we have introduced four new parameters, which each need their own prior! For these,we can be maximally uninformative. For $\tilde{\langle A \rangle}$, we can impose a uniform prior over a reasonable range for the size of *E. coli*,

$$
g(\tilde{\langle A \rangle}) = \begin{cases}\left(\tilde{\langle A \rangle}_\text{max} - \tilde{\langle A \rangle}_\text{min}\right)^{-1} & \tilde{\langle A \rangle}_\text{min} \leq \tilde{\langle A \rangle} \leq \tilde{\langle A \rangle}_\text{max}\\
0 & \text{otherwise}
\end{cases}.
$${#eq:area_hyperprior}


We know that the value of $\alpha$ is limited by the bit depth of our camera, though not necessarily an integer. We can assign a uniform prior for $\alpha$ as

$$
g(\tilde{\alpha}) = \begin{cases} \left(\tilde{\alpha}_\text{max} - \tilde{\alpha}_\text{min}\right)^{-1} & \tilde{\alpha}_\text{min} \leq \tilde{\alpha} \leq \tilde{\alpha}_\text{max}\\
0 & \text{otherwise}
\end{cases}.
$${#eq:alpha_hyperprior}

For $\tilde{\sigma_\alpha}$ and $\tilde{\sigma_A}$, we can use a maximally uninformative Jeffreys priors, 

$$
g(\tilde{\sigma_\alpha}) = {1 \over \tilde{\sigma_\alpha}},
$${#eq:sigma_alpha_hyperprior}

and 

$$
g(\tilde{\sigma_A}) = {1 \over \tilde{\sigma_A}}.
$${#eq:sigma_area_hyperprior}

Using [@eq:bayes] through [@eq:sigma_area_hyperprior], we can formulate the complete posterior probability distribution,
$$
\begin{aligned}
&g(\alpha, \tilde{\alpha}, \langle A \rangle, \tilde{\langle A \rangle}, \langle N_\text{tot} \rangle, \tilde{\sigma_\alpha}, \sigma_A, \tilde{\sigma_A}, \sigma_{I_A}\,\vert \, \sigma_N, \mu_N, D) = \\
& {1 \over \tilde{\langle A \rangle}_\text{max} - \tilde{\langle A \rangle}_\text{min}}{1 \over \tilde{\alpha}_\text{max} - \tilde{\alpha}_\text{min}}{1 \over \left(\tilde{\sigma_\alpha}\tilde{\sigma_A}\right)^{n+2}}\prod\limits_{i=1}^n{1 \over \left(\sigma_{A_i}^2\sigma_{{I_A}_i}^2\sigma_N^2\right)^{k_i / 2}}\times\\
&\prod\limits_{j=1}^{k_i}\exp\left[-{\left(A_{i,j}-\langle A \rangle_i\right)^2 \over 2\sigma_{A_i}^2} - {\left({I_{A_{i,j}} - {\alpha_i \langle N_\text{tot}\rangle_{i} \over \langle A\rangle_i}}\right)^2 \over 2\sigma_{{I_A}_i}^2}-{\left(\alpha_i - \tilde{\alpha}\right)^2 \over 2\tilde{\sigma_{\alpha}}^2}-{\left(\langle A \rangle_i - \tilde{\langle A \rangle}\right)^2 \over 2\tilde{\sigma_A}^2}-{\left(\langle N_\text{tot}\rangle_i - \mu_N\right)^2 \over 2\sigma_N^2}\right],
\end{aligned}
$${#eq:posterior}

where the normalization constants have been dropped for notational clarity.



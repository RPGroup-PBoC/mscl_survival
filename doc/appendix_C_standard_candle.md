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

In the previous section, we assumed maximally uninformative priors for the most-likely parameter values for $\alpha$ and $\langle N_\text{tot} \rangle$. While this is a fair approach to take for these parameters, we are not completely ignorant for how these values are distributed across biological replicates. A major assumption of our model is that the most-likely value of $\alpha$ and $\langle A \rangle$ for each biological replicate are comparable, so long as the experimental error between them is accounted for. In other words, we are assuming that the most-likely value for each replicate is drawn from a global distribution. While each replicate may have a unique value, they are all related to one another. With enough replicates, we would be able to adequately sample this distribution by treating each replicate independently. Unfortunately, proper sampling of this distribution requires an extensive amount of experimental work, making inferential approaches more attractive. 

In this approach, often called a multi-level or hierarchical model, is schematized in [@Fig:hierarchical_model]. Here, we use  an informative prior for $\alpha$ and $\langle A \rangle$ for each biological replicate. This informative prior probability distribution can be described by the summary statistics, often called hyper-parameters, which are then used as the most-likely values for the parameters. As this approach allows us to get a picture of the probability distribution of the hyper-parameters, we are able to report a point estimate and error that captures all known sources of the statistical error.

![**Schematic of hierarchical model structure**. The hyper-parameter probability distributions (top panel) are used as an informative prior for the most-likely parameter values for each biological replicate (middle panel). The single-cell measurements of cell area and areal intensity (bottom panel) are used as data in the evaluation of the likelihood](../figs/hierarchical_model.png){#Fig:hierarchical_model}

The choice for the functional form for the informative prior is often not obvious and can require other experimental approaches or back-of-the-envelope estimates to approximate. Each experiment in this work was carefully constructed to minimize the day-to-day variation. This involved adhering to well-controlled growth temperatures and media composition, harvesting of cells at comparable optical densities, and ensuring identical imaging parameters. As the experimental variation is minimized, we can use our knowledge of the underlying biological processes to guess at the approximate functional form. For similar reasons presented in the previous section, cell size is controlled by a myriad of independent processes. As each replicate is independent of another, it is reasonable to assume a normal distribution for the average cell area for each replicate. This normal distribution is described by a mean $\tilde{\langle A \rangle}$ and variance $\tilde{\sigma}_{\langle A \rangle}$. Therefore, the prior for $\langle A \rangle$ for $n$ biological replicates can be written as
$$
g(\langle A \rangle\, \vert\, \tilde{\langle A \rangle}, \tilde{\sigma}_{\langle A \rangle}) \propto {1 \over \tilde{\sigma}_{\langle A \rangle}^n}\prod\limits_{j=1}^{n}\exp\left[-{(\langle A \rangle_j - \tilde{\langle A \rangle})^2 \over 2 \tilde{\sigma}_{\langle A \rangle}^2}\right].
$${#Eq:hyperprior_area}
In a similar manner, we can assume that the calibration factor for each replicate is normally distributed with a mean $\tilde{\alpha}$ and variance $\tilde{\sigma}_\alpha$,
$$
g(\alpha\,\vert\,\tilde{\alpha}, \tilde{\sigma}_\alpha) \propto {1 \over \tilde{\sigma}_\alpha^n}\prod\limits_{j=1}^n \exp\left[-{(\alpha_j - \tilde{\alpha})^2 \over 2\tilde{\sigma}_\alpha^2}\right].
$$

With the inclusion of two more normally distributed parameters, we have introduced four new parameters, each of which needing their own prior. However, our knowledge of the reasonable values for the hyper-parameters has not changed from those described for a single replicate. We can therefore use the same maximally uninformative Jeffreys priors given in [@Eq:Jeffreys] for the variances and the uniform distributions given in [@Eq:alpha_uniform_prior] and [@Eq:area_uniform_prior] for the means.

Stitching all of this work together generates the full posterior probability distribution for the best-estimate of $\alpha$ and $\langle A \rangle$ shown in [@Eq:avg_ian] given $n$ replicates of $k$ single cell measurements, 
$$
\begin{aligned}
g(\tilde{\alpha}, \tilde{\sigma}_\alpha, \tilde{\langle A \rangle}, \tilde{\sigma}_{\langle A \rangle}, &\{\langle N_\text{tot} \rangle, \langle A \rangle, \alpha, \sigma_\rho\}\,\vert\, [\rho, A], \mu_N, \sigma_N) \propto\\
&{1 \over (\tilde{\alpha}_\text{max} - \tilde{\alpha}_\text{min})(\tilde{\langle A \rangle}_\text{max} - \tilde{\langle A \rangle}_\text{min})\sigma_N^n(\tilde{\sigma}_\alpha\tilde{\sigma}_{\langle A \rangle})^{n + 1}}\,\times\,\\
&\prod\limits_{j=1}^n\exp\left[-{(\langle N \rangle_j^{(i)} - \mu_N)^2 \over 2\sigma_N^2} - {(\alpha_j - \tilde{\alpha})^2 \over 2\tilde{\sigma}_\alpha^2} - {(\langle A \rangle_j - \tilde{\langle A \rangle})^2 \over 2\tilde{\sigma}_{\langle A \rangle}^2}\right]\,\times\,\\ 
&{1 \over (\sigma_{\rho_j}\sigma_{\langle A \rangle_j})^{k + 1}}\prod\limits_{i=1}^{k_j}\exp\left[-{(A_j^{(i)} - \langle A \rangle_j)^2 \over 2\sigma^{(i)2}_{\langle A \rangle_j}} - {\left(\rho^{(i)}_j - {\alpha_j \langle N_\text{tot}\rangle_j \over \langle A \rangle_j}\right)\over 2\sigma_{\rho_j}^{(i)2}}\right].
\end{aligned}
$${#Eq:cal_factor_posterior}

While [@Eq:cal_factor_posterior] is not analytically solvable, it can be easily sampled using Markov chain Monte Carlo (MCMC). The results of the MCMC sampling for $\tilde{\alpha}$ and $\tilde{\langle A \rangle}$ can be seen in [@Fig:alpha_area_sampling]. From this approach, we found the most-likely parameter values of $3300^{+700}_{-700}$ a.u. per MscL channel and $5.4^{+0.4}_{-0.5}$ $\mu$m$^2$ for $\tilde{\alpha}$ and $\tilde{\langle A \rangle}$, respectively. Here, we've reported the median value of the posterior distribution for each parameter with the upper and lower bound of the 95\% credible region as superscript and subscript, respectively. 

![**Posterior distributions for hyper-parameters and replicate parameters.** (A) The posterior probability distribution for $\tilde{\alpha}$ and $\tilde{\langle A \rangle}$. Probability increases from light to dark red. The replicate parameter (blue) and hyper-parameter (red) marginalized posterior probability distributions for $\alpha$ (B) and $\langle A \rangle$ (C).](../figs/FigSX1.png){#Fig:posterior_samples}
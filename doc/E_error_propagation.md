# ESTIMATING STATISTICAL ERROR FOR SURVIVAL PROBABILITY

In Fig. 5 of the main text, we presented the survival probability as a function of effective channel number using logistic regression in addition to the probabilities calculated by binning the data in two different ways. For the latter calculation, we included an estimate of the statistical error associated with the survival probability. As we pooled all experimental data together, it did not sense to separate them by biological replicate and compute the standard error of the mean. Rather, we estimated the error given the number of cells in that particular bin. Here, we derive an expression for this calculation.

Given a collection of single cell measurements, we are interested in calculating the survival probability of the population. Using Bayes' theorem, we can write an expression for the probability distribution for the probability of survival $p_s$ as 

$$
P(p_s\vert n_s, N) = {P(n_s\vert p_s, N) P(p_s) \over P(n \vert N)},
$${#eq:binned_bayes}
where $n_s$ is the number of survivors and $N$ is the total number of cells in the population. The prior probability distribution $P(p_s)$ represents our knowledge of the survival probability being ignorant of the data. As we have no reason to assume one cell is more likely to survive than another, we can assume that this probability is uniformly distributed between $0$  and $1$,

$$
P(p_s) = \begin{cases}
1\, & \, 0 \leq p_s \leq 1 \\
0\, & \text{otherwise}.
\end{cases}
$${#eq:ps_prior}
The likelihood of the data given a total number of cells $N$ and a probability of survival $p_s$ is Binomially distributed and can be expressed as

$$
P(n_s\vert p_s, N) = {(n_s + 1)! \over n_s!(N - n_s)!}p_s^{n_s}(1 - p)^{N - n_s}
$${#eq:ps_likelihood}
where the numerator of the binomial coefficient includes the possibility of zero cells surviving in the population. As we are neglecting the denominator of {@eq:binned_bayes} and assuming a constant uniform prior for $p_s$, {@eq:ps_likelihood} is equivalent to the posterior distribution $P(p_s \vert n_s, N)$.

To compute the most probable survival probability $p_s^*$, we can identify the value of $p_s$ at which the derivative of {@eq:ps_likelihood} vanishes. This is equivalent to finding the point at which the derivative of the log posterior goes to zero,

$$
p_s^* = {d\log P \over dp_s} = {n_s \over p_s} - {n - n_t\over 1 - p_s} = 0.
$${#eq:pstar_def}
Solving for $p_s$ yields the most probable value, 

$$
p_s^* = {n_s \over N}.
$${#eq:most_prob_ps}

Assuming that $N >> p$,  we can approximate {@eq:ps_likelihood} as a Gaussian distribution. Making this approximation, we can compute the variance using the fact

$$
\sigma^2 = -\left({d^2 \log P \over dp_s^2}\biggr\vert_{p=p_s^*}\right)^{-1}.
$${#eq:variance_def}
The second derivative can be computed as

$$
{d^2 \log P \over dp_s^2} = -{n_sp_s^2 - N(1 - 2p_s) \over p_s^2 (1 - p_s^2)}.
$${#eq:second_derivative}
Substituting in {@eq:most_prob_ps} and taking the negative reciprocal gives the desired result, 

$$
\sigma^2 = {n_s(N - n_s) \over N^3}.
$${#eq:variance}
Given pooled data, we can therefore report the most likely value for the survival probability and an estimate of the error as $p_s = p_s^* \pm \sigma.$




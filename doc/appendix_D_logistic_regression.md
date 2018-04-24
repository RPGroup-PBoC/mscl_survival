# 

## LOGISTIC REGRESSION ##

In this work, we were interested in computing the survival probability under a large hypo-osmotic shock as a function of MscL channel number. As the channel copy number distributions for each Shine-Dalgarno sequence mutant were broad and extensively overlapping, we chose to use a method of calculating survival probability through logistic regresion - a method that requires no binning of the data providing the least biased estimate of survival probability. Logistic regression is a technique that has been used in medical statistics since the late 1950's to describe diverse phenomena such as dose response curves, criminal recidivism, and survival probabilities for patients after treatment. It has also found much use in machine learning to tune a binary or categorical response given a continuous predictor signal. (**Find citations -- I know they are there**)

In this section, we derive the functional form of the logistic regression cost function and provide an intuitive interpretation for the seemingly arbitrary parameters. We conclude with a description of our estimation methods for the most-likely values of these parameters using Markov chain Monte Carlo.

### A Bayesian interpretation of $\beta_0$ and $\beta_1$ ###

The central challenge of this work is to estimate the probability of survival $p_s$ given only a measure of the total number of MscL channels in that cell. In other words, for a given measurement of $N_c$ channels, we want to know likelihood that cell would survive a given osmotic shock. The probability of observing a survival event $s$ given a measurement of $N_c$ channels can be stated using Bayes' theorem as
$$
g(s\,\vert\, N_c) = {f(N_c\,\vert\, s)g(s) \over f(N_c)}.
$${#eq:survival_bayes}
where $g$ and $f$ represent probability density functions over parameters and data respectively. The posterior distribution $g(s\,\vert\, N_c)$ is the quantity of interest and implicitly related to the probability of survival. The likelihood $g(N_c\,\vert\, s)$ tells us the probability of observing $N_c$ channels in this cell given that it survives. The quantity $g(s)$ captures all *a priori* knowledge we have regarding the probability of this cell surviving and the denominator $f(N_c)$ tells us the converse -- the probability of observing $N_c$ cells irrespective of the survival outcome.

Proper calculation of [@Eq:survival_bayes] requires that we have knowledge of $f(N_c)$, which is difficult to estimate. While we are able to give appropriate bounds on this term, such as a requirement of positivity and some knowledge of the maximum membrane packing density, it is not so obvious to determine the distribution between these bounds. Given this difficulty, it's easier to compute the odds of survival $\mathcal{O}(s\,\vert\, N_c)$, the probability of survival $s$ relative to death $d$,
$$
\mathcal{O}(s\,\vert\, N_c) = {g(s\,\vert\,N_c) \over g(d\,\vert\, N_c)} = {f(N_c\,\vert\, s)g(s) \over f(N_c\,\vert\,d)g(d)},
$${#Eq:odds_definition}
where the nuisance parameter $f(N_c)$ has not been cancelled out of our calculation.  

XXX log transform XX. Computing the log-odds gives 

$$
\log \mathcal{O}(s\,\vert\,N_c)= \log {g(s) \over g(d)} + \log {f(N_c \,\vert\, s )\over f(N_c\,\vert\, d)}.
$${#Eq:log_odds}



-----

### Deriving the functional form of logistic regression 

The functional form of the logistic regression model given by Eq. 1 in the main text includes two dimensionless and seemingly arbitrary coefficients, $\beta_0$ and $\beta_1$. Interpreting the probabilistic or physical meaning of these parameters is often an afterthought. However, there is a direct connection between Bayes' theorem (CITE DOWNEY) and these parameters which we will examine here in depth.

Suppose we are examining a single cell with an effective channel copy number of $n$ and we are interested in computing this its probability of survival. We can write an expression for this probability by using Bayes' theorem as

$$ P(s\vert n) = {P(n \vert s)P(s) \over P(n)},
$${#eq:bayes}
where $P(n\vert s)$ is the likelihood of having $n$ channels given an observation of survival, $P(s)$ is the probability of survival completely independent of the channel number measurement, and $P(n)$ is the probability of having $n$ channels irrespective of the cell fate. Calculating this probability requires that we have some knowledge of $P(n)$, which is a difficult quantity to estimate. We may know the bounds of the physically accessible values of $n$, but we have no knowledge of the distribution. For this reason, it's easier to consider the odds of survival with a given number of channels $n$ relative to the odds of survival at another channel copy number. For example, let's examine the odds of survival for a cell with $n = 0$ and another with $n = 1$. Using {@eq:bayes}, we can write the odds of survival for each copy number as

$$
\mathcal{O}(s \vert n=0 ) = {P(s \vert n=0) \over P(d \vert n=0)} = {P(n=0 \vert s)P(s) \over P(n=0 \vert d)P(d)},
$${#eq:explicit_odds_zero}
and

$$
\mathcal{O}(s \vert n=1) = {P(s \vert n=1) \over P(d \vert n=1)} = {P(n=1 \vert s)P(s) \over P(n=1 \vert d)P(d)}.
$${#eq:explicit_odds_one}
We see for each copy number, the ratio $P(s) / P(d)$ is always present. By computing the log-odds, we can separate two useful quantities,

$$
\log\mathcal{O}(s\vert n) = \log{P(s)\over P(d)} + \log{P([n = 0, n=1] \vert s)\over P([n=0, n=1] \vert d)},
$${#eq:log_odds}. 
The first term is the log-odds of survival independent of the number of channels. The second term 
is the log likelihood of observing $n$ channels given that the cell survives relative to when the cell dies. Rather than considering $n=0$ and $n=1$ separately, we can write the log likelihood ratio in one concise statement as

$$
 \log {P([n=0, n=1] \vert s)\over P([n=0, n=1]\vert d)} = \log {P(n=0\vert s) \over P(n=0\vert d)} + N_c\left[\log {P(n=1\vert s) \over P(n=1 \vert d)} - \log {P(n=0 \vert s) \over P(n=0 \vert d)}\right] = \log {P(N_c \vert s) \over P(N_c \vert d)},
$${#eq:single_llr_expression}
where $N_c \in [0, 1]$. The bracketed quantity is mathematically equivalent to the log-odds ratio of the two channel copy numbers $\mathcal{OR}(s\vert N_c)$.

We can stitch together these expressions to arrive at a single expression for the log-odds of survival given a channel copy number $n$. By combining {@eq:log_odds} with {@eq:single_llr_expression} we have

$$
\log\mathcal{O}(s\vert N_c) = \log \mathcal{O}(s)  + \log {P(n=0 \vert s) \over P(n=0 \vert s)} + N\left[\log {P(n=1 \vert s) \over P(n=1 \vert d)} - \log{P(n=0 \vert s) \over P(n=0 \vert d)}\right].
$${#eq:log_odds_step_one}
Finally, we can substitute for the second term in {@eq:log_odds_step_one} with {@eq:explicit_odds_zero} to arrive at

$$
\log \mathcal{O}(s \vert n=1) = \log\mathcal{O}(s\vert n=0) + N_c\mathcal{OR}(s\vert N).
$${#eq:bayes_logit}

In of the main text, Eq. 1 stated the model for logistic regression. Ignoring the use of $\log N_c$ as a predictor, this model is given as
$$
\log{p_s \over 1 - p_s} = \beta_0 + \beta_1 N_c,
$${#eq:main_text_logit}
which is of the same form as {@eq:bayes_logit}. By comparing these two expressions, we can make sense of the seemingly arbitrary coefficients. The intercept $\beta_0$ is the log-odds of survival with zero channels where as $\beta_1$ is the log-odds of survival with one channel relative to the log-odds of survival with zero. This can easily be expanded to cover any number of mechanosensitive channels, which is performed for the analysis presented in the main text.

## BAYESIAN PARAMETER ESTIMATION
  In this work, we chose to estimate the most-likely parameter values for
 $\beta_0$ and $\beta_1$ using a Bayesian definition of probability. We can
 compute the posterior probability distribution for these parameters given a collection of data $D$ which is given by Bayes' theorem,

$$
  P(\beta_0, \beta_1\vert D) = {P(D\vert\beta_0, \beta_1)P(\beta_0,\beta_1) \over P(D)}.
$${#eq:bayes_thm}
The data $D$ is composed of all single cell measurements of effective channel copy number, the statistical error in the channel copy number, and their survival ('True' or 'False'). The quantity $P(D\vert\beta_0, \beta_1)$ represents the likelihood of observing the data given the parameters. All prior knowledge of the coefficients for $\beta_0$ and $\beta_1$ independent of the observed data is captured by $P(\beta_0, \beta_1)$. Finally, the denominator $P(D)$ is the likelihood marginalized over the parameters. In the context of this work, this quantity serves simply as a normalization constant and an be ignored. 

To formulate the likelihood, we can abstract survival or death as a Bernoulli process with a probability of success $p_s$,

$$
f(s \vert p_s) = p_s^F(1 - p)^{1 - F},
$${#eq:bernoulli}
where $F$ is the fate of the cell with $1$ being survival and $0$ being death. Each measurement in the data set is treated independently, as is implied by {@eq:bernoulli} and Eq. 1 of the main text in which the probability is dependent on the log of the channel copy number $N_c$. The likelihood given in {@eq:bayes_thm} can be written for the entire data set as

$$
P(D\vert \beta_0, \beta_1) = \prod_{i=1}^k\left({1 \over e^{-\theta_i}}\right)^{F_i}\left(1 - {1 \over 1 + e^{-\theta_i}}\right)^{1 - F_i},
$${#eq:likelihood}
where $k$ is the total number of samples in the dataset, $s_i$ is the fate of the $i^{th}$ cel, and $\theta$ is defined as

$$
\theta = \log{p_s \over 1 - p_s} = \beta_0 + \beta_1 \log N_c.
$${#eq:theta_logit}

Each measurement of the channel copy number $N_c$ has an associated statistical error $\sigma_c$. To incorporate this uncertainty in our statistical model, we must include yet another prior for the "true" channel copy number for each measurement. As the copy number is dependent on many independent processes (e.g. translation, folding, membrane insertion, etc), it is reasonable to assume that the copy number is normally distributed. WE can therefore write a prior each single cell as

$$
P(n\vert N_c, \sigma_c) = {1 \over \sqrt{2\pi\sigma_c^2}}\exp\left[-{(n - N_c)^2 \over 2\sigma_c^2}\right],
$${#eq:nc_prior}.

An advantage of logistic regression is that there are no bounds on the values that $\beta_0$ and $\beta_1$ can take. This permits us to assume a flat distribution for the prior of both coefficients, allowing us to drop them from our expression. Putting all of these pieces together gives us the final posterior probability distribution,

$$
P(\beta_0,\beta_1, [n]\vert N_c, \sigma_c) = {1 \over 2\pi^{k/2}}\prod_{i=1}^k {1 \over \sigma_{c,i}} \exp\left[-{(n_i - N_{c,i})^2 \over 2\sigma_{c,i}^2}\right]\left({1 \over 1 + e^{-\beta_0 - \beta_1 \log n_i}}\right)^{F_i} \left(1 - {1 \over 1 + e^{-\beta_0 - \beta_1 \log n_i}}\right)^{1 - F_i}.
$${#eq:posterior}

Exact marginalization is not practical for this posterior. To find the most-likely parameter values, we used Markov chain Monte Carlo (MCMC) to sample directly form this distribution. To this end, we used the high efficiency "No U-turns" sampler (NUTS) packaged with the probabilistic programming language Stan. The posterior for the coefficients over the slow and fast shock rate experiments can be seen in [@fig:posterior]

![**Posterior probability distributions for logistic regression coefficients.** The density of the MCMC samples for the two coefficients of logistic regression for (A) slow and (B) fast shocks. The mode and 95% credible regions for each parameter are shown as the dot and crosshairs, respectively.](../figs/figSX1.png){#fig:posterior}


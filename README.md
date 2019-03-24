<div align="right">
Advanced Numerical Methods and Data Analysis - FS19-8,780
<br>
University of St. Gallen, 24.03.2019
<br>
<a href="https://github.com/peterlacour/Option-Pricing-using-COS-Method">Github Repository</a>

</div>

-------------



# Option Pricing using the Cosine Fast Fourier Transform Method


**Federico Baldi Lanfranchi** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; federico.baldilanfranchi@student.unisg.ch <br>
**Peter la Cour** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;peter.lacour@student.unisg.ch

## <div id="0"><a href="#0">Overview</a></div>


1. <a href="#2">Introduction</a>
2. <a href="#A2">Description of the COS-FFT Method</a>
	* <a href="#AA1">Carr-Madan FFT</a>
	* <a href="#AA2">COS Method</a>
3. <a href="#B2">Characteristic Functions</a>
    * <a href="#BB1">Black-Scholes Model </a>
	* <a href="#BB2">Heston Model</a>
	* <a href="#BB3">Variance Gamma Model</a>
	* <a href="#BB4">CGMY Model</a>
	* <a href="#BB5">Parameter Values</a>
4. <a href="#C2">Distributional Characteristics of the different Models</a>
5. <a href="#D2">Option Pricing with the COS-FFT Method</a>
6. <a href="#E2">Concluding Remarks</a>
7. <a href="#F2"> References </a>

## <div id="2"> <a href="#0">Introduction  </a> </div>
In the following, we compute prices for one day of European put and call options on the S&P 500 Index through the COS method, pioneered by Fang and Oosterle (2008). For this study, we consider options expiring on the March 20, 2020, whose maturity is close to one year at the time of writing, March 24, 2019. Option price data is retried from Yahoo, through the quantmod R package. We price both put and call options for every strike price quoted, also in cases where one of the two option types is not available for a particular strike. This approach is adopted since the related contingent claims can theoretically be replicated via the put-call parity relationship.

We compare estimated prices under four underlying diffusion assumptions, namely the Black-Scholes, Heston, Variance Gamma and CGMY models. In doing so, we employ calibrated parameters from some of the original papers. Consequently, we cannot make statements about the accuracy of the COS method with respect to other numerical option pricing techniques at this stage, since realised option prices should instead reflect current market parameters. We leave this analysis to successive work, where such parameters will be calibrated through numerical optimization. We thus take the Black-Scholes model, for which an analytical pricing formula is available, as a reference framework. We then comment on whether COS prices capture the intended deviation of other models from the BS setup, despite truncation of the domain and finite precision. We find that, in general, this seems to be the case.



## <div id="A2"> <a href="#0">Description of the COS-FFT Method</a> </div>

The COS mehod (Fang and Oosterlee, 2008) is a novel option pricing method based on a Fourier-cosine expansion of the density function. The algorithm can be applied to european plain vanilla options and to some options with early exercise rights. Cosine expansions offer an alternative to the industry-standard Fast Fourier Transform (FFT) approach developed in Carr and Madan (1999). Both  FFT and the COS method exploit knowledge of the characteristic function of the underlying process, which is often given, even when the density is unknown. However, the COS method can be significantly faster, while also circumventing some of the issues with the FFT approach, which are briefly highlighted in next section. FFT-based methods' for option pricing typically exhibit second order accuracy with computational complexity <img src="https://latex.codecogs.com/gif.latex?\inline&space;O(N&space;\log_2(N))" title="O(N \log_2(N))" /> (Fang and Oosterlee, 2008). The COS method instead can achieve exponential convergence, while keeping computational complexity linear. 



### <div id="BB1"> Carr-Madan FFT </div>

Since the characteristic function is equivalent to the Fourier transform of the density function in the risk-neutral domain the pricing problem can be solved in Fourier domain. FFT then relies on numerical integration through quadrature rules to solve the inverse Fourier integral and tranform the solution back to the time domain. The Fourier transform is taken with respect to the log strike price <img src="https://latex.codecogs.com/gif.latex?\inline&space;k&space;=&space;ln(K)" title="k = ln(K)" />, with <img src="https://latex.codecogs.com/gif.latex?\inline&space;x&space;=&space;ln(S_T)" title="x = ln(S_T)" />. The price for a European call option can then be expressed as:


<img src="https://latex.codecogs.com/gif.latex?C_T(k)&space;=&space;e^{-rT}&space;E^{Q}[(S_T-K)^&plus;]=&space;e^{-rT}&space;\int^\infty_k&space;(e^x&space;-&space;e^k)q(x)" title="C_T(k) = e^{-rT} E^{Q}[(S_T-K)^+]= e^{-rT} \int^\infty_k (e^x - e^k)q(x)" />


Where <img src="https://latex.codecogs.com/gif.latex?\inline&space;q(\cdot)" title="q(\cdot)" /> denotes the risk-neutral density function. When <img src="https://latex.codecogs.com/gif.latex?\inline&space;k" title="k" /> goes to <img src="https://latex.codecogs.com/gif.latex?\inline&space;-\infty" title="-\infty" />, the call price converges to <img src="https://latex.codecogs.com/gif.latex?\inline&space; S_0" title="S_0" />. As a result, the call pricing function is not square integrable, and its Fourier integral is not well defined. Since FFT evaluates the characteristic function at <img src="https://latex.codecogs.com/gif.latex?\inline&space;0" title="0" />, the introduction of a damping parameter <img src="https://latex.codecogs.com/gif.latex?\inline&space;e^{\alpha k}" title="e^{\alpha k}" /> is necessary to in order to shift the singularity away from the origin. Lewis (2001) and Lee (2004) show how this is equivalent to evaluating a contour integral in the complex plane, where the pole is shifted away from the real line on the immaginary axis. 

However, the Carr-Madan approach suffers from some disadvantages. First, quadrature rules employed for numerical integration can be inefficient due to the highly oscillatory nature of integrands. Further, for this same reason the choice of the damping parameter has profound implications in terms of efficiency. Higher values of <img src="https://latex.codecogs.com/gif.latex?\inline&space;\alpha" title="\alpha" /> can make integration less stable on the positive log-strike axis, so that the damping parameter must be chosen carefully. Oscillatory behaviour becomes more pronounced close to expiration, when the option approaches its intrisic value. Treating this issue requires a different adjustment factor based on the hyperbolic sine. 

Another issue with the FFT method is that the number of evaluation arguments <img src="https://latex.codecogs.com/gif.latex?\inline&space;N" title="N" /> must necessarily be a power of 2, which can be undesirable. Further, the grid size for the numerical integration is typically tied to the spacing between strike prices by the Nyquist relation:

<img src="https://latex.codecogs.com/gif.latex?\Delta&space;x&space;\cdot&space;\Delta&space;\omega&space;=&space;\frac{2&space;\pi}{N}" title="\Delta x \cdot \Delta \omega = \frac{2 \pi}{N}" />


where <img src="https://latex.codecogs.com/gif.latex?\inline&space;x" title="x" /> denotes the log price and <img src="https://latex.codecogs.com/gif.latex?\inline&space;w" title="w" /> denotes the argument of the characteristic function. Choosing a fine grid to improve integration accuracy increases the spacing between strike prices at which option prices are computed. Strikes are then pushed deeper in- and out-of-the-money, leading to a lower number of option prices laying in the relevant region.





### <div id="BB1"> COS Method </div>

The COS method does not rely on numerical integration to price contingent claims. Rather, it replaces the entire density with its cosine series expansion. The cosine expansion is defined for functions with support on <img src="https://latex.codecogs.com/gif.latex?\inline&space;[0,\pi]" title="[0,\pi]" />. Through change of variables, it can be generalised to functions with support on a generic finite interval <img src="https://latex.codecogs.com/gif.latex?\inline&space;[a,b]" title="[a,b]" />. 


<img src="https://latex.codecogs.com/gif.latex?\sum_{k=0}^{&plus;\infty}'&space;A_k&space;\cdot&space;\cos&space;\Big(k&space;\pi&space;(x-a)/(b-a)&space;\Big)" title="\sum_{k=0}^{+\infty}' A_k \cdot \cos \Big(k \pi (x-a)/(b-a) \Big)" />

<img src="https://latex.codecogs.com/gif.latex?A_k&space;=&space;\frac{2}{b-a}&space;\int^{b}_{a}&space;f(x)&space;\cos\Big(k&space;\pi&space;\frac{x-a}{b-a}\Big)dx" title="A_k = \frac{2}{b-a} \int^{b}_{a} f(x) \cos\Big(k \pi \frac{x-a}{b-a}\Big)dx" />


Where <img src="https://latex.codecogs.com/gif.latex?\inline&space;\sum&space;'" title="\sum '" /> indicates that the first summand is divided by 2. The density inside the equation for the <img src="https://latex.codecogs.com/gif.latex?\inline&space;A_k" title="A_k" /> terms can be substituted through an approximation of the characteristic function on the <img src="https://latex.codecogs.com/gif.latex?\inline&space;(a,b)" title="(a,b)" /> interval, which introduces a first error component, <i>i.e.</i> <img src="https://latex.codecogs.com/gif.latex?\inline&space;\varphi&space;\approx&space;\varphi_1&space;=&space;\int^b_a&space;e^{iux}&space;\cdot&space;f(x)&space;dx" title="\varphi \approx \varphi_1 = \int^b_a e^{iux} \cdot f(x) dx" />. This gives 

<img src="https://latex.codecogs.com/gif.latex?A_k&space;=&space;\frac{2}{b-a}&space;Re\Big\{\varphi_1\&space;\Big(&space;\frac{k&space;\pi}{b-a}\Big)&space;\cdot&space;\exp\Big(-i&space;\frac{k&space;a&space;\pi}{b-a}&space;\Big)&space;\Big\}" title="A_k = \frac{2}{b-a} Re\Big\{\varphi_1\ \Big( \frac{k \pi}{b-a}\Big) \cdot \exp\Big(-i \frac{k a \pi}{b-a} \Big) \Big\}" />

Which can be approximated by:

<img src="https://latex.codecogs.com/gif.latex?F_k&space;=&space;\frac{2}{b-a}&space;Re\Big\{\varphi&space;\Big(&space;\frac{k&space;\pi}{b-a}\Big)&space;\cdot&space;\exp\Big(-i&space;\frac{k&space;a&space;\pi}{b-a}&space;\Big)&space;\Big\}" title="F_k = \frac{2}{b-a} Re\Big\{\varphi \Big( \frac{k \pi}{b-a}\Big) \cdot \exp\Big(-i \frac{k a \pi}{b-a} \Big) \Big\}" />

We then substitute <img src="https://latex.codecogs.com/gif.latex?\inline&space;F_k" title="F_k" /> for <img src="https://latex.codecogs.com/gif.latex?\inline&space;A_k" title="A_k" />in the infinite series and truncate the sum, which introduces a second source of approximation error. The price of a contingent claim, which is given by an expectation in terms of the risk neutral density, can then be computed. In the case of a call option the equation reads:


<img src="https://latex.codecogs.com/gif.latex?C(x,T)&space;=&space;e^{-&space;r&space;\tau&space;}&space;\sum^{N-1}_{k=0}&space;'&space;Re\Big\{\varphi&space;\Big(&space;k&space;\pi&space;/&space;(b-a)&space;\Big)&space;\cdot&space;\exp\Big(-i&space;k&space;a&space;\pi/(b-a)&space;\Big)&space;\Big\}&space;V_k^{call}" title="C(x,T) = e^{- r \tau } \sum^{N-1}_{k=0} ' Re\Big\{\varphi \Big( k \pi / (b-a) \Big) \cdot \exp\Big(-i k a \pi/(b-a) \Big) \Big\} V_k^{call}" />


where <img src="https://latex.codecogs.com/gif.latex?\inline&space;V_k" title="V_k" /> call is the integral resulting from the inversion of summation and integral sign, which can be computed analytically. In particular:

<img src="https://latex.codecogs.com/gif.latex?V_k^{call}&space;=&space;\frac{2}{b-a}&space;K&space;(\chi_k(0,b)-\psi_k(0,b))" title="V_k^{call} = \frac{2}{b-a} K (\chi_k(0,b)-\psi_k(0,b))" />

Analytic forms for <img src="https://latex.codecogs.com/gif.latex?\inline&space;\chi_k" title="\chi_k" /> and <img src="https://latex.codecogs.com/gif.latex?\inline&space;\psi_k" title="\psi_k" /> are given in Fang and Oosterlee (2008).


<div align="right"><a href="#0">Back to top</a> </div>

## <div id="B2"> <a href="#0">Characteristic Functions</a> </div>

### <div id="BB1"> Black-Scholes Model </div>

The Black-Scholes model is arguably the most famous option pricing model. The model assumes that the underlying asset process can be described by a Geometric Brownian Motion. Thus, one of the main assumptions of the model is that the logarithmic returns are normally distributed. The stochastic differential equation takes the form:

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;dS_t&space;=&space;S_t&space;\left(&space;\mu&space;dt&space;&plus;&space;\sqrt{u_0}&space;\&space;\&space;dW_t&space;\right)&space;\notag&space;\end{align}" title="\begin{align} dS_t = S_t \left( \mu dt + \sqrt{u_t} \ \ dW_t \right) \notag \end{align}" />

where <img src="https://latex.codecogs.com/gif.latex?\mu" title="\mu" /> is the drift term and <img src="https://latex.codecogs.com/gif.latex?\inline&space;u_0" title="u_t" /> is the variance of the process.

From this, by making a change of measure using the risk neutral probability and after some manipulations, one gets the characteristic function of the Black Scholes model which follows a normal distribution: 

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\varphi_{BS}&space;=&space;exp&space;\left(&space;(&space;-&space;0.5&space;u_t&space;)&space;i&space;\Delta t&space;\omega&space;-&space;0.5&space;u_t&space;\Delta t&space;\omega^2&space;\right&space;)&space;\notag&space;\end{align}" title="\begin{align} \varphi_{BS} = exp \left( ( - 0.5 u_t ) i \Delta t \omega - 0.5 u_t \Delta t \omega^2 \right ) \notag \end{align}" />

Following Schmelzle (2010) page 39.

<details><summary>Click to see Matlab code</summary>
<p>

```Matlab
function phi_bs = bs_char_fn_v1(u_0, a, b, k, T)

%{
 This code computes the Characteristic Function for the Black-Scholes Model
 Notation follows Schmelzle (2010), eq. 5.10, p. 37

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_bs = bs_char_fn_v1(u_0, a, b, k, T)

 Inputs : u_0           - variance
        : T             - time to maturity
        : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - Vector of N evaluation intervals

Outputs : phi_bs        - characteristic function values [0:N-1] vector

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);

phi_bs = exp(( -0.5 * u_0 ) * 1i * T .* omega - 0.5 * u_0 * T * omega .^2);

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);

phi_bs = exp(( -0.5 * u_0 ) * 1i * T .* omega - 0.5 * u_0 * T * omega .^2);
```
</details>
</p>



### <div id="BB2"> Heston Model </div>

The Heston model incorporates a stochastic volatility term and can be described by this system of stochastic differential equations:

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;ds_t&space;&&space;=&space;(\mu&space;-0.5u_t)dt&space;&plus;&space;\sqrt{u_t}&space;\&space;\&space;dW_{1t}&space;\notag&space;\\&space;du_t&space;&=&space;\lambda(\bar{u}&space;-u_t)dt&space;&plus;&space;\eta&space;\sqrt{u_t}&space;\&space;\&space;dW_{2t}&space;\notag&space;\end{align}" title="\begin{align} ds_t & = (\mu -0.5u_t)dt + \sqrt{u_t} \ \ dW_{1t} \notag \\ du_t &= \lambda(\bar{u} -u_t)dt + \eta \sqrt{u_t} \ \ dW_{2t} \notag \end{align}" />

where <img src="https://latex.codecogs.com/gif.latex?\inline&space;\mu" title="\mu" /> is the drift term <img src="https://latex.codecogs.com/gif.latex?\inline&space;u_t" title="u_t" /> is the stochastic variance, <img src="https://latex.codecogs.com/gif.latex?\inline&space;u_0" title="u_0" /> is the initial volatility, <img src="https://latex.codecogs.com/gif.latex?\inline&space;\bar{u}" title="\bar{u}" /> is the long term variance, <img src="https://latex.codecogs.com/gif.latex?\inline&space;\lambda" title="\lambda" /> is the speed of mean reversion of the stochastic volatility, <img src="https://latex.codecogs.com/gif.latex?\inline&space;\eta" title="\eta" /> is the volatility of the volatility and <img src="https://latex.codecogs.com/gif.latex?\inline&space;W_{1t}" title="W_{1t}" /> and <img src="https://latex.codecogs.com/gif.latex?\inline&space;W_{2t}" title="W_{2t}" /> are two correlated Wiener processes.

The characteristic function of the log strike price can be written as:
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\begin{align}&space;\varphi_{hes}&space;=&space;&&space;exp\left(&space;i&space;\omega&space;\mu&space;\Delta&space;t&space;&plus;&space;\frac{u_t}{\eta^2}&space;\left(&space;\frac{1-&space;e^{-D&space;\Delta&space;t}}{1-Ge^{-D&space;\Delta&space;t}}&space;\right)&space;(\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;)&space;\right)&space;\times&space;\notag&space;\\&space;&&space;exp\left(&space;\frac{\lambda&space;\bar{u}}{\eta^2}&space;\left(&space;\Delta&space;t&space;*&space;(\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;-&space;D&space;)&space;-&space;2&space;log&space;\left(&space;\frac{1-Ge^{-D&space;\Delta&space;t}}{1-G}&space;\right&space;)&space;\right&space;)&space;\right&space;)&space;\notag&space;\end{align}" title="\begin{align} \varphi_{hes} = & exp\left( i \omega \mu \Delta t + \frac{u_t}{\eta^2} \left( \frac{1- e^{-D \Delta t}}{1-Ge^{-D \Delta t}} \right) (\lambda - i \rho \eta \omega ) \right) \times \notag \\ & exp\left( \frac{\lambda \bar{u}}{\eta^2} \left( \Delta t * (\lambda - i \rho \eta \omega - D ) - 2 log \left( \frac{1-Ge^{-D \Delta t}}{1-G} \right ) \right ) \right ) \notag \end{align}" />

For this project the parameters were taken from Fang and Osterle (2008) with <img src="https://latex.codecogs.com/gif.latex?\inline&space;\lambda" title="\lambda" /> = 1.5768, the long term variance <img src="https://latex.codecogs.com/gif.latex?\inline&space;\bar{u}" title="\bar{u}" /> = 0.0398, the initial variance <img src="https://latex.codecogs.com/gif.latex?\inline&space;u_0" title="u_0" /> = 0.0175, <img src="https://latex.codecogs.com/gif.latex?\inline&space;\eta" title="\eta" /> = 0.5751 and <img src="https://latex.codecogs.com/gif.latex?\inline&space;\rho" title="\rho" /> = -0.5711.

For more documentation on the Heston model, see Heston (1993) and for information on the Heston model with the Cosine Method see Fang and Oosterlee (2009).

<details><summary>Click to see Matlab code</summary>
<p>

```Matlab
function phi_hest = heston_char_fn_v2(mu, lambda, u_bar, u_0, eta, rho, a, b, k, T)

%{
 This code computes the Characteristic Function for the Heston Model
 Notation follows Fang and Oosterlee (2008), eq. 32, p. 8

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_hest = heston_char_fn_v2(mu, lambda, u_bar, u_0, eta, rho, a, b, k, T)


 Inputs : mu            - log price drift rate
        : lambda        - speed of mean reversion
        : u_bar         - mean (long run) volatility
        : u_0           - initial volatility
        : eta           - volatility of the volatility (vol of vol)
        : rho           - correlation between Wiener processes (W1 and W2)
        : T             - time to maturity
        : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - Vector of N evaluation intervals

Outputs : phi_hest      - characteristic function values [0:N-1] vector

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);


% D and G parameters: Fang and Oosterlee (2008) p. 8
D = sqrt( (lambda - 1i * rho * eta * omega) .^ 2 + (omega .^ 2 + 1i * omega) * eta ^ 2 );
G = (lambda - 1i * rho * eta * omega - D) ./ (lambda - 1i * rho * eta * omega + D);


% Characteristic function for the Heston Model: 
% Fang and Oosterlee (2008) p. 8
phi_hest = exp( 1i * omega * mu * T + ...
           (u_0 / eta ^ 2) * (1 - exp(-D * T)) ./ (1 - G .* exp(-D * T)) .* (lambda - 1i * rho * eta * omega - D) + ...
           (lambda * u_bar / eta ^ 2 * (T * (lambda - 1i * rho * eta * omega - D) - 2 * log( (1 - G .* exp(-D * T)) ./ (1 - G) ))));
      
  
%}
  
```
</details>
</p>

### <div id="BB3"> Variance Gamma Model </div>

The Variance Gamma model is obtained by evaluating a Black-Scholes style Geometric Brownian Motion at a random time change given by a gamma process with the purpose to control the skewness and kurtosis of the return distribution. (Madan et al. 1998) This creates an inifinite jump diffusion process to describe asset returns.

The characteristic function of the log stock price of the Variance Gamma model can be written as:

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\begin{align}&space;\varphi_{VG}&space;&&space;=&space;\frac{&space;exp\left(&space;\mu&space;&plus;&space;w_c&space;\right)^{i&space;\omega&space;}&space;}{&space;\left(&space;1&space;-&space;i&space;\theta&space;v&space;&plus;&space;0.5&space;u_t&space;v&space;\omega^2&space;\right)^{&space;\frac{\Delta&space;t}{v}}&space;}&space;\notag&space;\\&space;&&space;where&space;\&space;\&space;w_c&space;=&space;\frac{1}{v}&space;*&space;ln(1-&space;i&space;\theta&space;*&space;v&space;-0.5&space;u_t&space;v)&space;\notag&space;\end{align}" title="\begin{align} \varphi_{VG} & = \frac{ exp\left( \mu + w_c \right)^{i \omega } }{ \left( 1 - i \theta v + 0.5 u_t v \omega^2 \right)^{ \frac{\Delta t}{v}} } \notag \\ & where \ \ w_c = \frac{1}{v} * ln(1- i \theta * v -0.5 u_t v) \notag \end{align}" />

where <img src="https://latex.codecogs.com/gif.latex?\inline&space;w_c" title="w_c" /> is a drift correction term and <img src="https://latex.codecogs.com/gif.latex?\inline&space;\theta" title="\theta" /> is an additional drift term controlling the skewness which is set to -0.1436 following an estimation of Madan et al. (1998). <img src="https://latex.codecogs.com/gif.latex?\inline&space;v" title="v" /> is the variance rate of the gamma process controlling the kurtosis of the distribution and is set to 0.0403 following Carr et al. (2002).

<details><summary>Click to see Matlab code</summary>
<p>

```Matlab
function phi_vg = vg_char_fn(u_0, theta, a, b, k, T, v, mu)

%{
 This code computes the Characteristic Function for the Variance Gamma Process. 
 Notation follows Madan, Carr and Chang (1998), eq. 7, p. 8. 
 theta's average value in Madan, Carr, Chang (1998) is -0.1436, 
 theta = 0 gives the symmetric variance gamma,
 theta < 0 incorporates negative skewness in the model.
 
 Source: A.Itkin ”Pricing options with VG model using FFT”. 
                   The Variance Gamma and Related Financial Models. August 9, 2007 – p. 7

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_vg = vg_char_fn(u_0, v, theta, a, b, k, T)


 Inputs : u_0           - variance
        : theta         - drift term
        : T             - time to maturity
        : omega         - argument of the characteristic function
        : v             - variance rate
        : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - Vector of N evaluation intervals


Outputs : phi_vg        - characteristic function values [0:N-1] vector

%}

% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);

w = ( 1 / v ) * log(1 - theta * v - 0.5 * u_0 * v );

phi_vg = exp( ( mu * T + w  * T) ).^(1i * omega) ./ ( 1 - 1i * theta * v * omega + ( u_0 * v / 2 ) * omega.^2 ).^(T/v);
```
</details>
</p>

### <div id="BB4"> CGMY Model </div>

The CGMY model can be seen as an extension to the Variance Gamma model allowing for both allows for both diffusions and for jumps of both finite and infinite activity. (Carr et al., 2002) 

The characteristic function of the CGMY model can be written as:

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\varphi_{CGMY}&space;&&space;=&space;exp\left(&space;i&space;\omega&space;\mu&space;\Delta t&space;-&space;0.5&space;\omega^2&space;u_t&space;\Delta t&space;\right)&space;\cdot&space;exp\left(&space;\Delta t&space;C&space;\Gamma(-Y)\left[&space;(M&space;-&space;i&space;\omega&space;)^Y&space;-M^Y&space;&plus;&space;(G&plus;i&space;\omega)^Y&space;-&space;G^Y&space;\right]&space;\right&space;)&space;\notag&space;\end{align}" title="\begin{align} \varphi_{CGMY} & = exp\left( i \omega \mu \Delta t - 0.5 \omega^2 \sigma^2 \Delta t \right) \cdot exp\left( \Delta t C \Gamma(-Y)\left[ (M - i \omega )^Y -M^Y + (G+i \omega)^Y - G^Y \right] \right ) \notag \end{align}" />

where C can be described as a measure of the overall level of activity of the process, G controls the exponential decay on the right of the process and M controls the exponential decay on the left. For example, if G < M, the left tail of the distribution would be heavier
than the right tail. The parameter Y icharacterizes the fine structure of the jumps of the CGMY process and a more in depth description of it can be found in Carr et al. (2002).



For the purpose of this project the parameters C, G, M and Y were chosen according to Carr et al. (2002) with values of 24.79, 94.45, 95.79 and 0.2495 respectively.


<details><summary>Click to see Matlab code</summary>
<p>

```Matlab
function phi_cgmy = cgmy_char_fn(mu, u_0, C, G, Y, M, a, b, k, T)

%{

 This code computes the Characteristic Function for the CGMY model process. 
 Notation follows Fang and Oosterlee (2008) eq. 31, p. 8 and parameters are
 inspired by Carr et al. (2002).

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (22.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_cgmy = cgmy_char_fn(mu, u_0, C, G, Y, M, a, b, k, T)


 Inputs : mu            - r - q
        : T             - time to maturity
        : omega         - argument of the characteristic function
        : C             - measure of the overall level of activity of the
                          process, controls the kurtosis
        : G             - controls the exponential decay on the right tail 
                          of the process
        : M             - controls the exponential decay on the left tail 
                          of the process
        : Y             - characterizes the fine structure of the jumps of 
                          the CGMY process
        : u_0           - variance

Outputs : phi_cgmy      - characteristic function values [0:N-1] vector

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);


% Characteristic function
phi_cgmy = exp( 1i * omega * ( mu ) * T - 0.5 * omega.^2 * u_0 * T) .* ...
           exp( T * C * gamma( -Y ) * ( (M - 1i * omega).^Y - M^Y + (G + 1i * omega).^Y - G^Y ) );


```
</details>
</p>

### <div id="BB5"> Parameter Values </div>


Parameter values are given in the following:

<details><summary>Click to see Matlab code</summary>
<p>

```Matlab
%% Model Parameters

% Heston Model Parameters [defaults from Fang and Oosterlee (2008)]
lambda      =  1.5768;                          % Mean Reversion Speed
u_bar       =  0.0398;                          % Long-run Variance
u_0         =  0.0175;                          % Initial Variance
eta         =  0.5751;                          % Volatility of the Volatility Parameter
rho         = -0.5711;                          % Correlation between Wiener Processes


% VG Parameters [defaults from Madan et al. (1998) and Carr et al. (2002)]
theta      = -0.1436;
v          = 0.0403;


% CGMY Parameters [defaults from Carr et al. (2002)]
C          = 24.79;
G          = 94.45;
Y          = 0.2495;
M          = 95.79;

```
</details>
</p>

<div align="right"><a href="#0">Back to top</a> </div>


## <div id="C2"> <a href="#0">Distributional Characteristics of the different Models</a> </div>

The Figure below shows the distributional characteristics of the models considered in this project compared to a normal distribution. The distributions were obtained using the Fourier-cosine expansion using the given characteristic function with the assumed parameters.

Unsurprisingly, the fit of the Black-Scholes model is very close to the normal distribution as the distribution of the assumed underlying process is normal. For the Heston model the distribution has a higher kurtosis and negative skewness such that the left tail is fatter compared to a normal. 

Furthermore, the distribution of the assumed Variance Gamma model is close to normal (with a potential slight negative skew). This is mainly due to the choice of v and theta. Firstly, theta, the parameter controlling the skewness is chosen to be relatively small with -0.1436. Secondly, v, the parameter controlling the kurtosis, is chosen to be close to zero with 0.0403 and the Variance Gamma model approaches a geometric Brownian Motion as v approaches zero. (Daal and Madan, 2005)

Lastly, the CGMY distribution with the given parameters has a lower kurtosis and slight positive skew compared to the normal so that there is more probability mass in the tails of the distribution. The slight negative skew is due to M (95.79) being chosen to be larger than G (94.45) so that the left tail of the distribution is slightly heavier than the right tail. The kurtosis could be controlled via the parameter C which was chosen to be 24.79 following Carr et al. (2002).


<img src="Plots/DistributionPlots.png"
     alt="Call Option Plots"
     style="float: left; margin-right: 10px; padding-bottom: 30px;" />

<div align="right"><a href="#0">Back to top</a> </div>

## <div id="D2"> <a href="#0">Option Pricing with the COS-FFT Method</a> </div>

The figure below compares the call prices calculated using the Cos-FFT method compared to the analytical Black-Scholes solutions for 66 SPX options of the same expiration date of 20.03.2020 with strike prices ranging from 1275 to 3600. All computations are perfoed with N = 5000 evaluation points. The time to maturity is 1.0056, i.e. approximately one year. The price of the underlying is 2800.7, the risk-free rate is assumed to be approximately 2.39% (the current 3-month T-bill rate) and the continuous dividend yield is assumed to be 1.92% taken from [http://www.multpl.com/s-p-500-dividend-yield/table](http://www.multpl.com/s-p-500-dividend-yield/table). The option chain and parameters were obtained using the R code in `retrieve_option_data.R` and `project_functions.R`. In successive calibrations, a better estimate of the (expected) dividend yield may be inferred through put-call parity and regression.

From the figure it is clear that all four models have significant deviations from the call option price implied by the analytical Black-Scholes model. However, this says nothing about the accuracy of the pricing and it should not be expected to be precise pricing given that the parameters were chosen from what was computed in the original papers, rather than being calibrated with present market data. 

The deviation from the analytical Black-Scholes from the Black-Scholes calculated using the COS-FFT-Method is somewhat surprising. However, when considered in relative terms, compared to the option price, the distance between the two is not large. Furthermore, the first three models approach a call option price of zero as the strike price increases as would be expected. The CGMY model, however, does not seem to converge as quickly, and it exhibits relatively larger differences in pricing with respect to to other models. The most obvious reason for this phenomenon could be that the CGMY takes the highest number of unknown parameters. This may cause larger errors.

Despite the fact that the assumed parameters could increase errors, a relevant error component will be due to the truncation of the density function domains, as well as the truncation of the expansion. Truncation bounds are calculated following Fang and Oosterlee (2008), using the cumulants of the four models and the recommended scaling factors of 12 for the Heston model and 10 for the other three models.  

In addition to the truncation error, there is the possibility of a discretisation error having an influence on the pricing. The COS-FFT method provides only an approximation to the true values and thus, although more accurate than a pure FFT, can include such discretisation errors. A good example of the possibility of such an error is the comparison of the Black-Scholes distribution calculated using the COS-FFT and the true normal distribution shown in the Figure above.

In the end the pricing will only be as good as the approximation under the chosen assumptions. Analysing the performance of different characteristic functions using optimised parameters can therefore be a good approach to identify potential errors. The COS-FFT method makes this particularly simple and faster. However, to compare the true accuracy of the COS-FFT with the common FFT method, one needs the optimised parameters and ideally an evaluation sample of real option prices.

<img src="Plots/CallPlotsSPX.png"
     alt="Call Option Plots"
     style="float: left; margin-right: 10px; padding-bottom: 30px;" />



<div align="right"><a href="#0">Back to top</a> </div>


## <div id="E2"> <a href="#0">Concluding Remarks</a> </div>

In this project, we price one day of put and call options using the COS method. We report on how call option prices in particular relate to Black-Scholes benchmarks obtained analitically for four different models. In particular, we cover the Black-Scholes, Heston, Variance Gamma and CGMY models.

In summary, we find that option prices delivered by the COS method reflect the expected features of the assumed diffusions for the underlying price (and variance). Further, we note that the entire estimation process takes only around 0.025 seconds on our (average laptop) setup. However, the present implementation restricts us to comparing COS prices to a benchmark, since model parameters are taken as given, rather than being calibrated. Deviations from realised prices may thus be entirely an artifact of the model setup. Therefore, a comparison of the COS method to other pricing techniques could only be drawn in terms of speed, and it is thus not pursued here. 

Further additions could look to compare both speed and accuracy of the COS method with respect to other popular techniques. In particular, the performance of the four models considered under the COS-FFT algorithm could be contrasted with the FFT approach. Through a discretisation and simulation of the underlying price path, with techniques such as the Milstein scheme, COS and the Monte Carlo pricing could also be compared.

<div align="right"><a href="#0">Back to top</a> </div>


## <div id="E2"> <a href="#0">References</a> </div>

* Carr, P., Geman, H., Madan, D. and Yor, M. (2002). The Fine Structure of Asset Returns: An Empirical Investigation. The Journal of Business, 75(2), pp.305-333. 
* Carr, P. and Madan, D. (1999). Option valuation using the fast Fourier transform. The Journal of Computational Finance, 2(4), pp.61-73.
* Daal, E. and Madan, D. (2005). An Empirical Examination of the Variance‐Gamma Model for Foreign Currency Options. The Journal of Business, 78(6), pp.2121-2152.
* Fang, F. and Oosterlee, C. (2008). A Novel Pricing Method for European Options Based on Fourier-Cosine Series Expansions. SIAM Journal on Scientific Computing, 31(2), pp.826-848.
* Heston, S. (1993). A Closed-Form Solution for Options with Stochastic Volatility with Applications to Bond and Currency Options. Review of Financial Studies, 6(2), pp.327-343.
* Lee, R. W. (2004). Option Pricing by Transform Methods: Extensions, Uni$cation, and Error Control, Journal of Computational Finance 7(3), 51–86.
* Lewis, A. (2001). A Simple Option Formula for General Jump-Diffusion and other Exponential Lévy Processes, Envision Financial Systems and OptionCity.net, California.
* Madan, D., Carr, P. and Chang, E. (1998). The Variance Gamma Process and Option Pricing. Review of Finance, 2(1), pp.79-105.
* Schmelzle, M. (2010). Option Pricing Formulae using Fourier Transform: Theory and Application.

<div align="right"><a href="#0">Back to top</a> </div>

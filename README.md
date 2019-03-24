<div align="right">
Advanced Numerical Methods and Data Analysis - FS19-8,780
<br>
University of St. Gallen, 24.03.2019
<br>
</div>

-------------



# Option Pricing using the Cosine Fast Fourier Transform Method


**Federico Baldi Lanfranchi** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; federico.baldilanfranchi@student.unisg.ch <br>
**Peter la Cour** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;peter.lacour@student.unisg.ch

## <div id="0"><a href="#0">Overview</a></div>

1. <a href="#2">Introduction</a>
2. <a href="#A2">Description of the Cos-FFT Method</a>
3. <a href="#B2">Characteristic Functions</a>
   * <a href="#BB1">Black-Scholes Model </a>
	* <a href="#BB2">Heston Model</a>
	* <a href="#BB3">Variance Gamma Model</a>
	* <a href="#BB4">CGMY Model</a>
	* <a href="#BB5">Other Models </a>
4. <a href="#C2">Distributional Characteristics of the different Models</a>
5. <a href="#D2">Pricing Accuracy of the Cos-FFT Method</a>
6. <a href="#E2">Concluding Remarks</a>
7. <a href="#F2"> References </a>


## <div id="2"> <a href="#0">Introduction  </a> </div>








## <div id="A2"> <a href="#0">Description of the Cos-FFT Method</a> </div>


<div align="right"><a href="#0">Back to top</a> </div>

## <div id="B2"> <a href="#0">Characteristic Functions</a> </div>

### <div id="BB1"> Black-Scholes Model </div>

The Black-Scholes model is arguably the most famous option pricing model. The model assumes that the underlying asset process can be described by a Geometric Brownian Motion. One of the main assumptions of the model is that the logarithmic returns are normally distributed. The stochastic differential equation takes the form:

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\begin{align}&space;dS_t&space;=&space;S_t&space;\left(&space;\mu&space;dt&space;&plus;&space;\sqrt{u_0}&space;\&space;\&space;dW_t&space;\right)&space;\notag&space;\end{align}" title="\begin{align} dS_t = S_t \left( \mu dt + \sqrt{u_0} \ \ dW_t \right) \notag \end{align}" />

where <img src="https://latex.codecogs.com/gif.latex?\inline&space;\mu" title="\mu" /> is the drift term and <img src="https://latex.codecogs.com/gif.latex?\inline&space;u_0" title="u_0" /> is the variance of the process.

From this, by making a change of measure using the risk neutral probability and after some manipulations, one gets the characteristic function of the Black Scholes model which follows a normal distribution: 

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\varphi_{BS}&space;=&space;exp&space;\left(&space;(&space;\mu&space;-&space;0.5&space;u_0&space;)&space;i&space;T&space;\omega&space;-&space;0.5&space;u_0&space;T&space;\omega^2&space;\right&space;)&space;\notag&space;\end{align}" title="\begin{align} \varphi_{BS} = exp \left( ( \mu - 0.5 u_0 ) i T \omega - 0.5 u_0 T \omega^2 \right ) \notag \end{align}" />

### <div id="BB2"> Heston Model </div>

The Heston Model incorporates a stochastic volatility term and can be described by this system of stochastic differential equations:

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\begin{align}&space;ds_t&space;&&space;=&space;(&space;\mu&space;-&space;0.5&space;u_t&space;)&space;dt&space;&plus;&space;\sqrt{u_t}&space;\&space;\&space;dW_{1t}&space;\notag&space;\\&space;du_t&space;&&space;=&space;\lambda&space;(&space;\bar{u}&space;-&space;u_t&space;)&space;dt&space;&plus;&space;\eta&space;\sqrt{u_t}&space;\&space;\&space;dW_{2t}&space;\notag&space;\end{align}" title="\begin{align} ds_t & = ( \mu - 0.5 u_t ) dt + \sqrt{u_t} \ \ dW_{1t} \notag \\ du_t & = \lambda ( \bar{u} - u_t ) dt + \eta \sqrt{u_t} \ \ dW_{2t} \notag \end{align}" />

The characteristic function of the log strike price can be written as:

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\varphi_{hes}&space;&&space;=&space;exp\left(&space;i&space;\omega&space;\mu&space;T&space;&plus;&space;\frac{u_0}{\eta^2}&space;\left(&space;\frac{1-&space;e^{-DT}}{1-Ge^{-DT}}&space;\right)&space;(\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;)&space;\right)&space;\notag&space;\\&space;where&space;\&space;\&space;D&space;&&space;=&space;\sqrt{&space;(\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;)^2&space;&plus;&space;(&space;\omega^2&space;&plus;&space;i&space;\omega&space;)&space;\eta^2&space;}&space;\&space;\&space;and&space;\&space;\&space;G&space;=&space;\frac{\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;-&space;D&space;}{&space;\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;&plus;&space;D&space;}&space;\notag&space;\end{align}" title="\begin{align} \varphi_{hes} & = exp\left( i \omega \mu T + \frac{u_0}{\eta^2} \left( \frac{1- e^{-DT}}{1-Ge^{-DT}} \right) (\lambda - i \rho \eta \omega ) \right) \notag \\ where \ \ D & = \sqrt{ (\lambda - i \rho \eta \omega )^2 + ( \omega^2 + i \omega ) \eta^2 } \ \ and \ \ G = \frac{\lambda - i \rho \eta \omega - D }{ \lambda - i \rho \eta \omega + D } \notag \end{align}" />

For more documentation on the Heston model see Heston (1993) and for information on the Heston model with the Cosine Method see Fang and Osterle (2008).

### <div id="BB3"> Variance Gamma Model </div>

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\varphi_{VG}&space;&&space;=&space;\frac{&space;exp\left(&space;\mu&space;&plus;&space;w_c&space;\right)^{i&space;\omega&space;}&space;}{&space;\left(&space;1&space;-&space;i&space;\theta&space;v&space;&plus;&space;0.5&space;u_0&space;v&space;\omega^2&space;\right)&space;}&space;\notag&space;\\&space;where&space;\&space;\&space;w_c&space;&&space;=&space;\frac{1}{v}&space;*&space;ln(1-&space;i&space;\theta&space;*&space;v&space;-0.5&space;u_0&space;v)&space;\notag&space;\end{align}" title="\begin{align} \varphi_{VG} & = \frac{ exp\left( \mu + w_c \right)^{i \omega } }{ \left( 1 - i \theta v + 0.5 u_0 v \omega^2 \right) } \notag \\ where \ \ w_c & = \frac{1}{v} * ln(1- i \theta * v -0.5 u_0 v) \notag \end{align}" />



### <div id="BB4"> CGMY Model </div>

<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\varphi_{CGMY}&space;&&space;=&space;exp\left(&space;i&space;\omega&space;\mu&space;T&space;-&space;0.5&space;\omega^2&space;\sigma^2&space;T&space;\right)&space;\cdot&space;exp\left(&space;T&space;C&space;\Gamma(-Y)\left[&space;(M&space;-&space;i&space;\omega&space;)^Y&space;-M^Y&space;&plus;&space;(G&plus;i&space;\omega)^Y&space;-&space;G^Y&space;\right]&space;\right&space;)&space;\notag&space;\end{align}" title="\begin{align} \varphi_{CGMY} & = exp\left( i \omega \mu T - 0.5 \omega^2 \sigma^2 T \right) \cdot exp\left( T C \Gamma(-Y)\left[ (M - i \omega )^Y -M^Y + (G+i \omega)^Y - G^Y \right] \right ) \notag \end{align}" />


### <div id="BB5"> Other Models </div>




<div align="right"><a href="#0">Back to top</a> </div>



## <div id="C2"> <a href="#0">Distributional Characteristics of the different Models</a> </div>

Add Black Scholes, (FMLS?)

<img src="Plots/DistributionPlots.png"
     alt="Call Option Plots"
     style="float: left; margin-right: 10px; padding-bottom: 30px;" />

<div align="right"><a href="#0">Back to top</a> </div>

## <div id="D2"> <a href="#0">Pricing Accuracy of the Cos-FFT Method</a> </div>

Add Black Scholes and Puts?

<img src="Plots/CallPlots.png"
     alt="Call Option Plots"
     style="float: left; margin-right: 10px; padding-bottom: 30px;" />



<div align="right"><a href="#0">Back to top</a> </div>


## <div id="E2"> <a href="#0">Concluding Remarks</a> </div>

<div align="right"><a href="#0">Back to top</a> </div>


## <div id="E2"> <a href="#0">References</a> </div>

<div align="right"><a href="#0">Back to top</a> </div>
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



### <div id="BB2"> Heston Model </div>


<img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;\varphi_{hes}&space;&&space;=&space;exp\left(&space;i&space;\omega&space;\mu&space;T&space;&plus;&space;\frac{u_0}{\eta^2}&space;\left(&space;\frac{1-&space;e^{-DT}}{1-Ge^{-DT}}&space;\right)&space;(\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;)&space;\right)&space;\notag&space;\\&space;where&space;\&space;\&space;D&space;&&space;=&space;\sqrt{&space;(\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;)^2&space;&plus;&space;(&space;\omega^2&space;&plus;&space;i&space;\omega&space;)&space;\eta^2&space;}&space;\&space;\&space;and&space;\&space;\&space;G&space;=&space;\frac{\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;-&space;D&space;}{&space;\lambda&space;-&space;i&space;\rho&space;\eta&space;\omega&space;&plus;&space;D&space;}&space;\notag&space;\end{align}" title="\begin{align} \varphi_{hes} & = exp\left( i \omega \mu T + \frac{u_0}{\eta^2} \left( \frac{1- e^{-DT}}{1-Ge^{-DT}} \right) (\lambda - i \rho \eta \omega ) \right) \notag \\ where \ \ D & = \sqrt{ (\lambda - i \rho \eta \omega )^2 + ( \omega^2 + i \omega ) \eta^2 } \ \ and \ \ G = \frac{\lambda - i \rho \eta \omega - D }{ \lambda - i \rho \eta \omega + D } \notag \end{align}" />


### <div id="BB3"> Variance Gamma Model </div>





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
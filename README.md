## GVCM4GT

This repository contains R codes along with simulation results for "**Generalized Bayesian varying coefficient mixed models for group testing data**". Our model is try to estimate an individual-level regression model based on
group testing data that can capture the age-varying impact on
the Chlamydia risk with selection. To relate available information, we consider

$$
\text{logit}(\text{Pr}(\widetilde Y_i=1\mid \boldsymbol x_i, u_i))=\underbrace{\sum_{d=0}^p x_{id}\psi_d(u_i)}\_{\text{Age-varying Effects}} + \underbrace{\sum_{\ell=1}^L r_\ell(i)\gamma_\ell}\_{\text{Random Effect}} \quad\text{for }i=1,\ldots,N,
$$

where $\boldsymbol x_i=(x_{i0},x_{i1},\ldots,x_{ip})^\top$ and $\psi_d(u_i)=\delta_{1d}(\alpha_d+\delta_{2d}\beta_d(u_i))$ for binary indicators $\delta_{1d},\delta_{2d}$; see more details in the paper. In short, our stochastic search variable selection categorize each of covariates into three groups:

- $\delta_{1d}=0\longrightarrow$ insignificant effects.
- $\delta_{1d}=1$
	* $\delta_{2d}=0\longrightarrow$ age-independent effects.
	* $\delta_{2d}=1\longrightarrow$ age-varying effects.

To reproduce the results in the paper, we provide implementation details as follows. 

```sh
username@login001 ~$ git clone git@github.com:yizenglistat/GVCM4GT.git
username@login001 ~$ cd GVCM4GT
```

### Arguments

```r
# A demo example to run 500 repetitions in one machine.
task_id <- 1 						
nreps <- 500
Ns <- c(3000, 5000)
pool_sizes <- c(5, 10)
model_names <- c("m1", "m2")
testings <- c("AT", "DT", "IT")
N_test <- 600
sigma <- 0.5
```

- `task_id`
> The machine id. For example, 1,...,100 if running on the cluster. In this way, we will run 5 simulations independently on 100 nodes to have a total of 500 repetitions. 

- `nreps`
> The repetitions.

- `Ns`
> A vector of sample sizes.

- `pool_sizes`
> A vector of pool sizes.

- `model_names`
> A vector of model names. Different model names corresponds to different varying function sets.

- `testings`
> A vector of testing protocols such as AT (array testing), DT (Dorfman Testing) or IT (Individual Testing).

- `N_test`
> Number of knots values in inference for estimated varying functions. 

- `sigma`
> True random effect standard deviation

### Reproduce

After setting up the environment (`requirement.txt`) and arguments, one should be able to run the following code in `R` to reproduce simulation results in the paper.

```r
source('main.r')
```

After collecting `.RData` files under `output/`, one should be able to reproduce the results subsequently. The following demo figure and demo table show that $\textcolor{red}{\textbf{red}}$ means the $\textcolor{red}{\textbf{age-varying effects}}$, $\textcolor{blue}{\textbf{blue}}$ means the significant but $\textcolor{blue}{\textbf{age-independent effects}}$ and $\textbf{black}$ means the $\textbf{insignificant effects}$ can be both correctly identified and estimated; see details in the paper.

![figure](output/uniform_5000_m1_figure.png)


<table>
  
  <tr style="border-top: 2px solid black; border-bottom: 1.5px solid black;">
    <th rowspan="2" align="left">Parameter</th>
    <th rowspan="2" align="left">Summary</th>
    <th rowspan="2">IT</th>
    <th colspan="2">c=5</th>
    <th colspan="2">c=10</th>
  </tr>

<tr style="border-bottom: 1.5px solid black;">
    <td align="center" style="border-bottom: 1.5px solid black;">DT</td>
    <td align="center" style="border-bottom: 1.5px solid black;">AT</td>
    <td align="center" style="border-bottom: 1.5px solid black;">DT</td>
    <td align="center" style="border-bottom: 1.5px solid black;">AT</td>
  </tr>
  
  <tr>
    <td align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\psi_5(u)=0" alt="\psi_5(u)=0" /></td>
    <td>IP</td>
    <td>0.007</td>
    <td>0.007</td>
    <td>0.007</td>
    <td>0.007</td>
    <td>0.007</td>
  </tr>
  
  <tr>
    <td align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\psi_6(u)=0" alt="\psi_6(u)=0" /></td>
    <td style="border-bottom: 1.5px solid black;">IP</td>
    <td style="border-bottom: 1.5px solid black;">0.007</td>
    <td style="border-bottom: 1.5px solid black;">0.007</td>
    <td style="border-bottom: 1.5px solid black;">0.008</td>
    <td style="border-bottom: 1.5px solid black;">0.007</td>
    <td style="border-bottom: 1.5px solid black;">0.007</td>
  </tr>

   <tr>
    <td align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{\psi_1(u)=-1.0}" alt="\textcolor{blue}{\psi_1(u)=-1.0}" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{\text{IPF}}\color{black}{/\text{IPV}}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.994}\color{black}{/0.005}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.993}\color{black}{/0.007}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.993}\color{black}{/0.007}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.993}\color{black}{/0.007}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.993}\color{black}{/0.007}" alt="" /></td>
  </tr>
  
  <tr>
    <td align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{\psi_3(u)=-0.5}" alt="\textcolor{blue}{\psi_3(u)=-0.5}" /></td>
    <td style="border-bottom: 1.5px solid black;"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{\text{IPF}}\color{black}{/\text{IPV}}" alt="" /></td>
    <td style="border-bottom: 1.5px solid black;"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.993}\color{black}{/0.007}" alt="" /></td>
    <td style="border-bottom: 1.5px solid black;"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.993}\color{black}{/0.007}" alt="" /></td>
    <td style="border-bottom: 1.5px solid black;"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.994}\color{black}{/0.006}" alt="" /></td>
    <td style="border-bottom: 1.5px solid black;"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.987}\color{black}{/0.013}" alt="" /></td>
    <td style="border-bottom: 1.5px solid black;"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{blue}{0.989}\color{black}{/0.011}" alt="" /></td>
  </tr>

   <tr>
    <td align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{red}{\delta_{12}=\delta_{22}=1}" alt="\color{red}{\delta_{12}=\delta_{22}=1}" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{\text{IPF/}}\color{red}{\text{IPV}}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
  </tr>
  
  <tr style="border-bottom: 1.5px solid black;">
    <td align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{red}{\delta_{14}=\delta_{24}=1}" alt="\color{red}{\delta_{14}=\delta_{24}=1}" /></td>
   <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{\text{IPF/}}\color{red}{\text{IPV}}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
    <td><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\color{black}{0.000/}\color{red}{1.000}" alt="" /></td>
  </tr>

  <tr>
    <th rowspan="2" align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\psi_1(u)=\alpha_1=-1.0" alt="" /></th>
     <td>Bias(CP95)</td>
     <td>-0.008(0.962)</td>
     <td>-0.015(0.940)</td>
     <td>-0.004(0.942)</td>
     <td>-0.034(0.912)</td>
     <td>-0.013(0.932)</td>
  </tr>

  <tr align="left">
    <td>SSD(ESE)</td>
    <td>0.066(0.069)</td>
    <td>0.072(0.072)</td>
    <td>0.068(0.067)</td>
    <td>0.095(0.083)</td>
    <td>0.077(0.074)</td>
  </tr>

  <th rowspan="2" align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\psi_3(u)=\alpha_3=-0.5" alt="" /></th>
     <td>Bias(CP95)</td>
     <td>-0.007(0.942)</td>
     <td>-0.008(0.936)</td>
     <td>-0.002(0.938)</td>
     <td>-0.012(0.938)</td>
     <td>-0.001(0.938)</td>
  </tr>

  <tr align="left">
    <td>SSD(ESE)</td>
    <td>0.064(0.063)</td>
    <td>0.067(0.063)</td>
    <td>0.067(0.062)</td>
    <td>0.072(0.068)</td>
    <td>0.066(0.065)</td>
  </tr>

  <th rowspan="2" align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;\sigma=0.5\quad" alt="" /></th>
     <td>Bias(CP95)</td>
     <td>0.047(0.950)</td>
     <td>0.044(0.952)</td>
     <td>0.039(0.964)</td>
     <td>0.059(0.924)</td>
     <td>0.048(0.960)</td>
  </tr>

  <tr align="left" style="border-bottom: 1.5px solid black;">
    <td>SSD(ESE)</td>
    <td>0.062(0.079)</td>
    <td>0.062(0.080)</td>
    <td>0.058(0.078)</td>
    <td>0.067(0.084)</td>
    <td>0.062(0.081)</td>
  </tr>

  <th rowspan="2" align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;S_{e(1)}=0.95" alt="" /></th>
     <td>Bias(CP95)</td>
     <td></td>
     <td>-0.032(0.902)</td>
     <td>0.000(0.954)</td>
     <td>-0.038(0.914)</td>
     <td>0.002(0.964)</td>
  </tr>

  <tr align="left">
    <td>SSD(ESE)</td>
    <td></td>
    <td>0.081(0.053)</td>
    <td>0.024(0.021)</td>
    <td>0.084(0.062)</td>
    <td>0.041(0.034)</td>
  </tr>

  <th rowspan="2" align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;S_{e(2)}=0.98" alt="" /></th>
     <td>Bias(CP95)</td>
     <td></td>
     <td>-0.008(0.974)</td>
     <td>-0.001(0.988)</td>
     <td>-0.019(0.972)</td>
     <td>-0.006(0.996)</td>
  </tr>

  <tr align="left">
    <td>SSD(ESE)</td>
    <td></td>
    <td>0.021(0.024)</td>
    <td>0.016(0.017)</td>
    <td>0.050(0.046)</td>
    <td>0.021(0.033)</td>
  </tr>

  <th rowspan="2" align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;S_{p(1)}=0.98" alt="" /></th>
     <td>Bias(CP95)</td>
     <td></td>
     <td>0.002(0.990)</td>
     <td>0.000(0.920)</td>
     <td>-0.014(0.992)</td>
     <td>-0.010(0.990)</td>
  </tr>

  <tr align="left">
    <td>SSD(ESE)</td>
    <td></td>
    <td>0.011(0.014)</td>
    <td>0.012(0.011)</td>
    <td>0.032(0.052)</td>
    <td>0.027(0.033)</td>
  </tr>

  <th rowspan="2" align="left"><img src="https://latex.codecogs.com/png.latex?\dpi{300}&space;S_{p(2)}=0.99" alt="" /></th>
     <td>Bias(CP95)</td>
     <td></td>
     <td>0.000(0.966)</td>
     <td>-0.003(0.974)</td>
     <td>-0.003(0.922)</td>
     <td>-0.003(0.966)</td>
  </tr>

  <tr align="left">
    <td>SSD(ESE)</td>
    <td></td>
    <td>0.008(0.007)</td>
    <td>0.012(0.013)</td>
    <td>0.011(0.009)</td>
    <td>0.012(0.011)</td>
  </tr>

  <tr style="border-top: 1.5px solid black;">
    <td align="left">Cost</td>
    <td>AVGtest</td>
    <td>5000</td>
    <td>2943.15</td>
    <td>2971.33</td>
    <td>3567.84</td>
    <td>2943.73</td>
  </tr>

  <tr style="border-bottom: 2px solid black;">
    <td align="left">Savings</td>
    <td>Percent</td>
    <td>00.00%</td>
    <td>41.14%</td>
    <td>40.57%</td>
    <td>28.64%</td>
    <td>41.12%</td>
  </tr>

</table>

> IP: the inclusion probability of the any significant effect, i.e., $\alpha_d$ or $\beta_d(u)$. 

> IPF: the inclusion probability of the age-independent effect, i.e., $\alpha_d$ only.

> IPV: the inclusion probability of the age-varying effect, i.e., $\beta_d(u)$ only.

### Authors

* [Yizeng Li](https://yizengli.com)
* [Dewei Wang\*](https://sites.google.com/view/deweiwang)(the corresponding author)
* [Joshua M. Tebbs](https://people.stat.sc.edu/tebbs/)
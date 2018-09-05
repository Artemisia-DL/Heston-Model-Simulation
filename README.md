# Disclamer
This project was done while studing Heston Model for personal interest. Therefore the code leaves room to many improvements especially on the side of optimization and speed. The main aim was didatic, hence it may lack of elegance and precision, but I hope it maybe usefull to anyone and I encourage the few interested to improve it. I also tried to comment as clearly as possible the code, in order to mirror the Andersen paper (2007). 

I also thank TORBJÖRN ODELMAN whose Thesis ( Efficient Monte Carlo Simulation with Stochastic Volatility, 2009 ) inspired this project and from which I took and adapted the text below in order to better match my code. 

I strongly suggest the to read the Thesis above, accessible at this link: https://pdfs.semanticscholar.org/f155/24f29b9aebc3a241b8cdb3ed101162b66e57.pdf



# Heston-Model-Simulation
Heston Model Numerical Simulation with QE  Discretization Algorithm as in Andersen, L., (2007)

Heston Model Original Formulation:
   
    dSt = µ*St*dt + Vt^0.5 * St * dWt1      Security Price SDE
    dVt = κ(θ − Vt)dt + η*Vt^0.5 *dWt2      Variance SDE

For convenience during discretization the actual SDEs used where: 

    dXt = (r − Vt/2)dt + Vt^0.5 * dWt1     Security Return SDE computed with Ito's Lemma from St SDE
    dVt = κ(θ − Vt)dt + η*Vt^0.5 *dWt2     Variance SDE

# Quadratic Exponential Scheme
The Quadratic Exponential (QE) scheme is an algorithm specifically developed for square-root processes such as the Heston model. Due to the complicated and somewhat non-intuitive nature of this scheme, we will go through it rather briefly but yet as comprehensible as possible. For all the details the reader is refered to Andersen’s paper. The main idea, however, is to find distributions that behave similarly to the chi-square. With the knowledge that the conditional distribution of Vt+∆ given Vt is non-central chi-squared, Andersen makes use of the two following features:

1. With a moderate or high non-centrality parameter a non-central chi-square can be represented by a power function applied to a Gaussian variable. Hence, for sufficiently large values of Vt, one can write

        Vt+∆ = a(b + ZV )^2
                     
where a and b are constants depending on Vt and ∆, and ZV is a standard Gaussian random variable.

2. For low values of Vt, one can approximate the density of Vt+∆ by
        
        Vt+∆ = L^(−1)(u) =   0                          if 0 ≤ u < p
                             β^(-1)*ln((1−p)/(1-u) )    if p < u ≤ 1 
                     
Then, sampling is done by Vt+∆ = L−1(UV ), where u is a uniform random number.

First, to simplify the outline it is convenient to introduce the following variables:
    
    m = θ + (Vt − θ)e^(−κ∆)
    s2 = Vt*η^2 * e^(−κ∆)(1 − e^(−κ∆))/k  + θη^(2)/2k * (1 − e^(−κ∆))^2
    ψ = s^2/m^2

Further, with the aid of moment-matching techniques Andersen shows that:
    
    b2 = 2ψ^(−1) − 1 + 2(2ψ^(−1))^0.5 * (2ψ^(−1) − 1)^0.5
    a = m/(1 + b2)
    p = (ψ − 1)/(ψ + 1)
    β = 1 − p
    m = 2/m(1 + ψ)

With these constants defined, for some critical level ψc ∈ [1,2] ( set at 1.5 in Andersen )  we can formulate the QE algorithm for the variance process, which is what is found in the function *variancenew*

        1. Given Vt, compute m and s2 and ψ = s2/m2
        2. if ψ ≤ ψc :
        (a) Compute a and b from
        (b) Generate a standard Gaussian random number, ZV
        (c) Set Vt+∆ = a(b + ZV )2
        3. else, if ψ > ψc :
        (a) Compute β and p from (4.22)
        (b) Generate a uniform random number, UV
        (c) Set Vt+∆ = L−1(UV )

According to Andersen the exact choice of ψc is not of great importance. We later use ψc = 1.5 in our numerical tests. Because of the desirable correlation between the asset price and the variance, an Euler-like X-process leads to flaws, if combined with the QE scheme for the variance, see. Instead, one is suggested to look at the properties of the exact representation of X, derived in. The following scheme is proposed:

        Xt+∆ = Xt + r∆ + K0 + K1Vt + K2*Vt+∆ + (K3Vt + K4Vt+∆)^0.5 *Z 
        
where Z is a standard Gaussian random variable, not correlated with the random variables in the Variance-process, and
        
        K0 = −ρκθ∆/η
        K1 = γ1*∆*(κρ/η − 0.5)− ρ/η
        K2 = γ2∆(κρ/η−0.5)+ ρ/η
        K3 = γ1*∆*(1 − ρ^2)
        K4 = γ2*∆(1 − ρ^2), γ1,γ2 ∈ [0,1]   
        Where γ1,γ2 are the "averaging" coefficient to calculate the mean Variance between t and t+∆
        
# Simulation Engine

The simulation function runs simulations while periodically checking for convergence  accuracy level (more on this in the code notes).
Moreover, it saves ALL the paths generated and some statistic about each path ( such as min, max price which are useful in case of sigle barrier pricing ). Note that this feature makes the algorithm computationaly slow and it can be easily avoid for Barrier Pricing, nevertheless it is usuful during tests, hence is left. 

Once the accuracy level is reached or the usuer stop the simulation the fucntion also generates many plots to better visualize what went on during the simulation. Some of these plots are base on the paths data ( see above ) that were collected and saved. Again, not useful for the mere sake of the simualtion and easily optimizable. 



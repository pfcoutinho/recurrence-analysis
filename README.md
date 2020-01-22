# Recurrence Analysis

The files in this repository can be used to perform recurrence analysis of time series in MATLAB. With them you can generate the distance matrix (DM) and the recurrence plot (RP). Soon, you will also be able to carry out recurrence quantification analysis (RQA) and recurrence network analysis (RNA).

## Recurrence Plot

Here are some examples of distance matrices and recurrence plots, which are obtained from the distance matrices threshold after the application of the threshold. See how beautiful and informative these binary matrices are!

### Stochastic process

![Distance matrix and recurrence plot of a sequence drawn from the uniform distribution](imgs/img1.png)

![The sequence drawn from the uniform distribution](imgs/img1_seq.png)

**Figure 1** _Distance matrix and recurrence plot of a sequence drawn from the uniform distribution_ U(0,1)_. In this example, the embedding dimension = 1 and the threshold = 0.2._

### Sinusoidal function

![Distance matrix and recurrence plot of a sinusoidal function](imgs/img2.png)

![The sinusoidal function](imgs/img2_seq.png)

**Figure 2** _Distance matrix and recurrence plot of a sinusoidal function given by_ x(t) = sin(2\*pi\*0.05\*t)_. In this example, the embedding dimension = 2, the time delay = 3, and the threshold = 0.2._

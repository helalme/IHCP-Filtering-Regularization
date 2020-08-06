### Filtering and regularization of inverse heat conduction problems (IHCP)

#### Background:
IHCPs are severely ill-posed having the problem with uniqueness and stability. Whatever the techniques (regularization, optimization) we use, there
are always the concern of at least highly noisy solution. An algorithm was developed in 2012 to solve 3D IHCP where boundary condition is missing, more
specifically, temperature history at the source is missing. For example, outside surface temperature history of a pipe can be measured, but if we do not
have the possibility to access inside wall of the pipe, we would not know inner wall temperature history. In such case the developed algorithm can be applied.
The detail of the algorithm and a software module written in C++ can be found by following the link: <https://github.com/helalme/TAM>  

#### Problem definition:
Computed solution either from the developed algorithm or using a regularization (e.g. Tikhonov regularization) method gives very noisy solution. To have the
accurate estimate of the inner wall temperature history (missing BC) we need use filtering techniques such as moving average filter, Kalman filter etc. In this
repository, related Matlab scripts (algorithm for solving IHCP, adding noise/perturbation to data, Tikhonov regularization, Kalman filter) are provided along with
necessary input/output data files. 

#### Results:
A simulated dataset (inner wall and outer surface temperature history) is the starting point of this problem, where outer wall temperature history has been
perturbed with a very small range of [-0.05, 0.05], in order to realize high noised in the solution of inner wall temperature history. Then the noisy solution
has been gone through moving average (MA) filter and Kalman filer. We see that there are errors in the computed solution upto 26% (80°C). If we use MA as post filter,
the error drastically goes down from 80°C to less than 5°C. If we use MA as pre and post filter, it comes down to less than 4°C. Interestingly, if we use a
Kalman filter with only two iteration, the error goes down less than 3°C, on average it comes down to less than 0.5%. Results and graphs have been shown in the
PDF file *FilteringEffectOnSolutions.pdf*.   
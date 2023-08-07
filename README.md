# Mathematical modeling of polyA rhythm in mammalian cells
The original work is published in *"Critical role of deadenylation in regulating poly(A) rhythms and circadian gene expression"*
[X. Yao, S. Kojima and J. Chen, PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007842)


## ODE system and numerical simulation
We used Satelli's (AB) sampling scheme for the Sobol indices calculation in the paper. Sample size N = 100,000; 10 independent repeats are simualted
for average and std. Latin hyper cube method is used for parameter sampling from their assumed distributions.


## Analytical approximation of ODE system solution
Symbolic equations in MATLAB are used to approximate the solutions of the polyA tail length ODE system. Check SI of the paper for more info.

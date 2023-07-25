# twosampleDPM

R codes for two-sample testing with DPM-Poisson-Gamma model, accompanying the paper: **Nonparametric Bayes multiresolution testing for high-dimensional rare events** by Jyotishka Datta, Sayantan Banerjee , and David B. Dunson. 


Please see the R code `2gp_simulation_demo.r` for running the numerical example in section 3 of the paper with 35 points generated from the two non-homogeneous Poisson process. 

$$
\begin{align*}
\lambda_1(x) & = 2\exp \left(-\{(x-50)/10 \}^2 \right) +20\exp\left(-\{(x-10)/10\}^2 \right); x \in \mathbb{R}^+ \\
\lambda_2(x) & = 20\exp \left(-\{(x-50)/10 \}^2 \right) +2\exp\left(-\{(x-10)/10\}^2 \right); x \in \mathbb{R}^+,
\end{align*} 
$$







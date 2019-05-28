---
title: Initial results
format: beamer
...

Species abundance distributions
---------------------------
\includegraphics[width=\linewidth]{figs/species_abundance_curve_0.pdf}

Species abundance distributions
---------------------------
\includegraphics[width=\linewidth]{figs/species_abundance_curve_4.pdf}

Species abundance distributions with lognormal curve
---------------------------
\includegraphics[width=\linewidth]{figs/species_abundance_curve_5.pdf}

Shannon diversity
-----------------
\Large The Shannon diversity $S$ is defined $$S = -\sum_{i=1}^N p_i \ln(p_i).$$

\normalsize
$\implies$ maximum for $y_i = 1/N$ for all $i$

$\implies$ minimum for $y_i = 1$, $y_j = 0$ for $j \neq i$

(one of many metrics for "diversity" of a population)

Shannon diversity 
---------------------------
\includegraphics[width=\linewidth]{figs/shannon_div_0}

Shannon diversity with lognormal distribution
---------------------------
\includegraphics[width=\linewidth]{figs/shannon_div_1}

Unifrac\footnote{image from mothur.org/w/images/5/5b/UnweightedUniFracMeasure.jpg}
----------------
\includegraphics[width=\linewidth]{figs/UnweightedUniFracMeasure.jpg}

Unifrac distance matrix (weighted)
---------------------------
\centering
\includegraphics[width=.9\linewidth]{figs/unifrac_distance_matrix}

Unifrac distance matrix (unweighted)
---------------------------
\centering
\includegraphics[width=.9\linewidth]{figs/unweighted_distance_matrix}

$L_2$-norm distance matrix (weighted)
---------------------------
\centering
\includegraphics[width=.9\linewidth]{figs/l2_distance_matrix}

Future work?
------------
\Large
+ Estimate size of available microbiome pool (species abundance curve)
+ Coarse-grain at different taxonomic resolutions
+ Include metadata analyses
+ How can we address nestedness?

# Sampling from the Solution Space of Genome-Scale Metabolic Models

In this repo, we provide implementations on the flux sampling guidelines and best practices we describe on 
our chapter on the *“Flux Balance Analysis”* book, in the protocol series Methods in Molecular Biology, by Springer Nature.

This repository provides implementations of the flux sampling guidelines and best practices described in our chapter on _“Flux Balance Analysis”_, 
to be published in the [Methods in Molecular Biology protocol series](https://link.springer.com/series/7651) by Springer Nature.

We discuss flux sampling implementation under different scenarios and highlight some of its challenges. 
We apply sampling both within the cell, making use of the [`dingo`](https://github.com/geomScale/dingo) Python library [^1], and on the extracellular space, using the [MAMBO](./scripts/mambo.py) approach[^2]. 

For any trouble, feel free to [open an issue](https://github.com/hariszaf/metabolic_toy_model/issues) specifying that you are using the `sampling` branch of this repo.


## References
[^1]: Apostolos Chalkis, Vissarion Fisikopoulos, Elias Tsigaridas, Haris Zafeiropoulos, dingo: a Python package for metabolic flux sampling, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae037, https://doi.org/10.1093/bioadv/vbae037
[^2]: Garza, D.R., van Verk, M.C., Huynen, M.A. et al. Towards predicting the environmental metabolome from metagenomics with a mechanistic model. Nat Microbiol 3, 456–460 (2018). https://doi.org/10.1038/s41564-018-0124-8
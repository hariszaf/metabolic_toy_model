# Sampling from the Solution Space of Genome-Scale Metabolic Models


Hello friend. 

In this branch you will find the material for the tutorial on [`microbetag`](https://microbetag.readthedocs.io/)
in the framework of the [*Data integration with Microbial Networks and Community Models*](http://msysbiology.com/microbialdataintegration/) summer school. 

If you are reading this branch before the summer school and see anything funny, keep in mind that is a "living" branch, but also feel free to contact us, so we fix it before the school!

You may see [here](#contact) how to contact us! 

## Agenda

15:45 - 17:45


<table>
  <thead>
    <tr>
      <th style="text-align:center;">Time</th>
      <th style="text-align:center;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center" ><nobr>15:45-16:00<nobr></td>
      <td>Setup the scene</td>
    </tr>
    <tr>
      <td align="center" ><nobr>16:00-16:30<nobr></td>
      <td> Introduction in flux sampling: an (unbiased) constraint-based approach </td>
    </tr>
    <tr>
      <td align="center" >16:30-16:50</td>
      <td> <b>Hands-on (part A):</b> Sampling tasks using <code>dingo</code>, <code>OptGP</code> and <code>CHRR</code> implementations </td>
    </tr>
    <tr>
      <td colspan="2" align="center"><b>☕   BREAK (16:50-17:00)   ☕ </b></td>
    </tr>
    <tr>
      <td align="center" >17:00-17:15</td>
      <td> Leveraging flux sampling statistical power </td>
    </tr>
    <tr>
      <td align="center" >17:15-17:35</td>
      <td> <b>Hands-on (part B):</b> Statistical exploration of flux sampling data using <code>dingo-stats</code> </td>
    </tr>
    <tr>
      <td align="center" >16:45-17:00</td>
      <td>Q & A</td>
    </tr>
  </tbody>
</table>

## Material

- [slides]()
- [notebook](./flux_sampling.ipynb)


## How to work with this repo 

For the needs of this tutorial, make sure you have: 

 - a [GitHub](https://github.com/) account 


To run the notebook of this tutorial, you can either: 

  - [fire a GitHub codespace](https://github.com/codespaces/), checking out to the `ai-school25` branch, or
  - clone the repo on your IFB core account and jump into the `ai-school25` branch
  - get repo and rependencies to perform the tutorial locally (instructions for Linux systems)

For example:

    git clone https://github.com/hariszaf/metabolic_toy_model.git
    cd metabolic_toy_model
    git checkout kul2025


Make sure you have conda or miniconda available on your computing environment. 
GitHub codespace brings it by default, but if you prefer running this tutorial locally, and you don't have conda, 
then you may follow the instructions [here](https://www.anaconda.com/docs/getting-started/miniconda/install).

You can now jump to the [`microbetag_tutorial.ipynb`](./microbetag_tutorial.ipynb) notebook! 



## Contact

Feel free to [open an issue](https://github.com/hariszaf/metabolic_toy_model/issues) specifying the branch you are referring to. 

We also welcome the most your contributions! 
You can contribute to any of our events by following the guidelines [here](https://dev.to/javigong/how-to-contribute-to-an-open-source-project-on-github-1hbo). 
Just make sure to check out the appropriate branch before getting started.




























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
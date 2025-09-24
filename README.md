# `microbetag`: metabolic secrets behind microbial co-occurrence

Hello friend. 

In this branch you will find the material for the tutorial on [`microbetag`](https://microbetag.readthedocs.io/)
in the framework of the [*Data integration with Microbial Networks and Community Models*](http://msysbiology.com/microbialdataintegration/) summer school. 

If you are reading this branch before the summer school and see anything funny, keep in mind that is a "living" branch, but also feel free to contact us, so we fix it before the school!

You may see [here](#contact) how to contact us! 

## Agenda
<!-- 
|     Time      |                           Description                            |
|:-------------:|:----------------------------------------------------------------:|
|                         Part A (13:30 - 15:00)                                  |
| 13:30 - 13:50 | Setup the scene  |
| 13:50 - 14:00 | Reading the last page first |
| 13:00 - 14:20 | Build a co-occurrence network — a first example of running `microbetag` partially. | 
| 14:20 - 14:40 | Annotating nodes: literature and genome-derived phenotypic traits |
| 14:40 - 15:00 | Annotating edges (part A): pathway complmenetarity |
|                         Break (15:00 - 15:30)                                  |
|                         Part B (15:30 - 17:00)                                  |
| 15:30 - 16:00 | Annotating edges (part B): seed complementarity |
| 16:00 - 16:45 | Two real-world example cases |
 -->


<table>
  <thead>
    <tr>
      <th style="text-align:center;">Time</th>
      <th style="text-align:center;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center" ><nobr>13:30-13:50<nobr></td>
      <td>Setup the scene</td>
    </tr>
    <tr>
      <td align="center" ><nobr>13:50-14:00<nobr></td>
      <td>Reading the last page first</td>
    </tr>
    <tr>
      <td align="center" >14:00-14:20</td>
      <td>Build a co-occurrence network — a first example of running <code>microbetag</code> partially</td>
    </tr>
    <tr>
      <td align="center" >14:20-14:40</td>
      <td>Annotating nodes: literature and genome-derived phenotypic traits</td>
    </tr>
    <tr>
      <td align="center" >14:40-15:00</td>
      <td>Annotating edges (part A): pathway complmenetarity</td>
    </tr>
    <tr>
      <td colspan="2" align="center"><b>☕   BREAK (15:00-15:30)   ☕ </b></td>
    </tr>
    <tr>
      <td align="center" >15:30-16:00</td>
      <td>Annotating edges (part B): seed complementarity</td>
    </tr>
    <tr>
      <td align="center" >16:00-16:45</td>
      <td>Two real-world example cases and an on-the-fly one: which one you ll go for?</td>
    </tr>
    <tr>
      <td align="center" >16:45-17:00</td>
      <td>Q & A</td>
    </tr>
  </tbody>
</table>

## Material

- [slides](https://docs.google.com/presentation/d/15WvhB9Vff3fWYVFUMFNaGt8xiR1J4lhcKyq7AeR3YVU/edit?usp=sharing)
- [notebook](./microbetag_tutorial.ipynb)


## How to work with this repo 

For the needs of this tutorial, make sure you have: 

 - a [GitHub](https://github.com/) account 
 - [Cytoscape](https://cytoscape.org/download.html) ($\geq$ 3.9)
 - the [MGG](https://apps.cytoscape.org/apps/mgg) add-on


To run the notebook of this tutorial, you can either: 

  - [fire a GitHub codespace](https://github.com/codespaces/), checking out to the `kul2025` branch, ot
  - clone the repo locally and build a `conda` environment as described in the corresponding branch.

For example:

    git clone https://github.com/hariszaf/metabolic_toy_model.git
    cd metabolic_toy_model
    git checkout kul2025


Make sure you have conda or miniconda available on your computing environment. 
GitHub codespace brings it by default, but if you prefer running this tutorial locally, and you don't have conda, 
then you may follow the instructions [here](https://www.anaconda.com/docs/getting-started/miniconda/install).

You can now jump to the [`microbetag_tutorial.ipynb`](./microbetag_tutorial.ipynb) notebook! 

## About codespaces

Based on [GitHub's documentation](https://docs.github.com/en/codespaces/developing-in-a-codespace/deleting-a-codespace):

GitHub Codespaces are automatically deleted after they have been stopped and have remained inactive for a defined number of days. 
The retention period for each codespace is set when the codespace is created and does not change. The default retention period is 30 days. 



## Contact

Feel free to [open an issue](https://github.com/hariszaf/metabolic_toy_model/issues) specifying the branch you are referring to. 

We also welcome the most your contributions! 
You can contribute to any of our events by following the guidelines [here](https://dev.to/javigong/how-to-contribute-to-an-open-source-project-on-github-1hbo). 
Just make sure to check out the appropriate branch before getting started.

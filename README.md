# Metabolic modeling playground

## About this repo

This repository contains material designed to teach, explain, and showcase the fundamentals and methods of metabolic modeling.

We organize workshops and classes, each associated with its own dedicated branch. The table below lists the branches and their corresponding events held so far.

|                                 Branch name                                 |                          Description                      |
| :-------------------------------------------------------------------------: | :-------------------------------------------------------: |
| [`antony25`](https://github.com/hariszaf/metabolic_toy_model/tree/antony25) |                                  ["Metabolic models applied to microbiomes" workshop](https://metabolicmodelingantony2025.onrender.com/) @ INRAE/PROSE, Antony, France                                                                        |
|     [`duth`](https://github.com/hariszaf/metabolic_toy_model/tree/duth)     | ["Introduction to metabolic modeling" workshop](https://docs.google.com/presentation/d/1w0fhaz9G74UtEp7qEqdKYbYJpboj_SXjU-J2IxFrlhs/edit?usp=sharing) at Master in Biomedical Informatics, DUTH, Greece                                          |
| [`sampling`](https://github.com/hariszaf/metabolic_toy_model/tree/sampling) |                                                             ["Sampling from the Solution Space of Genome-Scale Metabolic Models" chapter]()                                                                            |


In the `main` branch, you'll find two key folders: [`scripts/`](./scripts/) and [`files`](./files/). These serve as the ground base of the repository, which is why they are included here.
They support the reconstruction of three human gut–related metabolic toy models and provide resources to work with two of the most widely used namespaces in metabolic modeling: [ModelSEED](https://github.com/ModelSEED/ModelSEEDDatabase) and [BiGG](http://bigg.ucsd.edu/).


## How to work with this repo 

No matter which event you wish to go for, you can either [fire a GitHub codespace](https://github.com/codespaces/) using the branch of your choice to be checked out on creation, 
or clone the repo locally and build a `conda` environment as described in the corresponding branch.

For example:

    ```bash
    git clone https://github.com/hariszaf/metabolic_toy_model.git
    cd metabolic_toy_model
    git checkout duth
    ```
Now, you could follow the instructions on [prep_env.ipynb](https://github.com/hariszaf/metabolic_toy_model/blob/duth/prep_env.ipynb) for how to build your local `conda` environment.

Since each branch has its own goals, it also comes with its own set of requirements.


## Contact

Feel free to [open an issue](https://github.com/hariszaf/metabolic_toy_model/issues) specifying the branch you are referring to. 

We also welcome the most your contributions! 
You can contribute to any of our events by following the guidelines [here](https://dev.to/javigong/how-to-contribute-to-an-open-source-project-on-github-1hbo). 
Just make sure to check out the appropriate branch before getting started.

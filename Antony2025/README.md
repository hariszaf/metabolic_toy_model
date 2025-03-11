# Metabolic models on microbiomes 

[Web-site](https://metabolicmodelingantony2025.onrender.com)

Organizers: 
* [Daniel Garza](https://danielriosgarza.github.io)
* [Meine Boer](https://www.nioz.nl/en/about/organisation/staff/meine-boer)
* [Haris Zafeiropoulos](https://hariszaf.github.io) 


## Program

### Day 1

#### Afternoon session

1. [Set up your coding environment for this workshop](./preparingYourEnvironment.ipynb)

2. [Introduction to `cobra`](./introductionToCOBRApy.ipynb)

3. [Reconstruction of a draft genome-scale metabolic model](./reconstructingDraftGSMMs.ipynb)



### Day 2

#### Morning session 

4. [Gap-filling a GEM with `DNNGIOR`](./gapfillingGSMMs.ipynb)

5. [Training `DNNGIOR` with your own reaction set](./DNNGIORtraining.ipynb)

6. [Assess the quality of your model](./modelQuality.ipynb)

#### Afternoon session


7. [Metabolic model analysis methods](./computationalMethods.ipynb)


#### Applications and project-based discussion

8. [Looking for *"partners"* in gut colonization](./gutColonization.ipynb) and [estimating the impact of methanol perturbation on anaerobic digestion microcosms](./methanolImpact.ipynb)


## Day 3 

Open talks
<!-- TODO (Haris Zafeiropoulos, 2025-03-11): 
Would be nice to have the presentations either on this repo or in a Google folder or something. 
-->


## Easy way to get all requirements 


In case you run the workshop locally (and not on Google Collab) you can get all the dependencies for your `conda` environment, simply by:

```bash 
# In case you have not already the repo locally
git clone https://github.com/hariszaf/metabolic_toy_model.git

cd metabolic_toy_model/Antony2025

conda create -n gsmmWorkshop python=3.11.11

conda activate gsmmWorkshop

pip install -r requirements.txt
```


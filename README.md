# Intro in metabolic modeling

## Topics to be discussed

We will cover some first topics regarding the metabolic modeling reconstruction and
some typical methods for their constraint-based analysis (CBA). 

Some concepts of CBA will be covered too. 



## How to work with this repo 

You can either clone this repo by: 
    ```bash
    git clone https://github.com/hariszaf/metabolic_toy_model.git
    ```

Or, download it by just clicking on the `Code > Download ZIP`.

Then, from within the root folder of the repo, initiate a `conda` environment and install the required libraries by running:

    conda create -n met_class python=3.10 --no-default-packages
    conda activate met_class

Please, make sure you have `cmake` and `swig`, depending your OS

- for mac
    brew install cmake swig

- for linux
  sudo apt-get install cmake swig

- for windows 
  

    python3 -m pip install -r requirements.txt

Once everything is set, you may run 

    jupyter notebook 

and jump into the `msc_class.ipynb` notebook!





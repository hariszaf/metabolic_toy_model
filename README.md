# Intro in metabolic modeling

## Topics to be discussed

We will cover some first topics regarding the metabolic modeling reconstruction and
some typical methods for their constraint-based analysis (CBA). 

Some concepts of CBA will be covered too. 

## How to work with this repo 

You can either clone this repo by: 

    ```bash
    git clone https://github.com/hariszaf/metabolic_toy_model.git
    cd metabolic_toy_model
    ```

Or, download it by just clicking on the `Code > Download ZIP`.

Then, from within the root folder of the repo, initiate a `conda` environment and install the required libraries by running:

    conda create -n met_class python=3.10 --no-default-packages
    conda activate met_class

Make sure you have `cmake`. If not already available, depending on your OS, run the following to get it:

- for **macOS**
  
    brew install cmake

- for **Linux**

  sudo apt-get install cmake

- In case of Windows I guess you are using a Linux a Windows Subsystem for Linux (WSL) or something similar, so it would be as the Linux case. 

Once you make sure `cmake` is there, you may run:

    python3 -m pip install -r requirements.txt


Once everything is set, you may run:

    jupyter notebook 

and jump into the `msc_class.ipynb` notebook!


For any trouble, feel free to contact [Haris](mailto:haris.zafeiropoulos@kuleuven.be) directly.



# panelScope

Multi-view gene panel characterization for spatially resolved omics

![](https://github.com/SydneyBioX/panelScope/blob/main/figure1.png)
We present panelScope, a framework based on a diverse collection of metrics to characterize a gene panel, allowing researchers to determine whether a panel is well-suited to their study’s objectives. We demonstrate the utility of panelScope by generating multiple-views of gene panels that describe their ability to capture cell types of interest, enrichment for biological pathways, or the amount of redundant information. In parallel, we leverage these metrics as loss functions in a genetic algorithm for panel design, where users can choose to weight each characterization category. Importantly, we have implemented this framework in an interactive web platform, which includes a library of pre-existing gene panels that users can compare to their own gene panels. Thus, by quantitatively summarizing a panel from multiple views, panelScope enables the design of panels that can capture diverse information relevant to one’s specific research questions. 



## Table of Contents

* [Installation](#Installation)
* [Usage](#Usage)
* [Citation](#Citation)


## Installation

To correctly use **panelScope** via your local device, we suggest first create a conda environment by:

~~~shell
conda create -n <env> python=3.9
conda activate <env>
~~~

After entering an entire new virtue environment, the recommended way to install depending packages are:

~~~shell
conda install pandas
conda install scanpy
~~~

## Usage

There are two main parts implemented within panelScope, the metrics-computation part and the optimization part.

The metrics-computation part is written by R, which describes numerical properties of a selected gene panel based on provided single cell dataset.

The optimization part is written by Python using an efficient version of evolution algorithm. As shown in the following figure, an iterative optimization process containing 4 steps is evolved. At the initialization step, a number of gene panels (default as 50) are randomly selected from the input dataset as initial population. Then, an evaluation-selection-generation loop would be repeated for a number of times (default as 5000). We use scoring functions related to feature diversity, pathway diversity, panel entropy, spatial specificity and variation recovery as our evaluation metrics. The algorithm would randomly pick better performed gene panels and use them to generate new panels to replace the old, not-well-performed ones.

![](https://github.com/SydneyBioX/panelScope/blob/main/figure2.png)

 It has 4 arguements, respectively are:

~~~
--dataset_path: should be a str. The absolute path of the reference single cell dataset within .rds format.

--panel_num: should be an int. This pre-defined number indicates the amount of genes within the final panel.

--search_space_path: should be a str. The absolute path of a json file which contains all candidate genes.

--objmode: should be a str. This str should be a choice in ["overall","cts","corr","pathway","spatial","tv"].
          These represent the objective functions: "Overall", "Panel entropy score", "Feature diversity score",
          "Pathway diversity score", "Moran’s I" and "transcriptional-variability-based functions". The users
          could optimize the gene panel by anyone of them. One step further, the users could also define their
          own functions and replace the R codes for more speficial panel selection.
~~~

To experience our algorithm, you can try with our demo data by:
~~~
python main.py --dataset_path ./demo_hca_10x2.rds --panel_num 100 --search_space_path ./search.json --objmode overall
~~~



## Citation

If you find our codes useful, please consider citing our work:

~~~bibtex
@article{panelScope,
  title={Multi-view gene panel characterization for spatially resolved omics},
  author={Daniel Kim+, Wenze Ding+, Akira Nguyen Shaw, Marni Torkel, Cameron J Turtle, Pengyi Yang, Jean Yang*
},
  journal={},
  year={2025},
}
~~~

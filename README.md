# Introduction

This package contains the txtlsim toolbox, which can be used to simulate the chemical reations that occur in the Transcription-Translation (TX-TL) cell-free gene expression system developed at the University of Minnesota (Vincent Noireaux) and Caltech (Richard Murray). The main features of this toolbox are its ability to track the loading of enzymatic machinery, the consumption of resources, the ease of setting up models, the automatic accounting of retroactivity effects, and the extensibility the reaction networks generated. 
<!-- 
The second toolbox is called mcmc_simbio. This is a concurrent Bayesian parameter inference toolbox for MATLAB Simbiology models. The Bayesian parameter inference is performed using a modification of Aslak Grinsted's MATLAB implementation of the affine invariant ensemble Metropolis-Hastings MCMC sampler (ref). We have added support for what we call 'concurrent' parameter inference, which refers to the capability to estimate a common set of parameters that get used simultaneously and in arbitrary combinations in multiple experiments/models. More information can be found below.  -->

# The Toolbox

## Getting the toolbox and running some simple examples

Clone the repository using 

```
git clone https://github.com/vipulsinghal02/txtlsim_buildacell.git
```

into a directory you wish to put the toolbox in. Alternatively, if you do not plan to version control the toolbox, you can simply download it as a zip file using the green button on the [main page](https://github.com/vipulsinghal02/txtlsim_buildacell). 

Lets call the directory where you cloned or downloaded the repository `trunk`, i.e., this is where directories like `core`, `components`, `examples` and `mcmc_simbio` are. Open MATLAB, and set the current working directory to `trunk`. Type in `txtl_init` into the MATLAB command line. This initializes the txtlsim toolbox. 

Check if the toolbox is installed properly by running the follwoing examples:

### txtlsim examples
Start with typing in `geneexpr` into the command window. This should run a constitutive gene expression example in the toolbox, and you should see a plot with three subplots (protein; mRNA and DNA; and resource usage) appear. There should not be any error messages in the command window. 

Next, run the `negautoreg` example in the command line. Again, you should see a plot of the species in the system, and no errors. This example simulates the negative autoregulation circuit. 

Type in `edit TXTL_tutorial` into the command line to open the tutorial file. You can also use the MATLAB publisher button to publish this file, and look at it in the MATLAB help file markup. We recommend running through the examples in this file, and exploring the reactions, species, etc set up in this file. Familiarity with the MATLAB Simbiology command line is helpful here. To learn more about Simbiology, go to the [Getting Started Using the Simbiology Command Line](https://www.mathworks.com/help/simbio/gs/simbiology-command-line-tutorial.html) page. 
<!-- 
### mcmc_simbio examples

Next, open and explore the mcmc_simbio estimation examples given in the files `proj_mcmc_tutorial`, `proj_mcmc_tutorial_II`, and `proj_mcmc_tutorial_III` in the `trunk\mcmc_simbio\proj\` directory. We strongly recommend you skim through the `mcmc_info.m` and `data_info` files (`trunk\mcmc_simbio\models_and_supporting_files\` or type `help mcmc_info` and `help data_info` into the MATLAB command line) to gain an understanding of some of the key functionalities of the parameter inference toolbox. Along with the three tutorial files, the `mcmc_info.m` and the `data_info.m` files provide an initial idea of the capabilities of the toolbox.  -->

## References

More information can be found in the following references: 

V. Singhal, Z. A. Tuza, Z. Z. Sun and R. M. Murray, (2020). A MATLAB Toolbox for Modeling the Behavior of Genetic Circuits in Cell-Free Transcription-Translation Systems. In preparation. 

Z. A. Tuza, V. Singhal, J. Kim and R. M. Murray, "An in silico modeling toolbox for rapid prototyping of circuits in a biomolecular “breadboard” system," 52nd IEEE Conference on Decision and Control, Firenze, 2013, pp. 1404-1410.
doi: 10.1109/CDC.2013.6760079


Vipul Singhal, 2018
California Institute of Technology # txtlsim_buildacell

QARTOD Testing
==============================
<!-- [![Build Status](https://github.com/@cooleyky/qartod_testing/workflows/Tests/badge.svg)](https://github.com/@cooleyky/qartod_testing/actions) -->
[![codecov](https://codecov.io/gh/@cooleyky/qartod_testing/branch/main/graph/badge.svg)](https://codecov.io/gh/@cooleyky/qartod_testing)
[![License:MIT](https://img.shields.io/badge/License-MIT-lightgray.svg?style=flt-square)](https://opensource.org/licenses/MIT)
[![pypi](https://img.shields.io/pypi/v/qartod_testing.svg)](https://pypi.org/project/qartod_testing)
<!-- [![conda-forge](https://img.shields.io/conda/dn/conda-forge/qartod_testing?label=conda-forge)](https://anaconda.org/conda-forge/qartod_testing) -->
[![Documentation Status](https://readthedocs.org/projects/qartod_testing/badge/?version=latest)](https://qartod_testing.readthedocs.io/en/latest/?badge=latest)


## Purpose
Scripts to check performance of QARTOD tests on OOI data in both production and development. This includes providing statistics illustrating the composition of QARTOD flags for a given dataset, running QC tests on datasets locally for comparison with expected results, and perform data deep dives in preparation for new QARTOD tests.

## Glossary
Coming soon!

## Project Organization
Since this repository hosts scripts and notebooks for a few different tasks all related to development of QARTOD tests and quantifying their performance, the organization of this repository follows some rules to help clarify which parts are used in each task.

### notebooks
All notebooks are numbered at the beginning of the file name to indicate the task and steps through which data was processed. The first digit will be the same for notebooks that are for the same task, and the second digit indicates the order that the notebooks should be run to arrive at the same results starting from step "1".

So far the numbering of notebooks for different tasks is as follows:
<ul><li><strong>0#: </strong> QARTOD test flag statistics and local test comparison for QC tests in production</li>
<li><strong>1#: </strong> CGSN SPKIR (downwelling spectral irradiance) data deep dive</li>
<li><strong>2#: </strong> Assessment of CGSN ADCP data for QARTOD planning</li>
<li><strong>4#: </strong> Find CGSN QARTOD tests not executed</li></ul>

### qartod_testing
This is the directory where the project source code lives. Each of the individual libraries within this directory is focused on a single task (or even half a task in the case of QC flag statistics and running a QC test locally).  

### data
As you use the notebooks, you will see that any data downloaded from OOINet and subsequently processed is saved in these folders. Git will not track the contents of these folders but they are here as placeholders.


## Setup
Working with the notebooks and modules contained in this repository requires installation of Git and Python on your machine of choice. The miniconda distribution ([link to installer download](https://docs.conda.io/en/latest/miniconda.html)) is sufficient to get started and progress through the rest of the setup below, which assumes that the environment will initially be setup using the `conda` library. On Windows, you will need to choose the option that adds Anaconda to the PATH environment variable during the miniconda installation. More details on how to do this are available from [this tutorial](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/setup-git-bash-conda/) from [Earth Lab](https://www.earthdatascience.org/).

### Create and activate qartod_test environment [^1]
From a Git Bash terminal in the `qartod_testing` working directory, you'll setup the python environment with `conda` and install the source code as a local development package by following these steps:

    # configure the OOI python environment
    conda env create -f environment.yml
    conda init # might be required for windows users if environment is not active
    conda activate qartod_test

    # you can check the active environment by running
    conda env list

    # add qartod_testing module to qartod_test environment path
    conda develop .

Next, we'll add two other OOI-specific modules from GitHub to the qartod_test environment. 

### Using OOINet and ooi-data-explorations libraries locally
To use the notebooks in this repository as-is, users will need to fork the [OOINet](https://github.com/reedan88/OOINet) and [ooi-data-explorations](https://github.com/oceanobservatories/ooi-data-explorations/tree/master/python) GitHub repositories and clone these to their local GitHub directory. This whole process can be done from the GitHub Desktop program using the included links to each repository and is explained [here](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop). 

The following commands will add the paths to ooi_data_explorations and OOINet modules to the qartod_test environment path:
    
    # Check that the qartod_test env is still starred
    conda env list
    
    # Navigate to the ooi-data-explorations directory
    cd path/to/ooi-data-explorations/python

    # install the ooi_data_explorations package as a local development package in qartod_test environment
    conda develop .

    # Navigate to the OOINet directory
    cd path/to/OOINet

    # install the ooinet package as a local development package in qartod_test environment
    conda develop .

Where `.` denotes the current working directory.

### Setup credentials for M2M data requests
In your home directory, you can save your M2M credentials in a text file named `.netrc`. The OOINet M2M module uses the netrc library to retrieve credentials from this file so that they don't have to be entered manually for every request. The text in this file should be formatted as follows:
    
       machine ooinet.oceanobservatories.org
           login <YOUR-OOI-API-USERNAME>
           password <YOUR-API-TOKEN>
           
If you have access to the development platform and intend to work with data from that server as well, you will need another similar entry with the server web address in place of ooinet.oceanobservatories.org. More detailed information on setting up this file can be found here: https://github.com/oceanobservatories/ooi-data-explorations/tree/master/python#access-credentials

### Working from a branch in your forked repository
Whether you're planning on just trying out the included notebooks in this repository or adding a feature that you will want to push to this project later, working within a uniquely-named branch on your fork will save you from headache and confusion down the road. This can be done from either the terminal or in the GitHub Desktop program. I've found that the most useful branch names are brief, descriptive names of the feature I want to add or change, or I will use an overall goal for working with this respository.

### Using your qartod_test environment in Jupyter Notebooks
Now that the environment is set up and you have a specific branch to work from, you'll need to choose a kernel when you run a notebook for the first time. You'll choose the kernel that is named the same as the qartod_test environment. This way, all modules that will be imported into the notebooks in this repository are already installed in the environment being used. If the `qartod_test` environment is not automagically added as a kernel option, you can add it from a terminal with the following command: `python -m ipykernel install --user --name=qartod_test`

### A tip for managing your environment
If you need to add multiple packages to the environment after you have already created it, using `conda install` can quickly result in a "broken" environment if conda cannot resolve incompatible packages. The first step is making sure that you add new dependencies to environment.yml. Next, you can either install individual packages with pip or update the environment from the edited environment file with conda via `conda env update -f environment.yml` once the qartod_test environment is activated.

--------
[^1]: These directions are modified from https://github.com/oceanobservatories/ooi-data-explorations/tree/master/python#obtaining-the-code-and-configuring-the-environment

<p><small>Project based on the <a target="_blank" href="https://github.com/jbusecke/cookiecutter-science-project">cookiecutter science project template</a>.</small></p>

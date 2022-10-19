# Pan-cancer analysis of beyondcell drug signatures

This repo hold the code of the aforementioned project. For collaborators: start
here to setup your dev. environment.

# Setup the environment for code and data sharing

## Get the code

You'll need to keep a local copy of the  preprocessing/analysis codebase: 
1. Install git.
1. Ask @SGMartin for the appropiate collaborator rights.
1. Clone the repo.

The default branch is called **main**. Feel free to push directly to this branch
unless stated otherwise.

## Get the data

Due to GitHub limitations when it comes to binary/big files (i.e. from single
cell experiments), we are hosting the required data files separately. To collaborate, 
you'll need to **sync** with us.

1. Download and install [**Syncthing**](https://syncthing.net/).
1. (Optional, but recommended) Set **Syncthing** as a [_daemon_](https://docs.syncthing.net/users/faq.html#how-do-i-run-syncthing-as-a-daemon-process-on-linux).
1. Ask @SGMartin to share the project folder. 
1. Verify that you are getting a local copy of the files in your chosen location.

## Install the required packages

1. Install [**conda**](https://docs.conda.io/en/latest/miniconda.html) in case you don't have it already. 
1. Locate the **envs** directory in your local copy of the repo. There you'll find 1 or more .yml files.
1. Each **yml** file is named after a certain step of the analysis pipeline. Run:

```
conda install -f [path_to_env_file.yml] -n  bc_meta
```

This create a new conda environment called **bc_meta** with the required dependencies
(i.e. R, Seurat, beyondcell etc.). You are now ready to collaborate with us.

## Set up the config

1. In your local copy of the repo, you'll find a template called **basics template.tsv**
inside the config directory. 
1. Fill every field and save the file as basics.tsv. 
1. Remember **NOT TO UPLOAD** your local copy of basics.tsv to _GitHub_. It 
should be in the _.gitignore_ file already, but please double check before commiting
changes to the repo.

## Useful tips

1. You may want to launch Rstudio and load the _Rproj._ file in your local copy
of the repo. This will set the default working directory to the repository.
1. We are (loosely) following these guidelines for 
coding conventions: [tidy guidelines](https://style.tidyverse.org/).

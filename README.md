# mscape-template

This repository is a template for creating new repositories containing code that 
will run on mSCAPE. It serves as a guide for code layout and files will need 
amending to fit the repo purpose.

There is guidance documentation including SOPs that contain further information on required repository
structure, development cycles, and making pull requests. Please read
this guidance documentation before using the template.

Below the dashed line is a README template to be edited.

---------------------
# Project or Repo Name

Brief description of project here

## Installation

Add installation instructions here. Ideally include commands to make
the process as easy as possible for users.

Clone repo and create environment:

`git clone git@github.com:ukhsa-collaboration/project-name.git`

`conda env create -n project-name `

`conda activate project-name`

Installation for users: 

`cd project-name`

`pip install .`

Installation for developers (installs code in editable mode): 

`cd project-name`

`pip install --editable '.[dev]'`

## Usage

Include command line arguments (e.g. the output displayed when using -h)
for reference. Example commands can also be helpful.

```
project-name --input <path> --output <path>
```

## Inputs (optional - useful if you have a CLI)

You may wish to use a table to list out the args:

| Argument | Required | Description |
| -------- | ------- | ------- |
| --input, -i | Yes | Input for command line use  |
| --output, -o | Yes | Output for command line use  |

## Outputs (optional - useful if your tool creates lots of output files)

Explain the output files in detail here.

## Other sections

Add other sections as appropriate for your repo. This may include
instructions on updating the repo, instructions on adding new
references, troubleshooting etc.

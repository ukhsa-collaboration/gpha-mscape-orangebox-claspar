############
Installation
############

****************************************
🔧 Installation for Command Line Use 🔧
****************************************

1) Set up environment - It is recommended that you install claspar into a
suitable environment. For example:

Create a new conda environment

.. code-block:: bash

    conda create -n claspar

    conda activate claspar


2a.) Either clone the repo and install from there:

.. code-block:: bash

    git clone https://github.com/ukhsa-collaboration/gpha-mscape-orangebox-claspar.git

    cd gpha-mscape-orangebox-claspar

    pip install .

*or:*

2b.) Install directly into the environment using pip:

.. code-block:: bash

    pip install git+https://github.com/ukhsa-collaboration/gpha-mscape-orangebox-claspar.git


***********************************
🔨 Installation for developers 🔨
***********************************
.. code-block:: bash

    git clone git@github.com:ukhsa-collaboration/gpha-mscape-orangebox-claspar.git

    cd gpha-mscape-orangebox-claspar

    pip install --editable '.[dev]'

    pre-commit install

**Be sure to install the pre-commit hooks in the repo if developing and pushing to the repo.**

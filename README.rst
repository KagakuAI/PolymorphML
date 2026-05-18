# PolymorphML
Molecular polymorphism prediction with machine learning

Installation
--------------------------------------------------------------------

Basic requirements:

.. code-block:: bash

    conda create -n polymorph python=3.10
    conda activate polymorph
    pip install rdkit molfeat scikit-learn xgboost

For Jupyter Notebook:

.. code-block:: bash

    conda activate polymorph
    pip install ipykernel
    python -m ipykernel install --user --name polymorph --display-name "polymorph"
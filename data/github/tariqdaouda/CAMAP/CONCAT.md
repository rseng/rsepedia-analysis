CAMAP
=====

Codon Arrangement MAP Predictor, predicting MHC-I Associated Peptides presentation from mRNA
--------------------------------------------------------------------------------------------

CAMAP's main aim is the acceleration of vaccine development through the prediction of epitope candidates on the surface of cells.


Short description:
------------------

After installation CAMAP can be easily called from the command line:

.. code:: shell

 camap

CAMAP's options can be listed with:

.. code:: shell

 camap --help

And here's an example command you can try:

.. code:: shell

 camap 162 -e 5 -n my_first_test -w 4 -o sgd --device cpu -subf pytorch-sgd_test

Note that specifying 5 epochs is ridiculous, you should try aiming for > 1000. For instance, for our 2020 paper we used 4000.

##################--------A brief introduction to iSumoK-PseAAC Server---------######################

Sumoylation is the post-translational modification that is involved in the adaption of cells and the functional properties of a large number of proteins. Sumoylation has key importance in subcellular concentration, transcriptional synchronization, chromatin remodeling, response to stress, and regulation of mitosis. Sumoylation is associated with developmental defects and many diseases in humans such as cancer and Huntington's, Alzheimer's, and Parkinson's diseases, Spin cerebellar ataxia 1, and amyotrophic lateral sclerosis. The covalent bonding of Sumoylation is essential to inheriting part of the operative characteristics of some other proteins. For that reason, the recognition of feasible Sumoylation sites has importance for research to find out the solution to many diseases.

############ Some remarks on used python code files ###########

The "Dataset" directory includes the positive and negative sequences files in fasta format.

The script in app.py is critical for controlling the overall functionality of the webserver. It includes all the necessary procedures that are used to interact with the user and perform operations requested upon input and navigations.

The script in extractFeatures.py includes all the important procedures that compute the features from the given protein sequences and make predictions based on the trained models using "iSumoK_Model.pkl" file. Furthermore, this python script also includes all the implementations of processes used to implement statistical moments.

passenger_wsgi.py python script file is the main application loader file for accessing the webserver.

The "templates" directory includes all the webpages in HTML that are used throughout during the application processing.

The "static" directory includes the datasets that are used in training and testing the computational model. It also includes the python packages for the current model's implementation. Most of these models are implemented and tested using Scikit-Learn library for python.

The "requirements.txt" file includes the list of all the python pips that are required to install and run the current project code.

For more information and queries kindly contact: ahmad.hassan@umt.edu.pk

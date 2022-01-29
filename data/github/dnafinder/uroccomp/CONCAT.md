# uroccomp
Compare two unpaired ROC curves. <br/>
The ROC graphs are a useful tecnique for organizing classifiers and
visualizing their performance. ROC graphs are commonly used in medical
decision making.
This function compares two unpaired ROC curves using my previously
submitted routine ROC
(http://www.mathworks.com/matlabcentral/fileexchange/19950)
If this file is absent, uroccomp will try to download it from FEX.

Syntax: uroccomp(x,y,alpha)

Input: x and y - These are the data matrix. 
                 The first column is the column of the data value;
                 The second column is the column of the tag: 
                 unhealthy (1) and healthy (0).
         alpha - significance level (default 0.05)

Output: The ROC plots;
        The z-test to compare Areas under the curves

  run uroccompdemo to see an example

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2009) uROCcomp: compare two unpaired ROC curves.
http://www.mathworks.com/matlabcentral/fileexchange/23020

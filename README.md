# compute_RkCC
Function to compute the K-category correlation coefficient implemented in Python. The K-category correlation is a measure of classification performance and may be considered a multiclass generalization of the [Matthews correlation coefficient](https://en.wikipedia.org/wiki/Matthews_correlation_coefficient).

## Comparing two K-category assignments by a K-category correlation coefficient

### Abstract


Predicted assignments of biological sequences are often evaluated by Matthews correlation coefficient. However, Matthews correlation coefficient applies only to cases where the assignments belong to two categories, and cases with more than two categories are often artificially forced into two categories by considering what belongs and what does not belong to one of the categories, leading to the loss of information. Here, an extended correlation coefficient that applies to K-categories is proposed, and this measure is shown to be highly applicable for evaluating prediction of RNA secondary structure in cases where some predicted pairs go into the category “unknown” due to lack of reliability in predicted pairs or unpaired residues. Hence, predicting base pairs of RNA secondary structure can be a three-category problem. The measure is further shown to be well in agreement with existing performance measures used for ranking protein secondary structure predictions. 

Paper is available at http://www.sciencedirect.com/science/article/pii/S1476927104000799

Paper author's server and software is available at http://rk.kvl.dk/






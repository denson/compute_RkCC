# -*- coding: utf-8 -*-
"""
Spyder Editor

Created on Fri Mar 31 14:35:34 2017

Function to compute the K-category correlation coefficient
 
@author: Denson Smith


Comparing two K-category assignments by a K-category correlation coefficient

Abstract


Predicted assignments of biological sequences are often evaluated by 
Matthews correlation coefficient. However, Matthews correlation coefficient applies 
only to cases where the assignments belong to two categories, and cases with more 
than two categories are often artificially forced into two categories by considering 
what belongs and what does not belong to one of the categories, leading to the 
loss of information. Here, an extended correlation coefficient that applies to 
K-categories is proposed, and this measure is shown to be highly applicable for 
evaluating prediction of RNA secondary structure in cases where some predicted 
pairs go into the category “unknown” due to lack of reliability in predicted pairs 
or unpaired residues. Hence, predicting base pairs of RNA secondary structure can 
be a three-category problem. The measure is further shown to be well in agreement 
with existing performance measures used for ranking protein secondary structure 
predictions. 

Server and software is available at http://rk.kvl.dk/

Paper is available at http://www.sciencedirect.com/science/article/pii/S1476927104000799


"""


__author__ = 'Denson Smith'

import numpy as np
 
 
from sklearn.metrics import classification_report ,confusion_matrix

def compute_RkCC(CM):     
    
    '''
    Function to compute the K-category correlation coefficient
    http://www.sciencedirect.com/science/article/pii/S1476927104000799
    
    http://rk.kvl.dk/suite/04022321447260711221/
    
    
    Parameters
    ----------
    CM : k X k confusion matrix of int
    
    n_samples : int
    
    
    Returns
    -------
    RkCC: float
    
    
    '''
    rows, cols = np.shape(CM)
    
    RkCC_numerator=0
    for k_ in range(cols):
        for l_ in range(cols):
            for m_ in range(cols):
    
                this_term = (CM[k_,k_] * CM[m_,l_]) - \
                    (CM[l_,k_] * CM[k_,m_])
    
                RkCC_numerator = RkCC_numerator + this_term
    
    RkCC_denominator_1=0           
    for k_ in range(cols):
        RkCC_den_1_part1=0
        for l_ in range(cols):
            RkCC_den_1_part1= RkCC_den_1_part1+CM[l_,k_]
    
        RkCC_den_1_part2=0
        for f_ in range(cols):
            if f_ != k_:
    
                for g_ in range(cols):
    
                    RkCC_den_1_part2= RkCC_den_1_part2+CM[g_,f_]
    
        RkCC_denominator_1=(RkCC_denominator_1+(RkCC_den_1_part1*RkCC_den_1_part2))
    
    
    
    RkCC_denominator_2=0           
    for k_ in range(cols):
        RkCC_den_2_part1=0
        for l_ in range(cols):
            RkCC_den_2_part1= RkCC_den_2_part1+CM[k_,l_]
    
        RkCC_den_2_part2=0
        for f_ in range(cols):
            if f_ != k_:
    
                for g_ in range(cols):
    
                    RkCC_den_2_part2= RkCC_den_2_part2+CM[f_,g_]
    
        RkCC_denominator_2=(RkCC_denominator_2+(RkCC_den_2_part1*RkCC_den_2_part2))
    
    RkCC = (RkCC_numerator)/(np.sqrt(RkCC_denominator_1)* np.sqrt(RkCC_denominator_2))
    
    return RkCC

def print_cm(cm, labels, hide_zeroes=False, hide_diagonal=False, hide_threshold=None):
    """pretty print for confusion matrixes"""
    columnwidth = max([len(x) for x in labels]+[5]) # 5 is value length
    empty_cell = " " * columnwidth
    # Print header
    print "    " + empty_cell,
    for label in labels: 
        print "%{0}s".format(columnwidth) % label,
    print
    # Print rows
    for i, label1 in enumerate(labels):
        print "    %{0}s".format(columnwidth) % label1,
        for j in range(len(labels)): 
            cell = "%{0}.1f".format(columnwidth) % cm[i, j]
            if hide_zeroes:
                cell = cell if float(cm[i, j]) != 0 else empty_cell
            if hide_diagonal:
                cell = cell if i != j else empty_cell
            if hide_threshold:
                cell = cell if cm[i, j] > hide_threshold else empty_cell
            print cell,
        print
        
        

if __name__ == "__main__":
    
    
    print('Using scikit-learn')
    
    y_true = [0, 1, 2, 2, 2]
    y_pred = [0, 0, 2, 2, 1]
    target_names = ['class 0', 'class 1', 'class 2']
     
    labels = [0,1,2]
    
    confusion_matrix = confusion_matrix(y_true,y_pred,labels = labels)
    
    print(confusion_matrix)
     
    print(classification_report(y_true, y_pred, target_names=target_names))
     
    
    RkCC = compute_RkCC(confusion_matrix)
     
    print('Rk correlation coefficient = %.4f' % RkCC)
    
    print('')
    print("Example from author's website.")
    confusion_matrix = np.array([[213, 21, 12], [89, 459, 90], [29, 39, 958]])
    
    
    RkCC = compute_RkCC(confusion_matrix)
    
    print(confusion_matrix)
    print('Rk correlation coefficient = %.4f' % RkCC)
        
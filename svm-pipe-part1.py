from sklearn.svm import SVC
import os, sys, random
from sklearn.externals import joblib
import numpy as np


# input: negatives.table and positives.table

raw_data_positives, raw_data_negatives = [], []

def read_raw_data_tablefile(fname):
    rawdata = []
    with open(fname, 'r') as dataflow:
        for line in dataflow.readlines()[1:]:
            line = line.replace(',', '.')
            toks = line[:-1].split('\t')
            if len(toks) >= 4 and line[0] == 'A':
                rawdata.append([float(toks[2]), float(toks[3])])
                #rawdata.append([float(toks[2]), float(toks[3]), float(toks[4] / 100.0)]) # 3D case
    return rawdata

raw_data_positives = read_raw_data_tablefile(sys.argv[1]) # 1 arg = positives table file name
raw_data_negatives = read_raw_data_tablefile(sys.argv[2]) # 2 arg = negatives table file name

alldata = raw_data_positives + raw_data_negatives
class_labels = [1 for x in range(len(raw_data_positives))] + [0 for x in range(len(raw_data_negatives))]

datavector = []
for i, j in zip(alldata, class_labels):
    datavector.append([i[0], i[1], j])
random.shuffle(datavector)

num_test_samples = -300  # use this number of samples as test set

Train_in = [[x[0], x[1]] for x in datavector[:num_test_samples]]
Train_out = [x[2] for x in datavector[:num_test_samples]]
Test_in = [[x[0], x[1]] for x in datavector[num_test_samples:]]
Test_out = [x[2] for x in datavector[num_test_samples:]]

clf = SVC(kernel='poly', degree=4, C=0.8) # you can try to adjust this (kernel='sigmoid', degree=4, C=0.7) also works well
clf.fit(Train_in, Train_out)
Scoring = clf.score(Test_in, Test_out)
print '{}'.format(Scoring)

# do not need to dump and load model! but note that due to random shuffling and splitting SVM can learn different frontiers!
'''
joblib.dump(clf, 'svm_model.clf')
clf = joblib.load('./svm_model.clf')
'''

ofile = open(sys.argv[3], 'w') # output table in a file
for x in range(0,101):
    xVal = float(x)/100
    for y in range(0,101):
        yVal = float(y)/100
        predVal = clf.predict([[xVal,yVal]])
        ofile.writelines("{0:.2f}\t{1:.2f}\t{2:.2f}\n".format(xVal, yVal, predVal[0]))
        #print("{0:.2f}".format(xVal)+"\t"+"{0:.2f}".format(yVal)+"\t"+"{0:.2f}".format(predVal[0]))

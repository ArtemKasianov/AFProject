from sklearn.externals import joblib
SVM_clf = joblib.load('./svm_model.clf')
for x in range(0,101):
    xVal = float(x)/100
    for y in range(0,101):
	yVal = float(y)/100
	
	predVal = SVM_clf.predict([[xVal,yVal]])
	print(str(xVal)+"\t"+str(yVal)+"\t"+str(predVal[0]))


#SVM_clf.predict([[0.9,0.1], [0.1,0.9], [0.9,0.9], [0.1,0.1]]) # -> [1 0 1 0]

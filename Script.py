#!/usr/bin/env python3


import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import glob, os
import re

from pandas import *
import pandas as pd
from scipy.stats import ttest_ind

from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import SGDClassifier
from sklearn import tree
from sklearn.neural_network import MLPClassifier
from sklearn import neighbors
from sklearn import svm

import numpy as np


def ChooseFolder():
    root = tk.Tk()
    root.withdraw()
    root.update()
    print("Select the Working Folder (PythonProject)")
    dirpath = filedialog.askdirectory(initialdir=os.getcwd(),
                                       title="Select the Working Folder")
    root.update()
    root.destroy()


    if not len(glob.glob(dirpath+"/*.txt"))==2 and not len(glob.glob(filename+"/*.csv"))==3:
        messagebox.showinfo("Visualizer error", "Wrong number of files, maybe you didn't choose the right directory?")
        sys.exit("Error, wrong file or directory")
    if not len(glob.glob(dirpath+"/temp/*"))==0:
        messagebox.showinfo("Visualizer error", "Folder not empty, cannot be used for the temp files")
        sys.exit("Error, wrong file or directory")
    if not len(glob.glob(dirpath+"/out/*"))==0:
        messagebox.showinfo("Visualizer error", "Folder not empty, cannot be used for the output")
        sys.exit("Error, wrong file or directory")
    return(dirpath)

def FiltList(tabfo):
    print(tabfo)

    try:
        muta = open(tabfo +'/mut.csv',"r").read().split('\n')
    except FileNotFoundError:
        sys.exit("Mutation file cannot be found, was it renamed or moved?")


    print("Preparing lists")


    try:
        m = open(tabfo +'/pheno_r_file.csv',"r").read().split('\n')
    except FileNotFoundError:
        sys.exit("Phenotype file cannot be found, was it renamed or moved?")


    phen=[]

    for n in range(len(m)-1):
        x=m[n].split(",")
        for i in range(len(x)):
            if x[3]=="POS" and x[4]=="NEG":
                if x[10]=="Ductal/NST" :
                    phen.append(x[1])

    phen=sorted(set(phen))

    lismut=["AKT1", "AKT2", "AKT3", "CDKN1A", "CDKN1B", "ERBB2", "ERBB3", "ERBB4", "FGFR1", "FGFR2", "FGFR3",
            "FGFR4", "GRB2", "HIF1A", "HRAS", "IGF1R", "IRS1", "KIT", "MTOR", "NRAS", "PDGFRA", "PDGFRB",
            "PIK3CB", "PIK3R1", "PRKCA", "PRKCB", "PRKCG", "PTEN", "RICTOR", "RPS6KB1", "RPTOR", "SOS1",
            "TSC1", "TSC2"]

    lmu=[]

    for n in range (1, len(muta)-1):
        x = muta[n].split(',')
        x[9]=x[9]
        if x[9] in lismut  :
            lmu.append(x[1])

    lmu=sorted(set(lmu))

    Hmu=[]
    blmu=[]

    for n in range (1, len(muta)-1):
        x = muta[n].split(',')
        x[9]=x[9]
        if x[9]=="PIK3CA":
            if x[12] == "1047" or x[12] == "545" or x[12] == "542":
                Hmu.append(x[1])
            else:
                blmu.append(x[1])

    Hmu=sorted(set(Hmu))

    blmu=sorted(set(blmu))


    try:
        p = open(tabfo +"/CN.csv","r").read().split('\n')
    except FileNotFoundError:
        sys.exit("CopyNumber file cannot be found, was it renamed or moved?")

    li=[]
    lik=[]

    pa=p[0].split(',')

    for i in range(1, len(p)-1):
        x=p[i].split(',')
        if x[7] in lismut:
            for n in range(7,len(x)):
                if x[n]=="-2" or x[n]=="2":
                    li.append(pa[n])
                else: lik.append(pa[n])

    lif=sorted(set(li))

    return(phen, lmu, Hmu, blmu, lif)


def Filt(er, mu, hm, bl, cn, file, tabfo):

    try:
        di = open(file,"r").read().split('\n')
    except FileNotFoundError:
        sys.exit("Gene expression file cannot be found, was it renamed or moved?")

    fname=re.split(r'[\.-\/]',file)[-2]

    print("Preparing the " + fname + " dataset")

    der=[""]*len(di)

    paz=di[0].split(',')


    c=[0]*len(paz)

    for i in range(1, len(paz)):
        if  paz[i] not in er:
            c[i]= 1


    print("Number of samples for each ER status:\n ER + ",c.count(0)-1,"ER - ", c.count(1))


    print("Samples processed")

    for n in range(len(di)-1):
        if n % 10000 == 0 : print(str(n) + "/" + str(len(di)))
        x=di[n].split(",")
        der[n]=x[0]
        for i in range(1, len(c)):
            if c[i] == 0 :
                der[n]=der[n]+ "," + x[i]


    di= der
    der=[""]*len(di)

    paz=di[0].split(',')


    c=[0]*len(paz)


    for i in range(1, len(paz)):
        if  paz[i] in hm :
            c[i]= 2
        elif paz[i] in mu or paz[i] in cn :
            c[i]= 1


    print("Number of samples for each subset:\n Gray ", c.count(0)-1,"Purple ", c.count(1),"Red ", c.count(2))

    print("Samples processed")

    for n in range(len(di)-1):
        if n % 10000 == 0 : print(str(n) + "/" + str(len(di)))
        x=di[n].split(",")
        der[n]=x[0]
        for i in range(1, len(c)):
            if c[i] != 1 :
                der[n]=der[n]+ "," + x[i]



    di= der
    der=[""]*len(di)

    paz=di[0].split(',')

    c=[0]*len(paz)


    for i in range(len(paz)):
        if  paz[i] in bl and paz[i] not in hm:
            c[i]= 1


    print("Number of samples for each subset:\n Gray and Red ", c.count(0)-1,"Blue ", c.count(1))



    print("Samples processed")

    for n in range(len(di)-1):
        if n % 10000 == 0 : print(str(n) + "/" + str(len(di)))
        x=di[n].split(",")
        der[n]=x[0]
        for i in range( len(c)):
            if c[i] == 0 and i!=0 :
                der[n]=der[n]+ "," + x[i]
            elif c[i] == 1 or i==0:
                if i == 0:
                    with open(tabfo + '/out/data_blu_'+ fname+'.csv', 'a') as out:
                        out.write(x[i])
                else:
                    with open(tabfo + '/out/data_blu_'+ fname+'.csv', 'a') as out:
                        out.write("," + x[i] )
        with open(tabfo + '/out/data_blu_'+ fname+'.csv', 'a') as out:
            out.write("\n")

    di= der
    der=[""]*len(di)

    paz=di[0].split(',')


    c=[1]*len(paz)

    for i in range(1, len(paz)):
        if  paz[i] in hm:
            c[i]= 2


    print("Number of samples for each subset:\n Gray ",c.count(1)-1,"Red ", c.count(2))


    for n in range(len(di)):
        with open(tabfo + '/out/Vdataset'+ fname+'.csv', 'a') as out:
            out.write(di[n] + "\n")

    filelist = glob.glob(tabfo + "/temp/*.csv")
    for f in filelist:
        os.remove(f)



    print(fname + " dataset complete")
    input("Press enter to continue")

    return(c[1:len(c)])




def SampleSelection (tabfo):
    classes=[]
    for file in glob.glob(tabfo+"/*.txt"):

        phen, lmu, Hmu, blmu, lif = FiltList(tabfo)
        os.system('clear')
        classes.append(Filt(phen, lmu, Hmu, blmu, lif, file, tabfo))

    return(classes)

def TopFeatures(tabfo, ngen):
    print("Script is trying to obtain the top "+ str(ngen) +" features." )
    mu = pd.read_csv(tabfo + '/out/VdatasetTraining.csv', sep= ",", index_col=0)

    mu = mu.transpose()

    indWT = [i for i, x in enumerate(clasf[0]) if x == 1]
    indMU = [i-1 for i, x in enumerate(clasf[0]) if x == 2]

    tt=ttest_ind(mu.iloc[indWT], mu.iloc[indMU])


    mu = mu.transpose()

    mu = mu.assign(pvalue=pd.Series(tt[1]).values)
    mu = mu.sort_values('pvalue')


    mu.head(n=5)
    mu.drop(mu.columns[[-1]], axis=1, inplace=True)
    mu=mu.head(n=ngen)
    return(mu, list(mu.index))



def NumGenes(prompt="Please enter how many features you want to use: "):

    while True:
        num_str = input(prompt).strip()
        if all(c in '0123456789' for c in num_str):
            break
        else:
            sys.exit("Please write an integer")
    return int(num_str)


def TestML(mu, clasf):
    dataset_array = mu.values
    get=np.transpose(dataset_array)

    target=clasf[0]


    clf1 = LogisticRegression(penalty="l2", solver="liblinear", max_iter=200)
    clf2 = RandomForestClassifier(n_jobs=100, random_state=None, n_estimators=1000, max_features=10)
    clf3 = GaussianNB()
    clf4 = SGDClassifier(loss="log", penalty="l2", max_iter=1000, tol=None)
    clf5 = tree.DecisionTreeClassifier(max_features=6, max_depth=4)
    clf6 = MLPClassifier(solver='adam', alpha=1e-5,
                        hidden_layer_sizes=(5, 2), max_iter= 1000)
    clf7 = neighbors.KNeighborsClassifier(15)
    clf8 = svm.SVC(kernel='linear')

    eclf = VotingClassifier(estimators=[('lr', clf1), ('rf', clf2), ('gnb', clf3), ('SGD', clf4), ('dt', clf5),
                                        ('mlp', clf6), ('nn', clf7), ('svc', clf8)], voting='hard')
    print("Testing different Machine Learning Methods:\n")

    for clf, label in zip([clf1, clf2, clf3, clf4, clf5, clf6, clf7, clf8, eclf],
                          ['Logistic Regression', 'Random Forest', 'Naive Bayes', 'Ensemble',
                           'Stochastic Gradient Descent', 'Decision Tree', 'Multilayer Perceptron','Nearest Neighbor', 'Support Vector Machines']):
        scores = cross_val_score(clf, get, target, cv=10, scoring='accuracy')


        print("Accuracy: %0.2f (+/- %0.2f) [%s]" % (scores.mean(), scores.std(), label))




def ModelBuild(mu, clasf):
    print("The script is building the model")

    dataset_array = mu.values
    get=np.transpose(dataset_array)

    target=clasf[0]

    clf = LogisticRegression(penalty="l2", solver="liblinear", max_iter=200)

    print("The script is testing the model on the training data")

    clfe = clf.fit(get, target)

    targg=clfe.predict(get[:, :])

    cv=cross_val_score(clf, get, target, cv=10, scoring='accuracy')


    print("10 Fold Cross Validation Accuracies (Training Set):\n"+"\n".join(map(str,cv)))

    return(clfe)

def ValidationSet(glist, clasf, clfe):
    muv = pd.read_csv(tabfo + '/out/VdatasetValidation.csv', sep= ",", index_col=0)

    muv = muv.transpose()

    muv = muv[glist]

    valid_array = muv.values
    getv = valid_array

    targetv=clasf[1]

    targgv=clfe.predict(getv[:, :])

    cvv=cross_val_score(clfe, getv, targetv, cv=10, scoring='accuracy')

    print("10 Fold Cross Validation Accuracies (Validation Set):\n"+"\n".join(map(str,cvv)))

def FindBlue(blc, glist, clfe, file):
    blt = pd.read_csv(file, sep= ",", index_col=0)

    blt = blt.transpose()

    blt = blt[glist]

    bt_array = blt.values
    getbt = bt_array


    targbt=clfe.predict(getbt[:, :])

    btMU = [i-1 for i, x in enumerate(targbt) if x == 2]


    btn=list(blt.index[btMU])

    try:
        muta = open(tabfo +'/mut.csv',"r").read().split('\n')
    except FileNotFoundError:
        sys.exit("Mutation file cannot be found, was it renamed or moved?")



    for n in range (1, len(muta)-1):
        x = muta[n].split(',')
        x[9]=x[9]
        #print(x[1])
        if x[9]=="PIK3CA" and x[1] in btn:
            blc.append(x[12])

    return(blc)



tabfo=ChooseFolder()

print("The working folder is " + tabfo)

clasf= SampleSelection(tabfo)

ngen = NumGenes()
mu, glist = TopFeatures(tabfo, ngen)

TestML(mu, clasf)

np.random.seed(None)

clfe = ModelBuild(mu, clasf)

ValidationSet(glist, clasf, clfe)

blc=[]

for file in glob.glob(tabfo+"/out/data_*.csv"):
    FindBlue(blc, glist, clfe, file)


print("List of mutations in the blue subset classified similiar to the hotspots ones:\n" + "\n".join(map(str,list(set(blc)))))

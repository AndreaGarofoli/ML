
## from sklearn import datasets
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import SGDClassifier
from sklearn import tree
from sklearn.neural_network import MLPClassifier
from sklearn import neighbors


ge = open("/Users/andrea/Documents/BI/BI/B/PhD/bla.txt","r").read().split('\n')

get=np.transpose(np.loadtxt(ge))

target=[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

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

x=0
y=0

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(6, 6), sharey=True)

for clf, label in zip([clf1, clf2, clf3, clf4, clf5, clf6, clf7, eclf], 
                      ['Logistic Regression', 'Random Forest', 'naive Bayes', 'Ensemble', 
                       'SGDClassifier', 'Decision Tree', 'MLPClassifier','NN', 'SVC']):
    scores = cross_val_score(clf, get, target, cv=10, scoring='accuracy')
    axes[y, x].boxplot(scores)
    #axes[x, y].set_title(label, fontsize=fs)
    if x is not 3: x+=1
    else: 
        x=0
        y+=1
    
    print("Accuracy: %0.2f (+/- %0.2f) [%s]" % (scores.mean(), scores.std(), label))


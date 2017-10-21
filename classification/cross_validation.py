# -*- coding=utf-8 -*-
import numpy as np
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import accuracy_score

class mySVC:
    @staticmethod
    def do(train_data, train_label, test_data, test_label=None):
        clf = SVC()
        clf.fit(train_data, train_label)
        predicts = clf.predict(test_data)
        acc = None

        print test_label
        print predicts
        if test_label is not None:
            acc = accuracy_score(test_label, predicts)
            print acc
        return predicts, acc

class myLinearSVC:
    @staticmethod
    def do(train_data, train_label, test_data, test_label=None):
        svm = LinearSVC()
        svm.fit(train_data, train_label)
        predicts = svm.predict(test_data)
        acc = None
        print test_label
        print predicts
        if test_label is not None:
            acc = accuracy_score(test_label, predicts)
            print acc
        return predicts, acc

class myXgboost:
    @staticmethod
    def do(training_data, train_label, test_data, test_label=None):
        import xgboost as xgb
        param = {}
        # use softmax multi-class classification
        param['objective'] = 'multi:softmax'
        # scale weight of positive examples
        param['eta'] = 0.1
        param['max_depth'] = 6
        param['silent'] = 1
        param['nthread'] = 4
        param['num_class'] = np.max(train_label) + 1
        plst = param.items()
        dTrain = xgb.DMatrix(training_data, label=train_label)
        num_round = 100
        bst = xgb.train(plst, dTrain, num_round)
        dTest = xgb.DMatrix(test_data)
        yPred = bst.predict(dTest)
        print test_label
        print yPred
        acc = None
        if test_label is not None:
            acc = accuracy_score(test_label, yPred)
            print acc
        return yPred, acc

# coding=gbk
# k折交叉验证
from sklearn.model_selection import KFold

k_nums = 2
class KCrossValidation:

    @staticmethod
    def do(data, label, method, method_name):
        average_score = 0.0
        kf = KFold(n_splits=k_nums, shuffle=True)
        for train_index, test_index in kf.split(data, label):
            train_data, test_data = data[train_index], data[test_index]
            train_label, test_label = label[train_index], label[test_index]
            predicted_res, score = method(train_data, train_label, test_data, test_label)
            average_score += score
        average_score /= k_nums
        print 'function name is ', method_name, 'average score is ', average_score

class RunTest:
    @staticmethod
    def do(train_feature, train_label, test_feature, method, method_name):
        predict_res, acc = method(train_feature, train_label, test_feature, None)
        return predict_res, acc


def do_pca(features, n_components=100):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components=n_components, copy=True)

    features_pca = pca.fit_transform(features)
    return features_pca


if __name__ == '__main__':
    features = np.load(
        '/home/give/PycharmProjects/LiRads/Washout/histgram_splited_features.npy'
    )
    for i in range(len(features)):
        features[i, :] = features[i, :] / np.sum(features[i, :])
    print np.shape(features)
    labels = np.load(
        '/home/give/PycharmProjects/LiRads/Washout/histgram_labels.npy'
    )
    print np.shape(labels)
    features = do_pca(features, n_components=20)
    # print labels
    KCrossValidation.do(features, labels, myLinearSVC.do, 'svm')
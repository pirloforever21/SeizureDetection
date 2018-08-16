'''
This program is dedicated for Albany Medical College, Neuroscience & Experimental Therapeutics Department, Dr. Yunfei Huang's Laboratory.
Do not distribute without permission. For Copy Right information, please visit https://github.com/pirloforever21/SeizureDetection.

This program is designed to perform seizure detection on genetic modified mice, and produce figures of potential seizure.
This program take .txt/.csv files exported from Datawave and save .pdf results to the home folder.
For detail information, please read the user's manual.

Author: Zhenhuan(Neyo) Yang
Date: 08/15/2018
Version: 1.0
'''

from PyQt5 import QtWidgets
from PyQt5 import QtCore
import sys
import csv
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import datetime

class PrettyWidget(QtWidgets.QWidget):

    def __init__(self):
        super(PrettyWidget, self).__init__()
        self.initUI()

    def initUI(self):
        self.setGeometry(600,300, 1000, 600)
        self.setWindowTitle('Seizure Detection v1.0')

        # Grid Layout
        grid = QtWidgets.QGridLayout()
        self.setLayout(grid)

        self.lbl0 = QtWidgets.QLabel(self, alignment=QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.lbl0.setText('Built and maintained by Zhenhuan(Neyo) Yang')
        self.lbl0.adjustSize()
        self.lbl0.setOpenExternalLinks(True)
        grid.addWidget(self.lbl0, 0, 0)

        urlLink = "<a href=\"https://github.com/pirloforever21/SeizureDetection\">Visit My Github for Copy Right and More</a>"
        self.lbl1 = QtWidgets.QLabel(self, alignment=QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        #self.lbl1.setText('Built and maintained by Zhenhuan(Neyo) Yang\n' +'urlLink')
        self.lbl1.setText(urlLink)
        self.lbl1.adjustSize()
        self.lbl1.setOpenExternalLinks(True)
        grid.addWidget(self.lbl1, 1, 0)

        grid.setRowStretch(2, 1)  # <---- Sets the stretch factor of row to stretch
        # Label indicator
        self.lbl2 = QtWidgets.QLabel(self, alignment=QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lbl2.setText('Select data to start...')
        self.lbl2.adjustSize()
        grid.addWidget(self.lbl2, 6, 1)

        # Import data Button
        btn1 = QtWidgets.QPushButton('Select Data', self)
        btn1.resize(btn1.sizeHint())
        btn1.clicked.connect(self.getData)
        grid.addWidget(btn1, 3, 0)

        self.datapath = ''

        # Import names Button
        btn2 = QtWidgets.QPushButton('Select Labels', self)
        btn2.resize(btn2.sizeHint())
        btn2.clicked.connect(self.getNames)
        grid.addWidget(btn2, 3, 1)

        self.namespath = ''

        # Run Button
        btn3 = QtWidgets.QPushButton('Run', self)
        btn3.resize(btn3.sizeHint())
        btn3.clicked.connect(self.Run)
        grid.addWidget(btn3, 4, 0)

        self.names=[]
        self.fig = []

        # Save Button
        btn4 = QtWidgets.QPushButton('Save',self)
        btn4.resize(btn4.sizeHint())
        btn4.clicked.connect(self.Save)
        grid.addWidget(btn4, 4, 1)

        grid.setRowStretch(5, 1)  # <---- Sets the stretch factor of row row to stretch .

        self.show()

    def getData(self):
        self.datapath, _ = QtWidgets.QFileDialog.getOpenFileName(self,'Text File','','*.txt')
        if self.datapath == '':
            pass
        else:
            self.lbl2.setText('Data selected!')
            self.lbl2.adjustSize()

    def getNames(self):
        self.namespath,_ = QtWidgets.QFileDialog.getOpenFileName(self,'Text File','','*.txt')
        if self.namespath == '':
            pass
        else:
            self.lbl2.setText('Labels selected!')
            self.lbl2.adjustSize()

    def Run(self):
        df = pd.read_csv(str(self.datapath),header = None).dropna(axis='columns')
        m, n = df.shape  # m-by-n matrix

        names = []
        with open(str(self.namespath)) as csvfile:
            headreader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for name in headreader:
                names.extend(name)

        names = [name.strip() for name in names]
        names = [name for name in names if name]

        if len(names) == n:
            df.columns = names
            self.names = names
        else:
            names = []
            for column in range(n):
                names.extend(['Channel' + str(column)])
            df.columns = names
            self.names = names

        def spike(eeg):
            epoch = round(len(eeg) / 1000)  # change millisecond into second
            index = []
            collect = []
            for i in range(round(epoch)):
                power, freqs = mlab.psd(eeg.iloc[(1000 * i):(1000 * (i + 1))], NFFT=256, Fs=256)
                collect.append(max(power))  # Collect all the maximum power
            l = len(collect)
            average = sum(collect) / l
            index = [i for i in range(l) if collect[i] > 3 * average]  # Collect abnormal maximum power
            return index

        def clustering(index):
            clusters = [[index[0]]]
            k = 0
            for i in range(len(index) - 1):
                if index[i + 1] - index[i] < 20:
                    clusters[k].extend([index[i + 1]])
                else:
                    clusters.append([index[i + 1]])
                    k = k + 1
            return clusters

        def prune(clusters):
            pruned = []
            l = len(clusters)
            for i in range(l):
                if len(clusters[i]) > 10:
                    pruned.append(clusters[i])
            k = len(pruned)
            return pruned

        self.fig = ['fig' + str(column) for column in range(n)]

        for column in range(n):
            index = spike(df[names[column]])
            clusters = clustering(index)
            pruned_clusters = prune(clusters)

            k = len(pruned_clusters)
            '''
            with open(str(names[column]) + '.txt', 'w') as f:
                for cluster in range(k):
                    a = pruned_clusters[cluster][0] - 60
                    b = pruned_clusters[cluster][-1] + 60
                    start = str(datetime.timedelta(seconds=a))
                    end = str(datetime.timedelta(seconds=b))

                    f.write(str(start)+'-')
                    f.write(str(end)+'\n')
            '''
            self.fig[column] = plt.figure(figsize=(16,k*2))
            ax = ['ax' + str(cluster) for cluster in range(k)]
            for cluster in range(k):
                a = pruned_clusters[cluster][0] - 60
                b = pruned_clusters[cluster][-1] + 60
                start = str(datetime.timedelta(seconds=a))
                end = str(datetime.timedelta(seconds=b))

                ax[cluster] = self.fig[column].add_subplot(k, 1, cluster + 1)
                ax[cluster].plot(df[names[column]].iloc[(a * 1000):(b * 1000)], 'blue')
                ax[cluster].set_xticklabels([])
                ax[cluster].set_yticklabels([])
                ax[cluster].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', length=0)
                ax[cluster].spines["top"].set_visible(False)
                ax[cluster].spines["bottom"].set_visible(False)
                ax[cluster].spines["right"].set_visible(False)
                ax[cluster].spines["left"].set_visible(False)
                ax[cluster].set_xlabel(start + '-' + end, fontsize=12, rotation=0)
                ax[cluster].xaxis.set_label_position("bottom")

        self.lbl2.setText('Done!')
        self.lbl2.adjustSize()

    def Save(self):
        number = len(self.names)
        for column in range(number):
            self.fig[column].savefig(str(self.names[column]) + '.pdf')

        self.lbl2.setText('Saved!')
        self.lbl2.adjustSize()

def main():
    app = QtWidgets.QApplication(sys.argv)
    w = PrettyWidget()
    app.exec_()


if __name__ == '__main__':
    main()
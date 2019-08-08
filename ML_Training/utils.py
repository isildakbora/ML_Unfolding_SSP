import keras.callbacks
from IPython.display import clear_output
import matplotlib.pyplot as plt
import numpy as np
import ROOT

class PlotLearning(keras.callbacks.Callback, ):
    plt.style.use('seaborn')
    
    def __init__(self, monitor=['loss','mean_squared_error']):
        self.monitor = monitor
        
    def on_train_begin(self, logs={}):
        self.i = 0
        self.x = []
        self.monitored_1 = []
        self.val_monitored_1 = []
        self.monitored_2 = []
        self.val_monitored_2 = []
        self.fig = plt.figure()
        
        self.logs = []

    def on_epoch_end(self, epoch, logs={}):
        
        self.logs.append(logs)
        self.x.append(self.i)
        self.monitored_1.append(logs.get(self.monitor[0]))
        self.val_monitored_1.append(logs.get('val_'+self.monitor[0]))
        self.monitored_2.append(logs.get(self.monitor[1]))
        self.val_monitored_2.append(logs.get('val_'+self.monitor[1]))
        self.i += 1
        f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(10,5))
        
        clear_output(wait=True)
        
        ax1.set_yscale('log')
        ax1.plot(self.x, self.monitored_1, label = self.monitor[0])
        ax1.plot(self.x, self.val_monitored_1, label= 'val_' + self.monitor[0])
        ax1.legend()
        
        ax2.plot(self.x, self.monitored_2, label = self.monitor[1])
        ax2.plot(self.x, self.val_monitored_2, label = 'val_' + self.monitor[1])
        ax2.legend()
        
        plt.show();
        
plot = PlotLearning()

def array2TH1F(array, bins, name = 'historam', title = 'histogram'):  
    if (len(bins)==0):
        # if binning is not specified, use numpy's auto bin finder
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges
        bins  = np.histogram_bin_edges(array, 'auto')
    else:
        bins = np.array(bins, dtype=float)
        
    histo = ROOT.TH1F(name, title, len(bins)-1, bins)
        
    for x in array:
        histo.Fill(x)
        
    return histo

import keras, numpy
import matplotlib.pyplot as plt
from IPython.display import clear_output

class PlotLearning(keras.callbacks.Callback):
    plt.style.use('seaborn')
    
    def on_train_begin(self, logs={}):
        self.i = 0
        self.x = []
        self.losses = []
        self.val_losses = []
        self.acc = []
        self.val_acc = []
        self.fig = plt.figure()
        
        self.logs = []

    def on_epoch_end(self, epoch, logs={}):
        
        self.logs.append(logs)
        self.x.append(self.i)
        self.losses.append(logs.get('loss'))
        self.val_losses.append(logs.get('val_loss'))
        self.acc.append(logs.get('mean_squared_error'))
        self.val_acc.append(logs.get('val_mean_squared_error'))
        self.i += 1
        f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(10,5))
        
        clear_output(wait=True)
        
        ax1.set_yscale('log')
        ax1.plot(self.x, self.losses, label='loss')
        ax1.plot(self.x, self.val_losses, label='val_loss')
        ax1.legend()
        
        ax2.plot(self.x, 1/numpy.array(self.acc), label='1 / mse')
        ax2.plot(self.x, 1/numpy.array(self.val_acc), label='1 / validation mse')
        ax2.legend()
        
        plt.show();
        
plot = PlotLearning()
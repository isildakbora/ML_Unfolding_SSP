#!/usr/bin/env python3.7

from ROOT import TMVA, TFile, TTree, TCut, TMath
from subprocess import call
from os.path import isfile

from keras.models import Sequential
from keras.layers.core import Dense, Activation, Dropout
from keras.regularizers import l2
from keras.optimizers import SGD, Adadelta
from keras.initializers import glorot_uniform, VarianceScaling

# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

output = TFile.Open('TMVA.root', 'RECREATE')
factory = TMVA.Factory('TMVARegression', output, '!V:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=Regression')

# Load data
data = TFile.Open("ML_Unfolding_tree_100k.root")
tree = data.Get('outTree')

dataloader = TMVA.DataLoader('dataset')

# Define the input variables that shall be used for the MVA training
dataloader.AddVariable("pt1", "pt1", "GeV", 'F');
dataloader.AddVariable("pt2", "pt2", "GeV", 'F');
dataloader.AddVariable("eta1", "eta1", "units", 'F');
dataloader.AddVariable("eta2", "eta2", "units", 'F');
dataloader.AddVariable("phi1", "phi1", "units", 'F');
dataloader.AddVariable("phi2", "phi2", "units", 'F');
dataloader.AddVariable("iso1", "iso1", "units", 'F');
dataloader.AddVariable("iso2", "iso2", "units", 'F');

dataloader.AddVariable("met", "met", "units", 'F');
dataloader.AddVariable("met_eta", "met_eta", "units", 'F');
dataloader.AddVariable("met_phi", "met_phi", "units", 'F');
dataloader.AddVariable("scalar_ht", "scalar_ht", "units", 'F');

dataloader.AddVariable("electron_size", "electron_size", "units", 'I');
dataloader.AddVariable("muon_size", "muon_size", "units", 'I');
dataloader.AddVariable("photon_size", "photon_size", "units", 'I');
dataloader.AddVariable("jet_size", "jet_size", "units", 'I');

dataloader.AddVariable("mass_Z_reco", "mass_Z_reco", "GeV", 'F');

# Add the variable carrying the regression target
dataloader.AddTarget("mass_Z_gen");

dataloader.AddRegressionTree(tree, 1.0)

N_Train =  TMath.FloorNint(0.5 * tree.GetEntriesFast()); # 80% of the sample
N_Test  =  tree.GetEntriesFast() - N_Train; # remaining 20% of the sample
dataloader.PrepareTrainingAndTestTree(TCut(''), 'nTrain_Regression=' + str(N_Train) + ':SplitMode=Random:NormMode=NumEvents:!V')

# Generate model

# Define model
n_input = 17
model = Sequential()
model.add(Dense(128, activation='tanh', input_dim=n_input))
model.add(Dense(64, activation='tanh', input_dim=n_input))
model.add(Dense(32, activation='tanh', input_dim=n_input))
model.add(Dense(64, activation='tanh', input_dim=n_input))
model.add(Dense(128, activation='tanh', input_dim=n_input))
model.add(Dense(1, activation='linear', input_dim=n_input))

# Set loss and optimizer
weight_initializer = VarianceScaling(scale=1.0, mode='fan_in', distribution='normal', seed=None)
#weight_initializer=glorot_uniform(seed=None)
model.compile(loss="mean_squared_error", optimizer=Adadelta(lr=1e-02), initializer=weight_initializer)

# Store model to file
model.save('model.h5')
model.summary()

# Book methods
factory.BookMethod(dataloader, TMVA.Types.kPyKeras, 'PyKeras','H:!V:VarTransform=D,G:FilenameModel=model.h5:NumEpochs=100:BatchSize=32')
factory.BookMethod(dataloader, TMVA.Types.kBDT, 'BDTG',"!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=-1:PruneMethod=NoPruning:PruneStrength=50:MaxDepth=5")

# Run TMVA
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()	
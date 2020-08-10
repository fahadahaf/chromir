import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.autograd import Variable
from torch.autograd import Function # import Function to create custom activations
from torch.nn.parameter import Parameter # import Parameter to create custom activations with learnable parameters
from torch import optim # import optimizers for demonstrations


class Basset(nn.Module):
    def __init__(self, paramsdict):
        super(Basset, self).__init__()
        self.paramdict = paramsdict

        self.layer1  = nn.Sequential(
            nn.Conv1d(in_channels=self.paramdict['inputchannels'],
                      out_channels=self.paramdict['Conv1numfilters'],
                      kernel_size=self.paramdict['Conv1filtersize'],
                      padding=self.paramdict['Conv1padding'],
                      bias=False), #if using batchnorm, no need to use bias in a CNN
            nn.BatchNorm1d(num_features=self.paramdict['Conv1numfilters']),
            nn.ReLU() if self.paramdict['Conv1usesoftplus']==False else nn.Softplus(),
            nn.MaxPool1d(kernel_size=self.paramdict['Conv1maxpoolsize']))
        self.dropout1 = nn.Dropout(p=0.2)

        self.layer2 = nn.Sequential(
            nn.Conv1d(in_channels=self.paramdict['Conv1numfilters'],
                      out_channels=self.paramdict['Conv2numfilters'],
                      kernel_size=self.paramdict['Conv2filtersize'],
                      padding=self.paramdict['Conv2padding']),
            nn.BatchNorm1d(num_features=self.paramdict['Conv2numfilters']),
            nn.ReLU())
            #nn.MaxPool1d(kernel_size=3))
        self.dropout2 = nn.Dropout(p=0.2)

        self.layer3 = nn.Sequential(
            nn.Conv1d(in_channels=self.paramdict['Conv2numfilters'],
                      out_channels=self.paramdict['Conv3numfilters'],
                      kernel_size=self.paramdict['Conv3filtersize'],
                      padding=self.paramdict['Conv3padding']),
            nn.BatchNorm1d(num_features=self.paramdict['Conv3numfilters']),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=self.paramdict['Conv3maxpoolsize']))
        self.dropout3 = nn.Dropout(p=0.2)

        self.fc1 = nn.Linear(in_features=self.paramdict['Fc1inputsize'], out_features=self.paramdict['Fc1outputsize'])
        self.relu4 = nn.ReLU()
        self.dropout4 = nn.Dropout(p=0.4)

        self.fc2 = nn.Linear(in_features=self.paramdict['Fc1outputsize'], out_features=self.paramdict['Fc2outputsize'])
        self.relu5 = nn.ReLU()
        self.dropout5 = nn.Dropout(p=0.4)

        self.fc3 = nn.Linear(in_features=self.paramdict['Fc2outputsize'], out_features=self.paramdict['numtargets'])


    def forward(self, inputs):
        
        
        output = self.layer1(inputs)
        output = self.dropout1(output)

        output = self.layer2(output)
        output = self.dropout2(output)

        output = self.layer3(output)
        output = self.dropout3(output)

        output = output.reshape(output.size(0), -1)

        output = self.fc1(output)
        output = self.relu4(output)
        output = self.dropout4(output)

        output = self.fc2(output)
        output = self.relu5(output)
        output = self.dropout5(output)

        output = self.fc3(output)

        assert not torch.isnan(output).any()

        return output

class Basset_embd(Basset):
    def __init__(self):
        return
    
class Generalized:
    def __init__(self):
        return

class Generalized_embd(Generalized):
    def __init__(self):
        return
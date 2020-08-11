import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.autograd import Variable
from torch.autograd import Function # import Function to create custom activations
from torch.nn.parameter import Parameter # import Parameter to create custom activations with learnable parameters
from torch import optim # import optimizers for demonstrations


class Basset(nn.Module):
    def __init__(self, wvmodel=None, params):
        super(Basset, self).__init__()
		self.CNN1filters = params['CNN1_filters']
		self.CNN1filterSize = params['CNN1_filtersize']
		self.CNN1poolSize = params['CNN1_poolsize']
        self.CNN1padding = params['CNN1_padding']
        self.CNN1useSoftplus = params['CNN1_usesoftplus']
        self.CNN2filters = params['CNN2_filters']
		self.CNN2filterSize = params['CNN2_filtersize']
		self.CNN2poolSize = params['CNN2_poolsize']
        self.CNN2padding = params['CNN2_padding']
        self.CNN3filters = params['CNN3_filters']
		self.CNN3filterSize = params['CNN3_filtersize']
		self.CNN3poolSize = params['CNN3_poolsize']
        self.CNN3padding = params['CNN3_padding']
        self.FC1inputSize = params['FC1_inputsize']
        self.FC1outputSize = params['FC1_outputsize']
        self.FC2outputSize = params['FC2_outputsize']
		self.numClasses = params['num_classes']
        self.useEmbeddings = params['use_embeddings']
        if not self.useEmbeddings:
			self.numInputChannels = params['input_channels'] #number of channels, one hot encoding
		else:
            self.embSize = params['embd_size']
            weights = torch.FloatTensor(wvmodel.wv.vectors)
		    self.embedding = nn.Embedding.from_pretrained(weights, freeze=False)
			self.numInputChannels = self.embSize

        self.layer1  = nn.Sequential(
            nn.Conv1d(in_channels=self.numInputChannel,
                      out_channels=self.CNN1filters,
                      kernel_size=self.CNN1filterSize,
                      padding=self.CNN1padding,
                      bias=False), #if using batchnorm, no need to use bias in a CNN
            nn.BatchNorm1d(num_features=self.CNN1filters),
            nn.ReLU() if self.CNN1useSoftplus==False else nn.Softplus(),
            nn.MaxPool1d(kernel_size=self.CNN1poolSize))
        self.dropout1 = nn.Dropout(p=0.2)

        self.layer2 = nn.Sequential(
            nn.Conv1d(in_channels=self.CNN1filters,
                      out_channels=self.CNN2filters,
                      kernel_size=self.CNN2filterSize,
                      padding=self.CNN2padding,
                      bias=False),
            nn.BatchNorm1d(num_features=self.CNN2filters),
            nn.ReLU())
        self.dropout2 = nn.Dropout(p=0.2)

        self.layer3 = nn.Sequential(
            nn.Conv1d(in_channels=self.CNN2filters,
                      out_channels=self.CNN3filters,
                      kernel_size=self.CNN3filterSize,
                      padding=self.CNN3padding,
                      bias=False),
            nn.BatchNorm1d(num_features=self.CNN3filters),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=self.CNN3poolSize))
        self.dropout3 = nn.Dropout(p=0.2)

        self.fc1 = nn.Linear(in_features=self.FC1inputSize, out_features=self.FC1outputSize)
        self.relu4 = nn.ReLU()
        self.dropout4 = nn.Dropout(p=0.4)

        self.fc2 = nn.Linear(in_features=self.FC1outputSize, out_features=self.FC2outputSize)
        self.relu5 = nn.ReLU()
        self.dropout5 = nn.Dropout(p=0.4)

        self.fc3 = nn.Linear(in_features=self.FC2outputSize, out_features=self.numClasses)

    def forward(self, inputs):
        if self.useEmbeddings:
            output = self.embedding(inputs)
            output = output.permute(0,2,1)
        else:
		    output = inputs
        output = self.layer1(output)
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


class PositionalEncoding(nn.Module):
    # Taken from: https://nlp.seas.harvard.edu/2018/04/03/attention.html
    "Implement the PE function."
    def __init__(self, d_model, dropout, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)
       
        # Compute the positional encodings once in log space.
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0., max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0., d_model, 2) *
                             -(math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)
       
    def forward(self, x):
        x = x + Variable(self.pe[:, :x.size(1)],
                         requires_grad=False)
        return self.dropout(x)


class AttentionNet: #for the model that uses CNN, RNN (optionally), and MH attention
    def __init__(self, params, wvmodel=None):
        super(AttentionNet, self).__init__()
		self.numMultiHeads = params['num_multiheads']
		self.SingleHeadSize = params['singlehead_size']#SingleHeadSize
		self.MultiHeadSize = params['multihead_size']#MultiHeadSize
		self.usepooling = params['use_pooling']
		self.pooling_val = params['pooling_val']
		self.readout_strategy = params['readout_strategy']
		self.kmerSize = params['embd_kmersize']
		self.useRNN = params['use_RNN']
		self.useCNN = params['use_CNN']
		self.usePE = params['use_posEnc']
		self.useCNNpool = params['use_CNNpool']
		self.RNN_hiddenSize = params['RNN_hiddensize']
		self.numCNNfilters = params['CNN_filters']
		self.filterSize = params['CNN_filtersize']
		self.CNNpoolSize = params['CNN_poolsize']
        self.CNNpadding = params['CNN_padding']
		self.numClasses = params['num_classes']
        self.useEmbeddings = params['use_embeddings']
		if not self.useEmbeddings:
			self.numInputChannels = params['input_channels'] #number of channels, one hot encoding
		else:
            self.embSize = params['embd_size']
            weights = torch.FloatTensor(wvmodel.wv.vectors)
		    self.embedding = nn.Embedding.from_pretrained(weights, freeze=False)
			self.numInputChannels = self.embSize
		
		if self.usePE:
			self.pe = PositionalEncoding(d_model = self.numInputChannels, dropout=0.1)
		
		if self.useCNN and self.useCNNpool:
			self.layer1  = nn.Sequential(nn.Conv1d(in_channels=self.numInputChannels, out_channels=self.numCNNfilters,
										 kernel_size=self.filterSize, padding=self.CNNpadding, bias=False),nn.BatchNorm1d(num_features=self.numCNNfilters),
										 nn.ReLU(),nn.MaxPool1d(kernel_size=self.CNNpoolSize))
			self.dropout1 = nn.Dropout(p=0.2)
        
		if self.useCNN and self.useCNNpool == False:
			self.layer1  = nn.Sequential(nn.Conv1d(in_channels=self.numInputChannels, out_channels=self.numCNNfilters,
										 kernel_size=self.filterSize, padding=self.CNNpadding, bias=False),
										 nn.BatchNorm1d(num_features=self.numCNNfilters),nn.ReLU())
			self.dropout1 = nn.Dropout(p=0.2)
		
		if self.useRNN:
			self.RNN = nn.LSTM(self.numInputChannels if self.useCNN==False else self.numCNNfilters, self.RNN_hiddenSize, num_layers=2, bidirectional=True)
			self.dropoutRNN = nn.Dropout(p=0.4)
			self.Q = nn.ModuleList([nn.Linear(in_features=2*self.RNN_hiddenSize, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
			self.K = nn.ModuleList([nn.Linear(in_features=2*self.RNN_hiddenSize, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
			self.V = nn.ModuleList([nn.Linear(in_features=2*self.RNN_hiddenSize, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
		
		if self.useRNN == False and self.useCNN == False:
			self.Q = nn.ModuleList([nn.Linear(in_features=self.numInputChannels, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
			self.K = nn.ModuleList([nn.Linear(in_features=self.numInputChannels, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
			self.V = nn.ModuleList([nn.Linear(in_features=self.numInputChannels, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
		
		if self.useRNN == False and self.useCNN == True:
			self.Q = nn.ModuleList([nn.Linear(in_features=self.numCNNfilters, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
			self.K = nn.ModuleList([nn.Linear(in_features=self.numCNNfilters, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
			self.V = nn.ModuleList([nn.Linear(in_features=self.numCNNfilters, out_features=self.SingleHeadSize) for i in range(0,self.numMultiHeads)])
			
		self.RELU = nn.ModuleList([nn.ReLU() for i in range(0,self.numMultiHeads)])
		self.MultiHeadLinear = nn.Linear(in_features=self.SingleHeadSize*self.numMultiHeads, out_features=self.MultiHeadSize)#50
		self.MHReLU = nn.ReLU()
		
		self.fc3 = nn.Linear(in_features=self.MultiHeadSize, out_features=self.numClasses)
		
	def attention(self, query, key, value, mask=None, dropout=0.0):
        #based on: https://nlp.seas.harvard.edu/2018/04/03/attention.html
		d_k = query.size(-1)
		scores = torch.matmul(query, key.transpose(-2, -1)) / math.sqrt(d_k)
		p_attn = F.softmax(scores, dim = -1)
		p_attn = F.dropout(p_attn, p=dropout,training=self.training)
		return torch.matmul(p_attn, value), p_attn
	
	def forward(self, inputs):
        if self.useEmbeddings:
            output = self.embedding(inputs)
            output = output.permute(0,2,1)
        else:
		    output = inputs

		if self.usePE:
			output = self.pe(output)
		
		if self.useCNN:
			output = self.layer1(output)
			output = self.dropout1(output)
			output = output.permute(0,2,1)
		
		if self.useRNN:
			output, _ = self.RNN(output)
			F_RNN = output[:,:,:self.RNN_hiddenSize]
			R_RNN = output[:,:,self.RNN_hiddenSize:] 
			output = torch.cat((F_RNN,R_RNN),2)
			output = self.dropoutRNN(output)
		
		attn_concat = torch.Tensor([]).to(device)
		for i in range(0,self.numMultiHeads):
			query, key, value = self.Q[i](output), self.K[i](output), self.V[i](output)
			attnOut,p_attn = self.attention(query, key, value, dropout=0.2)
			attnOut = self.RELU[i](attnOut)
			if self.usepooling:
				attnOut = self.MAXPOOL[i](attnOut.permute(0,2,1)).permute(0,2,1)
			attn_concat = torch.cat((attn_concat,attnOut),dim=2)
		
		output = self.MultiHeadLinear(attn_concat)
		output = self.MHReLU(output)

		if self.readout_strategy == 'normalize':
			output = output.sum(axis=1)
			output = (output-output.mean())/output.std()
	
		output = self.fc3(output)	
		assert not torch.isnan(output).any()
		return output
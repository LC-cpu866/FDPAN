clear; clc;

data = xlsread('ringnorm\ringnorm.xlsx');
labels = xlsread('ringnorm\label.xlsx');

nclust = length(unique(labels));
numAnchor = 8;
k = 85;


tic;
cl = FDPAN(data, nclust, numAnchor, k);
rt = toc;

[F1, ACC, ARI, AMI, NMI] = evaluation_index(labels, cl);
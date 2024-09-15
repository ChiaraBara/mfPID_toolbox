clear all; close all; clc;

load('data_application.mat');
addpath([pwd,'\functions\']);

lab_serie = {'RESP','RR','SAP'};

%%% parameters
k = 5;
lags{1} = [0];
lags{2} = [0 1];
lags{3} = [0 1];

%%% embedding matrix contruction
B = mfPID_B_lags(data,1,{2,3},lags);

%%% PID measures
iy = [1:length(lags{1})];
ix1 = [length(lags{1})+1:length(lags{1})+length(lags{2})];
ix2 = [length(lags{1})+length(lags{2})+1:length(lags{1})+length(lags{2})+length(lags{3})];

out = mfPID_2sources_mixed_mex(B,iy,ix1,ix2,k);      
I = out.I;
I1 = out.I1;      
I2 = out.I2;
U1 = out.U1;
U2 = out.U2;
S = out.S;
R = out.R;
pY = out.pY;
Ispec = out.Ispec;
        

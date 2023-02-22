%% Neuron time series data source: https://osf.io/na4f9
%% Neuron connectivity data source: https://www.wormatlas.org/neuronalwiring.html#Connectivitydata

clear

%% Load the downloaded dataset and choose a particular time-series
load WT_NoStim
x=WT_NoStim(5);

%% List of neurons to use and their IDs in the dataset
names=["AVAL";"AVAR";"SMDVL";"SMDVR";"RIVL";"RIVR";"RIML";"RIMR";"SMDDL";"SMDDR";"AIBL";"AIBR";"RIBL";"RIBR";"AVBL";"AVBR"]';
IDs=[46,45,42,44,71,79,93,83,90,113,74,68,77,76,84,72];

%% Use only the time series from selected neurons for training the reservoir and subsequent analysis
x=x.deltaFOverF_bc;
data=x(:,IDs)';

%% Plot Time Series Data
figure

imagesc(data)
colorbar
yticks(1:1:16)
yticklabels(names)
xticks(1000:1000:3000)
xticklabels((1000:1000:3000)*0.3575)%duration between frames is 0.3575s from fps data
xlabel('Time (seconds)')
ylabel('Neuron Name')
set(gca,'FontSize',15)


%% Train a Reservoir Computer
%% Data for Reservoir training
measurements1 = data(:,1:end-1);% + z;
measurements2 = data(:, 2:end);

%% Build and train a reservoir computer
resparams=struct();
[num_inputs,~] = size(measurements1);
resparams.radius = 0.9; % spectral radius
resparams.degree = 5; % connection degree
approx_res_size = 3000; % reservoir size
resparams.N = floor(approx_res_size/num_inputs)*num_inputs; % actual reservoir size divisible by number of inputs
resparams.sigma = 0.4; % input weight scaling
resparams.bias=0;%bias
resparams.leakage=1;%leakage rate
resparams.train_length = size(data,2)-10; % number of points used to train
resparams.num_inputs = num_inputs; 
resparams.predict_length = 2000; % number of predictions after training
resparams.predict_length_max=resparams.predict_length;
resparams.beta = 10^-3;%10^(-10+ibeta); %regularization parameter

[xx, w_out, A, win,r] = train_reservoir(resparams, measurements1,measurements2);%Train and save w_out matrix for future use

%% Connectivity Estimation with RC
av_length=resparams.train_length-1;%number of time-steps to average connection strengths over

J=zeros(size(win,2));%This matrix stores inferred connection strengths
B=A+win*w_out;
for it=resparams.train_length-av_length:resparams.train_length
        
    xx=r(:,it);
    A2=B*xx;
        
    mat1=zeros(size(w_out));

    for i1=1:resparams.N
        mat1(:,i1)=w_out(:,i1)*(sech(A2(i1)))^2;
    end
    J=J+abs(mat1*(win))/(av_length);
end

%% Averaging left-right pair connections (note that we do not include self-connection strengths in analysis)
J=(J(1:2:end,1:2:end)+J(2:2:end,2:2:end)+J(1:2:end,2:2:end)+J(2:2:end,1:2:end))/4;
J=J-diag(diag(J));%set the diagonals to zero
% J is an 8-by-8 matrix now

%% Ground truth connection
load J0;%8-by-8 ground truth connectivity
J0=J0-diag(diag(J0));%set the diagonals to zero

%% Plot the adjacency matrices side-by-side
figure
subplot(1,2,1)
imagesc(J0)
pbaspect([1 1 1])
colorbar
colormap gray
xticks(1:1:8)
yticks(1:1:8)
xticklabels(["AVA","SMDV","RIV","RIM","SMDD","AIB","RIB","AVB"])
yticklabels(["AVA","SMDV","RIV","RIM","SMDD","AIB","RIB","AVB"])
xtickangle(45)
title("Ground Truth Connections")
set(gca,'FontSize',15)

subplot(1,2,2)
imagesc(J)
pbaspect([1 1 1])
colorbar
colormap gray
xticks(1:1:8)
yticks(1:1:8)
xticklabels(["AVA","SMDV","RIV","RIM","SMDD","AIB","RIB","AVB"])
yticklabels(["AVA","SMDV","RIV","RIM","SMDD","AIB","RIB","AVB"])
xtickangle(45)
title("Inferred Connections")
set(gca,'FontSize',15)
%% Load link scores and non-link scores
J=J(~eye(size(J)));%load off-diagonals of inferred matrix
J0=J0(~eye(size(J0)));%load off-diagonals of ground truth matrix

J_link=J(J0~=0);%store link scores
J_nonlink=J(J0==0);%store link scores

%% Plot link and non-link score histograms
figure
h = histogram([J_link;J_nonlink],20,'Normalization','count', 'DisplayStyle', 'stairs','LineWidth',3,'EdgeColor','k','EdgeAlpha',0.5,'LineStyle','-.');
hold on
histogram(J_link,h.BinEdges,'Normalization','count', 'DisplayStyle', 'stairs','LineWidth',3,'EdgeColor','r','EdgeAlpha',0.5,'LineStyle','-');
histogram(J_nonlink,h.BinEdges,'Normalization','count', 'DisplayStyle', 'stairs','LineWidth',3,'EdgeColor','b','EdgeAlpha',0.5,'LineStyle','-.');

legend('All scores','Link scores','Non-link scores')
xlabel('Reservoir computer score')
ylabel('Bin counts')
set(gca,'Fontsize',20)
hold off

%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, wout, A, win,states] = train_reservoir(resparams, data1,data2)

A = generate_reservoir(resparams.N, resparams.radius, resparams.degree);
q = resparams.N/resparams.num_inputs;
win = zeros(resparams.N, resparams.num_inputs);
for i=1:resparams.num_inputs
    rng(i+21,'twister');%random seed for reproducibility
    ip = resparams.sigma*(-1 + 2*rand(q,1));
    win((i-1)*q+1:i*q,i) = ip;
end
states = reservoir_layer(A, win, data1, resparams);
wout = train(resparams, states, data2(:,1:resparams.train_length));
x = states(:,end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = generate_reservoir(size, radius, degree)
   rng(10,'twister');%random seed for reproducibility
sparsity = degree/size;
while 1
A = sprand(size, size, sparsity);
e = max(abs(eigs(A)));

if (isnan(e)==0)%Avoid NaN in the largest eigenvalue, in case convergence issues arise
    break;
end

end
A = (A./e).*radius;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function states = reservoir_layer(A, win, input, resparams)

states = zeros(resparams.N, resparams.train_length);
for i = 1:resparams.train_length-1
    states(:,i+1) = (1-resparams.leakage)*states(:,i)+resparams.leakage*tanh(A*states(:,i) + win*input(:,i)+resparams.bias);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_out = train(params, states, data)

beta = params.beta;
rng(18,'twister');%random seed for reproducibility
idenmat = beta*speye(params.N);
w_out = data*transpose(states)*pinv(states*transpose(states)+idenmat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = load('dataStats.mat');
data = data.DataStats;

%% short latency responses
SLMeans = data(:,1);
SLstds = data(:,2);
SLstderr =  data(:,3);
LLMeans = (data(:,4));
LLstds = (data(:,5));
LLstderr =  data(:,6);
MSL= mean(nonzeros(SLMeans));
MLL = mean(nonzeros(LLMeans));
cellType = data(:,14);  % 0-OFF, 1-ON, 2-ON-OFF, 10000 Unknown
OFF = find(cellType == 0);
ON = find(cellType == 1);
ONOFF = find(cellType == 2);
unknown = find(cellType == 10000);

STDEfactor = 5;
scatterError(SLMeans(OFF), LLMeans(OFF), STDEfactor*SLstderr(OFF), STDEfactor*LLstderr(OFF));hold on;
scatterError(SLMeans(ON), LLMeans(ON), STDEfactor*SLstderr(ON), STDEfactor*LLstderr(ON));
scatterError(SLMeans(ONOFF), LLMeans(ONOFF), STDEfactor*SLstderr(ONOFF), STDEfactor*LLstderr(ONOFF));
scatterError(SLMeans(unknown), LLMeans(unknown), STDEfactor*SLstderr(unknown), STDEfactor*LLstderr(unknown));
axis([0 25 0 25]);
xx = [0 25]; yy=[0 25];
plot(xx,yy,'k--');
scatter(MSL, MLL, '^k')

%% Excitatory/inhibitory components
NumberExcite = data(:,7);
NumberInhib = data(:,8);
hist([NumberExcite NumberInhib])

%% ratio PC1/PC2
ratio = data(:,9);
Histcenters = [1 2 3 4 5 10 20 30];
h = hist(ratio, Histcenters);
bar(h)

%% RMSE improvement pc1 vs pc1,pc2
RMSE1 = nonzeros(data(:,10));
RMSE2 = nonzeros(data(:,11));
scatter(RMSE1, RMSE2);
xx = [0 0.2]; yy=[0 0.2];hold on;
plot(xx,yy,'k--');

%% Thresholds among cell types
positiveThreshold = data(:,12);
negativeThreshold = abs(data(:,13));
cellType = data(:,14);  % 0-OFF, 1-ON, 2-ON-OFF, 10000 Unknown

OFF = find(cellType == 0);
ON = find(cellType == 1);
ONOFF = find(cellType == 2);
unknown = find(cellType == 10000);

scatter(positiveThreshold(OFF), negativeThreshold(OFF), '*');hold all;
scatter(positiveThreshold(ON), negativeThreshold(ON), '^');
scatter(positiveThreshold(ONOFF), negativeThreshold(ONOFF), 'o');
scatter(positiveThreshold(unknown), negativeThreshold(unknown), '.');
plot([0 600],[0 600],'k--')

% stat p-value for differences between positive negative threshold among
% cell types
[h,p,ci,stats] = ttest2(positiveThreshold(OFF),negativeThreshold(OFF));p
[h,p,ci,stats] = ttest2(positiveThreshold(ON),negativeThreshold(ON));p
[h,p,ci,stats] = ttest2(positiveThreshold(ONOFF),negativeThreshold(ONOFF));p

% stat p-value for differences between between cell types
[h,p,ci,stats] = ttest2(positiveThreshold(OFF),positiveThreshold(ON)); p
[h,p,ci,stats] = ttest2(positiveThreshold(OFF),positiveThreshold(ONOFF)); p
[h,p,ci,stats] = ttest2(positiveThreshold(ON),positiveThreshold(ONOFF)); p



%% Comparison between single electrode stimulation and STA_P
Threshold_single = data(:,15);
Threshold_twoElec = data(:,16);
Threshold_threeElec = data(:,17);

Threshold_STAP = data(:,18);

for i=1:length(Threshold_single)
    if(i==14)
        continue;
    end
    a1 = Threshold_single(i);
    a2 = Threshold_twoElec(i);
    a3 = Threshold_threeElec(i);
    temp =[a1 a2 a3];
    bb = min(temp);
    numOfelec(i) = find( temp == bb );
    b(i) = bb;
end
b(14)=[];
numOfelec(14) = [];
Threshold_STAP(14) = [];

oneElec = find(numOfelec == 1);
TwoElec  = find(numOfelec == 2);
threeElec  = find(numOfelec == 3);

scatter(Threshold_STAP(oneElec), b(oneElec), '*');hold all;
scatter(Threshold_STAP(TwoElec), b(TwoElec), '^')
scatter(Threshold_STAP(threeElec), b(threeElec), 'o')
plot([0 600],[0 600],'k--')


%% Correlation coefficient between STAP STAN
CorrCoeff = data(:,19);
x = -1:0.1:0.1;
hist(CorrCoeff,x)


%% Electrode influence on cell: array 0,0 point is electrode 1

electrodesX = 0.875*[0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
electrodesY = [0 1 2 3 0.5 1.5 2.5 3.5 0 1 2 3 0.5 1.5 2.5 3.5 0 1 2 3];

CellX = [2.785 0.875 1.53125 1.3125 2.0125,...
    2.2 1.4 1.31 2.9 2,...
    3 2 2 -.2 1.6,...
    1.6 0.6 1.85 2.3 1.5];
CellY = [2.5 3 2.75 1.2 1.3,...
    1 2.25 1.75 1.5 2,...
    1.75 2.1 2.2 1.5 1.5,...
    2.7 2 1.8 2.25 1.75];

parameterFiles = {'2014Apr25\param_C1.mat', '2014May07\param_C2.mat', '2014May08\param_C3.mat',...
    '2014May11\param_C2.mat', '2014May17\param_C2.mat', '2014May21\param_C1.mat',...
    '2014May21\param_C2.mat','2014Sep01\param_C1.mat','2014Sep04\param_C2.mat',...
    '2014Sep08\param_C1.mat','2014Sep08\param_C2.mat','2014Sep10\param_C1.mat',...
    '2014Sep10\param_C2.mat','2014Sep17\param_C3.mat','2014Sep20\param_C1.mat',...
    '2014Sep20\param_C2.mat','2014Sep30\param_C1.mat','2014Sep30\param_C2.mat',...
    '2014Oct02\param_C1.mat','2014Oct02\param_C2.mat'};
clc;
cell2calc = [2 5 7:13 15:20];
rad_influence=[];
min_dist=[];

for i=cell2calc
    cell = [CellX(i) CellY(i)];
    cell = repmat(cell, 20, 1);
    dist_cell_elec = (electrodesX - cell(:,1)').^2 + (electrodesY - cell(:,2)').^2;
    dist_cell_elec=sqrt(dist_cell_elec);
    file = (['response\' parameterFiles{i}]);
    load(file);
    STA_P = abs(param.STA_P);
    STA_N = abs(param.STA_N);
    
    a = (STA_P).*dist_cell_elec;
    b = sum(a)/sum(STA_P);
    rad_influenceP(i) = b;%sum(a*min(dist_cell_elec));
    min_distP(i) = min(dist_cell_elec);
    
    a = (STA_N).*dist_cell_elec;
    b = sum(a)/sum(STA_N);
    rad_influenceN(i) = b;%sum(a*min(dist_cell_elec));
    min_distN(i) = min(dist_cell_elec);
    
end

mean(nonzeros(rad_influenceP));
mean(nonzeros(rad_influenceN));
min_distP;


figure;
scatter(nonzeros(rad_influenceP), nonzeros(rad_influenceN));hold on;
plot([0.5 1.8],[0.5 1.8],'--k')

[h,p,ci,stats] = ttest2(nonzeros(rad_influenceP),nonzeros(rad_influenceN))




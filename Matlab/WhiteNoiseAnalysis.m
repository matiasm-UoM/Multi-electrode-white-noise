function WhiteNoiseAnalysis( )
%GENERATERESPONSE Summary of this function goes here
%   Detailed explanation goes here

clc;
clear all;
close all;


%GUI variables
hFigure=[];
cmap =[];
haxes1=[];
haxes2=[];
DataFolder=[];     
dataNonEvoked=[];
LoadData=[];
dataEvoked=[];
ValidationData=[];
ValidationDataNotEvoked=[];
WNstd=[];
folders=[];
nbinsLabel=[];
TotalStim=[];
nbins=15;
cellLoc=[];
d1=[];
d2=[];
surfFigure=[];
StimPulseLength = 1.05e-3; %1.05 ms
workfolder = [];
LoadedCell=[];
FirstEig=[];
NextEig=[];
delay1 = [];
delay2 = [];

init();
makeGUI();



    %%---------------------------------------------------------------------
    function init()
        folders = dir('response');
        
        for i=3:length(folders)
            currentfilename = folders(i).name;
            fld{i-2} = currentfilename; 
        end
        folders = fld;
        fld=[];
        
    end
    
    
    %%---------------------------------------------------------------------
    function makeGUI()
        hFigure = figure('deleteFcn', @figureCloseCallback,...
            'NumberTitle', 'off', 'units', 'normalized',...
            'Name', 'Name ', 'Position',[0.05 0.3 0.6 0.6]);
        set(hFigure, 'Toolbar', 'figure');
        
        uipanel('Title','Parameters','FontSize',8,...
             'Position', [0.05 0.45 0.22 0.5]);
        DataFolder = uicontrol('Style', 'popupmenu', 'String', folders,...
            'units', 'normalized', 'position', [0.06 0.83 0.2 0.04],'FontSize',8, 'Value', 1);
        LoadData=uicontrol('Style', 'PushButton', 'String', 'Load Data',...
            'units', 'normalized', 'position', [0.06 0.73 0.20 0.04],'FontSize',10,...
            'callback', @LoadDataFiles);
%         delaystr = ['0-' num2str(delay1*1000) ' ms|' num2str(delay1*1000) '-' num2str(delay2*1000) ' ms|',...
%             num2str(delay2*1000) '-' num2str(delay3*1000) ' ms|' num2str(delay3*1000) '-' num2str(delay4*1000) ' ms'];
%         Delay = uicontrol('Style', 'popupmenu', 'String', delaystr,...
%             'units', 'normalized', 'position', [0.06 0.78 0.2 0.04],'FontSize',8, 'Value', 1);
        d1 = uicontrol('Style', 'edit', 'String', '0',...
            'units', 'normalized', 'position', [0.06 0.78 0.05 0.04],'FontSize',8,...
            'callback', @SetDelay1);
        d2 = uicontrol('Style', 'edit', 'String', '5',...
            'units', 'normalized', 'position', [0.11 0.78 0.05 0.04],'FontSize',8,...
            'callback', @SetDelay2);   
        uicontrol('Style', 'text', 'String', 'ms',...
            'units', 'normalized', 'position', [0.16 0.78 0.05 0.04],'FontSize',8);
        
        
        haxes1=axes('units', 'normalized', 'position', [0.06 0.1 0.20 0.3]);
        
        haxes2=axes('units', 'normalized', 'position', [0.3 0.1 0.6 0.85]);
        
    end

    %%---------------------------------------------------------------------
    function SetDelay1(gcbo, eventdata)
        delay1 = get(d1, 'String');
        delay1=str2double(delay1);
    end

    %%---------------------------------------------------------------------
    function SetDelay2(gcbo, eventdata)
        delay2 = get(d2, 'String');
        delay2=str2double(delay2);
    end
    %%---------------------------------------------------------------------
    function LoadDataFiles(gcbo, eventdata);
        clc;
        delete(figure(2));
        delete(figure(3));
        delete(figure(4));
        delete(figure(5));
        SetDelay1(gcbo, eventdata);
        SetDelay2(gcbo, eventdata);
        delay1 = delay1/1000 + StimPulseLength;
        delay2 = delay2/1000 + StimPulseLength;
        
         folder = get(DataFolder,'Value');
         str = folders{folder};
         
         file = ['response\' str];
         C = findFiles(file, 'Evoked_M');
         C=C/2;
         if(C>1)
             [i indices] = findFiles(file, 'Evoked_M');

             inputText = ['Cell to load (' mat2str(indices) ' cells available):'];
             C2Load = input(inputText);
         else
             [i indices] = findFiles(file, 'Evoked_M');
             C2Load = indices;
         end
           
         workfolder = file;
         LoadedCell = C2Load;
         file1 = [file '\Evoked_M' num2str(C2Load) '.txt'];
         file2 = [file '\NonEvoked_M' num2str(C2Load) '.txt']; 

         wavegratingfile1 = [file '\EvokedGrating_M' num2str(C2Load) '.txt'];
         wavegratingfile2 = [file '\NonEvokedGrating_M' num2str(C2Load) '.txt'];
             
         
         cmap = load('cmap.mat');
         cmap = cmap.cmap;
         
         dataEvoked=[];
         dataNonEvoked=[];
         TotalStim=[];
         [h1, dataEvoked] = hdrload(file1);

         cellLoc= h1(2,:)
         cellLoc = cellLoc(9:end);
         cellLoc = [str2double(cellLoc(1:3)) str2double(cellLoc(5:7)) ];
         if(isnan(cellLoc))
             cellLoc=[];
         end
         
         
         %assignin('base', 'dataEvoked', dataEvoked);
         [h2, dataNonEvoked] = hdrload(file2);
         %Remove any shorting stimuli
         [rw, cl] = find(dataNonEvoked(:,:) == -10000);
         dataNonEvoked(rw, :)=[];
         [rw, cl] = find(dataEvoked(:,:) == -10000);
         dataEvoked(rw, :)=[];
         
         SpikeTimes = dataEvoked(:, 21:end);
         dataEvoked(:,21:end)=[];
         [DEvoked, datanotevoked]=FindDelayEvokedStimuli(dataEvoked, dataNonEvoked, SpikeTimes);
         dataEvoked = DEvoked;
         dataNonEvoked = datanotevoked;
         
         %remove column of zero values
%          id=[];
%          for i=1:20
%              if(dataEvoked(1,i)==0)
%                  id(end+1)=i;
%              end
%          end
%          dataEvoked(:,id)=[];
%          dataNonEvoked(:,id)=[];
         assignin('base', 'dataEvoked', dataEvoked);
         assignin('base', 'dataNonEvoked', dataNonEvoked);

         
         
         %ONLY USE Percent2use % of the data to analyse, remainder to verify
         Pecent2Use = 0.80;
         
         if(exist(wavegratingfile1) ==2)
             useGrating = input('Use wave grating as validation data?: ', 's');
             strcontains = strfind(['y', 'e', 's'],useGrating);
             
             if(isempty(strcontains))
                 [A B] = size(dataEvoked);
                 NoStimTriggeredSpikes = A
                 Ap = floor(A*Pecent2Use);
                 A2=floor(linspace(1, NoStimTriggeredSpikes, NoStimTriggeredSpikes-Ap));
                 ValidationData = dataEvoked(A2, :);
                 ValidationDatalength=length(ValidationData)
                 dataEvoked(A2,:) = [];

                 [A B] = size(dataNonEvoked);
                 Ap = floor(A*Pecent2Use);
                 A2=floor(linspace(1, A, A-Ap));
                 ValidationDataNotEvoked = dataNonEvoked(A2, :);
                 dataNonEvoked(A2,:) = [];

                 TotalStim = [dataEvoked; dataNonEvoked];
%                  ValidationDataNotEvoked = [ValidationDataNotEvoked; ValidationData];
                 
             else
                 % USE WAVE GRATING DATA TO VALIDATE MODEL
                 [A B] = size(dataEvoked);
                 NoStimTriggeredSpikes = A
                 
                 [h1, wg1] = hdrload(wavegratingfile1);
                 [h2, wg2] = hdrload(wavegratingfile2);
                 
                 ValidationData = wg1;
                 size(ValidationData)
                 ValidationDataNotEvoked = wg2;
%                  delay1 = 10e-3;
                 SpikeTimes = ValidationData(:, 21:end);
                 ValidationData(:,21:end)=[];
                 [DEvoked, datanotevoked]=FindDelayEvokedStimuli(ValidationData, ValidationDataNotEvoked, SpikeTimes);
                 ValidationData = DEvoked;
                 ValidationDataNotEvoked = datanotevoked;
                 
                 TotalStim = [dataEvoked; dataNonEvoked];
%                  ValidationDataNotEvoked = [ValidationDataNotEvoked; ValidationData];
                 
             end
             
         else
             
                [A B] = size(dataEvoked);
                 NoStimTriggeredSpikes = A
                 Ap = floor(A*Pecent2Use);
                 A2=floor(linspace(1, NoStimTriggeredSpikes, NoStimTriggeredSpikes-Ap));
                 ValidationData = dataEvoked(A2, :);
                 ValidationDatalength=length(ValidationData)
                 dataEvoked(A2,:) = [];

                 [A B] = size(dataNonEvoked);
                 Ap = floor(A*Pecent2Use);
                 A2=floor(linspace(1, A, A-Ap));
                 ValidationDataNotEvoked = dataNonEvoked(A2, :);
                 dataNonEvoked(A2,:) = [];

                 TotalStim = [dataEvoked; dataNonEvoked];
%              ValidationDataNotEvoked = [ValidationDataNotEvoked; ValidationData];
         end
         
         
         uicontrol('Style', 'Text', 'String', '# of bins',...
            'units', 'normalized', 'position', [0.06 0.68 0.2 0.04],'FontSize',10,'Visible','on');
        nbinsLabel=uicontrol('Style', 'edit', 'String', '15',...
            'units', 'normalized', 'position', [0.06 0.63 0.20 0.04],'FontSize',10);
         CompCov=uicontrol('Style', 'PushButton', 'String', 'Compute parameters',...
            'units', 'normalized', 'position', [0.06 0.58 0.20 0.04],'FontSize',10,...
            'callback', @ComputeParameters);
         
         axes(haxes1);
         cla;
         axes(haxes2);
         cla;
         
    end


        %%---------------------------------------------------------------------
    function [DEvoked, datanotevoked]=FindDelayEvokedStimuli(DEvoked, datanotevoked, SpikeTimes)
        
            [A B] = find(SpikeTimes(:,1)>delay1 & SpikeTimes(:,1)<=delay2);
            assignin('base','A',A);
            assignin('base','B',B);
            aa=[];index=1;
            for ii=1:length(A)
                if(B(ii)~=1)
                    aa(index)=ii;
                    index=index+1;
                end
            end
            A(aa)=[];
            
            temp1 = DEvoked;
            temp1(A,:) =[]; %stim not evoking spike within delay1
            temp2 = DEvoked(A,:); %stim evoking a spike within delay1
            DEvoked = temp2;
            datanotevoked = [datanotevoked;temp1];
         
    
    
    end


        %%---------------------------------------------------------------------
    function ComputeParameters(gcbo, eventdata);
%         clc;
        axes(haxes1);
        cla;
        axes(haxes2);
        cla;
        nbins = str2double(get(nbinsLabel, 'String'));
        
        EvokedCov = cov(dataEvoked);
        [ST_EigVect, ST_EigVal, ST_EXPLAINED] = pcacov(EvokedCov);
         % compute change in covariance
        
        covT = cov(TotalStim);
        delCov = EvokedCov - covT;
        [EigVect, EigVal, EXPLAINED] = pcacov(delCov);
        PC1 = EXPLAINED(1)
        PC2 = EXPLAINED(2)
        PCA1PCA2ratio = EXPLAINED(1)/EXPLAINED(2)
        
        
        %Divide stim evoked into positive, negative projection against PC1
%         figure(33);
%         plot(EigVect(:,1));
        eigVectnumber =1;
        arbitraryEigVect = 2;
        FirstEig=eigVectnumber;
        NextEig=arbitraryEigVect;
%         EigVect(:,1) = zeros(20,1);
%         EigVect(6,1) = 1;
%         EigVect(:,2) = zeros(20,1);
%         EigVect(11,2) = 1;
        
        assignin('base','EigVect',EigVect);
        mxE1 = max(EigVect(:,eigVectnumber));
        mnE1 = min(EigVect(:,eigVectnumber));
        if(abs(mxE1)<abs(mnE1))
            thisIsTrue = 1
            EigVect = -EigVect;
        end
       
        evoked = dataEvoked*EigVect(:,eigVectnumber);
        evokedP = evoked>=0;
        evokedP = dataEvoked(evokedP,:);
        evokedN = evoked<0;
        evokedN = dataEvoked(evokedN,:);

        %Divide total stim into positive, negative projection against PC1
        Total = TotalStim*EigVect(:,eigVectnumber);
        TotalP = Total>=0;
        TotalP = TotalStim(TotalP,:);
        TotalN = Total<0;
        TotalN = TotalStim(TotalN,:);

        %project onto PC1, PC2 plane for visualisation
        Proj_N_PC1 = evokedN*EigVect(:,eigVectnumber);
        Proj_N_PC2 = evokedN*EigVect(:,arbitraryEigVect);
        Proj_P_PC1 = evokedP*EigVect(:,eigVectnumber);
        Proj_P_PC2 = evokedP*EigVect(:,arbitraryEigVect);

        Proj_T_PC1 = TotalStim*EigVect(:,eigVectnumber);
        Proj_T_PC2 = TotalStim*EigVect(:,arbitraryEigVect);

        axes(haxes2);
        figure(10);
        scatter(Proj_T_PC1, Proj_T_PC2, 'k.');hold on;
        scatter(Proj_P_PC1, Proj_P_PC2,'r+', 'LineWidth', 2);
        scatter(Proj_N_PC1, Proj_N_PC2,'b*', 'LineWidth', 2);
        set(gca,'Color',[1 1 1]);
        axis([-500 500 -500 500]);
        
        
        % Find STA for positive and negative
            %ensemble mean
        MeanT = mean(TotalStim);
        Proj_MeanT_P1 = MeanT*EigVect(:,eigVectnumber);
        Proj_MeanT_P2 = MeanT*EigVect(:,arbitraryEigVect);

        STA_P = mean(evokedP);
        Proj_STAP_P1 = STA_P*EigVect(:,eigVectnumber);
        Proj_STAP_P2 = STA_P*EigVect(:,arbitraryEigVect);

        STA_N = mean(evokedN);
        Proj_STAN_P1 = STA_N*EigVect(:,eigVectnumber);
        Proj_STAN_P2 = STA_N*EigVect(:,arbitraryEigVect);
        
        
        scatter(Proj_STAP_P1, Proj_STAP_P2, 500, 'd', 'LineWidth', 1.5,...
            'MarkerEdgeColor','k', 'MarkerFaceColor','w');
        scatter(Proj_STAN_P1, Proj_STAN_P2, 500, 'o','LineWidth', 1.5,...
            'MarkerEdgeColor','k', 'MarkerFaceColor','w');
        scatter(Proj_MeanT_P1, Proj_MeanT_P2, 'filled', 'w');
%         plot([Proj_STAN_P1, Proj_MeanT_P1], [Proj_STAN_P2, Proj_MeanT_P2], 'b');
%         plot([Proj_STAP_P1, Proj_MeanT_P1], [Proj_STAP_P2, Proj_MeanT_P2], 'b');
        xlabel('Principal component 1');
        ylabel('Arbitrary component');

        STA_vector_P = (STA_P - MeanT);
        STA_vector_P = STA_P./norm(STA_P);
        STA_vector_N = (STA_N - MeanT);
        STA_vector_N = -STA_N./norm(STA_N);

        
        %project evokedP onto STA_vector_P and bin
        Proj_P_STAP = evokedP*STA_vector_P';
        Proj_P_STAP = sort(Proj_P_STAP);
        Proj_P_Total = TotalP*STA_vector_P';
        Proj_P_Total= sort(Proj_P_Total);
        
        Proj_N_STAN = evokedN*STA_vector_N';
        Proj_N_STAN = sort(Proj_N_STAN);
        Proj_N_Total = TotalN*STA_vector_N';
        Proj_N_Total= sort(Proj_N_Total);
        
        
              
        [xx, Pspike]=findPspike([Proj_N_STAN; Proj_P_STAP],[Proj_N_Total; Proj_P_Total], 2*nbins,'');
        NanPspk = find(isnan(Pspike));
        Pspike(NanPspk) =[];
        xx(NanPspk) =[];
        axes(haxes1);
         figure(11);
         plot(xx, Pspike, 'ob');hold on;
         % Set up fittype and options.
        ft = fittype( 'c/(1+exp(-a*(x+b))) +  f-(f/(1+exp(-d*(x+e)))) + s', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( ft );
        opts.Display = 'Off';
        opts.Lower = [0 -1000 0 0 -10 0 0];
        opts.StartPoint = [0.1 -150 0.8 0.1 300 0.8 0.1];
        opts.Upper = [1 10 1 1 1000 1 1];
        [fitresult, gof] = fit(xx(:),Pspike(:), ft, opts );
        gof
        aP=fitresult.a;
        bP=fitresult.b;
        cP=fitresult.c;
        aN=fitresult.d;
        bN=fitresult.e;
        cN=fitresult.f;
        spont = fitresult.s;
%         const = fitresult.g
        disp(['aP =' num2str(aP) ';']);
        disp(['bP =' num2str(bP) ';']);
        disp(['cP =' num2str(cP) ';']);
        disp(['aN =' num2str(aN) ';']);
        disp(['bN =' num2str(bN) ';']);
        disp(['cN =' num2str(cN) ';']);
        disp(['spont =' num2str(spont) ';']);
        
        param.STA_P = STA_vector_P; param.STA_N = STA_vector_N;
        param.aP = aP; param.bP = bP; param.cP = cP;
        param.aN = aN; param.bN = bN; param.cN = cN; 
        param.spont = spont;
%         param.const=const;
        assignin('base', 'param', param);
        %save model parameters
        fname = [workfolder '\param_C' num2str(LoadedCell) '.mat'];
        save(fname, 'param');
        
        xlim([-500 500])
        plot(fitresult,'--b' );
        legend off;
        xlabel('Current (uA)', 'FontSize', 20);
        ylabel('Spike probability', 'FontSize', 20);
        ylim([0 1.1]);
        set(gca,'FontSize', 20);
        
        
        % significance
        std_errorP = std(evokedP)/sqrt(length(evokedP));
        std_errorN = std(evokedN)/sqrt(length(evokedN));
        electrodes = 1:20;
        %make a vector containing only significant electrodes
        k = 10;
        for jj=1:length(electrodes)
            l1 = k*std_errorP(jj);
            l2 = abs(STA_P(jj));
            if(l1>l2)
                signifSTAP(jj) = 0;
            else
                signifSTAP(jj) = STA_P(jj);
            end
            l1 = k*std_errorN(jj);
            l2 = abs(STA_N(jj));
            if(l1>l2)
                signifSTAN(jj) = 0;
            else
                signifSTAN(jj) = STA_N(jj);
            end
        end
        assignin('base', 'STA_N', STA_N);
        assignin('base', 'STA_P', STA_P);
        
        %Calculate correlation coefficient of STA_P, STA_N
        A = [STA_P(:) STA_N(:)];
        R = corrcoef(A)
        
        
        % Create surface probability plot
        Generate2dSpikePlot(EigVect);
        
        
       
%         % Information content
%         C = bsxfun(@minus,evokedP,STA_P); 
%         covP = cov(C);
%         C = bsxfun(@minus,evokedN,STA_N);
%         covN = cov(C);
%         CalcInformation(covT, covP, covN, STA_vector_P, STA_vector_P,EigVect, evokedP, evokedN);  
        
        
        
        
        
        % --- Validate model ---%
        %Divide validation data into positive, negative projection against PC1
        
        ValidationP = (ValidationData*STA_vector_P')>=0;
        ValidationP = ValidationData(ValidationP,:);
        ValidationP = ValidationP*STA_vector_P';
        ValidationP = sort(ValidationP);
        
        ValidationN = (ValidationData*STA_vector_N')<0;
        ValidationN = ValidationData(ValidationN,:);
        ValidationN = (ValidationN*STA_vector_N');
        ValidationN = sort(ValidationN);
         
        VNotEvokedP = (ValidationDataNotEvoked*STA_vector_P')>=0;
        VNotEvokedP = ValidationDataNotEvoked(VNotEvokedP,:);
        VNotEvokedP = VNotEvokedP*STA_vector_P';
        VNotEvokedP = sort(VNotEvokedP);
        VNotEvokedN = (ValidationDataNotEvoked*STA_vector_N')<0;
        VNotEvokedN = ValidationDataNotEvoked(VNotEvokedN,:);
        VNotEvokedN = VNotEvokedN*STA_vector_N';
        VNotEvokedN = sort(VNotEvokedN);

        
%         [xx1, PspikeV_P]=findPspike(ValidationP,VNotEvokedP, nbins, 'v');
%         [xx2, PspikeV_N]=findPspike(ValidationN,VNotEvokedN, nbins, 'v');
        [xx, PspikeV]=findPspike([ValidationN;ValidationP],[VNotEvokedN;VNotEvokedP], 2*nbins, 'v');
        
%         PspikeV=[PspikeV_N PspikeV_P];
        NanPspk = find(isnan(PspikeV));
        PspikeV(NanPspk) = [];
        xx(NanPspk) =[];

        
         %find the predicted output
        validationData = [ValidationN;ValidationP;VNotEvokedP;VNotEvokedN];
        Predicted = (cP./(1+exp(-aP.*(validationData+bP))) ) + (cN - (cN./(1+exp(-aN.*(validationData+bN)))) ) + spont;
       
        
        % Set up fittype and options.
        ft = fittype( 'c/(1+exp(-a*(x+b))) +  f-(f/(1+exp(-d*(x+e))))  + s', 'independent', 'x', 'dependent', 'y' );
%         ft = fittype( 'c/(1+exp(-a*(x+b))) ', 'independent', 'x', 'dependent', 'y' );
%         opts = fitoptions( ft );
%         opts.Display = 'Off';
%         opts.Lower = [0 -1000 0 ];
%         opts.StartPoint = [0.01 -150 0.8];
%         opts.Upper = [1 10 1];
        opts.Display = 'Off';
        opts.Lower = [0 -1000 0 0 -100 0 0];
        opts.StartPoint = [0.1 -150 0.8 0.1 300 0.8 0.1];
        opts.Upper = [1 100 1 1 1000 1 1];
        [fitresultV, gofV] = fit(xx(1:end-1)',PspikeV(1:end-1)', ft, opts );
        gofV
        a1=fitresultV.a
        b1=fitresultV.b
        c1=fitresultV.c
        a2=fitresultV.d
        b2=fitresultV.e
        c2=fitresultV.f
        sp2 = fitresultV.s
%         constV = fitresultV.g
        

         
        VNE = [VNotEvokedN;VNotEvokedP];
%         VTotal = [validationData; VNE];
        VTotal = [ValidationData;ValidationDataNotEvoked];
        VTotalp = [ValidationData;ValidationDataNotEvoked]*STA_vector_P';
        VTotaln = [ValidationData;ValidationDataNotEvoked]*STA_vector_N';
        Vpredicted = (cP./(1+exp(-aP.*(VTotalp+bP))) ) + (cN - (cN./(1+exp(-aN.*(VTotaln+bN)))) ) +spont;
        
%         assignin('base', 'VTotalp', VTotalp);
%         assignin('base', 'VTotaln', VTotaln);
%         aP2 =0.018947;
%         bP2 =-138.0917;
%         cP2 =0.99732;
%         aN2 =0.020273;
%         bN2 =174.1758;
%         cN2 =1;
%         dt2 = load('VTotal2.mat');
%         VTotal2n = dt2.VTotaln;
%         VTotal2p = dt2.VTotalp;
%         Predicted2 =(cP2./(1+exp(-aP2.*(VTotal2p+bP2))) ) + (cN2 - (cN2./(1+exp(-aN2.*(VTotal2n+bN2)))) ) ;      
%         Vpredicted=(1/PCA1PCA2ratio)*Vpredicted+(1-1/PCA1PCA2ratio)*Predicted2;
        
        EvokedV = ones(1, length(ValidationData));
        NEvokedV = zeros(1, length(ValidationDataNotEvoked));
        VT = [EvokedV NEvokedV];
        vbins=10;%nbins;
        step = 1/vbins;
        for ii=1:vbins
            x1 = 0 + (ii-1)*step;
            x2 = ii*step;
            evoked = (Vpredicted >= x1 & Vpredicted < x2);
            ev1 = VT(evoked);
            lenEv1 = length(ev1);
            
            
            PspikeVT(ii) = mean(ev1);
            xxVT(ii) = 0.5*(x2+x1);
            
            if(lenEv1<30)
                PspikeSTE(ii) = std(ev1);
            else
                PspikeSTE(ii) = std(ev1)/sqrt(lenEv1);
            end

        end
        
        
        
        % plot figures
        plotfigures=1;
        if(plotfigures==1)
            figure(4);
            plot(electrodes,evokedP); hold all;
            errorbar(electrodes, STA_P, std_errorP,'k', 'LineWidth', 3);
            title('STA P', 'FontSize', 20);
            xlim([0 21]);
            xlabel('Electrode number');
            ylabel('Current amplitude (uA)');
            figure(5);
            plot(electrodes,evokedN); hold all;
            errorbar(electrodes, STA_N, std_errorN,'k', 'LineWidth', 3);
            title('STA N','FontSize', 20);
            xlim([0 21]);
            xlabel('Electrode number');
            ylabel('Current amplitude (uA)');

            %Plot STA_N and STA_P
            figure(2);
            PlotResponseStimulus_disc(signifSTAP);
            colormap(cmap);h=colorbar;
            title(h, 'Pulse amplitude (uA)','FontSize', 20, 'Rotation', -90);
            hold on;
            if(~isempty(cellLoc))
                scatter(cellLoc(1)*1.5, cellLoc(2), 40, 'filled', 'g');
            end
            mxP = max(STA_P);
            mnP = min(STA_P);
            if(abs(mnP)>abs(mxP))
                caxis([-abs(mnP) abs(mnP)+1]);
            elseif(abs(mxP)>abs(mnP))
                caxis([-abs(mxP) abs(mxP)]+1);
            end
            title('Positive stimulus selectivity', 'FontSize', 20);
            set(gca,'FontSize', 20);
            caxis([-200 200]);
            figure(3);
            PlotResponseStimulus_disc(signifSTAN);colormap(cmap);
            h=colorbar;
            title(h, 'Pulse amplitude (uA)','FontSize', 20, 'Rotation', -90);
            hold on;
            if(~isempty(cellLoc))
                scatter(cellLoc(1)*1.5, cellLoc(2), 40, 'filled', 'g');
            end
            mxN = max(STA_N);
            mnN = min(STA_N);
            if(abs(mnN)>abs(mxN))
                caxis([-abs(mnN) abs(mnN)+1]);
            elseif(abs(mxN)>abs(mnN))
                caxis([-abs(mxN) abs(mxN)+1]);
            end
            title('Negative stimulus selectivity', 'FontSize', 20);
            set(gca,'FontSize', 20);
            caxis([-200 200]);
            
            figure(6);
%             subplot(2, 1, 1);
%             h = bar(ST_EXPLAINED);
%             set(h, 'FaceColor', [0.5 0.5 0.5]);
%             title('Amount of variance explained by principal components');
%             subplot(2,1,2);
            h = bar(EXPLAINED);
            set(h, 'FaceColor', [0.5 0.5 0.5]);
            title('Change in covariance expalined by each component');
            xlabel('Principal component');
            ylabel('% change');
            
            figure();
%             scatter(xx, PspikeV, 'filled','bo');hold on;
            xR = -600:600;
            validationModel = feval(fitresultV, xR);
            
            
%             scatter(xx, PspikeV,'ob');
            plot(validationData,Predicted, '*k', 'LineWidth', 3);hold on;
            h = plot(xR,validationModel, '--b' );
            set(h, 'LineWidth', 3);
            legend('Model predicted','Validation fit' );
            xlabel('Amplitude','FontSize', 20);
            ylabel('Spike probability','FontSize', 20);
            set(gca,'FontSize', 20);
            ylim([0 1.05]);
            
            figure();
            scatter(Vpredicted, VT);
            figure(19);
    %         scatter(xxVT, PspikeVT, 'k');
            hold all;
%             plot(xxVT, PspikeVT, 'LineWidth', 3);
            errorbar(xxVT, PspikeVT, PspikeSTE, 'LineWidth', 3);
            plot(0:0.1:1,0:0.1:1, 'k--') 
            xlabel('Predicted spike probability');
            ylabel('Actual spike probability');

        end
        
            
        MSE = (PspikeVT - (xxVT)).^2;
        Nanmse = isnan(MSE);
        MSE(Nanmse) =[];
        MEANMSE = mean(MSE)
        RMSE = sqrt(MEANMSE)
        
        
    end
    %%---------------------------------------------------------------------
    function Generate2dSpikePlot(EigVect)
        
        projEvokedPC1 = dataEvoked*EigVect(:,FirstEig);
        projEvokedPC2 = dataEvoked*EigVect(:,NextEig);
        
        projTotalPC1 = TotalStim*EigVect(:,FirstEig);
        projTotalPC2 = TotalStim*EigVect(:,NextEig);
        
        
        step1=(max(projEvokedPC1)-min(projEvokedPC1) )/(2*nbins+1);
        step2=(max(projEvokedPC2)-min(projEvokedPC2) )/(2*nbins+1);
        id=1;
        surfFigure = figure(20);clf;
        axes('units', 'normalized', 'position', [0.05 0.05 0.75 0.75]);colormap((gray));
        grid off;
%         scatter(projEvokedPC1, projEvokedPC2, 'r+');hold on;
%         scatter(projTotalPC1, projTotalPC2,'k.', 'LineWidth', 2);
        for i=1:2*nbins+2
            for j=1:2*nbins+2
                
                x1=min(projEvokedPC1)+(i-1)*step1;
                x2=min(projEvokedPC1)+(i)*step1;
                y1=min(projEvokedPC2)+(j-1)*step2;
                y2=min(projEvokedPC2)+(j)*step2;
                 
%                 plot([x1 x1], [y1 y2], 'k--'); hold on;
%                 plot([x1 x2], [y1 y1], 'k--')
                
                spikesInBinX = find(projEvokedPC1>=x1 & projEvokedPC1<x2);
                ycomponent = projEvokedPC2(spikesInBinX);
                spikesInBox = find(ycomponent>=y1 & ycomponent<y2);
                spikesInBox = length(spikesInBox);
                
                TotalInBinX = find(projTotalPC1>=x1 & projTotalPC1<x2);
                ycomponent = projTotalPC2(TotalInBinX);
                TotalInBox = find(ycomponent>=y1 & ycomponent<y2);
                TotalInBox = length(TotalInBox);
                
                pspike(j,i) = spikesInBox/TotalInBox;
                X(j,i) = (x1 + x2)/2;
                Y(j,i) = (y1 +y2)/2;
                
                if(isnan(pspike(j,i)))
                    pspike(j,i) = NaN;
                else
                    x(id) = (x1 + x2)/2;
                    y(id) = (y1 +y2)/2;
                    z(id) = pspike(j,i);
                    id=id+1;
                    patch([x1+step1/2 x1+step1/2 x2+step1/2 x2+step1/2],...
                        [y1+step2/2 y2+step2/2 y2+step2/2 y1+step2/2], pspike(j,i));
                end
                
            end
        end
        
        xlims = [min(projEvokedPC1) max(projEvokedPC1)]
        xlim(xlims);ylim([min(projEvokedPC2) max(projEvokedPC2)])
        xlims = [min(min(X)) max(max(X))]
        xlim(xlims);ylim([min(min(Y)) max(max(Y))])
%         set(gca, 'XTick', []);
%         set(gca, 'YTick', []);
        view(0, 90);
        xlabel('PC1');
        ylabel('Arbitrary component');
        
        assignin('base', 'X', X);
        assignin('base', 'Y', Y);
        assignin('base', 'pspike', pspike);
        
        % NOTE: the bar plot on top and sides may appear to extend past the
        % surface plot. This is because of the low sampling in these
        % regions, the surf plot cannot plot properly. 
        
        figure(surfFigure);
        evokedT = sort(projEvokedPC1); TotalT = sort(projTotalPC1);
        axes('position', [0.05 0.82 0.75 0.16]); 
        x1 = linspace(min(projEvokedPC1), max(projEvokedPC1),2*nbins+2);
        y1 = hist(projEvokedPC1, x1 );
        y2 = hist(projTotalPC1, x1 );
        bar(x1,y2/(max(y2)),'FaceColor',[0.7 0.7 0.7]); hold on;
        bar(x1,y1/(max(y2)));xlim(xlims);    
         
        
        axes('position', [0.82 0.05 0.16 0.75]); 
        x1 = linspace(min(projEvokedPC2), max(projEvokedPC2),2*nbins+2);
        y1= hist(projEvokedPC2,x1); 
        y2 = hist(projTotalPC2, x1);
        bar(x1,y2/max(y2),'FaceColor',[0.7 0.7 0.7]); hold on;
        bar(x1,y1/max(y2)); xlim([min(min(Y)) max(max(Y))]); set(gca, 'view', [90 90]);
        
    end

%%---------------------------------------------------------------------
    function CalcInformation(covT,covP, covN, staP, staN, EigVect, evokedP, evokedN)
        assignin('base', 'staP', staP);
        assignin('base', 'staN', staN);
        assignin('base', 'evokedP', evokedP);
        lenSTA = length(staP);
        [EigVectP eigVal explained] = pcacov(covP);
        EigVectN = pcacov(covN);
        
        for i=1:1
            subspace = EigVectP(:,1:i);
            total0 = TotalStim*subspace;
            evoked1 = evokedP*subspace;
            
            mean0 = mean(total0);
            mean1 = mean(evoked1);
            C = bsxfun(@minus,total0,mean0); 
            cov0 = cov(C);
            C = bsxfun(@minus,evoked1,mean1); 
            cov1 = cov(C);
            
            BBT = cov1^-1*(cov1 + mean1'*mean1)*cov1
%             
%             DQP_P(i) = 0.5*( log(det(cov0)/det(cov1)) - i + trace(cov0^-1*cov1) +...
%                 (mean0-mean1)*cov0^1*(mean0-mean1)');
            
            
%             a = trace(cov1^-1 * cov0);
%             b = (mean1 - mean0)*cov1^-1*(mean1-mean0)';
%             c = i;
%             d = log(det(cov0)/det(cov1));
%             DQP_P(i) = (0.5*(a +b-c-d))/log(2);
            
            
%             a = trace(cov1+mean1*mean1');
%             d = log(det(cov1));
%             DQP_P(i) = (0.5*(a -d + mean1*mean1' - i))/log(2);
%             
            
%             subspace = EigVectP(:,1:i);
%             a = trace(subspace'*(covP+(staP'*staP))*subspace);
%             b = log(det(subspace'*covP*subspace));
%             DQP_P(i) = 0.5*(a - b -i);
            
%             subspace = EigVectN(:,1:i);
%             a = trace(subspace'*(covN+staN'*staN)*subspace);
%             b = log(det(subspace'*covN*subspace));
%             DQP_N(i) = 0.5*(a - b -i);
        end
%         assignin('base', 'DQP_P',DQP_P);
%         figure;
%         subplot(2,1,1);
%         bar(real(DQP_P))
%         subplot(2,1,2);
%         bar(real(DQP_N))
    end

    %%---------------------------------------------------------------------
    function [xx, Pspike]=findPspike(Projection,Proj_Total, bins, validation)
        
        %         step = (max(Proj_N_STAN)-min(Proj_N_STAN) )/nbins;
%         xn = linspace(min(Proj_N_STAN),max(Proj_N_STAN),nbins);
%         for ii=1:nbins
% 
%             if(ii==1)
%                 x1 = max(Proj_N_STAN);
%                 x2 = max(Proj_N_STAN)-ii*step;
%             elseif(ii==nbins)
%                 x1 = max(Proj_N_STAN)-(ii-1)*step;
%                 x2 = max(Proj_N_STAN)-ii*step-0.1;
%             else
%                 x1 = max(Proj_N_STAN)-(ii-1)*step;
%                 x2 = max(Proj_N_STAN)-ii*step;
%             end
% 
%             pS = length(find(Proj_N_STAN>x2 & Proj_N_STAN<=x1));
%             if(isempty(pS))
%                 pS=0;
%             end
%             pT = length(find(Proj_N_Total>x2 & Proj_N_Total<=x1));
%             PspikeN(nbins+1-ii) = pS/pT;
% 
%         end
        
        step = floor(size(Projection)/bins);
        step = step(1);
        xx=[];
        b1 = 1;
        b2 = step;
        index=1;
        while(1)
            xT=Projection(b1:b2);
            x1 = find(Projection>=xT(1) & Projection<xT(end));
            x2=find(Proj_Total>=xT(1) & Proj_Total<xT(end));
            
            if(isempty(x1))
               x1 = find(Projection>=xT(1) & Projection<=xT(end)); 
               x2=find(Proj_Total>=xT(1) & Proj_Total<=xT(end));
               if(index==1)
                    x2=find(Proj_Total<=xT(end));
               end
            else
                if(index==1)
                    x2=find(Proj_Total<xT(end));
               end
            end
            b1=x1(end)+1;
            b2=b1+step;
            x1=Projection(x1);
            x2=Proj_Total(x2);
            
            pS = length(x1);
            pT = length(x2);
            xx(index)=x1(end);
            if(strcmp(validation,'v'))
                Pspike(index) = pS/(pT+pS);
            else
                Pspike(index) = pS/pT;
            end
            
            if( ( b1>length(Projection) ) || (b2>=length(Projection)) )
                
                if(b1<=length(Projection))
                    xT=Projection(b1:end);
                    x1 = find(Projection>=xT(1));
                    x1=Projection(x1);
                    x2=find(Proj_Total>=xT(1));
                    pS = length(x1);
                    pT = length(x2);
                    xx(end+1)=x1(end);
                    if(strcmp(validation,'v'))
                        Pspike(end+1) = pS/(pT+pS);
                    else
                        Pspike(end+1) = pS/pT;
                    end
                    
                end
                
                break
            end
            
            index=index+1;
        end
        
    
    end

    %%---------------------------------------------------------------------
    function figureCloseCallback(gcbo, eventdata);
        delete(hFigure);
        
    end

  %%---------------------------------------------------   
   function SaveData(Data2Save, type)
       val=get(WNstd, 'String');
       
       Date = date;
       Day = sprintf('%02.0f', day(Date));
       [mNum, mstr] = month(date)
       Month = mstr;
       Year = sprintf('%02.0f', year(Date));
       
       Todaysdate = [Year Month Day];
       Folder = [WorkFolder '\' Todaysdate];
       
       if(exist(Folder, 'file') == 7)
       else
           mkdir(Folder);
       end
       
       if (type==1)
           saveString = [Todaysdate '_Response'];
           index = findFiles(Folder, saveString);
           index = sprintf('%03.0f', index);
           saveString = [Folder '\' saveString '_std' val '_' index '.txt'];
       elseif(type==2)
           saveString = [Todaysdate '_Threshold'];
           index = findFiles(Folder, saveString);
           index = sprintf('%03.0f', index);
           saveString = [Folder '\' saveString '_std' val '_' index '.txt'];
       
       end
       
       save(saveString, 'Data2Save', '-ASCII');
       
       
     
   end


  %%---------------------------------------------------   
   function [index indices] = findFiles(folder, strStem)
       if nargout < 1
            indices=[];
       else
           indices=[];
       end
        
       a1 = folder;
       files=dir(a1);
       l = length(files);
       nfiles = length(files);
       index = 0;
       
       for i=3:nfiles
           currentfilename = files(i).name;
           str = strStem;
           %indices(i-2)=str2double(currentfilename(end-4));
           if(strfind(currentfilename, str))
               indices(index+1)=str2double(currentfilename(end-4));
               index = index +1;
               
           end
       end
       
       indices=sort(indices);
       indices= unique(indices);
       
   end

end


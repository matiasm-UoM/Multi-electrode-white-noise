function WhiteNoiseAnalysis( )
%GENERATERESPONSE Summary of this function goes here
%   Detailed explanation goes here

clc;
% clear all;
% close all;


%GUI variables
hFigure=[];
cmap =[];
haxes1=[];
haxes2=[];
DataFolder=[];     
dataNonEvoked=[];
LoadData=[];
dataEvoked=[];
WNstd=[];
folders=[];
nbinsLabel=[];
TotalStim=[];
cellLoc=[];
d1=[];
d2=[];
SpikeLogicVect=[];
StimPulseLength = 1.05e-3; %1.05 ms

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
         C = findFiles(file, 'AllStim_M');
         
         if(C>1)
%              [i indices] = findFiles(file, 'Evoked_M');
             [i indices] = findFiles(file, 'AllStim_M');

             inputText = ['Cell to load (' mat2str(indices) ' cells available):'];
             C2Load = input(inputText);
         else
%              [i indices] = findFiles(file, 'Evoked_M');
             [i indices] = findFiles(file, 'AllStim_M');
             C2Load = indices;
         end


         file1 = [file '\AllStim_M' num2str(C2Load) '.txt'];
             
         
         cmap = load('cmap.mat');
         cmap = cmap.cmap;
         
         TotalStim=[];
         [h1, TotalStim] = hdrload(file1);

         cellLoc= h1(2,:)
         cellLoc = cellLoc(9:end);
         cellLoc = [str2double(cellLoc(1:3)) str2double(cellLoc(5:7)) ];
         if(isnan(cellLoc))
             cellLoc=[];
         end
         
         

         %Remove any shorting stimuli
         [rw, cl] = find(TotalStim(:,:) == -10000);
         TotalStim(rw, :)=[];
         assignin('base', 'TotalStim', TotalStim);
         TotalStim(all(TotalStim == 0,2),:) = [];   % remove any rows of zeros
         SpikeTimes = TotalStim(:, 21:end);
         TotalStim(:,21:end)=[];
         [SpikeLogicVect]=FindDelayEvokedStimuli(TotalStim, SpikeTimes);
         disp(['Number of spikes: ' num2str(length(find(SpikeLogicVect)))]);
         
         
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
    function [SpikeLogicVect]=FindDelayEvokedStimuli(DEvoked, SpikeTimes)
        
            SpikeLogicVect = zeros(1, length(SpikeTimes));
            [A] = find(SpikeTimes(:,1)>delay1 & SpikeTimes(:,1)<=delay2);
            SpikeLogicVect(A) = 1;
            SpikeLogicVect=logical(SpikeLogicVect);
            assignin('base','SpikeLogicVect',SpikeLogicVect);
    end

 
         %%---------------------------------------------------------------------
    function ComputeParameters(gcbo, eventdata);
        clc;
        axes(haxes1);
        cla;
        axes(haxes2);
        cla;
        
        tstim  = TotalStim;
        tstim = tstim(:,any(tstim));    %remove any columns of zeros 
        TotalStim = tstim;
        LengthEvokedStim = length(find(SpikeLogicVect))
        LengthTotalStim = length(SpikeLogicVect)
        
        NumberOfPerm = 1000;        
        HypothesisNotMet=1;
        ExcitatorySignificantComponents = 0;
        InhibitorySignificantComponents = 0;

        
        Evoked = tstim(SpikeLogicVect,:);
        EvokedCov = cov(Evoked) - cov(tstim);
        [Eig, EigValST_orig, explained] = pcacov(EvokedCov);
        CovTotal = cov(tstim);
        StimVariance = mean(diag(CovTotal));
        disp(['Stimulus variance: ' num2str(StimVariance)]);
        
        test = TotalStim*EigValST_orig(:,1);
        figure;
        hist(test,31)
            
        while(HypothesisNotMet)
            
%             assignin('base', 'tstim', tstim);
            Evoked = tstim(SpikeLogicVect,:);
            EvokedCov = cov(Evoked)- cov(tstim);
            [Eig, EigValST, explained] = pcacov(EvokedCov); 
            [aCov bCov] = size(EvokedCov);
        
            VariancePerm=[];
            circShift = randperm(LengthTotalStim,NumberOfPerm);
            for i=1:NumberOfPerm
                randPerm = circshift(SpikeLogicVect,circShift(i),2);
                Evoked = tstim(randPerm,:);
                EvokedCov = cov(Evoked)- cov(tstim);
                [EigVect, EigVal, EXPLAINED] = pcacov(EvokedCov);

                VariancePerm(end+1,:) = EigVal;
            end

            MeanEigenval = mean(VariancePerm);
            stdEigenval = std(VariancePerm);
            assignin('base', 'VariancePerm',VariancePerm)
            
%             axes(haxes2);
            figure;
            bounds = 2; 
            xx = 1:length(MeanEigenval);
            X = [xx, fliplr(xx)];
            Y = [MeanEigenval+bounds*stdEigenval fliplr(MeanEigenval-bounds*stdEigenval)];
            fill(X,Y,'r');hold on;
            plot(xx, MeanEigenval,'k');
            plot(xx, MeanEigenval+bounds*stdEigenval,'--k',xx, MeanEigenval- bounds*stdEigenval, '--k')
            scatter(xx, EigValST, 'k', 'filled' );
            xlabel('Component');
            ylabel('Variance');
            
            
            
            first = 0;
            last = 0;
            numberOfCompFound = ExcitatorySignificantComponents+InhibitorySignificantComponents;
            excitatory = abs(EigValST(1)-MeanEigenval(1));
            inhibitory = abs(EigValST(end-numberOfCompFound)-MeanEigenval(end-numberOfCompFound));
            if(excitatory>inhibitory)
                first = 1;
                excitatory=excitatory;
                componentNumber = ExcitatorySignificantComponents + 1
            else
                last = 1;
                inhibitory=inhibitory;
                componentNumber = length(EigValST) - InhibitorySignificantComponents
            end
            lim1 = MeanEigenval(1)+bounds*stdEigenval(1);
            lim2 = MeanEigenval(1)-bounds*stdEigenval(1);
            lim3 = MeanEigenval(end-numberOfCompFound)+bounds*stdEigenval(end-numberOfCompFound);
            lim4 = MeanEigenval(end-numberOfCompFound)-bounds*stdEigenval(end-numberOfCompFound);
            
            figure(21)
            if(numberOfCompFound==0)
                scatter(1:length(EigValST_orig), EigValST_orig, 'k', 'filled' );hold all;
                distance = abs(EigValST_orig-MeanEigenval');
                assignin('base','EigValST_orig',EigValST_orig);
                assignin('base','MeanEigenval',MeanEigenval);
                distance=sort(distance,'descend');
                disp(['first largest: ' num2str(distance(1))]);
                disp(['second largest: ' num2str(distance(2))]);
                
            end
            scatter(xx+ExcitatorySignificantComponents, EigValST);
            
            drawnow;
            
            if(first)
                disp('largest component is excitatory');
            end
            if(last)
                disp('largest component is inhibitory');
            end
            
            [a b] = size(tstim);
            [a b] = size(Eig);
            if(first)
                if(EigValST(1)>lim1 || EigValST(1) < lim2)
                
                ExcitatorySignificantComponents = ExcitatorySignificantComponents+1;
                HypothesisNotMet=1;
                %remove current principal component from tstim
                tstim = tstim*Eig(:,2:end);
                tstim = tstim*Eig(:,2:end)';
                else
                    HypothesisNotMet=0;
                end
            elseif(last)
                if(EigValST(end-numberOfCompFound)>lim3 || EigValST(end-numberOfCompFound) < lim4)
                
                InhibitorySignificantComponents = InhibitorySignificantComponents+1;
                HypothesisNotMet=1;
                %remove current principal component from tstim
                tstim = tstim*Eig(:,1:end-(numberOfCompFound+1));
                tstim = tstim*Eig(:,1:end-(numberOfCompFound+1))';
                else
                    HypothesisNotMet=0;
                end
            end
            
            [a b] = size(tstim)
            [a b] = size(Eig)
%             if(numberOfCompFound>0)
%                  HypothesisNotMet=0;
%             end
        end
       
     disp(['Number of excitatory significant components = ' num2str(ExcitatorySignificantComponents)]);
     disp(['Number of inhibitory significant components = ' num2str(InhibitorySignificantComponents)]);
        
        
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


%% Matlab Spike Sorter 
%This is an open-code free Matlab tool which allows user to make offline
%manual and semi-automatic sorting of spikes. This pretends to be a start
%point to the development of a full functionally Matlab tool for spike sorting. 
%Developed by Joaquin Mansilla Yulan and Dr. Ing. Sergio Lew
%Facultad de Ingeniería - Universidad de Buenos Aires
%Githu


%% INITIALIZATION FUNCTIONS

function varargout = SpikeSort(varargin)
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
        'gui_Singleton',  gui_Singleton, ...
        'gui_OpeningFcn', @SpikeSort_OpeningFcn, ...
        'gui_OutputFcn',  @SpikeSort_OutputFcn, ...
        'gui_LayoutFcn',  [] , ...
        'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
        % End initialization code - DO NOT EDIT
    end

% --- Executes just before SpikeSort is made visible.
%Declare all global variables
function SpikeSort_OpeningFcn(hObject, eventdata, handles, varargin)

    % Define Global Variables
    %Parameters:
    global Channels; % Channels is the number of channels
    global Fs %Sampling Frecuency 
    global DataBlock % Longest Raw Data package length to process at once
    global Scale
    global SL   %Sample Length (Bits)
    global norm %Normalization 1:Yes 0:NO
    global GcaColor %Color of figures background
    global WfsColor %Color of Wfs signals
    global SelColor %Color of selected signals
    global RawColor %Color of raw data signals
    global BlockTime %Waiting time for plotting figures
    global RecordPathName %Path name of data to be recorded
    global Ls % ms length of each waveform
    global ps % Point Size
    global DensityRes %Grid resolution of density plot (3D histogram)
    global ISI %InterSpike Interval

    %Raw Data variables:
    global PathName %Path name of the loaded data
    global RwDPlot %Raw Data Spikes plot object
    global ThPlot %Threshold plot object
    global fileType %Type of loaded file: 1- Binary File, 2- 1 Channel.mat File, 3- Saved File
    global matdata %Raw data from 1 channel -mat data
    global data % Raw Data 
    global MaxDataLen; % Raw Data package length to process at once
    global threshold % Amplitude threshold to detect Spikes
    global dataView % Raw Data displayed
    global slider_step; % step of slidebar (Raw Data plot)
    global normScale % Multiplying factor to normalize data


    %Waveform variables:
    global WFmatrix % WaveForm Matrix
    global SpikeTime % Time Position of each Spike (Waveform)
    global SelectedUnitsIDX % Selection vector
    global showRange; % samples to show
    global MouseSelection  % ploted points of pointer selection
    global MouseSelectionI %Index of pointer selection

    %Features Variables:
    global FeaturePoints % 2D or 3D Feature points
    global FeaTime  % One Feature (to plot over time)
    global slice %PLot object of Slice 1 and 2

    %Cluster Variables:
    global UnitStatus   % Units - Clusters vector
    global Unitn %Number of Units (clusters)
    

    %Program Variables:
    global ListC % List of Clustering Algorithms
    global ListF % List of Features Algorithms
    global handles2Spk %Handles to SpikeSorter GUI
    global handles2Par %Handles to Parameters GUI

    %Set default parameters
    Set_default_parameters(handles);
    handles2Spk=handles;
    handles.output = hObject;
    guidata(hObject, handles);
    set(handles.WFaxes,'xtick',[])
    set(handles.WFaxes,'ytick',[])
    set(handles.WFaxes,'Color',GcaColor);
    set(handles.Featureaxes,'xtick',[]);
    set(handles.Featureaxes,'ytick',[]);
    set(handles.Featureaxes,'Color',GcaColor);
    set(handles.FeatureTimeaxes,'xtick',[])
    set(handles.FeatureTimeaxes,'ytick',[])
    set(handles.FeatureTimeaxes,'Color',GcaColor);
    set(handles.Signalaxes,'Color',GcaColor);
    ListC=search_fnc_from_file('clustering_algorithm.m');
    ListF=search_fnc_from_file('Features.m');
    set(handles.Clus_Alg,'String',ListC);
    set(handles.Feat1,'String',ListF);
    set(handles.Feat2,'String',ListF);
    set(handles.Feat3,'String',strvcat('---',ListF));
    set(handles.FeatTime,'String',ListF);
    set(handles.Feat2,'Value',2);
%     axes(handles.logoUBA)
%     imshow('logo.png')
%     axis off
%     axis image

% --- Outputs from this function are returned to the command line.
function varargout = SpikeSort_OutputFcn(hObject, eventdata, handles)
    varargout{1} = handles.output;









%% FILE, MEMORY AND GENERAL FUNCTIONS


%Open Parameters GUI
function Parameter_Configuration_Callback(hObject, eventdata, handles)
    global handles2Par
    handles2Par=guidata(Parameters);

    
% LOAD FILE

%Loads a Binary File Raw Data.
function File_Callback(hObject, eventdata, handles)
    [auxFileName,auxPathName,FilterIndex] = uigetfile('*.*'); % Get the file location and name of the selected file
    if auxFileName~=0

        % Clear all variables and plots
        Reset_Variables(handles);
        Reset_plots(handles);

        global SL
        global DataBlock
        global fileType
        global data 
        global FileName  
        global PathName
        global MaxDataLen; 
        global slider_step; 
        global Channels

        % File Info
        fileType=1;
        FileName=auxFileName;
        PathName=auxPathName;
        fileInfo = dir([PathName FileName]);
        fileSize = fileInfo.bytes;
        Set_RecordPath;

        % Load File
        MaxDataLen=min(DataBlock,fileSize/(Channels*2*5));  % The longest Raw Data channel package to process at once will be DataBlock *SL bits
        fid=fopen([PathName FileName],'rb'); %opens a file in binary read mode.
        SL_str=num2str(SL);
        eval(['data=fread(fid,[Channels MaxDataLen],''int' SL_str ''');']); %Reads binary data from the opened file
        
        % Set Slider Step
        slider_step=min(1,1/fix((fileSize/(SL/8*MaxDataLen*Channels)-1))); %Set the slider step of the slide bar (raw data plot)
        set(handles.slider1,'SliderStep',[slider_step/10 slider_step]);
        
        %Scale or Normalize
        Get_normfactor(handles);
        Data_scaling(handles)

        %Quantization
        data=eval(['int' SL_str '(data);']);
        
        % Filter
        Channel_and_Filter(handles);

        % Plot Raw Data
        plot_channel(handles);
    end

%Load 1 Channel
%Loads One Channel Raw Data contained in a .mat variable file.
% --------------------------------------------------------------------
function Channel_Raw_Callback(hObject, eventdata, handles)
    [auxFileName,auxPathName,FilterIndex] = uigetfile('*.*'); % Get the file location and name of the selected file
    if auxFileName~=0

        % Clear all variables and plots
        Reset_Variables(handles);
        Reset_plots(handles);

        global SL
        global DataBlock
        global fileType
        global data
        global matdata
        global FileName  
        global PathName
        global MaxDataLen; 
        global slider_step; 
        global Channels
        global norm
        global Scale


        % Save File Info
        fileType=2;
        FileName=auxFileName;
        PathName=auxPathName;
        Set_RecordPath;

        % Load File
        set(handles.Channels,'String','1');
        Channels=1;
        aux=load([PathName FileName]); %Reads raw data from the opened .mat file
        aux2=fields(aux);
        matdata=aux.(cell2mat(aux2));
        SL_str=num2str(SL);
        %Scale or Normalize
        if norm==1
            matdata=2^(SL-1)*(matdata/(max(max(matdata),abs(min(matdata)))));
        else
            matdata=matdata/Scale;
        end
        %Quantization
        matdata=eval(['int' SL_str '(matdata)']);
        MaxDataLen=min(DataBlock,length(matdata)/5);  % The longest Raw Data channel package to process at once will be DataBlock *SL bits
        data=matdata(1:MaxDataLen);

        % Filter
        Channel_and_Filter(handles);

        % Set Slider Step
        slider_step=min(1,1/fix(length(matdata)/MaxDataLen-1)); %Set the slider step of the slide bar (raw data plot)
        set(handles.slider1,'SliderStep',[slider_step/10 slider_step]);

        % Plot Raw Data
        plot_channel(handles);
    end

% Load Recorded Units
%Loads all recorded results from a recorded sort.
function Recorded_Units_Callback(hObject, eventdata, handles)
    [auxFileName,auxPathName,FilterIndex] = uigetfile('*.*'); % Get the file location and name of the selected file
    if auxFileName~=0 % if loads any file 

        %Clear all variables and plots
        Reset_Variables(handles);
        Reset_plots(handles);

        % Load Recordes DAta
        load([auxPathName auxFileName]); % Warning: Resets global Workspace also
        if (max(SpikeTime)>40000) %from old Spike Sorter
            SpikeTime=SpikeTime/(Fs*60);
        end
        get_Features(handles);
        plot_Wfs(handles);

        %Set other global variables
        global Unitn
        global FileName  
        global PathName
        global fileType

        Unitn=max(UnitStatus);

        % Save File Info
        FileName=auxFileName;
        PathName=auxPathName;
        fileType=3;
        Set_RecordPath;
    end

% Save Spikes and Spike Time
% Saves Waveforms data, time position of each waveform and clusters (Units)
function Save_Callback(hObject, eventdata, handles)
    global FileName
    global PathName
    global WFmatrix
    global UnitStatus
    global SpikeTime
    global fileType
    global BlockTime
    global RecordPathName

    % Block objects and start error timer
    time=start_timer(handles,BlockTime);
    s = warning('off');

    %Save Data
    if fileType==3  %if file is recorded data
        fold=[RecordPathName FileName];
    else
        mkdir(RecordPathName);
        Channel = get(handles.Channel,'Value');
        fold=[RecordPathName FileName '_C' num2str(Channel) '.mat']; %Set output adress
    end
        if exist(fold, 'file')==2  %Overwrite check
            ovw=questdlg('Saved data already exists. Do you want to overwrite?' ,'Overwrite data','Yes','No','No');
            if strcmp(ovw,'Yes')
                save(fold,'WFmatrix','SpikeTime','UnitStatus');
                msgbox(['Data saved: ' fold]);
            end
        else
            save(fold,'WFmatrix','SpikeTime','UnitStatus');
            msgbox(['Data saved: ' fold]);
        end
        
    % Unblock objects and end error timer
    warning(s)
    end_timer(time,1,handles);

%Set Record Directory Path
function Set_RecordPath
    global RecordPathName
    global fileType
    global PathName
    global FileName
    if fileType==3  %if file is recorded data
            RecordPathName=PathName;
        else
            RecordPathName=[PathName 'Saves_' FileName  '\']; 
    end

%Save mean Spikes  
%Save the mean signal (centroid) of each cluster.
function Spike_mean_Callback(hObject, eventdata, handles)
    global UnitStatus
    global WFmatrix
    
    auxcd=cd;
    path=which('SpikeSort.m');
    path=path(1:(end-11));
    cd(path);
    path=[path 'Spikes.mat']; %Set output adress
    for i=(1:max(UnitStatus));
        Spikes(i,:)=mean(WFmatrix((find(UnitStatus==i)),:));
        figure(i)
        plot(Spikes(i,:))
    end
    if (exist('Spikes.mat','file')==2)
            k=[];
            aux=Spikes;
            load('Spikes.mat');
            for i=1:(length(aux(:,1)))
                for j=1:(length(Spikes(:,1)))
                    if (isequal(aux(i,:),Spikes(j,:)))
                        k=[k;i];
                    end
                end
            end
            aux(k,:)=[];
            Spikes=[Spikes;aux];
            save(path,'Spikes')
    else
        save(path,'Spikes')
    end
    cd(auxcd);
    
%Get normazlization factor of all raw data
function Get_normfactor(handles)
    global norm
    global normScale
    global slider_step
    global data
    global SL
    
    if norm==1
        max_abs=0;
        ndataBlocks=1/slider_step;
        for i=0:ndataBlocks
            slider=i*slider_step;
            %Get position 
            seek_on_data(slider)
            aux=max(max(max(data)),abs(min(min(data))));
            if aux>max_abs
                max_abs=aux;
            end
        end
        normScale=2^(SL-1)/max_abs;
    end
    
%Raw Data Scaling
%Scale or Normalize Raw Data on Data according to norm parameter value
function Data_scaling(handles)
    global data
    global norm
    global Scale
    global normScale
    if norm==1
        data=data*normScale;
    else
        data=data/Scale; 
    end
    
%Reset Variables
%Reset all variables to initial values
function Reset_Variables(handles)
    global RwDPlot 
    global ThPlot 
    global matdata 
    global data
    global threshold 
    global dataView 
    global MaxDataLen
    global slider_step 
    global WFmatrix 
    global FeaturePoints 
    global FeaTime
    global SelectedUnitsIDX 
    global UnitStatus 
    global SpikeTime 
    global showRange 
    global Unitn 
    global slice 
    global MouseSelection 
    global MouseSelectionI 
    global ListC
    global ListF
    global handles2Spk
    global fileType
    global FileName
    global PathName

    fileType=0;
    FileName=[];
    PathName=[];
    RwDPlot=[]; 
    ThPlot=[]; 
    matdata=[];
    handles2Spk=handles;
    ListC=search_fnc_from_file('clustering_algorithm.m');
    ListF=search_fnc_from_file('Features.m');
    MaxDataLen=0;
    data=[];
    threshold=0;
    dataView=[];
    FeaturePoints=[];
    FeaTime=[];
    SelectedUnitsIDX=[];
    UnitStatus=[];
    SpikeTime=[];
    WFmatrix=[];
    showRange=1;
    Unitn=0;
    slice=zeros(1,2);
    set(handles.Slice_1,'String','0');
    set(handles.Slice_2,'String','0');
    MouseSelection=0;
    MouseSelectionI=0;
    slider_step=1;
    set(handles.Channel,'Value',1);
    set(handles.Clusters,'Value',1);
    set(handles.Fil,'Value',1);
    set(handles.Th,'String',0);
    set(handles.slider1,'Max',1);
    set(handles.slider1,'Min',0);

%Set parameters to default values
function Set_default_parameters(handles)
    global Channels
    global Fs 
    global DataBlock 
    global Scale
    global SL   
    global norm 
    global GcaColor 
    global WfsColor
    global SelColor
    global RawColor
    global BlockTime 
    global RecordPathName 
    global Ls 
    global ps 
    global DensityRes 
    global ISI

    Channels=25;
    Fs=30030;
    DataBlock=801000;
    Scale=2.048;
    SL=16;
    norm=0;
    GcaColor=[0 0.15 0.35];
    WfsColor=[0.35 0.6 0.55];
    SelColor=[1 1 1];
    RawColor=[0.8 0.3 0.3];
    BlockTime=7;
    RecordPathName=[];
    Ls=4;
    ps=25;
    DensityRes=[100 100];
    ISI=2;

% Free Memory
function FreeMem_Callback(hObject, eventdata, handles)
    clear
    clear global

% Fix 
% Unblock all objects 
function Fix_Callback(hObject, eventdata, handles)
    unblock_objects(handles);

% Logo_Button
function Logo_Button_Callback(hObject, eventdata, handles)
    cdata=get(hObject,'cdata');
    msgbox(['Matlab Spike Sorter was developed by:' sprintf('\n') ...
    'Joaquin Mansilla Yulan' sprintf('\n') 'Dr. Ing. Sergio Lew' ...
    sprintf('\n') 'Instituto de Ingeniería Biomédica - Universidad' ...
    ' de Buenos Aires'],'Spike Sorter','custom',cdata);


    
%TIMER FUNCTIONS

%Initiate block timer
%Block all objects while timer on 
%This function is usually used when SpikeSorter is plotting
function [time]=start_timer(handles,sec)
    time=timer('TimerFcn',{@end_timer,handles},'StartDelay',sec);
    block_objects(handles);
    start(time);

%Ends block timer
%Unblock all objects and stop timer
%This function is usually used after plotting
function end_timer(obj,event,handles)
    s = warning('off');
    delete(obj);
    warning(s);
    unblock_objects(handles);

%Unblock
%Unblocks buttons after plotting
function unblock_objects(handles)
    aux=get(handles.SelectWFs,'Enable');
    set(handles.SpikeSort, 'pointer', 'arrow') 
    if strcmp(aux,'off');
        set(handles.SelectWFs,'Enable','on');
        set(handles.SelWfsPol,'Enable','on');
        set(handles.PlotColors,'Enable','on');
        set(handles.showReduced,'Enable','on');
        set(handles.Sort,'Enable','on');
        set(handles.Clusters,'Enable','on');
        set(handles.Right_1,'Enable','on');
        set(handles.Slice_1,'Enable','on');
        set(handles.Left_1,'Enable','on');
        set(handles.Right_2,'Enable','on');
        set(handles.Slice_2,'Enable','on');
        set(handles.Left_2,'Enable','on');
        set(handles.AlignMax,'Enable','on');
        set(handles.AlignMin,'Enable','on');
        set(handles.SelectData,'Enable','on');
        set(handles.SelFtPol,'Enable','on');
        set(handles.Add_Unit,'Enable','on');
        set(handles.ClearUnits,'Enable','on');
        set(handles.PlotUnitsSel,'Enable','on');
        set(handles.DeselectData,'Enable','on');
        set(handles.DeleteData,'Enable','on');
        set(handles.PlotFeatures,'Enable','on');
        set(handles.ShowSpikes,'Enable','on');
        set(handles.FeatTime,'Enable','on');
        set(handles.SelRectFeaTime,'Enable','on');
        set(handles.SelFeaTimePol,'Enable','on');
    end

%Block
%Blocks buttons while plotting
function block_objects(handles)
    aux=get(handles.SelectWFs,'Enable');
    set(handles.SpikeSort, 'pointer', 'watch')
    if strcmp(aux,'on');
        set(handles.SelectWFs,'Enable','off');
        set(handles.SelWfsPol,'Enable','off');
        set(handles.PlotColors,'Enable','off');
        set(handles.showReduced,'Enable','off');
        set(handles.Sort,'Enable','off');
        set(handles.Clusters,'Enable','off');
        set(handles.Right_1,'Enable','off');
        set(handles.Slice_1,'Enable','off');
        set(handles.Left_1,'Enable','off');
        set(handles.Right_2,'Enable','off');
        set(handles.Slice_2,'Enable','off');
        set(handles.Left_2,'Enable','off');
        set(handles.AlignMax,'Enable','off');
        set(handles.AlignMin,'Enable','off');
        set(handles.SelectData,'Enable','off');
        set(handles.SelFtPol,'Enable','off');
        set(handles.Add_Unit,'Enable','off');
        set(handles.ClearUnits,'Enable','off');
        set(handles.PlotUnitsSel,'Enable','off');
        set(handles.DeselectData,'Enable','off');
        set(handles.DeleteData,'Enable','off');
        set(handles.PlotFeatures,'Enable','off');
        set(handles.ShowSpikes,'Enable','off');
        set(handles.FeatTime,'Enable','off');
        set(handles.SelRectFeaTime,'Enable','off');
        set(handles.SelFeaTimePol,'Enable','off');
    end

%Search name functions
%Create a list of function names in a script delimited by %##%
%Used to create a list of feature alghoritms available
%Used to create a list of clustering alghoritms available
function [List]=search_fnc_from_file(file)
    List=[];
    f=fopen(file,'r');
    f=fread(f,'char');
    f=char(f');
    i=strfind(f,'HERE');
    Istart=strfind(f(i:end),'%##')+i+2;
    Ifinal=strfind(f(i:end),'##%')+i-2;
    for j=1:length(Ifinal)
        if j==1
            List(j,:)=f(Istart(j):Ifinal(j));
        else
            auxlist=f(Istart(j):Ifinal(j));
            List=strvcat(List,auxlist);
        end      
    end












%% RAW DATA FUNCTIONS


%Change channel
%Callback to channel selection. Plots the selected channel after filtering.
function Channel_Callback(hObject, eventdata, handles)
    global SelColor
    global SelectedUnitsIDX
    
    % Reset threshold to 0
    Reset_threshold(handles);
    %Change channel and Filter
    Channel_and_Filter(handles);
    %Plot raw data
    plot_channel(handles);
    %Plot threshold
    plot_Th(handles);
    %Plot selected Spikes on raw data
    plot_sel_rawdata(handles,SelectedUnitsIDX,SelColor);

%Reset threshold to 0 value
function Reset_threshold(handles)
    global threshold
    threshold=0;
    set(handles.Th,'String',num2str(round(threshold)));

%Change Channels parameter from SpikeSort
function Channels_Callback(hObject, eventdata, handles)
    global Channels
    Channels=str2double(get(handles.Channels,'String'));

% Slider movement.
%Callback to slidebar movement.  
function slider1_Callback(hObject, eventdata, handles)
    global SelColor
    global SelectedUnitsIDX 

    %Slide into Raw Data
    slider=get(hObject,'Value');
    move_raw_data(handles,slider);
    %Filter
    Channel_and_Filter(handles);
    %Plot raw data
    plot_channel(handles);
    %Plot threshold
    plot_Th(handles);
    %Plot selected Spikes on raw data
    plot_sel_rawdata(handles,SelectedUnitsIDX,SelColor);
    
%Filter
% Changes the channel and applies a High-pass filter to actual channel signal.
%Also, creates dataView from data and channel & filter parameters
function Channel_and_Filter(handles)
    global data
    global dataView
    global Fs
    Channel = get(handles.Channel,'Value'); 
    aux=get(handles.Fil,'String');
    f=str2num(cell2mat(aux(get(handles.Fil,'Value'))));
    %Apply High-Pass filter at f cutoff frequency
    if f>0
        h= firls(300,[0 0.8*f/Fs 1.2*f/Fs 1],[0 0 1 1]);
        dataView=conv(double(data(Channel,:)),h,'same');
        dataView=int16(dataView);
    else
        dataView=data(Channel,:);
    end

%Filter current signal
function Fil_Callback(hObject, eventdata, handles)
    %Filter
    Channel_and_Filter(handles);
    %Plot raw data
    plot_channel(handles);
    %Plot threshold
    plot_Th(handles);

%Move threshold up
function Up_Callback(hObject, eventdata, handles)
    global threshold
    global dataView
    %Set new Threshold
    threshold = threshold + max(dataView)/40;
    set(handles.Th,'String',num2str(round(threshold)));
    %Plot Threshold
    plot_Th(handles);

%Move threshold down
function Down_Callback(~, eventdata, handles)
    global threshold
    global dataView
    %Set new Threshold
    threshold = threshold - max(dataView)/40;
    set(handles.Th,'String',num2str(round(threshold)));
    %Plot Threshold
    plot_Th(handles);

%Set threshold
function Th_Callback(hObject, eventdata, handles)
    %Set new Threshold
    global threshold
    threshold=str2num(get(handles.Th,'String'));
    %Plot Threshold
    plot_Th(handles);

%Plot Raw Data
function plot_channel(handles)
    global dataView
    global slider_step
    global Fs
    global GcaColor
    global RawColor
    global norm

    if isempty(dataView)==0
        %Time axis
        W=10000;
        p=(get(handles.slider1,'Value')/slider_step)*length(dataView)/Fs; % get inital position of time axis
        t=(1/Fs+p):1/Fs:(length(dataView)/Fs+p); %Time axis
        %Plot data
        axes(handles.Signalaxes);
        hold off
        plot(t,dataView,'Color',RawColor)
        %Set graphic features
        set(gca,'Color',GcaColor);
        xlim([p (length(dataView)/Fs+p)])
        ylim([mean(double(dataView(1:W)))-15*std(double(dataView(1:W))),mean(double(dataView(1:W)))+15*std(double(dataView(1:W)))])
        xlabel('Time (S)','FontSize', 9);
        if norm==0
            ylabel('Voltage (uV)','FontSize', 9);
        else
            ylabel('Amplitude','FontSize', 9);
        end
    end

%Plot Threshold
function plot_Th(handles)
    global ThPlot
    global threshold
    %Deletes old threshold plot and plot new one
    axes(handles.Signalaxes);
    hold on
    delete(ThPlot);
    lim=axis;
    x=[lim(1) lim(2)];
    ThPlot=plot(x,[threshold threshold],'Color',[0 0.8 0]);

%Plot selected Spikes on RawData
function plot_sel_rawdata(handles,IDX,SColor)
    global SpikeTime
    global MaxDataLen
    global dataView
    global slider_step
    global Ls
    global Fs
    global RwDPlot
    
   if not(isempty(IDX) & isempty(SpikeTime)) % Plots selected waveforms
       L=round((Ls*Fs/(2*1000)));
       delete(RwDPlot);
       axes(handles.Signalaxes);
       hold on
       pack=get(handles.slider1,'Value')/slider_step;
       pii=round(pack*MaxDataLen);
       pf=round(MaxDataLen+pii);
       auxSpikeTime=SpikeTime(IDX);
       auxSpikeTime=round(auxSpikeTime(auxSpikeTime>(pii/(Fs*60)) & auxSpikeTime<(pf/(Fs*60)))*60*Fs);
       t=[];
       %spike=[];
       dataux=[];
       for i=1:length(auxSpikeTime)
           %spike=[spike auxSpikeTime(i)-pii-L:auxSpikeTime(i)-pii+(L-1)];
           %t=(spike+pii)/F;
           dataux=[dataux 0 dataView(auxSpikeTime(i)-pii-L:auxSpikeTime(i)-pii+(L-1)) 0];
           taux=(auxSpikeTime(i)-L:auxSpikeTime(i)+(L-1))/Fs;
           t=[t (auxSpikeTime(i)-(L-1))/Fs taux (auxSpikeTime(i)+L)/Fs];
       end
       RwDPlot=plot(t,dataux,'Color',SColor);
   end
    
%Asign to Data a new portion of Raw Data according to slider
function move_raw_data(handles,slider)
    global SL
    global matdata
    global fileType
    global slider_step;
    global FileName
    global PathName
    global MaxDataLen;
    global data
    global Channels
    
    switch fileType 
        case 1 %Binary File
                %Moves data to slider position
                seek_on_data(slider)
                %Scale or Normalize
                Data_scaling(handles)
                %Quantization
                SL_str=num2str(SL);
                data=eval(['int' SL_str '(data);']);
                
        case 2 % One channel .mat file
            ii=round((slider/slider_step)*MaxDataLen);
            data=matdata(ii+1:ii+MaxDataLen);
    end

%Asign to Data a new portion of binary Raw Data file according to slider (double type)
%Subfunction of move_raw_data
function seek_on_data(slider)
    global SL
    global slider_step;
    global FileName
    global PathName
    global MaxDataLen;
    global data
    global Channels
    
    seek=round(SL/8*Channels*MaxDataLen*(slider/slider_step));
    rem=mod(seek,(Channels*SL/8));
    seek=seek-rem;
    fid=fopen([PathName FileName],'rb');
    if (fseek(fid,seek,-1)==0)
        %Load File
        SL_str=num2str(SL);
        eval(['data=fread(fid,[Channels MaxDataLen],''int' SL_str ''');']); %Reads binary data from the opened file
    end
    fclose(fid);




    
    


















%% WAVEFORMS FUNCTIONS


%WFs local detection
function PrevWf_Callback(hObject, eventdata, handles)
    global dataView
    global WFmatrix
    global SpikeTime
    global Ls
    global UnitStatus
    global Fs

    %Detects Spikes under or over threshold on Dataview
    L=round((Ls*Fs/(2*1000)));
    WFmatrix=[];
    SpikeTime=[];
    [idxpos]=Spike_position;

    %Create Spike Waveforms
    for n=1:length(idxpos)
        WFmatrix(n,:)=dataView(idxpos(n)-L:idxpos(n)+L-1);
        UnitStatus(n,:)=0;
        SpikeTime(n,:)=idxpos(n);
    end
    SpikeTime=SpikeTime/(Fs*60); %Time ocurrence of Spikes in minutes
    %Plot Waveforms
    plot_Wfs(handles);
    set(handles.WFsN,'String',[num2str(size(WFmatrix,1)) ' WFs']);
    %Get Feature points
    get_Features(handles);

%WFdetection
%Detects Waveforms (Spikes) over or under the threshold. 
function WFdetection_Callback(hObject, eventdata, handles)
    global dataView
    global MaxDataLen;
    global threshold
    global WFmatrix
    global SpikeTime
    global Ls
    global UnitStatus
    global slider_step
    global Fs

    %Reset Variables
    WFmatrix=[];
    SpikeTime=[];
    UnitStatus=[];
    auxWFmatrix=[];
    %Obtain Parameters
    L=round((Ls*Fs/(2*1000)));
    spkidx=0;
    ndataBlocks=1/slider_step; %Warning: Last datablock is lost if the loaded file is not a multiple of 40 MB
    aux=get(handles.Fil,'String');
    f=str2double(cell2mat(aux(get(handles.Fil,'Value'))));
    if (f==0 || threshold==0)
        errordlg('No Filter or Threshold found','Error');
    else
        for i=0:ndataBlocks
            slider=i*slider_step;
            % get data block
            move_raw_data(handles,slider);
            %%Filter
            Channel_and_Filter(handles);
            % Detect Spikes position
            [idxpos]=Spike_position;
            %Creation of waveforms
            if isempty(idxpos)==0
                for n=1:length(idxpos)
                    spkidx=spkidx+1;
                    auxWFmatrix(n,:)=dataView(idxpos(n)-L:idxpos(n)+L-1);
                    SpikeTime(spkidx,:)=idxpos(n)+i*MaxDataLen;     
                end
                WFmatrix=[WFmatrix;auxWFmatrix];
            end
                set(handles.percentDetection,'String',[num2str(i/ndataBlocks*100) '% completed'])
                pause(0.1)
                auxWFmatrix=[];
        end
        UnitStatus=zeros(spkidx,1);
        SpikeTime=SpikeTime/(Fs*60);

        %Plot Wfs, Raw Data and Features
        set(handles.percentDetection,'String','')
        set(handles.WFsN,'String',[num2str(size(WFmatrix,1)) ' WFs']);
        %Get Features
        get_Features(handles);
        %Reset Dataview 
        slider=get(handles.slider1,'Value');
        move_raw_data(handles,slider); 
        plot_erase(handles)
    end


%Detects Spike position on Dataview signal
function [idxpos]=Spike_position

    global Ls
    global threshold
    global dataView
    global Fs
    global ISI

    L=round((Ls*Fs/(2*1000)));
    idx=zeros(1,length(dataView));
    %Threshold Detection
    %Detects points surpassing the threshold
    if(threshold>0)
        idx(dataView > threshold)=1;
    else
        idx(dataView < threshold)=1;
    end
    %Detect Spikes position 
    idxpos=find(diff(idx)==1); %Detect local maximum or minimum 
    aux_idxpos=idxpos(end);
    idxpos=idxpos(diff(idxpos)>(ISI*(Fs/1000))); %%%% delete spikes ISI<2 ms;
    idxpos=[idxpos aux_idxpos]; %idxpos correction due to deletion of last idx (diff function use) 

    %Delete Spikes on edges shorter than 2*L
    if isempty(idxpos)==0
        if (idxpos(end)+L > length(dataView))
            idxpos=idxpos(1:end-1);
        end
        if (idxpos(1)-L < 1)
            idxpos=idxpos(2:end);
        end
    end

%Align Max.
%Ailgns waveforms according to their maximum value near to center.
function AlignMax_Callback(hObject, eventdata, handles)
    Align(handles,'max');

%Align Min
%Aligns waveforms according to their minimum value near to center.
function AlignMin_Callback(~, ~, handles)
    Align(handles,'min');

%Align to max or min
%Align Waveforms to one center by taking their maximum or minimum at
%reference point. It will only allign waveforms whose min or max is near to the
%center
function Align(handles,cond)
     global WFmatrix
     global Ls
     global slice
     global Fs
     L=round((Ls*Fs/(2*1000)));
     slice=zeros(1,2);
    eval(['[~,x]=' cond '(WFmatrix,[],2);']);
    ux= unique(x);
    if (length(ux)~=1)
       [N,X]= hist(x,ux);
       [~,auxx]=max(N);
       center=X(auxx);
       NotalgIndex=find(x~=center);
       for i=1:length(NotalgIndex);
           d=center-x(NotalgIndex(i));
           if (abs(d)<=round(L/5))
               aux=WFmatrix(NotalgIndex(i),:)';
               aux=circshift(aux,d);
               WFmatrix(NotalgIndex(i),:)=aux;
           end
       end
    end
    get_Features(handles);
    plot_erase(handles);


%Plot Waveforms
function plot_Wfs(handles)
    global BlockTime
    global WFmatrix
    global norm
    %Start Timer and block other functions
    time=start_timer(handles,BlockTime);
    drawnow;

    [objplot,h]=Wfs_plot_cfg(handles);
    plot2img(h,objplot);
    xlabel('Time (mS)','FontSize', 9);
    if norm==0
        ylabel('Voltage (uV)','FontSize', 9);
    else
        ylabel('Amplitude','FontSize', 9);
    end
    set(handles.WFsN,'String',[num2str(size(WFmatrix,1)) ' WFs'])

    %End Timer and unblock other functions
    end_timer(time,1,handles);

%PlotColors.
%Show the Waveforms plot with different colors for each waveform.
function PlotColors_Callback(~, ~, handles)
    global BlockTime
    global norm
    time=start_timer(handles,BlockTime);

    color = get(handles.PlotColors,'Value');
    if(color==1)
        [objplot,h]=Wfs_color_cfg(handles);
    else
        [objplot,h]=Wfs_plot_cfg(handles);
    end
    plot2img(h,objplot);
    xlabel('Time (mS)','FontSize', 9);
    if norm==0
        ylabel('Voltage (uV)','FontSize', 9);
    else
        ylabel('Amplitude','FontSize', 9);
    end

    end_timer(time,1,handles);
    set(handles.PlotColors,'Value',0);

%Plot selected FeaturePoints
function plot_sel_WFs(handles,IDX,SColor)  
    global WFmatrix
    global Fs
    
    axes(handles.WFaxes)
    hold on
    t=(1:length(WFmatrix(1,:)))*1000/Fs;
    objplot=plot(t,WFmatrix(IDX,:)','Color',SColor);
    h=gca;
    plot2img(h,objplot);

%Plot data on Waveforms axes as plot objects
function [objplot,h]=Wfs_plot_cfg(handles)
    global WFmatrix
    global showRange
    global Fs
    global GcaColor
    global WfsColor
    global slice
    slice=zeros(1,2);
    axes(handles.WFaxes)
    hold off
    t=(1:length(WFmatrix(1,:)))*1000/Fs;
    objplot=plot(t,WFmatrix(1:showRange:end,:)','Color',WfsColor);
    xlim([t(1),t(end)])
    set(gca,'Color',GcaColor);
    h=gca;

%Plot data on Waveforms axes as plot objects colored
function [objplot,h]=Wfs_color_cfg(handles)
    global WFmatrix
    global showRange
    global Fs
    global GcaColor
    axes(handles.WFaxes)
    hold off
    t=(1:length(WFmatrix(1,:)))*1000/Fs;
    objplot=plot(t,WFmatrix(1:showRange:end,:)');
    xlim([t(1),t(end)])
    set(gca,'Color',GcaColor);
    h=gca;





















%% FEATURES FUNCTIONS


%All Featureaxes
%Callback to "Plot Features". Plot Feature Points.
function PlotFeatures_Callback(hObject, eventdata, handles)
    get_Features(handles);
    plot_all_Features(handles);

%Get all Feature Points
%Calls the "Features.m" file to get the feature Points.
function get_Features(handles)
    global WFmatrix
    global FeaturePoints
    global SelectedUnitsIDX
    global FeaTime
    SelectedUnitsIDX=[];
    FeaturePoints=zeros(length(WFmatrix(:,1)),3);

    Feat1=get(handles.Feat1,'Value');
    Feat2=get(handles.Feat2,'Value');
    Feat3=get(handles.Feat3,'Value')-1;
    Ftime=get(handles.FeatTime,'Value');

    FeaturePoints(:,1)=Features(Feat1);
    FeaturePoints(:,2)=Features(Feat2);
    if Feat3==0
        FeaturePoints=FeaturePoints(:,1:2);
    else
        FeaturePoints(:,3)=Features(Feat3);
    end
    FeaTime=Features(Ftime);

%Feature Time
% Callback to Feature selection of Feature vs Time plit
% This is a reduced version of the previous function 
function FeatTime_Callback(hObject, eventdata, handles)
    global FeaTime
    Ftime=get(handles.FeatTime,'Value');
    FeaTime=Features(Ftime);
    plot_FeatureTime(handles)


%SLICES - PLOT ON WAVEFORM AXES
 
%Slice 1

%Callback to input value of the Slice X. Position the Slice on the selected
%time
function Slice_1_Callback(hObject, eventdata, handles)
    change_slice(handles,1)

%Callback to left slice button. Move the slice X to the left one sample.
function Left_1_Callback(hObject, eventdata, handles)
    move_slice_left(handles,1)

%Callback to right slice button. Move the slice X to the right one sample.
function Right_1_Callback(hObject, eventdata, handles)
    move_slice_right(handles,1)

    
%Slice_2

%Callback to input value of the Slice Y. Position the Slice on the selected
%time
function Slice_2_Callback(hObject, eventdata, handles)
    change_slice(handles,2)

%Callback to left slice button. Move the slice Y to the left one sample.
function Left_2_Callback(hObject, eventdata, handles)
    move_slice_left(handles,2)
    
%Callback to right slice button. Move the slice Y to the right one sample.
function Right_2_Callback(hObject, eventdata, handles)
    move_slice_right(handles,2)

%Change Slice position to input value
function change_slice(handles,n)
    global slice
    global Fs

    if slice(n)~=0
        delete(slice(n))
    end
    eval(['pos=str2num(get(handles.Slice_' num2str(n) ',''String''));']);
    if pos>(128*1000/Fs)
        pos=roundn(128*1000/Fs,-2);
    end
    if pos<0
        pos=0;
    end
    plot_slice(handles,n,pos)

%Move Slice to left
function move_slice_left(handles,n)
    global slice
    global Fs

    if slice(n)~=0
        delete(slice(n))
    end
    eval(['pos=str2num(get(handles.Slice_' num2str(n) ',''String''));']);
    pos=roundn((pos-(1000/Fs)),-2);
    if pos<0
        pos=0;
    end
    plot_slice(handles,n,pos)
    
%Move Slice to right
function move_slice_right(handles,n)
    global slice
    global Fs

    if slice(n)~=0
        delete(slice(n))
    end
    eval(['pos=str2num(get(handles.Slice_' num2str(n) ',''String''));']);
    pos=roundn((pos+(1000/Fs)),-2);
    if pos>(128*1000/Fs)
        pos=roundn(128*1000/Fs,-2);
    end
    plot_slice(handles,n,pos)

%Plot slice
function plot_slice(handles,n,pos)
    global slice
    color=['gr' 'r'];
    axes(handles.WFaxes)
    lim=axis;
    hold on
    slice(n)=plot([pos pos],[lim(3) lim(4)],color(n));
    eval(['set(handles.Slice_' num2str(n) ',''String'',num2str(pos));' ]);

    
%SLICES - PLOT FEATURES
    
%Plot all Feature Points
function plot_all_Features(handles)
    global BlockTime
 
    time=start_timer(handles,BlockTime);
    drawnow;
    
    [h,objplot]=plot_Features_cfg(handles);
    plot2img(h,objplot);
    plot_FeatureTime(handles)
    set(handles.Featureaxes,'ButtonDownFcn',@(hObject,eventdata)SpikeSort('Featureaxes_ButtonDownFcn',hObject,eventdata,guidata(hObject)))

    end_timer(time,1,handles);
    
%Plot Feature Time axes
function plot_FeatureTime(handles)
    global BlockTime
    
    time=start_timer(handles,BlockTime);
    drawnow;
    
    [h,objplot]=plot_FeatureTime_cfg(handles);
    plot2img(h,objplot);
    xlabel('Time (Min)','FontSize', 7);

    end_timer(time,1,handles);
    
%Plot Features Points
function [h,objplot]=plot_Features_cfg(handles);
    global FeaturePoints
    global showRange
    global ps
    global GcaColor
    global WfsColor

    aux=size(FeaturePoints);
    axes(handles.Featureaxes);
    h=gca;
    hold off
    if get(handles.Density,'Value')
        color=plot_density;
        objplot=scatter(FeaturePoints(1:showRange:end,1),FeaturePoints(1:showRange:end,2),ps,color(1:showRange:end,:,:,:),'.');

    elseif (aux(2)==3)
            color=plot_3D_Feature;
            objplot=scatter(FeaturePoints(1:showRange:end,1),FeaturePoints(1:showRange:end,2),ps,color(1:showRange:end,:,:,:),'.');
        else
            objplot=scatter(FeaturePoints(1:showRange:end,1),FeaturePoints(1:showRange:end,2),ps,WfsColor,'.');  
    end
    set(gca,'Color',GcaColor);    

%Plot Feature vs Time Points
function [h,objplot]=plot_FeatureTime_cfg(handles)
    global FeaTime
    global SpikeTime
    global showRange
    global ps
    global GcaColor
    global WfsColor

    axes(handles.FeatureTimeaxes);
    hold off
    objplot=scatter(SpikeTime(1:showRange:end),FeaTime(1:showRange:end),ps,WfsColor,'.');
    xlim([SpikeTime(1) SpikeTime(end)])
    set(gca,'Color',GcaColor);
    h=gca;

%Plot selected FeaturePoints
function plot_sel_FeaPoints(handles,IDX,SColor)
    global FeaturePoints
    global ps
    axes(handles.Featureaxes);
    hold on
    objplot=scatter(FeaturePoints(IDX,1),FeaturePoints(IDX,2),ps+10,SColor,'.');
    h=gca;
    plot2img(h,objplot);

%Plot selected Feature_Time
function plot_sel_FeaTime(handles,IDX,SColor)
    global FeaTime
    global SpikeTime
    global ps
    axes(handles.FeatureTimeaxes)
    hold on
    objplot=scatter(SpikeTime(IDX),FeaTime(IDX),ps+10,SColor,'.');
    h=gca;
    plot2img(h,objplot);
    
%Color
%Creates a color scale according to dot density on features plot
function Color=plot_density
    global FeaturePoints
    global DensityRes
    auxFeaturePoints=FeaturePoints(:,1:2);
    [N,C]=hist3(auxFeaturePoints,[DensityRes(1) DensityRes(2)]);
    x=cell2mat(C(1));
    y=cell2mat(C(2));

    lx=(x(2)-x(1))/2;
    ly=(y(2)-y(1))/2;
    maxn=max(max(N));
    I=[];
    for i=1:length(x)
        for j=1:length(y)
            I=find(FeaturePoints(:,1)>(x(i)-lx) & FeaturePoints(:,1)<(x(i)+lx) & FeaturePoints(:,2)>=(y(j)-ly) & FeaturePoints(:,2)<=(y(j)+ly));
             if (isempty(I)==0)
                 for k=1:length(I)
    %                 Color(I(k),:)=[0 0.8 0.8]*(N(i,j)/maxn)+[0 0.2 0.2];
                      Color(I(k),:)=[1 0.2 0]*(N(i,j)/maxn)+[0 0.3 1]*((1-N(i,j)/maxn));
                end
            end
        end
    end
    
    %Creates a color scale according to Feature (3) value
function Color=plot_3D_Feature
    global FeaturePoints 
    aux=length(FeaturePoints(:,3));
    Color=zeros(aux,3);
    maxf=max(FeaturePoints(:,3));
    minf=min(FeaturePoints(:,3));
    normFeaturePoints=(FeaturePoints(:,3)-minf)/(maxf-minf);
    for i=1:aux(1);
        Color(i,:)=[0 0.85 0]*normFeaturePoints(i)+[0 0.15 0];
    end














%% UNITS - CLUSTERING FUNCTIONS


%Sort
%Sort Spikes and create clusters with Automatic Algorithms.
%It calls the "clustering_algorithm.m" file to get the points according to the algorithm selected.
function Sort_Callback(hObject, eventdata, handles)
    global FeaturePoints
    global WFmatrix
    global UnitStatus   % N:cluster index
    global Unitn
    global ListC
    global BlockTime

    time=start_timer(handles,BlockTime);

    i=get(handles.Clusters,'Value');
    sizeW=size(WFmatrix);
    UnitStatus=zeros(sizeW(1),1);

    aux=get(handles.Clusters,'String');
    if strcmp(aux(i),'0')==0
        K=str2double(aux{i});
        Value=get(handles.Clus_Alg,'Value');
        [IDX, ~] = clustering_algorithm(FeaturePoints, K,Value,ListC);
        Unitn=max(IDX);
        for n=1:Unitn
            UnitStatus(IDX==n)=n;
            plot_unit(handles,n)
        end
    end

    end_timer(time,1,handles);

%Add_Unit
%Callback to "Add Unit" button. Creates a new cluster (Unit) of points
%according to selection.
function Add_Unit_Callback(hObject, eventdata, handles)
    global SelectedUnitsIDX
    global UnitStatus   % N:cluster index
    global Unitn
    global BlockTime

    if SelectedUnitsIDX~=0
        time=start_timer(handles,BlockTime);

        Unitn=Unitn+1;
        if Unitn==1
            plot_erase(handles)
        end
        UnitStatus(SelectedUnitsIDX)=Unitn;
        SelectedUnitsIDX=[];
        plot_unit(handles,Unitn)

        end_timer(time,1,handles);
    end

%ClearUnits
%Erase the clustering done.
function ClearUnits_Callback(hObject, eventdata, handles)
    global UnitStatus   % N:cluster index
    global Unitn

    Unitn=0;
    UnitStatus(1:end)=0;

    for i=1:5
        j=num2str(i);
        eval(['axes(handles.Unit' j ')']);
        cla
        eval(['set(handles.Unit' j ',''Visible'',''off'')']);
        eval(['set(handles.Un' j ',''Visible'',''off'')']);
        eval(['axes(handles.Hist' j ')']);
        cla
        eval(['set(handles.Hist' j ',''Visible'',''off'')']);
        eval(['set(handles.length' j ',''Visible'',''off'')']);
    end

    plot_erase(handles);
    
%Show Spikes
%Show spikes of each Unit separetly and also a shoot frequency histogram of each one. 
function ShowSpikes_Callback(hObject, eventdata, handles)
    global Unitn
    if Unitn>=1
        for n=1:Unitn
            plot_unit_Spikes(handles,n);
        end
    end

%Plot last created Unit
%Part of Callback to "Add unit"
function plot_unit(handles,n)
    global UnitStatus
    color=[0 0 1; 0 1 0 ; 1 0 0 ; 0 1 1; 1 0 1; 1 1 0];
    plot_sel_WFs(handles,UnitStatus==n,color(n,:));
    plot_sel_FeaPoints(handles,UnitStatus==n,color(n,:));
    plot_sel_FeaTime(handles,UnitStatus==n,color(n,:));

%Plot Spikes of each Cluster
function plot_unit_Spikes(handles,n)              
    global SpikeTime
    global WFmatrix
    global UnitStatus
    global GcaColor
    global WfsColor

    %Plot Unit Waveforms
    color=[0 0 1; 0 1 0 ; 1 0 0 ; 0 1 1; 1 0 1; 1 1 0];
    naux=num2str(n);
    x=find(UnitStatus==n);
    aux=diff(SpikeTime(x))*60*1000;
    l=num2str(length(aux));
    aux=aux(aux(:)<1000);
    eval(['axes(handles.Unit' naux ')']);
    hold off
    plot(WFmatrix(x,:)','Color',color(n,:))
    xlim([1 length(WFmatrix(1,:))])
    set(gca,'Color',GcaColor);
    eval(['set(handles.Unit' naux ',''xtick'',[])']);
    eval(['set(handles.Unit' naux ',''ytick'',[])']);
    eval(['set(handles.Un' naux ',''Visible'',''on'')']);
    eval(['set(handles.length' naux ',''String'',l)']);
    eval(['set(handles.length' naux ',''Visible'',''on'')']);

    %Create Shot Rate histogram
    eval(['axes(handles.Hist' naux ')']);
    hold off
    [n,x]=hist(aux,500);
    bar(x,n,'facecolor',WfsColor);
    hold on
    set(gca,'Color',GcaColor);
    xlabel('Time (mS)','FontSize', 9)

%Plot selected Units
%Part of Callback to "Plot Units"
function plot_selected_units(handles)
    global Unitn
    global BlockTime

    time=start_timer(handles,BlockTime);
    drawnow;

    checked=zeros(1,5);
    for i=1:5
        eval(['checked(' num2str(i) ')=get(handles.U' num2str(i) ',''Value'');']);
    end

    for i=1:Unitn
        if checked(i)==1 
            plot_unit(handles,i)
        end
    end
    end_timer(time,1,handles);


%Plot Selected units Callback
function PlotUnitsSel_Callback(hObject, eventdata, handles)
    plot_selected_units(handles);












%% SELECTION FUNCTIONS


%Rectangle Selection Calllbacks

% WFs 
%Detects Waveforms inside the rectangle.
function SelectWFs_Callback(hObject, eventdata, handles)
    global WFmatrix
    global SelectedUnitsIDX
    global Fs
    global BlockTime

    time=start_timer(handles,BlockTime);

    aux=size(WFmatrix);
    rect = getrect(handles.WFaxes); % [xmin ymin width height]
     rect(1)=rect(1)*Fs/1000;
     rect(3)=rect(3)*Fs/1000;
    if rect(1)<1
        rect(3)=rect(3)-(1-rect(1));
        rect(1)=1;

    end
    if (rect(1)+rect(3))>129
        rect(3)=aux(2)-rect(1);
    end
    rect(1)=round(rect(1));
    rect(3)=round(rect(3));
    minWFY=min(WFmatrix(:,rect(1):rect(1)+rect(3)),[],2);
    maxWFY=max(WFmatrix(:,rect(1):rect(1)+rect(3)),[],2);
    SelectedUnitsIDX=[SelectedUnitsIDX;find((minWFY>=rect(2) & minWFY<=rect(2)+rect(4)) | (maxWFY>=rect(2) & maxWFY<=rect(2)+rect(4)))];
    SelectedUnitsIDX=unique(SelectedUnitsIDX);
    plot_selected(handles);
    end_timer(time,1,handles);

%Featureaxes
%Detects Feature points inside the rectangle.
function SelectData_Callback(hObject, eventdata, handles)
    global FeaturePoints
    global SelectedUnitsIDX
    global BlockTime

    time=start_timer(handles,BlockTime);

    rect = getrect(handles.Featureaxes); 
    SelectedUnitsIDX=[SelectedUnitsIDX;find(FeaturePoints(:,1)>rect(1) & FeaturePoints(:,2)>rect(2) ...
        & FeaturePoints(:,1)<rect(1)+rect(3) & FeaturePoints(:,2)<rect(2)+rect(4))];
    plot_selected(handles);

    end_timer(time,1,handles);

%FeatureTime
%Detects points inside the rectangle of the Feature vs Time plot.
function SelRectFeaTime_Callback(hObject, eventdata, handles)
    global SpikeTime
    global FeaTime
    global SelectedUnitsIDX
    global BlockTime

    time=start_timer(handles,BlockTime); 

    rect = getrect(handles.FeatureTimeaxes); % [xmin ymin width height]
    SelectedUnitsIDX=[SelectedUnitsIDX;find(SpikeTime(:)>rect(1) & FeaTime(:)>rect(2) ...
        & SpikeTime(:)<rect(1)+rect(3) & FeaTime(:)<rect(2)+rect(4))];
    plot_selected(handles);

    end_timer(time,1,handles);


%Selection  Polygon Calllbacks

%Wfs
%Detects Waveforms inside the polygon.
function SelWfsPol_Callback(hObject, eventdata, handles)
    global Fs

    block_objects(handles); 
    set(handles.SpikeSort, 'pointer', 'arrow')

    axes(handles.WFaxes)
    h=imfreehand;
    aux=getPosition(h);
    x=(aux(:,1)')*Fs/1000;
    y=aux(:,2)';
    delete(h);
    detectSelectionWfs(x,y);
    plot_selected(handles);
    unblock_objects(handles); 

%Featureaxes
%Detects Feature points inside the polygon.
function SelFtPol_Callback(hObject, eventdata, handles)
    block_objects(handles); 
    set(handles.SpikeSort, 'pointer', 'arrow')
    
    axes(handles.Featureaxes)
    h=imfreehand;
    a=getPosition(h);
    x=a(:,1)';
    y=a(:,2)';
    delete(h);
    detectSelectionFeatures(x,y)
    plot_selected(handles);
    
    unblock_objects(handles); 



%FeaTIme
%Detects points inside the polygon of the Feature vs Time plot
function SelFeaTimePol_Callback(hObject, eventdata, handles)

    block_objects(handles); 
    set(handles.SpikeSort, 'pointer', 'arrow')
    
    axes(handles.FeatureTimeaxes)
    h=imfreehand;
    a=getPosition(h);
    x=(a(:,1)');
    y=a(:,2)';
    delete(h);
    detectSelectionFeaTime(x,y)
    plot_selected(handles);
    
    unblock_objects(handles); 


%Detection of selection (subfunctions of previous functions)

%Detection Wfs Polygon
function detectSelectionWfs(x,y)
    global WFmatrix   
    global SelectedUnitsIDX
    CantWfs=size(WFmatrix,1); %cantidad de Wfs
    aux=size(WFmatrix);
    minsc=round(min(x));
    maxsc=round(max(x));
    if minsc<1
        minsc=1;
    end
    if maxsc>aux(2)
        maxsc=aux(2);
    end
    auxWFs=double(WFmatrix(:,minsc:maxsc));
    j=minsc:maxsc;
    for i=1:CantWfs
        aux=auxWFs(i,:);
        in=inpolygon(j,aux,x,y);
        in=mean(in);
        if in>0
        SelectedUnitsIDX=[SelectedUnitsIDX;i];
        end
    end
    SelectedUnitsIDX=unique(SelectedUnitsIDX);

%Detection Featureaxes Polygon
function detectSelectionFeatures(x,y)  
    global SelectedUnitsIDX
    global FeaturePoints
    CantF=size(FeaturePoints,1); 
    for i=1:CantF
        in=inpolygon(FeaturePoints(i,1),FeaturePoints(i,2),x,y);
        if in==1
        SelectedUnitsIDX=[SelectedUnitsIDX;i];
        end
    end
        SelectedUnitsIDX=unique(SelectedUnitsIDX);

%Detection FeaTime Polygon
function detectSelectionFeaTime(x,y)  
    global SelectedUnitsIDX
    global FeaTime
    global SpikeTime
    CantF=size(FeaTime,1); 
    for i=1:CantF
        in=inpolygon(SpikeTime(i),FeaTime(i),x,y);
        if in==1
        SelectedUnitsIDX=[SelectedUnitsIDX;i];
        end
    end
        SelectedUnitsIDX=unique(SelectedUnitsIDX);

    
%Only one Feature Selection
%Callback to move the pointer over the Features points plot after clicking
%on it. Highlights the Spike corrresponding to the nearest feature point to
%the pointer.
function SpikeSort_WindowButtonMotionFcn(hObject, ~, handles)
    global FeaturePoints
    global MouseSelection
    global MouseSelectionI
    set(handles.Featureaxes,'ButtonDownFcn','');
    set(hObject, 'WindowButtonMotionFcn','');
    PanelPos=get(handles.Feat_Panel,'Position');% Position of the Feature Panel window on the screen
    FeaPos=get(handles.Featureaxes,'Position'); % Position of the Feature plot window on the panel
    SpkPos=get(hObject,'Position'); % Position of Spike sorter window on the screen in pixels
    PanelPos([1 3])=PanelPos([1 3])*SpkPos(3);
    PanelPos([2 4])=PanelPos([2 4])*SpkPos(4);
    Mouse=get(hObject,'CurrentPoint'); % Position of the pointer on the screen

    if (Mouse(1)>(FeaPos(1)*PanelPos(3)+PanelPos(1)) && Mouse(1)<(PanelPos(1)+(FeaPos(1)+FeaPos(3))*PanelPos(3)) ...
            && Mouse(2)>(PanelPos(2)+FeaPos(2)*PanelPos(4)) && Mouse(2)<(PanelPos(2)+(FeaPos(2)+FeaPos(4))*PanelPos(4))) % If the pointer is inside the Features plot window
        axes(handles.Featureaxes);
        lim=axis;
        nx=lim(2)-lim(1);
        ny=lim(4)-lim(3);
        Point=get(handles.Featureaxes,'CurrentPoint');
        %dist=((FeaturePoints(:,1)-Point(1,1))/nx).^2 + ((FeaturePoints(:,2)-Point(1,2))/ny).^2;
        aux1=find(abs(FeaturePoints(:,1)-Point(1,1))<nx/50);
        aux2= abs(FeaturePoints(aux1,2)-Point(1,2))<ny/50;
        aux1=aux1(aux2);
        if isempty(aux1)
            if MouseSelection~=0
                delete(MouseSelection(1));
                delete(MouseSelection(2));
                delete(MouseSelection(3));
                MouseSelection=0;
                MouseSelectionI=0;
            end
        else  
        dist=((FeaturePoints(aux1,1)-Point(1,1))/nx).^2 + ((FeaturePoints(aux1,2)-Point(1,2))/ny).^2;
        [~,I]=min(dist);
        I=aux1(I);
        if (I~=MouseSelectionI)
        MouseSelectionI=I;
            if MouseSelection==0
                plot_MouseSelection(I,handles);
            else
                delete(MouseSelection(1));
                delete(MouseSelection(2));
                delete(MouseSelection(3));
                plot_MouseSelection(I,handles);
            end
        end
    end
    path=@(hObject,eventdata)SpikeSort('SpikeSort_WindowButtonMotionFcn',hObject,eventdata,guidata(hObject));
    set(handles.SpikeSort, 'WindowButtonMotionFcn',path);
    else
        path=@(hObject,eventdata)SpikeSort('Featureaxes_ButtonDownFcn',hObject,eventdata,guidata(hObject));
        set(handles.Featureaxes,'ButtonDownFcn',path);
        if MouseSelection~=0
            delete(MouseSelection(1));
            delete(MouseSelection(2));
            delete(MouseSelection(3));
            MouseSelection=0;
            MouseSelectionI=0;
        end
    end

%Plot unique point selection
function plot_MouseSelection(I,handles)
    global MouseSelection
    global WFmatrix
    global SpikeTime
    global FeaTime
    global FeaturePoints
    global ps
    global Fs
    global SelColor

    axes(handles.WFaxes)
    hold on
    t=(1:length(WFmatrix(1,:)))*1000/Fs;
    MouseSelection(1)=plot(t,WFmatrix(I,:)','Color',SelColor);

    axes(handles.Featureaxes);
    hold on
    MouseSelection(2)=scatter(FeaturePoints(I,1),FeaturePoints(I,2),ps+10,SelColor,'.');

    axes(handles.FeatureTimeaxes)
    hold on
    MouseSelection(3)=scatter(SpikeTime(I),FeaTime(I),ps+10,SelColor,'.');
        
%Activate one Feature Selection by a click on Feature axes
function Featureaxes_ButtonDownFcn(hObject, eventdata, handles)
    path=@(hObject,eventdata)SpikeSort('SpikeSort_WindowButtonMotionFcn',hObject,eventdata,guidata(hObject));
    set(handles.SpikeSort, 'WindowButtonMotionFcn',path)
     
%Remove Selection

%Callback to "Deselect". Deselect all the selected data.
function DeselectData_Callback(hObject, eventdata, handles)
    global SelectedUnitsIDX
    SelectedUnitsIDX=[];
    plot_erase(handles);

% --- Delete Selected Data
%Erase all the selected data.
function DeleteData_Callback(hObject, eventdata, handles)
    global FeaturePoints
    global WFmatrix
    global SelectedUnitsIDX
    global SpikeTime
    global UnitStatus
    global FeaTime
    WFmatrix=WFmatrix(setxor(1:size(WFmatrix,1),SelectedUnitsIDX),:);
    SpikeTime=SpikeTime(setxor(1:length(SpikeTime),SelectedUnitsIDX));
    FeaturePoints=FeaturePoints(setxor(1:size(FeaturePoints,1),SelectedUnitsIDX),:);
    FeaTime=FeaTime(setxor(1:length(FeaTime),SelectedUnitsIDX));
    UnitStatus=UnitStatus(setxor(1:length(UnitStatus),SelectedUnitsIDX));
    SelectedUnitsIDX=[];
    plot_erase(handles);














%% GENERAL PLOT FUNCTIONS


%Plot all the figures without selection.
function plot_erase(handles)
    global BlockTime
    global fileType

    time=start_timer(handles,BlockTime);
    drawnow;

    %Waveforms
    plot_Wfs(handles)
    %Features
    plot_all_Features(handles)
    
    %Raw Data
    if fileType~=3
        Channel_and_Filter(handles);
        plot_channel(handles);
    end

    end_timer(time,1,handles);

%Plot selected data
function plot_selected(handles)
    global SelectedUnitsIDX
    global BlockTime
    global SelColor
    global fileType
  
    time=start_timer(handles, BlockTime);
    
    if not(isempty(SelectedUnitsIDX)) 
        drawnow;
        plot_sel_FeaPoints(handles,SelectedUnitsIDX,SelColor)
        plot_sel_FeaTime(handles,SelectedUnitsIDX,SelColor)
        plot_sel_WFs(handles,SelectedUnitsIDX,SelColor)
        %Raw Data
        if fileType~=3
            plot_sel_rawdata(handles,SelectedUnitsIDX,SelColor)
        end 
    end
    
    end_timer(time,1,handles);


%Converts the plot figure into an image figure to free up memory space
%The plot objects are deleted and it shows a capture of the plot
function plot2img(h,objplot)
    axes(h)
    lim=axis;
    %Get pixel position of axes
    Pos=round(getpixelposition(h));
    g=getframe(h,[0 0 Pos(3) Pos(4)]); %takes a capture of axes Method 1
    % g=getframe(h); Method 2
    s=size(g.cdata); %Method 1
    [g.cdata,xc,yc]=Im_correction(g.cdata,s); %corrects major dissaligments Method 1
    % g=getframe(h,[xc yc Pos(3) Pos(4)]); Method 2
    %Create Pixel Grid
    xstep=((lim(2)-lim(1))/(s(2)-1));
    ystep=-((lim(4)-lim(3))/(s(1)-1));
    x=lim(1):xstep:lim(2);
    y=lim(4):ystep:lim(3);
    delete(objplot); %Free memory space, delete Wfs plot objects
    hold off
    ylim([lim(3) lim(4)]); %Set limits
    xlim([lim(1) lim(2)]); %Set limits
    hold on
    im=image(x,y,g.cdata); %plots the capture 
    set(im,'HitTest','off')        

%Show Reduced.
%Shows only 10% of the Spikes and Feature Points. 
%This function is helpfull when plotting times are too long.
function showReduced_Callback(hObject, eventdata, handles)
    global showRange;
    reduced = get(handles.showReduced,'Value');
    if(reduced==1)
        showRange=10;
    else
        showRange=1;
    end
        
%Reset all_plots
function Reset_plots(handles)
    axes(handles.WFaxes)
    cla
    axes(handles.Featureaxes)
    cla
    set(handles.Featureaxes,'ButtonDownFcn',@(hObject,eventdata)SpikeSort('Featureaxes_ButtonDownFcn',hObject,eventdata,guidata(hObject)))
    axes(handles.FeatureTimeaxes)
    cla
    for i=1:5
        j=num2str(i);
        eval(['axes(handles.Unit' j ')']);
        cla
        eval(['set(handles.Unit' j ',''Visible'',''off'')']);
        eval(['set(handles.Un' j ',''Visible'',''off'')']);
        eval(['axes(handles.Hist' j ')']);
        cla
        eval(['set(handles.Hist' j ',''Visible'',''off'')']);
        eval(['set(handles.length' j ',''Visible'',''off'')']);
    end

%Corrects capture of getframe from outside axe acquisitions
function [cdata,xc,yc]=Im_correction(cdata,s)
    global GcaColor
    xc=0;
    yc=0;
    x_found=0;
    y_found=0;
    nx=1;
    ny=1;
    Gcacolor255=zeros(1,1,3);
    for k=1:3
        Gcacolor255(1,1,k)=round(GcaColor(k)*255);
    end
    %find pixel top row not in axes
    while isequal(cdata(nx,10,:),Gcacolor255)==0
            x_found=nx;
            nx=nx+1;
    end
    if x_found==0
        %find pixel bot row not in axes
        while isequal(cdata(end-(nx-1),10,:),Gcacolor255)==0
            x_found=nx;
            nx=nx+1;
        end
        %Delete bot rows
        if x_found>0
            cdata=cdata(1:s(1)-x_found,:,:);
            cdata_aux=zeros(x_found,s(2),s(3));
            for i=1:x_found
                for j=1:s(2)
                        cdata_aux(i,j,:)=Gcacolor255;
                end
            end
            %Add top rows
            cdata=[cdata_aux; cdata];
            xc=x_found;
        end
    %Delete top rows
    else
            cdata=cdata((x_found+1):s(1),:,:);
            cdata_aux=zeros(x_found,s(2),s(3));
            for i=1:x_found
                for j=1:s(2)
                    cdata_aux(i,j,:)=Gcacolor255;
                end
            end
            %Add bot rows
            cdata=[cdata;cdata_aux];
            xc=-x_found;
    end

    %find pixel left column not in axes
    while isequal(cdata(10,ny,:),Gcacolor255)==0
            y_found=ny;
            ny=ny+1;
    end

    if y_found==0
        %find pixel right column not in axes
        while isequal(cdata(10,end-(ny-1),:),Gcacolor255)==0
                y_found=ny;
                ny=ny+1;
        end
        %Delete right columns
        if y_found>0
            cdata=cdata(:,1:s(2)-y_found,:);
            cdata_aux=zeros(s(1),y_found,s(3));
            for j=y_found
                for i=1:s(1)
                    cdata_aux(i,j,:)=Gcacolor255;
                end
            end
            %Add left columns
            cdata=[cdata_aux cdata];
            yc=y_found;
        end
    %Delete left columns
    else
            cdata=cdata(:,(y_found+1):s(2),:);
            cdata_aux=zeros(s(1),y_found,s(3));
            for j=y_found
                for i=1:s(1)
                    cdata_aux(i,j,:)=Gcacolor255;
                end
            end
            %Add right columns
            cdata=[cdata cdata_aux];
            yc=-y_found;
    end
    
%Refresh all axes (Wfs, Features and FeaturevsTime)
function Refresh_axes_Callback(hObject, eventdata, handles)
    plot_erase(handles)
    plot_selected(handles)  
        
        

















%% Objects Creation DO NOT EDIT

function Channel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Feat1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Feat2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function FeatTime_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Scale_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function Clusters_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider2_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function Signalaxes_CreateFcn(hObject, eventdata, handles)
function popupmenu5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Th_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Slice_1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Slice_2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Fil_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Channels_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Frequency_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Feat3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function PointSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Clus_Alg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Length_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [calfilename_local] = SpectroSemantic_Neuro_model_glm_savio(MatfilePath,ValidKth_i, SWITCH, ParamModel,Cellname)
% [calfilename_local] = SpectroSemantic_Neuro_model_glm_savio(MatfilePath,ValidKth_i, SWITCH, ParamModel,Cellname)
% [OptimalFreqCutOff] = SpectroSemantic_Neuro_model_glm_savio(MatfilePath, SWITCH, ParamModel,Cellname)
% [PG_Index,FanoFactor_Index, Wins] = SpectroSemantic_Neuro_model_glm_savio(MatfilePath, SWITCH, ParamModel,Cellname)
%% Get the environment to figure out on which machine/cluster we are
fprintf(1,'The environment is: %s\n',getenv('HOSTNAME'))

if ~isempty(strfind(getenv('HOSTNAME'),'ln')) || ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))%savio Cluster
    Savio=1;
    fprintf(1, 'We are on savio!\n')
    addpath(genpath('/global/home/users/jelie/CODE/SingleUnitModels'));
    addpath(genpath('/global/home/users/jelie/CODE/GeneralCode'));
    addpath(genpath('/global/home/users/jelie/CODE/tlab/src/slurmbot/matlab'));
elseif ismac()
    Savio=0;
    Me = 1;
    fprintf(1, 'We are on my Mac Book Pro!\n')
    addpath(genpath('/Users/elie/Documents/CODE/SingleUnitModels'));
    addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'));
    addpath(genpath('/Users/elie/Documents/CODE/STRFLab/trunk'));
else %we are on strfinator or a cluster machine
    Savio = 0;
    Me = 0;
    fprintf(1, 'Hum.. We must be on one of the lab machine\n')
    addpath(genpath('/auto/fhome/julie/Code/SingleUnitModels'));
    addpath(genpath('/auto/fhome/julie/Code/GeneralCode'));
    addpath(genpath('/auto/fhome/julie/Code/strflab'));
end
%% Start a timer for the function
TimerVal=tic;

if nargin<1
    MatfilePath = '/auto/tdrive/julie/k6/julie/matfile/FirstVoc1s_Site3_L1250R1650_e13_s0_ss1.VariousKNeigh.mat';
end

%% Deal with input parameters
if nargin<3
    SWITCH = struct();
end
if ~isfield(SWITCH,'FanoFactor') || isempty(SWITCH.FanoFactor)
    SWITCH.FanoFactor=0;
end
if ~isfield(SWITCH,'BestBin') || isempty(SWITCH.BestBin)
    SWITCH.BestBin=0;
end
if ~isfield(SWITCH,'Models') || isempty(SWITCH.Models)
    SWITCH.Models=0;
end
if ~isfield(SWITCH,'InfoCal') || isempty(SWITCH.InfoCal)
    SWITCH.InfoCal=1;%Set to 1 if you want to calculate information on spike trains and change the name of the output file so they indicate "Info"
end

if ~isfield(SWITCH,'GaussWin') || isempty(SWITCH.GaussWin)
    SWITCH.GaussWin=0;%Set to 1 if you want to add the data about the gaussian window sizes used to convolve spike trains to the output file
end

if nargin<4
    ParamModel = struct();
end
if ~isfield(ParamModel,'LINK') || isempty(ParamModel.LINK)
    ParamModel.LINK='log'; %'identity'
end
if ~isfield(ParamModel,'DISTR') || isempty(ParamModel.DISTR)
    ParamModel.DISTR='poisson';%'normal'
end
if ~isfield(ParamModel,'NeuroRes') || isempty(ParamModel.NeuroRes)
    ParamModel.NeuroRes = 'count_gaussfiltered';
end
if  ~isfield(ParamModel,'MinWin') || isempty(ParamModel.MinWin)
    ParamModel.MinWin = 1; % end point of the first analysis window (spectrogram and neural response)
end
if ~isfield(ParamModel,'MaxWin') || isempty(ParamModel.MaxWin)
    ParamModel.MaxWin = 600; %end point of the last anaysis window for...
    ... neural response and end point of the largest analysis window for...
        ... spectrogram
end
if ~isfield(ParamModel,'MaxWin_cumInfo') || isempty(ParamModel.MaxWin_cumInfo)
    ParamModel.MaxWin_cumInfo = 600; %end point of the last anaysis window for...
    ... the calculation of cumulative information
end
if ~isfield(ParamModel,'Increment') || isempty(ParamModel.Increment)
    ParamModel.Increment = 10; %increase the size of the spectro window with a Xms pace
end
if ~isfield(ParamModel,'NeuroBin') || isempty(ParamModel.NeuroBin)
    ParamModel.NeuroBin = 10; % size of the window (ms) within which the neural response is analyzed
                               % The end of the window of analysis is
                               % determined by the Increment and ResDelay (see below).
end
if ~isfield(ParamModel,'ResDelay') || isempty(ParamModel.ResDelay)
    ParamModel.ResDelay = 0; % Delay in ms between the end of the...
    ... spectrogram window and the end of the neural response window
end

% Number of bootstraps
if ~isfield(ParamModel, 'NbBoot_Info') || isempty(ParamModel.NbBoot_Info)
    ParamModel.NbBoot_Info = 100;
end

if ~isfield(ParamModel, 'NbBoot_CumInfo') || isempty(ParamModel.NbBoot_CumInfo)
    ParamModel.NbBoot_CumInfo = 20;
end


if nargin<5
    [path,Cellname,ext]=fileparts(MatfilePath);
end

%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L1100R1450_e21_s0_ss1.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L2000R1600_e27_s1_ss1.mat')
if Savio %savio Cluster
    Dir_local='/global/scratch/jelie/MatFiles/';
    Res=loadfromTdrive_savio(MatfilePath, Dir_local);
elseif Me
    Dir_local='/Users/elie/Documents/CODE/data/matfile/FirstVoc1sMat/';
    if ~exist('ext','var')
        [~,Cellname,ext]=fileparts(MatfilePath);
    end
    Res = load([Dir_local Cellname ext]);
else
    Res = load(MatfilePath);
end



%% Get ready saving files and directories
%OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatSavio');

if Savio
    OutputDir_local='/global/scratch/jelie/MatFiles/ModMatInfo';
    OutputDirEx_local='/global/home/users/jelie/JobExecutableFiles';
    OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatInfo');
elseif Me
    if SWITCH.InfoCal || SWITCH.BestBin || SWITCH.FanoFactor || SWITCH.GaussWin
        OutputDir_local='/users/elie/Documents/CODE/data/matfile/ModMatInfo';
    else
        OutputDir_local='/users/elie/Documents/CODE/data/matfile/ModMatAcOnly';
    end
    OutputDir_final=OutputDir_local;
else
    if SWITCH.InfoCal || SWITCH.BestBin || SWITCH.FanoFactor || SWITCH.GaussWin
        OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatInfo');
    else
        OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatAcOnly');
    end
    OutputDir_local=OutputDir_final;
end
if SWITCH.InfoCal || SWITCH.BestBin || SWITCH.FanoFactor || SWITCH.GaussWin
    calfilename_local=fullfile(OutputDir_local,['InfoPoissonGF_' Res.Site '.mat']);
    calfilename_final=fullfile(OutputDir_final,['InfoPoissonGF_' Res.Site '.mat']);
else
    calfilename_local=fullfile(OutputDir_local,['Models_GLMPoisson_' Res.Site '.mat']);
    calfilename_final=fullfile(OutputDir_final,['Models_GLMPoisson_' Res.Site '.mat']);
end
outfilename_local=fullfile(OutputDir_local,['slurm_out*' Res.Site '*.txt']);
outfilename_final=fullfile(OutputDir_final,['slurm_out*' Res.Site '*.txt']);

if SWITCH.Models
    PrevData=0;
    if Savio
        try
            DoneCalc=loadfromTdrive_savio(calfilename_final, OutputDir_local,1);
            fprintf('Found some data for this unit\n')
                if isfield(DoneCalc, 'Deviance') && isfield(DoneCalc, 'LL') && isfield(DoneCalc, 'LambdaChoice') && isfield(DoneCalc, 'Model') && isfield(DoneCalc, 'PropVal') && isfield(DoneCalc, 'Data') && isfield(DoneCalc, 'Wins') && ~isempty(DoneCalc.Model.MeanSpectroStim{1})
                    PrevData = 1;
                else
                    frpintf('Data are not complete enough to be used\n')
                    system(['rm ' calfilename_local]);
                    clear 'DoneCalc'
                    PrevData = 0;
                end


        catch err
            fprintf('No Previous Data available or complete, working from the first window\nThe error is %s\n',err.identifier, err.message);
            PrevData = 0;
        end
    %         if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
    %             fprintf('No previous Data working from the first window\n');
    %             PrevData = 0;
    %         elseif strcmp(err.identifier, 'MATLAB:load:cantReadFile')
    %             fprintf('Previous Data File corrupted working from the first window\n');
    %             PrevData = 0;
    %         else
    %             fprintf('Error loading previous Data: %s\nworking from the first window\n',err.identifier);
    %             PrevData = 0;
    %         end
    else
        try
            DoneCalc=load(calfilename_final);
            fprintf('Found some data for this unit\n')
                if isfield(DoneCalc, 'Deviance') && isfield(DoneCalc, 'LL') && isfield(DoneCalc, 'LambdaChoice') && isfield(DoneCalc, 'Model') && isfield(DoneCalc, 'PropVal') && isfield(DoneCalc, 'Data') && isfield(DoneCalc, 'Wins') && ~isempty(DoneCalc.Data.MeanSpectroStim{1})
                    PrevData = 1;
                else
                    fprintf('Data are not complete enough to be used\n')
                    system(['rm ' calfilename_local]);
                    clear 'DoneCalc'
                    PrevData = 0;
                end


        catch err
            fprintf('No Previous Data available or complete, working from the first window\nThe error is %s\n',err.identifier, err.message);
            PrevData = 0;
        end
    end

    % save the path for now if no previous file
    if ~ PrevData
        save(calfilename_local,'MatfilePath', '-append')
    end
elseif SWITCH.InfoCal || SWITCH.BestBin || SWITCH.FanoFactor || SWITCH.GaussWin
    PrevData=0;
    if Savio
        try
            DoneCalc=loadfromTdrive_savio(calfilename_final, OutputDir_local,1);
            fprintf('Found some data for this unit\n')
            PrevData = 1;

        catch err
            fprintf('No Previous Data available or complete, working from the first window\nThe error is %s\n',err.identifier, err.message);
            PrevData = 0;
        end
    %         if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
    %             fprintf('No previous Data working from the first window\n');
    %             PrevData = 0;
    %         elseif strcmp(err.identifier, 'MATLAB:load:cantReadFile')
    %             fprintf('Previous Data File corrupted working from the first window\n');
    %             PrevData = 0;
    %         else
    %             fprintf('Error loading previous Data: %s\nworking from the first window\n',err.identifier);
    %             PrevData = 0;
    %         end
    else
        try
            DoneCalc=load(calfilename_final);
            fprintf('Found some data for this unit\n')
                PrevData = 1;
        catch err
            fprintf('No Previous Data available or complete, working from the first window\nThe error is %s\n',err.identifier, err.message);
            PrevData = 0;
        end
    end
end
    
%% Get the data ready
if SWITCH.Models % For models we use vocalization sections, only the first element of each vocalization sequence
    % Select first sections
    Firsts = find(Res.Voc_orders == 1);
else
    Firsts = 1:length(Res.VocType);
end
% Need to get rid of mlnoise sections and whine sections when they
% exist. I construct a vector of indices of the right sections
DataSel=zeros(1,length(Firsts));

nvoc=0;
voctype=Res.VocType;
for ii=1:length(Firsts)
    dd = Firsts(ii);
    if strcmp(voctype{dd}, 'Ag')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Be')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'DC')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Di')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'LT')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Ne')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Te')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Th')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'song')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
        %     elseif strcmp(voctype{dd}, 'Wh')
        %         nvoc=nvoc+1;
        %         DataSel(nvoc)=dd;
    end
end
DataSel=DataSel(1:nvoc);


%% Select the spectrograms of the selected stims
if SWITCH.Models
    Spectro.spec = Res.Spectro(DataSel);
    Spectro.to = Res.Spectroto(DataSel);
    Spectro.fo = Res.Spectrofo(DataSel);

    if ~PrevData
        % save the Stimuli Spectrograms for now if not done previously
        save(calfilename_local,'Spectro')
    end
end

%% Extract Emitter ID
if SWITCH.Models
    Ename = cell(length(DataSel),1);
    Esex = cell(length(DataSel),1);
    Eage = cell(length(DataSel),1);
    Erelated = cell(length(DataSel),1);
    for ii=1:length(DataSel)
        dd=DataSel(ii);
        [Path,File,Ext] = fileparts(Res.Original_wavfiles{dd});
        Ename{ii} = File(1:10);
        Esex{ii} = File(12);
        Eage{ii} = File(13);
        Erelated{ii} = File(14);
    end
    Emitter.Ename = Ename;
    Emitter.Esex = Esex;
    Emitter.Eage = Eage;
    Emitter.Erelated = Erelated;
end

%% Compute coherence to determine the optimal window size, the scale at which neural response information is maximized
if SWITCH.BestBin
    ParamModel.Response_samprate = Res.Response_samprate;
    
    %Data processing
%     ParamModel.NeuroRes = 'count';
%     [HalfTrain1, HalfTrain2, NumTrials]=organiz_data4coherence(Res.Trials(DataSel),Res.PSTH(DataSel),ParamModel);
    
    % compute coherence
%     [CoherenceStruct]=compute_coherence_mean(HalfTrain1, HalfTrain2,Res.Response_samprate);
%     OptimalFreqCutOff.CoherenceRaw = CoherenceStruct.freqCutoff;
    
%     ParamModel.NeuroRes = 'count_gaussfiltered';
    %Data processing
%     [HalfTrain1, HalfTrain2, NumTrials]=organiz_data4coherence(Res.Trials_GaussFiltered(DataSel),Res.PSTH_GaussFiltered(DataSel),ParamModel);
%     [CoherenceStruct]=compute_coherence_mean(HalfTrain1, HalfTrain2,Res.Response_samprate);
    % Compute coherence on gaussian filtered spike trains
%     OptimalFreqCutOff.CoherenceGaussFilt = CoherenceStruct.freqCutoff;
    
    % Calculate the frequency of the gaussian filtered PSTH below which 99%
    % of the spectrum power density is contained
   
    
    OptimalFreqCutOff.Thresh = 80:99;
    for kk=1:5
        % First retrieve the PSTH calculated with the same # of nearest
        % neighbour Ntrial/d where d=1:Ntrials
        SignalTot_local = nan(nvoc,size(Res.PSTH_GaussFiltered{1},2));
        if kk<5 % Treating Neigh = NT/2, NT/3, NT/4, NT/5
            for vv=1:nvoc
                NT = size(Res.PSTH_GaussFiltered{vv},1);
                SignalTot_local(vv,:) = Res.PSTH_GaussFiltered{vv}(NT-kk,:);
            end
        else % Treating Neigh =1 = NT/NT
            for vv=1:nvoc
                SignalTot_local(vv,:) = Res.PSTH_GaussFiltered{vv}(1,:);
            end
            SignalTot_1dim=reshape(SignalTot_local',[size(SignalTot_local,1)*size(SignalTot_local,2),1]);
        end
        
        % Calculate the power spectrum of the signal
        Window = 0.2*Res.Response_samprate;
        [Pxx,F] = pwelch(SignalTot_1dim, Window, [],[],10000);
        Pxx_Perc = 100*cumsum(Pxx / sum(Pxx));
        
        % Identify the frequency cut-off corresponding to the list of
        % thresholds
        OptimalFreqCutOff.(sprintf('PowerSpectrumDensityKth%d',kk)) = nan(1,length(80:99));
        for Thresh = OptimalFreqCutOff.Thresh
            IndMax=find(Pxx_Perc > Thresh);
            OptimalFreqCutOff.(sprintf('PowerSpectrumDensityKth%d',kk))(Thresh) = F(IndMax(1));
        end
    end
    
    
    % According to this code 10ms is the best size for 97% of cells see
    % fig BestPSTHBin.fig
    % At that window size, the values of the FanoFactor over cells is very
    % close to 1. see fig PoissonFanoFactor.fig
    if PrevData
        save(calfilename_local,'MatfilePath', 'OptimalFreqCutOff', '-append');
    else
        %save(calfilename_local,'MatfilePath', 'OptimalFreqCutOff');
    end
end

%% Estimate Poisson assumption for data at the choosen bining
if SWITCH.FanoFactor
    [PG_Index,FanoFactor_Index, Wins] = PoissonGaussianNeuralResponses(Res.Trials(DataSel),ParamModel,SWITCH,Cellname);
    if PrevData
        save(calfilename_local,'PG_Index', 'FanoFactor_Index','Wins','-append');
    else
        save(calfilename_local,'PG_Index', 'FanoFactor_Index','Wins');
    end
    
end

%% Calculate information about stimuli along time
if SWITCH.InfoCal
    ParamModel.MarkovParameters_Cum_Info = [];% supressing the calculation of Markov approximation for the cumulative information
    ParamModel.ExactHist = [];% supressing the exact calculation of the cumulative information
    ParamModel.Response_samprate = Res.Response_samprate;
    
    Kth_Neigh_local =unique(cell2mat(Res.Kth_Neigh'));
    EffectifK = nan(1,length(Kth_Neigh_local));
    ValidKth = 1;
    for kk=1:length(Kth_Neigh_local)
        EffectifK(kk) = sum(cell2mat(Res.Kth_Neigh') == Kth_Neigh_local(kk));
        if kk>1 && (EffectifK(kk) == EffectifK(kk-1))
              ValidKth = ValidKth + 1;
        end
    end
      
    KthNeigh = Res.Kth_Neigh(DataSel);
    PSTH_All = Res.PSTH_GaussFiltered(DataSel);
    JK_All = Res.JackKnife_GaussFiltered(DataSel);
    %for kk=1:ValidKth
        PSTH_GaussFilteredK = cell(length(DataSel),1);
        JK_GaussFilteredK = cell(length(DataSel),1);
        
        % Find out the number of trials per stimulus and feed-in PSTH and
        % PSTH_JackKnife
        Ntrials_perstim = nan(length(DataSel),1);
        for ss=1:length(DataSel)
            Ntrials_perstim(ss) = length(Res.Trials{DataSel(ss)});
            IndKth = find(KthNeigh{ss}==ValidKth_i);
            PSTH_GaussFilteredK{ss} = PSTH_All{ss}(IndKth,:);
            JK_GaussFilteredK{ss} = JK_All{ss}{IndKth};
        end
        ParamModel.Mean_Ntrials_perstim = [mean(Ntrials_perstim) mean(Ntrials_perstim - 1)];
        % Calculate information
        Calfilename_localKth = sprintf('%s_Kth%d_%s',calfilename_local(1:(end-4)),ValidKth_i,calfilename_local((end-4):end));
        [ParamModel, Data_local.(sprintf('Kth%d',ValidKth_i)), InputData_local.(sprintf('Kth%d',ValidKth_i)), Wins]=info_cuminfo_callsemantic(PSTH_GaussFilteredK,JK_GaussFilteredK,Res.VocType(DataSel), ParamModel, Calfilename_localKth);
        
   % end
   if PrevData 
       fprintf(1,'appending to the file\n');
        Data.(sprintf('Kth%d',ValidKth_i)) = Data_local.(sprintf('Kth%d',ValidKth_i));
        InputData.(sprintf('Kth%d',ValidKth_i)) = InputData_local.(sprintf('Kth%d',ValidKth_i));
        save(Calfilename_localKth,'Data', 'InputData','Wins','ParamModel');
   else
       Data.(sprintf('Kth%d',ValidKth_i)) = Data_local.(sprintf('Kth%d',ValidKth_i));
       InputData.(sprintf('Kth%d',ValidKth_i)) = InputData_local.(sprintf('Kth%d',ValidKth_i));
       save(Calfilename_localKth,'Data', 'InputData','Wins','ParamModel');
   end
    
    ElapsedTime = toc(TimerVal);
    Days = floor(ElapsedTime/(60*60*24));
    ETRem = ElapsedTime - Days*60*60*24;
    Hours = floor(ETRem/(60*60));
    ETRem = ETRem - Hours*60*60;
    Minutes = floor(ETRem/60);
    ETRem = ETRem-60*Minutes;
    fprintf(1,'Code run for %d days %dh %dmin and %dsec\n',Days,Hours,Minutes,ETRem);
    
end

%return

%% Add the size of the Gaussian used to convolve spike trains as output in the file
if SWITCH.GaussWin
    HwidthSpikes = Res.HwidthSpikes(DataSel);
    save(calfilename_local,'HwidthSpikes', '-append');
end
    
%% Inject the data in the models
if SWITCH.Models
    if PrevData
        [LambdaChoice, Deviance, LL, Model, ParamModel, Data, PropVal, Wins] = GrowingModelsRidgeglmLambdaMod( Spectro, Res.VocType(DataSel), Res.PSTH_GaussFiltered(DataSel),Res.Trials_GaussFiltered(DataSel),Emitter, ParamModel,calfilename_local, DoneCalc);
    else
        [LambdaChoice, Deviance, LL, Model, ParamModel, Data, PropVal, Wins] = GrowingModelsRidgeglmLambdaMod( Spectro, Res.VocType(DataSel), Res.PSTH_GaussFiltered(DataSel),Res.Trials_GaussFiltered(DataSel),Emitter, ParamModel,calfilename_local);
    end
    ElapsedTime = toc(TimerVal);
    Days = floor(ElapsedTime/(60*60*24));
    ETRem = ElapsedTime - Days*60*60*24;
    Hours = floor(ETRem/(60*60));
    ETRem = ETRem - Hours*60*60;
    Minutes = floor(ETRem/60);
    ETRem = ETRem-60*Minutes;
    fprintf(1,'Corrected Code run for %d days %dh %dmin and %dsec\n',Days,Hours,Minutes,ETRem);
    fprintf(1,'Good calculation of AcSemAc\n');
    fprintf(1,'Threshold for coordinate descent corrected: relying on L2Norm of the parameters''vector\n');
    save(calfilename_local,'MatfilePath', 'LambdaChoice', 'Deviance','LL','Model','ParamModel','Data','PropVal','Wins','ElapsedTime','-append');
end

% if Savio
%     [Status1]=transfertoTdrive_savio(calfilename_local,calfilename_final);
%     [Status2]=transfertoTdrive_savio(outfilename_local,[OutputDir_final '/']);
%     if ~(Status1 || Status2)
%         system(['mv ' OutputDirEx_local '/JobToDoSavio/Ex*' Res.Site '*.txt ' OutputDirEx_local '/JobDoneSavio/'])
%     end
%     fprintf(1,'Ready to quit');
%     quit
% end

end


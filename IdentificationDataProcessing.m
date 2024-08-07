%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processing of the gait data
%
% By: Huawei Wang
% Date: 04/05/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


subjMass = [69.3, 97.8, 58.4, 75.2, 65.1, 56.0, 80.8, 61.5, 81, 61.6, 62.4, 69];
data_trials = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54", "run_63", "run_81", "run_99"];
subj = 11;
mass = subjMass(subj);
phase_list = [2];

for tPhase = phase_list

    for data_trial = data_trials

        dataFile = sprintf('D:/HuaweiWang/PortableSystem_DataAnalysis/Processed_data/Subj%02d/Subj%02d_%s.mat', ...
                           subj, subj, data_trial);
        osimModelFile = sprintf('D:/HuaweiWang/PortableSystem_DataAnalysis/Processed_data/Subj%02d/OS/UTmodel/gait2392_simbody_subj%02d_scaled_1.osim', ...
                           subj, subj);
%         MuscleParaFile = sprintf('D:/HuaweiWang/MuscleOptResults/MuscleParameterOptResults_linear/Subj%02d/mus_par.mat', subj);
        muscleNames = ["soleus_r", "lat_gas_r", "med_gas_r", "tib_ant_r"];
        coordNames = ["ankle_angle_r"];
        savingPath = sprintf('D:/HuaweiWang/ReflexIdParaData/Subj%02d/', subj);

        idDataPreparison(dataFile, osimModelFile, ....
            muscleNames, coordNames, mass, tPhase, savingPath)

    end
end

function idDataPreparison(dataFile, osimModelFile, ....
    muscleNames, coordNames, mass, tPhase, savingPath)

    if ~exist(savingPath, 'dir')
        mkdir(savingPath)
    end
    
    dataNames = strsplit(dataFile, '/');
    if exist(sprintf('%s/%s', savingPath, dataNames{end}), 'file')
        load(sprintf('%s/%s', savingPath, dataNames{end}));
    end

    % general parameters
    M = length(muscleNames);
    J = length(coordNames);

    % load muscle parameters
    % get muscle parameters
    [lce_opt0, lt_slack0, theta0, Fmax0] = ...
        getOsimMuscleParameter(osimModelFile, muscleNames);
       
    mus_par0 = [lce_opt0, lt_slack0, theta0, Fmax0/mass];

    % load propcessed experimental data
    processedData = importdata(dataFile);

    % get muscle activation
    muscleIndex = [];
    for musclename = muscleNames
        muscleIndex = [muscleIndex, find(processedData.EMG.DataLabel == musclename)];
    end
        mus_act = processedData.Resample.Sych.Average.EMG.ave_r(:, muscleIndex+1);

    % get joint torques
    jointIndex = [];
    for coordname = coordNames
        for coori = 1:length(processedData.Resample.Sych.IDTrqDataLabel)
            if contains(processedData.Resample.Sych.IDTrqDataLabel{coori}, coordname)
                jointIndex = [jointIndex, coori];
            end
        end
    end

    torque = processedData.Resample.Sych.Average.IDTrqData.ave_r(:, jointIndex)/mass;

    % get muscle length and moment arms
    ikData.data = processedData.Resample.Sych.Average.IKAngData.ave_r;
    ikData.colheaders = processedData.Resample.Sych.IKAngDataLabel;

    [lmt, ma] = getOsimMuscleLengthMA(osimModelFile, ikData, ...
        muscleNames(1:M), coordNames);

    % get the time interval
    hs = mean(processedData.Resample.Sych.Average.hsMatrix_right(:, 2) - processedData.Resample.Sych.Average.hsMatrix_right(:, 1))/100/100;

    % get ankle joint angle index
    for coori = 1:length(processedData.Resample.Sych.IKAngDataLabel)
        if contains(processedData.Resample.Sych.IKAngDataLabel{coori}, 'ankle_angle_r')
            ankle_ik_id = coori;
        end
    end

    % get Fy index
    for coori = 1:length(processedData.Resample.Sych.ForcePlateGRFDataLabel)
        if contains(processedData.Resample.Sych.ForcePlateGRFDataLabel{coori}, 'ground_force_vy')
            fy_id = coori;
        end
    end
    
    ankle_ik = processedData.Resample.Sych.Average.IKAngData.ave_r(:, ankle_ik_id);
    fy = processedData.Resample.Sych.Average.ForcePlateGRFData.ave_r(:, fy_id)/mass;
    
    if tPhase < 10
        [phase, fig, l1, l2] = getPhase(ankle_ik, fy);

        numPhase = round(max(phase));
    
    else
        N = length(ankle_ik);
        phase = round(linspace(0.500001, tPhase + 0.499999, N))';
        numPhase = round(max(phase));
        
        fig = figure();
        l1 = plot(ankle_ik, ...
                '-o', 'linewidth', 2);
        hold on
        l2 = plot(fy, ...
                '-', 'linewidth', 2);
        hold on
        title('PhaseSeparation')
        hold on
        for np = 1:numPhase
            phase_node = find(phase == np);
            plot([phase_node(end), phase_node(end)], ...
                [min(min(ankle_ik), min(fy)), max(max(ankle_ik), max(fy))], ...
            'b-', 'linewidth', 2)
            hold on
        end
    end
    
    dataNames = strsplit(dataFile, '/');
    load(sprintf('%s/%s', savingPath, dataNames{end}))
    
    idParaData.angle = ankle_ik;
    idParaData.mus_act = mus_act;
    idParaData.mus_par0 = mus_par0;
    idParaData.torque = torque;
    idParaData.lmt = lmt;
    idParaData.ma = ma;
    idParaData.hs = hs;
    idParaData.mass = mass;
    idParaData.(sprintf('phase%d', numPhase)) = phase;
    idParaData.mus_names = muscleNames;
    idParaData.coord_names = coordNames;

    dataNames = strsplit(dataFile, '/');
    save(sprintf('%s/%s', savingPath, dataNames{end}), 'idParaData');
    
    % save figure and close
    legend([l1, l2], ["ankle angle", "Fy"])
    savefig(fig, sprintf('%s/%s_phase%d.fig', savingPath, dataNames{end}(1:end-4), numPhase));
    close(fig)

end

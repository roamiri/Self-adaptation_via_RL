%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in dense urban femtocell network using 
%   greedy policy, i.e, using maximum power for all femtocells with random 
%   addition of femtocells to the network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FBS_out = PA_greedy( FBS_in, MBS, mue, Npower, fbsCount, femtocellPermutation, NumRealization, saveNum, kk)

%% Initialization
% clear all;
clc;
% format short
% format compact
%% Parameters
Pmin = 5;                                                                                                                                                                                                                                                                                                                                                                           %dBm
Pmax = 15; %dBm
%StepSize = (Pmax-Pmin)/Npower; % dB

%% Minimum Rate Requirements for N MUE users
q_mue = 4.0; q_fue=0.50;
FBS = FBS_in(1:fbsCount);
%% Initialize the new Agent (FBS)
%     
% sumQ = sumQ * 0.0;
% for j=1:size(FBS,2)-1
%     fbs = FBS{j};
%     sumQ = sumQ + fbs.Q; 
% end

%% Calc channel coefficients
    fbsNum = size(FBS,2);
    G = zeros(fbsNum+1, fbsNum+1); % Matrix Containing small scale fading coefficients
    L = zeros(fbsNum+1, fbsNum+1); % Matrix Containing large scale fading coefficients
    [G, L] = measure_channel_3GPP(FBS,MBS,mue,NumRealization);
    %% Main Loop
    for j=1:size(FBS,2)
        fbs = FBS{j};
%         fbs.P_index = index;
        fbs.P = Pmax;
        FBS{j} = fbs;
    end
    % calc FUEs and MUEs capacity
    SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, -174);
    mue.SINR = SINR_MUE_4(G, L, FBS, MBS, mue, -174);
    mue.C = log2(1+mue.SINR);
    for j=1:size(FBS,2)
        fbs = FBS{j};
        fbs.C_FUE = log2(1+SINR_FUE_Vec(j));
        FBS{j}=fbs;
    end
    answer.mue = mue;
    answer.FBS = FBS;
    for j=1:size(FBS,2)
        c_fue(1,j) = FBS{1,j}.C_FUE;
    end
    sum_CFUE = 0.0;
    for i=1:size(FBS,2)
        sum_CFUE = sum_CFUE + FBS{i}.C_FUE;
    end
    answer.C_FUE = c_fue;
    answer.sum_CFUE = sum_CFUE;
    answer.time = 0.0;
    QFinal = answer;
    save(sprintf('oct19/T1/pro_greedy_%d_%d.mat', fbsCount, saveNum),'QFinal');
    FBS_out = FBS;
end

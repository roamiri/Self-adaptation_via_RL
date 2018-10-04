%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in dense urban femtocell network using 
%   exhaustive search with random addition of femtocells to the network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FBS_out = PA_exhaustive_search( FBS_in, MBS, mue, Npower, fbsCount, femtocellPermutation, NumRealization, saveNum, kk)

%% Initialization
% clear all;
clc;
% format short
% format compact
%% Parameters
Pmin = 5;                                                                                                                                                                                                                                                                                                                                                                           %dBm
Pmax = 15; %dBm
actions = linspace(Pmin, Pmax, Npower);
%StepSize = (Pmax-Pmin)/Npower; % dB

%% Minimum Rate Requirements for N MUE users
q_mue = 4.0; q_fue=0.50;
FBS = FBS_in(1:fbsCount);
%% Initialize the new Agent (FBS)

%% Calc channel coefficients
%fbsNum = size(FBS,2);
%G = zeros(fbsNum+1, fbsNum+1); % Matrix Containing small scale fading coefficients
%L = zeros(fbsNum+1, fbsNum+1); % Matrix Containing large scale fading coefficients
[G, L] = measure_channel_3GPP(FBS,MBS,mue,NumRealization);
%% Main Loop
K = size(FBS,2); % Number of FBSs
R0max = 0;
Rmax = 0; % maximum sum rate for FUEs
Rate_array = zeros(1,K);
best_actions = zeros(1,K);
all_actions = permn(actions, K);
tStart = tic;
for i = 1:size(all_actions,1)
    p_ar = all_actions(i,:);
    
    % calc FUEs and MUEs transmission rate
    SINR_FUE_Vec = SINR_FUE_3(G, L, K, p_ar, MBS.P, -174);
    r0 = log2(1+SINR_MUE_5(G, L, K, p_ar, MBS.P, -174));

    if r0>=q_mue
        rate_array = log2(1+SINR_FUE_Vec);
        rsum = sum(rate_array);
        if(rsum>Rmax)
            R0max = r0;
            Rmax = rsum;
            Rate_array=rate_array;
            best_actions = p_ar;
        end
    end
end

final.r0 = R0max;
final.rsum=Rmax;
final.r = Rate_array;
final.p = best_actions;
final.time = toc(tStart);
save(sprintf('oct4/ex/pro_ex_%d_%d.mat', fbsCount, saveNum),'final');
FBS_out = FBS;
end

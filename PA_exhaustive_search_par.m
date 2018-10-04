%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in dense urban femtocell network using 
%   exhaustive search with random addition of femtocells to the network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FBS_out = PA_exhaustive_search_par( FBS_in, MBS, mue, Npower, fbsCount, femtocellPermutation, NumRealization, saveNum, kk)

if fbsCount<8
    FBS_out = FBS_in;
    return;
end
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

all_actions = permn(actions, 7);
tStart = tic;

Ans.r0=0; Ans.rsum=0; Ans.r=[]; Ans.p=[];
contain = zeros(11,18); %[R0max, Rmax, Rate_array, best_actions]

parfor j=1:11
    p8 = actions(j);
    fprintf('iter = %d \n', j);
    contain(j,:) = best_PA(G, L, K, MBS.P, all_actions, p8);
end

dummy = contain(1,2);
best = 1;
for j=1:11
    if rsum<contain(j,2)
        best =j;
    end
end

final.r0 = contain(best,1);
final.rsum=contain(best,2);
final.r = contain(best,3:10);
final.p = contain(best,11:18);
final.time = toc(tStart);
save(sprintf('oct4/ex/pro_ex_%d_%d.mat', fbsCount, saveNum),'final');
FBS_out = FBS;
end

%%
function [R0max, Rmax, Rate_array, best_actions] = best_PA(G, L, K, MBS_P, all_actions, p8)
    R0max = 0;
    Rmax = 0; % maximum sum rate for FUEs
    Rate_array = zeros(1,K);
    best_actions = zeros(1,K);
    for i = 1:size(all_actions,1)
        p_ar = [all_actions(i,:) p8];
    
        % calc FUEs and MUEs transmission rate
        SINR_FUE_Vec = SINR_FUE_3(G, L, K, p_ar, MBS_P, -174);
        r0 = log2(1+SINR_MUE_5(G, L, K, p_ar, MBS_P, -174));

        if r0>=4.0
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
end
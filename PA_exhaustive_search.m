%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in dense urban femtocell network using 
%   exhaustive search with random addition of femtocells to the network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FBS_out = PA_exhaustive_search( FBS_in, MBS, mue, Npower, fbsCount, femtocellPermutation, NumRealization, saveNum, kk)

if fbsCount<9
    FBS_out = FBS_in;
    return;
end
%% Initialization
% clear all;

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
all_actions = permn(actions, 7);
tStart = tic;
iterations = size(all_actions,1);

% p3_ar = permn(actions, 3);
p2_ar = permn(actions, 2);
% iterations1 = size(p2_ar,1);
for j=1:121
    fprintf('mini= %d ', j);
    for i = 1:iterations
        p_ar = [all_actions(i,:) p2_ar(j,:)];

        % calc FUEs and MUEs transmission rate
        SINR_FUE_Vec = SINR_FUE_3(G, L, K, p_ar, MBS.P, -174);
        r0 = log2(1+SINR_MUE_5(G, L, K, p_ar, MBS.P, -174));

        if r0>=q_mue
            rate_array = log2(1+SINR_FUE_Vec);
            if sum(rate_array>=q_fue)==fbsCount
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
end
    final.r0 = R0max;
    final.rsum=Rmax;
    final.r = Rate_array;
    final.p = best_actions;
    final.time = 0;
    save(sprintf('oct13/ex/pro_ex_%d_%d.mat', fbsCount, saveNum),'final');
% end
FBS_out = FBS;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in femtocell network using 
%   Reinforcement Learning with random adding of femtocells to the network
%   Independent Learning
%
function FBS_out = PA_beginQ( FBS_in, Npower, fbsCount,femtocellPermutation, NumRealization, saveNum, CL)

%% Initialization
% clear all;
clc;
% format short
% format compact
%% Parameters
Pmin = -20;                                                                                                                                                                                                                                                                                                                                                                           %dBm
Pmax = 25; %dBm
%StepSize = (Pmax-Pmin)/Npower; % dB

%% Minimum Rate Requirements for N MUE users
N = 3;
q_mue = 1.00; q_fue=1.0;
%% Q-Learning variables
% Actions
actions = linspace(Pmin, Pmax, Npower);

% States
states = allcomb(0:1, 0:3 , 0:3); % states = (I, dMUE , dBS)

% Q-Table
% Q = zeros(size(states,1) , size(actions , 2));
Q_init = ones(size(states,1) , Npower) * 0.0;
Q1 = ones(size(states,1) , Npower) * inf;
sumQ = ones(size(states,1) , Npower) * 0.0;
% meanQ = ones(size(states,1) , Npower) * 0.0;

alpha = 0.5; gamma = 0.9; epsilon = 0.1 ; Iterations = 50000;
%% Generate the UEs
 mue = UE(204, 207);
% mue(1) = UE(150, 150);
% mue(1) = UE(-200, 0);
% selectedMUE = mue(mueNumber);
MBS = BaseStation(0 , 0 , 50);
%%
% FBS = cell(1,2);
FBS = FBS_in(1:fbsCount);
%% Initialize the new Agent (FBS)
%     
% sumQ = sumQ * 0.0;
% for j=1:size(FBS,2)-1
%     fbs = FBS{j};
%     sumQ = sumQ + fbs.Q; 
% end
j=size(FBS,2);
fbs = FBS{j};
fbs = fbs.getDistanceStatus;
fbs = fbs.setQTable(Q_init);
FBS{j} = fbs;

%% Calc channel coefficients
    fbsNum = size(FBS,2);
    G = zeros(fbsNum+1, fbsNum+1); % Matrix Containing small scale fading coefficients
    L = zeros(fbsNum+1, fbsNum+1); % Matrix Containing large scale fading coefficients
    [G, L] = measure_channel(FBS,MBS,mue,NumRealization);
    %% Main Loop
%     fprintf('Loop for %d number of FBS :\t', fbsCount);
%      textprogressbar(sprintf('calculating outputs:'));
    count = 0;
    errorVector = zeros(1,Iterations);
    dth = 25; %meter
    extra_time = 0.0;
    for episode = 1:Iterations
%          textprogressbar((episode/Iterations)*100);
        sumQ = sumQ * 0.0;
        for j=1:size(FBS,2)
            fbs = FBS{j};
            sumQ = sumQ + fbs.Q; 
        end
        %Choosing action for stable FBSs
        for j=1:size(FBS,2)-1
            fbs = FBS{j};
            kk = fbs.s_index;
            if CL == 1 
                [M, index] = max(sumQ(kk,:));     % CL method
            else                                    
                [M, index] = max(fbs.Q(kk,:));   %IL method
            end
            fbs.P_index = index;
            fbs.P = actions(index);
            FBS{j} = fbs;
        end
        j = size(FBS,2);
        if (episode/Iterations)*100 < 80
            % Action selection with epsilon=0.1
                fbs = FBS{j};
                if rand<epsilon
                      index = floor(rand*Npower+1);
                      fbs.P_index = index;
                      fbs.P = actions(index);
                else
                    kk = fbs.s_index;
                    if CL == 1 
                        [M, index] = max(sumQ(kk,:));     % CL method
                    else                                    
                        [M, index] = max(fbs.Q(kk,:));   %IL method
                    end
                      fbs.P_index = index;
                      fbs.P = actions(index);
                end
                FBS{j} = fbs;
        else
                fbs = FBS{j};
                kk = fbs.s_index;
                
                if CL == 1 
                    [M, index] = max(sumQ(kk,:));     % CL method
                else                                    
                    [M, index] = max(fbs.Q(kk,:));   %IL method
                end
                fbs.P_index = index;
                fbs.P = actions(index);
                FBS{j} = fbs;
        end 
        % calc FUEs and MUEs capacity
        SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, -120);
        mue.SINR = SINR_MUE_4(G, L, FBS, MBS, mue, -120);
        mue.C = log2(1+mue.SINR);
        
        for j=1:size(FBS,2)
            fbs = FBS{j};
%             fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            fbs.C_FUE = log2(1+SINR_FUE_Vec(j));
            if fbs.C_FUE <= q_fue
                if (fbs.s_index>16), fbs.s_new = fbs.s_index-16; else, fbs.s_new = fbs.s_index; end
            else
                if (fbs.s_index>16), fbs.s_new = fbs.s_index; else, fbs.s_new = fbs.s_index+16; end
            end
            FBS{j}=fbs;
        end
        for j=1:size(FBS,2)
            fbs = FBS{j};
            qMax=max(fbs.Q,[],2);
            
            % CALCULATING REWARD
            beta = fbs.dMUE/dth;
            R = beta*fbs.C_FUE*(mue.C).^2 -(fbs.C_FUE-q_fue).^2 - (1/beta)*(mue.C-q_mue)^2;
            if j == size(FBS,2)
                d_reward = fbs.dr(episode) + (gamma^0) * R;
                fbs.dr = [fbs.dr d_reward];
            end
            kk = fbs.s_index;
            nextState = fbs.s_new;
            jjj = fbs.P_index;
            fbs.Q(kk,jjj) = fbs.Q(kk,jjj) + alpha*(R+gamma*qMax(nextState)-fbs.Q(kk,jjj));
            fbs.s_index = nextState;
            FBS{j}=fbs;
        end
        
        % break if convergence: small deviation on q for 1000 consecutive
%         errorVector(episode) =  sum(sum(abs(Q1-sumQ)));
%         if sum(sum(abs(Q1-sumQ)))<0.001 && sum(sum(sumQ >0))
%             if count>1000
% %                 episode;  % report last episode
%                 break % for
%             else
%                 count=count+1; % set counter if deviation of q is small
%             end
%         else
%             Q1=sumQ;
%             count=0;  % reset counter when deviation of q from previous q is large
%         end
    end
%     Q = sumQ;
    answer.mue = mue;
%     answer.Q = sumQ;
%     answer.Error = errorVector;
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
    answer.episode = episode;
    QFinal = answer;
    save(sprintf('Jun14/reward/pro_IL_%d_%d.mat', fbsCount, saveNum),'QFinal');
    FBS_out = FBS;
end

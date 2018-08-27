%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in femtocell network using 
%   Reinforcement Learning with random adding of femtocells to the network
%   Independent Learning
%
function FBS_out = PA_learn_rate( FBS_in, MBS, mue, Npower, fbsCount, femtocellPermutation, NumRealization, saveNum, CL)

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
N = 3;
q_mue = 4.0; q_fue=0.50;
%% Q-Learning variables
% Actions
actions = linspace(Pmin, Pmax, Npower);

% States
% states = allcomb(0:1, 0:3 , 0:3); % states = (I_MUE, I_FUE, dMUE , dBS)
states = allcomb(0:1, 0:1, 0:3 , 0:3);
% Q-Table
% Q = zeros(size(states,1) , size(actions , 2));
Q_init = ones(size(states,1) , Npower) * 0.0;
a_init = ones(size(states,1) , Npower) * 0;
Q1 = ones(size(states,1) , Npower) * inf;
sumQ = ones(size(states,1) , Npower) * 0.0;
% meanQ = ones(size(states,1) , Npower) * 0.0;

alpha = 0.5; gamma = 0.9; epsilon = 0.1 ; Iterations = 75000;
%% Generate the UEs
%  mue = UE(204, 207);
% mue(1) = UE(150, 150);
% mue(1) = UE(-200, 0);
% selectedMUE = mue(mueNumber);
% MBS = BaseStation(0 , 0 , 50);
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
fbs = fbs.setLearningRate(a_init);
FBS{j} = fbs;

%% Calc channel coefficients
    fbsNum = size(FBS,2);
    G = zeros(fbsNum+1, fbsNum+1); % Matrix Containing small scale fading coefficients
    L = zeros(fbsNum+1, fbsNum+1); % Matrix Containing large scale fading coefficients
    [G, L] = measure_channel_3GPP(FBS,MBS,mue,NumRealization);
    %% Main Loop
%     fprintf('Loop for %d number of FBS :\t', fbsCount);
%      textprogressbar(sprintf('calculating outputs:'));
    count = 0;
    errorVector = zeros(1,Iterations);
    dth = 25; %meter
    tt = tic;
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
                [M, index] = max(sumQ(kk,:));    %CL method
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
                kk = fbs.s_index;
                index = -1;
                if rand<epsilon
                      index = floor(rand*Npower+1);
                      fbs.P_index = index;
                      fbs.P = actions(index);
                else
                    if CL == 1 
                        [M, index] = max(sumQ(kk,:));    %CL method
                    else                                    
                        [M, index] = max(fbs.Q(kk,:));   %IL method
                    end
                      fbs.P_index = index;
                      fbs.P = actions(index);
                end
                fbs.alpha(kk,index) = 1 + fbs.alpha(kk,index);
                FBS{j} = fbs;
        else
                fbs = FBS{j};
                kk = fbs.s_index;
                if CL == 1 
                    [M, index] = max(sumQ(kk,:));    % CL method
                else                                    
                    [M, index] = max(fbs.Q(kk,:));   %IL method
                end
                fbs.P_index = index;
                fbs.P = actions(index);
                fbs.alpha(kk,index) = 1 + fbs.alpha(kk,index);
                FBS{j} = fbs;
        end 
        % calc FUEs and MUEs capacity
        SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, -174);
        mue.SINR = SINR_MUE_4(G, L, FBS, MBS, mue, -174);
        mue.C = log2(1+mue.SINR);
        
        for j=1:size(FBS,2)
            fbs = FBS{j};
%             fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            fbs.C_FUE = log2(1+SINR_FUE_Vec(j));
            if mue.C < q_mue, I_mue = 0; else, I_mue = 1; end
            if fbs.C_FUE < q_fue, I_fue = 0; else, I_fue=1; end
           fbs.s_new = 16*I_fue + fbs.index; % state set X_2
%            fbs.s_new = 16*I_mue + fbs.index; % state set X_3
%             fbs.s_new = fbs.index; % state set X_1
%               fbs.s_new = 32*I_mue+16*I_fue + fbs.index; % state set X_4
%             if mue.C <= q_mue
%                 if (fbs.s_index>16), fbs.s_new = fbs.s_index-16; else, fbs.s_new = fbs.s_index; end
%             else
%                 if (fbs.s_index>16), fbs.s_new = fbs.s_index; else, fbs.s_new = fbs.s_index+16; end
%             end
            FBS{j}=fbs;
        end
        for j=1:size(FBS,2)
            fbs = FBS{j};
            qMax=max(fbs.Q,[],2);
            
            % CALCULATING REWARD
            beta = fbs.dMUE/dth;
%             if mue.C < q_mue
%                 R = beta* fbs.C_FUE - (100/beta);
%             else
%                 R = beta* fbs.C_FUE - (1/beta)*(mue.C-q_mue)^2;
%             end
%              R = beta*fbs.C_FUE*(mue.C).^2 -(fbs.C_FUE-q_fue).^2 - (1/beta)*(mue.C-q_mue)^2;
             R = (fbs.C_FUE-q_fue).^3 + (1/beta)*(mue.C-q_mue)^3;
%              R = -(fbs.C_FUE-q_fue).^2 - (1/beta)*(mue.C-q_mue)^2;
%             R = fbs.C_FUE -(fbs.C_FUE-q_fue).^2;
%             if j == size(FBS,2)
%                 d_reward = fbs.dr(episode) + (gamma^episode) * R;
%                 fbs.dr = [fbs.dr d_reward];
%             end
            kk = fbs.s_index;
            nextState = fbs.s_new;
            jjj = fbs.P_index;
            alpha = 1/((fbs.alpha(kk,index))^0.77);
            fbs.Q(kk,jjj) = fbs.Q(kk,jjj) + alpha*(R+gamma*qMax(nextState)-fbs.Q(kk,jjj));
            fbs.s_index = nextState;
            FBS{j}=fbs;
        end
        
        % break if convergence: small deviation on q for 1000 consecutive
        fbs = FBS{size(FBS,2)};
        sumQ = fbs.Q;
        errorVector(episode) =  sum(sum(abs(Q1-sumQ)));
%         delta = 1 * errorVector(episode);
        if errorVector(episode)<0.01
            if count>1000
%                 episode;  % report last episode
                break % for
            else
                count=count+1; % set counter if deviation of q is small
            end
        else
            Q1=sumQ;
            count=0;  % reset counter when deviation of q from previous q is large
        end
    end
%     Q = sumQ;
    answer.mue = mue;
%     answer.Q = sumQ;
%     FBS{size(FBS,2)}.Error = errorVector;
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
    answer.time = toc(tt);
    QFinal = answer;
    save(sprintf('Aug26/T1/pro_x2_CL_%d_%d.mat', fbsCount, saveNum),'QFinal');
    FBS_out = FBS;
end

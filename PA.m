%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in femtocell network using 
%   Reinforcement Learning with random adding of femtocells to the network
%   Independent Learning
%
function PA_1( Npower, fbsCount,femtocellPermutation, NumRealization, saveNum, CL)

%% Initialization
% clear all;
clc;
% format short
% format compact
total = tic;
%% Parameters
Pmin = -20;                                                                                                                                                                                                                                                                                                                                                                           %dBm
Pmax = 25; %dBm
%StepSize = (Pmax-Pmin)/Npower; % dB

dth = 25;
Kp = 100; % penalty constant for MUE capacity threshold
Gmue = 1.37; % bps/Hz

K = 1000;
PBS = 50 ; %dBm
sinr_th = 1.64;%10^(2/10); % I am not sure if it is 2 or 20!!!!!
gamma_th = log2(1+sinr_th);
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
% sumQ = ones(size(states,1) , Npower) * 0.0;
% meanQ = ones(size(states,1) , Npower) * 0.0;

alpha = 0.5; gamma = 0.9; epsilon = 0.1 ; Iterations = 50000;
%% Generate the UEs
 mue = UE(204, 207);
% mue(1) = UE(150, 150);
% mue(1) = UE(-200, 0);
% selectedMUE = mue(mueNumber);
MBS = BaseStation(0 , 0 , 50);
%%
%Generate fbsCount=16 FBSs, FemtoStation is the agent of RL algorithm
FBS_Max = cell(1,16);
for i=1:3
%     if i<= fbsCount
        FBS_Max{i} = FemtoStation_3S(180+(i-1)*35,150, MBS, mue, 10);
%     end
end

for i=1:3
%     if i+3<= fbsCount
        FBS_Max{i+3} = FemtoStation_3S(165+(i-1)*30,180, MBS, mue, 10);
%     end
end

for i=1:4
%     if i+6<= fbsCount
        FBS_Max{i+6} = FemtoStation_3S(150+(i-1)*35,200, MBS, mue, 10);
%     end
end

for i=1:3
%     if i+10<= fbsCount
        FBS_Max{i+10} = FemtoStation_3S(160+(i-1)*35,240, MBS, mue, 10);
%     end
end

for i=1:3
%     if i+13<= fbsCount
        FBS_Max{i+13} = FemtoStation_3S(150+(i-1)*35,280, MBS, mue, 10);
%     end
end
%%
% 
FBS = cell(1,fbsCount);

for i=1:fbsCount
    FBS{i} = FBS_Max{femtocellPermutation(i)};
end

    %% Initialize Agents (FBSs)
%     permutedPowers = randperm(Npower,size(FBS,2));
    
    for j=1:size(FBS,2)
        fbs = FBS{j};
%         fbs = fbs.setPower(actions(permutedPowers(j)));
%         fbs.state(1,1) = 0;
        fbs = fbs.getDistanceStatus;
        fbs = fbs.setQTable(Q_init);
        FBS{j} = fbs;
    end
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
        
        if (episode/Iterations)*100 < 80
            % Action selection with epsilon=0.1
            for j=1:size(FBS,2)
                fbs = FBS{j};
                if rand<epsilon
%                     fbs = fbs.setPower(actions(floor(rand*Npower+1)));
                      index = floor(rand*Npower+1);
                      fbs.P_index = index;
                      fbs.P = actions(index);
                else
                    a = tic;
                    kk = fbs.s_index;
                    if CL == 1 
                        [M, index] = max(sumQ(kk,:));     % CL method
                    else                                    
                        [M, index] = max(fbs.Q(kk,:));   %IL method
                    end
%                     fbs = fbs.setPower(actions(index));
                      fbs.P_index = index;
                      fbs.P = actions(index);
                      
                end
                FBS{j} = fbs;
            end
        else
            for j=1:size(FBS,2)
                fbs = FBS{j};
                kk = fbs.s_index;
                
                if CL == 1 
                    [M, index] = max(sumQ(kk,:));     % CL method
                else                                    
                    [M, index] = max(fbs.Q(kk,:));   %IL method
                end
%                 fbs = fbs.setPower(actions(index));
                fbs.P_index = index;
                fbs.P = actions(index);
                FBS{j} = fbs;
            end
        end 
        % calc FUEs and MUEs capacity
        SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, -120);
        mue.SINR = SINR_MUE_4(G, L, FBS, MBS, mue, -120);
%             MUE = MUE.setCapacity(log2(1+MUE.SINR));
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
    save(sprintf('Jun13/R_1/pro_CL_%d_%d.mat', fbsCount, saveNum),'QFinal');
%     FBS_out = FBS;
end

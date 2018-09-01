%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           function for running PA from 1 to 40 femtocells
%                    Dense urban dual strip model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fullRun40(femtocellPermutation,saveNum)

%% Generate the UEs
mue = UE(25, 500);
MBS = BaseStation(0 , 0 , 33);

%Generate fbsCount=16 FBSs, FemtoStation is the agent of RL algorithm
FBS_Max = cell(1,10);
for i=1:5
    FBS_Max{i} = FemtoStation_dual_strip((i-1)*10+5,500-15, MBS, mue, 5, 0);
end

for i=1:5
    FBS_Max{i+5} = FemtoStation_dual_strip((i-1)*10+5,500+15, MBS, mue, 5, 1);
end

% for i=1:10
%     FBS_Max{i+20} = FemtoStation_dual_strip(155+(i-1)*10,185, MBS, mue, 5, 1);
% end
% 
% for i=1:10
%     FBS_Max{i+30} = FemtoStation_dual_strip(155+(i-1)*10,195, MBS, mue, 5, 1);
% end

FBS_in = cell(1,10);
for i=1:10
    FBS_in{i} = FBS_Max{femtocellPermutation(i)};
    FBS_in = PA_learn_rate(FBS_in, MBS, mue, 11, i,femtocellPermutation,1e2, saveNum, 0);
end

end

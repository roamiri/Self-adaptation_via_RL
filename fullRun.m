%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           function for running PA from 1 to 16 femtocells
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fullRun(femtocellPermutation,saveNum)

%% Generate the UEs
 mue = UE(204, 207);
MBS = BaseStation(0 , 0 , 50);

%Generate fbsCount=16 FBSs, FemtoStation is the agent of RL algorithm
FBS_Max = cell(1,16);
for i=1:3
        FBS_Max{i} = FemtoStation_3S(180+(i-1)*35,150, MBS, mue, 10);
end

for i=1:3
        FBS_Max{i+3} = FemtoStation_3S(165+(i-1)*30,180, MBS, mue, 10);
end

for i=1:4
        FBS_Max{i+6} = FemtoStation_3S(150+(i-1)*35,200, MBS, mue, 10);
end

for i=1:3
        FBS_Max{i+10} = FemtoStation_3S(160+(i-1)*35,240, MBS, mue, 10);
end

for i=1:3
        FBS_Max{i+13} = FemtoStation_3S(150+(i-1)*35,280, MBS, mue, 10);
end

FBS_in = cell(1,16);
for i=1:16
    FBS_in{i} = FBS_Max{femtocellPermutation(i)};
    FBS_in = PA_learn_rate(FBS_in, 32, i,femtocellPermutation,1e3, saveNum, 1);
end
end

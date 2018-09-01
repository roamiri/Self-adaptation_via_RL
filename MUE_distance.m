%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           function for running PA from 1 to 40 femtocells
%                    Dense urban dual strip model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MUE_distance(saveNum)

%% Generate the UEs
MBS = BaseStation(0 , 0 , 33);

femtocellPermutation = [7     2     8     6     9     1     4     3     5    10];
FBS_in = cell(1,10);
FBS_Max = cell(1,10);
for kk=1:10
    x = 15 * kk;
    mue = UE(25, 365+x);
    for i=1:5
        FBS_Max{i} = FemtoStation_dual_strip((i-1)*10+5,350+x, MBS, mue, 5, 0);
    end

    for i=1:5
        FBS_Max{i+5} = FemtoStation_dual_strip((i-1)*10+5,380+x, MBS, mue, 5, 1);
    end
    FBS_in = FBS_Max;
    PA_greedy(FBS_in, MBS, mue, 11, 5,femtocellPermutation,1e2, saveNum, kk);
end

end

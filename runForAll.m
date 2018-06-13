%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           function for running PA_IL_CL2 from 1 to 16 femtocells
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runForAll(femtocellPermutation,saveNum)

for i=1:16
    PA(32, i,femtocellPermutation,1e3, saveNum, 1);
end
end

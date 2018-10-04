%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Main Loop Runner in parallel:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parallel_R_pro(pref_poolSize)
parpool(pref_poolSize)
permutationsMat = zeros(500,10);

for i=1:100
    permutationsMat(i,:) = randperm(10,10);
end

parfor_progress(500);
 parfor i=1:100
    fullRun40(permutationsMat(i,:),i);
    pause(rand);
    parfor_progress;
 end
 parfor_progress(0); % Clean up
end

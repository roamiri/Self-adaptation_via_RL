
clear;
clc;

%%

D_R_Mat = cell(1,16);

% for j=1:10
    s = sprintf('Jun14/reward/pro_IL_16_1.mat');
    filename = strcat(s);
    if exist(s)
        load(filename);
        for i=1:16
            D_R_Mat{i} = QFinal.FBS{i}.dr;
        end
    end
% end
%%
for i=1:16
    figure;
    plot(D_R_Mat{i});
end
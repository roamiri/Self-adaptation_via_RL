
clear;
clc;

%%
MUE_C = [];    
min_FUE = [];
sum_FUE = [];
max_MUE = [];
max_FUE = [];
C_FUE_Mat = cell(1,16);
for i=1:16
    fprintf('FBS num = %d\t', i);
    maxmue = 0.;
    maxfue = 0.;
    mue_C = 0.;
    minfue = 0.;
    sumfue = 0.;
    c_fue_vec = zeros(1,i);
    Cnt = 0;
    lowCnt = 0;

    
    for j=1:100
        s = sprintf('Jun13/SO_1/pro_IL_%d_%d.mat',i,j);
        filename = strcat(s);
        if exist(s)
            load(filename);
%                 mue_C  = QFinal.mue.C;
%                 cc = sum(C(40000:size(C,2)))/(-40000+size(C,2)+1);
                mue_C = mue_C + QFinal.mue.C;
                sumfue = sumfue + QFinal.sum_CFUE;
                c_fue_vec = c_fue_vec + QFinal.C_FUE;
                Cnt = Cnt+1;
        end
    end
    fprintf('Total Cnt = %d\n',Cnt);
    max_MUE = [max_MUE maxmue];
    max_FUE = [max_FUE maxfue];
    MUE_C = [MUE_C mue_C/Cnt]; 
    sum_FUE = [sum_FUE sumfue/Cnt];
    C_FUE_Mat{i} = c_fue_vec./Cnt;
    min_FUE = [min_FUE min(C_FUE_Mat{i})];
end
%%
MUE_C_1 = [];    
min_FUE_1 = [];
sum_FUE_1 = [];
max_MUE_1 = [];
max_FUE_1 = [];
C_FUE_Mat_1 = cell(1,16);
for i=1:16
    fprintf('FBS num = %d\t', i);
    maxmue = 0.;
    maxfue = 0.;
    mue_C = 0.;
    minfue = 0.;
    sumfue = 0.;
    c_fue_vec = zeros(1,i);
    Cnt = 0;
    lowCnt = 0;
    
    for j=1:100
%         s = sprintf('Rref_1/R3_%d_%d.mat',i,j);
        s = sprintf('Jun14/learn_rate/pro_IL_1_%d_%d.mat',i,j);
        filename = strcat(s);
        if exist(s)
            load(filename);
%                 C = QFinal.mue.C_profile;
%                 cc = sum(C(40000:size(C,2)))/(-40000+size(C,2)+1);
                mue_C = mue_C + QFinal.mue.C;
                sumfue = sumfue + QFinal.sum_CFUE;
                c_fue_vec = c_fue_vec + QFinal.C_FUE;
                Cnt = Cnt+1;
        end
    end
    fprintf('Total Cnt = %d\n',Cnt);
    max_MUE_1 = [max_MUE_1 maxmue];
    max_FUE_1 = [max_FUE_1 maxfue];
    MUE_C_1 = [MUE_C_1 mue_C/Cnt]; 
    sum_FUE_1 = [sum_FUE_1 sumfue/Cnt];
    C_FUE_Mat_1{i} = c_fue_vec./Cnt;
    min_FUE_1 = [min_FUE_1 min(C_FUE_Mat_1{i})];
end
%%
MUE_C_2 = [];    
min_FUE_2 = [];
sum_FUE_2 = [];
max_MUE_2 = [];
max_FUE_2 = [];
C_FUE_Mat_2 = cell(1,16);
for i=1:16
    fprintf('FBS num = %d\t', i);
    maxmue = 0.;
    maxfue = 0.;
    mue_C = 0.;
    minfue = 0.;
    sumfue = 0.;
    c_fue_vec = zeros(1,i);
    Cnt = 0;
    lowCnt = 0;
    
    for j=1:100
%         s = sprintf('Rref_1/R3_%d_%d.mat',i,j);
        s = sprintf('Jun14/learn_rate/pro_IL_77_%d_%d.mat',i,j);
        filename = strcat(s);
        if exist(s)
            load(filename);
%                 C = QFinal.mue.C_profile;
%                 cc = sum(C(40000:size(C,2)))/(-40000+size(C,2)+1);
                mue_C = mue_C + QFinal.mue.C;
                sumfue = sumfue + QFinal.sum_CFUE;
                c_fue_vec = c_fue_vec + QFinal.C_FUE;
                Cnt = Cnt+1;
        end
    end
    fprintf('Total Cnt = %d\n',Cnt);
    max_MUE_2 = [max_MUE_2 maxmue];
    max_FUE_2 = [max_FUE_2 maxfue];
    MUE_C_2 = [MUE_C_2 mue_C/Cnt]; 
    sum_FUE_2 = [sum_FUE_2 sumfue/Cnt];
    C_FUE_Mat_2{i} = c_fue_vec./Cnt;
    min_FUE_2 = [min_FUE_2 min(C_FUE_Mat_2{i})];
end
%%
MUE_C_3 = [];    
min_FUE_3 = [];
sum_FUE_3 = [];
max_MUE_3 = [];
max_FUE_3 = [];
C_FUE_Mat_3 = cell(1,16);
for i=1:16
    fprintf('FBS num = %d\t', i);
    maxmue = 0.;
    maxfue = 0.;
    mue_C = 0.;
    minfue = 0.;
    sumfue = 0.;
    c_fue_vec = zeros(1,i);
    Cnt = 0;
    lowCnt = 0;
    
    for j=1:100
%         s = sprintf('Rref_1/R3_%d_%d.mat',i,j);
        s = sprintf('Jun15/ILCL/pro_CL_77_%d_%d.mat',i,j);
        filename = strcat(s);
        if exist(s)
            load(filename);
%                 C = QFinal.mue.C_profile;
%                 cc = sum(C(40000:size(C,2)))/(-40000+size(C,2)+1);
                mue_C = mue_C + QFinal.mue.C;
                sumfue = sumfue + QFinal.sum_CFUE;
                c_fue_vec = c_fue_vec + QFinal.C_FUE;
                Cnt = Cnt+1;
        end
    end
    fprintf('Total Cnt = %d\n',Cnt);
    max_MUE_3 = [max_MUE_3 maxmue];
    max_FUE_3 = [max_FUE_3 maxfue];
    MUE_C_3 = [MUE_C_3 mue_C/Cnt]; 
    sum_FUE_3 = [sum_FUE_3 sumfue/Cnt];
    C_FUE_Mat_3{i} = c_fue_vec./Cnt;
    min_FUE_3 = [min_FUE_3 min(C_FUE_Mat_3{i})];
end
%%
figure;
hold on;
grid on;
box on;
plot(ones(1,16)*1.0, '--k', 'LineWidth',1);
plot(MUE_C, '--*r', 'LineWidth',1,'MarkerSize',10);
plot(MUE_C_1, '--*b', 'LineWidth',1,'MarkerSize',10);
plot(MUE_C_2, '--*g', 'LineWidth',1,'MarkerSize',10);
% plot(MUE_C_4, '--*k', 'LineWidth',1,'MarkerSize',10);
title('MUE capacity','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
xlim([2 16]);
ylim([0 7]);
legend({'threshold','proposed RF','[9]'},'FontSize',14, 'FontWeight','bold');
%%
figure;
hold on;
grid on;
box on;
plot( ones(1,16)*1.0, '--k', 'LineWidth',1);
for i=1:16
    vec = C_FUE_Mat{i};
    vec_ref = C_FUE_Mat_1{i};
    vec_ilq = C_FUE_Mat_2{i};
%     vec_4 = C_FUE_Mat_4{i};
    for j=1:size(vec,2)
        plot(i,vec(j), '*r', 'LineWidth',1,'MarkerSize',10);
        plot(i,vec_ref(j), '*b', 'LineWidth',1,'MarkerSize',10);
        plot(i,vec_ilq(j), '*g', 'LineWidth',1,'MarkerSize',10);
%         plot(i,vec_4(j), '*k', 'LineWidth',1,'MarkerSize',10);
    end
end
plot(min_FUE, '--r', 'LineWidth',1,'MarkerSize',10);
plot(min_FUE_1, '--b', 'LineWidth',1,'MarkerSize',10);
plot(min_FUE_2, '--g', 'LineWidth',1,'MarkerSize',10);
% plot(min_FUE_4, '--k', 'LineWidth',1,'MarkerSize',10);
title('FUEs capacity','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
xlim([2 16]);
ylim([0 3.5]);
legend({'threshold','proposed RF','[9]'},'FontSize',14, 'FontWeight','bold');
%%
figure;
hold on;
grid on;
box on;
% plot( ones(1,16)*2.0, '--k', 'LineWidth',1 );
plot(sum_FUE, '--*r', 'LineWidth',1,'MarkerSize',10);
plot(sum_FUE_1, '--*b', 'LineWidth',1,'MarkerSize',10);
plot(sum_FUE_2, '--*g', 'LineWidth',1,'MarkerSize',10);
% plot(sum_FUE_4, '--*k', 'LineWidth',1,'MarkerSize',10);
title('SUM capacity of FUEs','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
xlim([2 16]);
ylim([0 16]);
legend({'proposed RF','[9]'},'FontSize',14, 'FontWeight','bold');


%%%
% figure;
% hold on;
% grid on;
% box on;
% plot( ones(2,16)*1.0, '--k', 'LineWidth',1 );
% plot(max_MUE, '--*r', 'LineWidth',1,'MarkerSize',10);
% title('max MUE capacity in high interference','FontSize',14, 'FontWeight','bold');
% xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
% ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
% legend({'threshold','RF1','RF2'},'FontSize',14, 'FontWeight','bold');
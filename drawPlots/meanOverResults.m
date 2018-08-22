
clear;
clc;

%%
MUE_C = [];    
min_FUE = [];
sum_FUE = [];
mean_FUE = [];
max_FUE = [];
C_FUE_Mat = cell(1,40);
for i=1:10
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
%         s = sprintf('Jun14/learn_rate/pro_IL_77_%d_%d.mat',i,j);
        s = sprintf('Aug21/T3/pro_IL_77_%d_%d.mat',i,j);
%         s = sprintf('Aug16/IL/pro_IL_77_%d_%d.mat',i,j);
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
   
    MUE_C = [MUE_C mue_C/Cnt]; 
    sum_FUE = [sum_FUE sumfue/Cnt];
    C_FUE_Mat{i} = c_fue_vec./Cnt;
    mean_FUE = [mean_FUE mean(C_FUE_Mat{i})];
    max_FUE = [max_FUE max(C_FUE_Mat{i})];
    min_FUE = [min_FUE min(C_FUE_Mat{i})];
end
%%
MUE_C_1 = [];    
min_FUE_1 = [];
sum_FUE_1 = [];
mean_FUE_1 = [];
max_FUE_1 = [];
C_FUE_Mat_1 = cell(1,16);
for i=1:10
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
        s = sprintf('Aug21/T4/pro_IL_77_%d_%d.mat',i,j);
%         s = sprintf('July10/ILCL/pro_CL_77_%d_%d.mat',i,j);
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
    
    MUE_C_1 = [MUE_C_1 mue_C/Cnt]; 
    sum_FUE_1 = [sum_FUE_1 sumfue/Cnt];
    C_FUE_Mat_1{i} = c_fue_vec./Cnt;
    min_FUE_1 = [min_FUE_1 min(C_FUE_Mat_1{i})];
    mean_FUE_1 = [mean_FUE_1 mean(C_FUE_Mat_1{i})];
    max_FUE_1 = [max_FUE_1 max(C_FUE_Mat_1{i})];
end
%%
MUE_C_2 = [];    
min_FUE_2 = [];
sum_FUE_2 = [];
mean_FUE_2 = [];
max_FUE_2 = [];
C_FUE_Mat_2 = cell(1,16);
for i=1:10
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
        s = sprintf('Aug21/T4/pro_CL_77_%d_%d.mat',i,j);
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
    
    MUE_C_2 = [MUE_C_2 mue_C/Cnt]; 
    sum_FUE_2 = [sum_FUE_2 sumfue/Cnt];
    C_FUE_Mat_2{i} = c_fue_vec./Cnt;
    min_FUE_2 = [min_FUE_2 min(C_FUE_Mat_2{i})];
    mean_FUE_2 = [mean_FUE_2 mean(C_FUE_Mat_2{i})];
    max_FUE_2 = [max_FUE_2 max(C_FUE_Mat_2{i})];
end
%%
MUE_C_3 = [];    
min_FUE_3 = [];
sum_FUE_3 = [];
mean_FUE_3 = [];
max_FUE_3 = [];
C_FUE_Mat_3 = cell(1,16);
for i=1:10
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
        s = sprintf('Aug20/Rings/IL/pro_IL_77_%d_%d.mat',i,j);
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
    
    MUE_C_3 = [MUE_C_3 mue_C/Cnt]; 
    sum_FUE_3 = [sum_FUE_3 sumfue/Cnt];
    C_FUE_Mat_3{i} = c_fue_vec./Cnt;
    min_FUE_3 = [min_FUE_3 min(C_FUE_Mat_3{i})];
    mean_FUE_3 = [mean_FUE_3 mean(C_FUE_Mat_3{i})];
    max_FUE_3 = [max_FUE_3 max(C_FUE_Mat_3{i})];
end
%%
figure;
hold on;
grid on;
box on;
plot(ones(1,10)*4.0, '--k', 'LineWidth',1);
plot(MUE_C, '--or', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
plot(MUE_C_1, '--ob', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
plot(MUE_C_2, '--og', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','g', 'MarkerEdgeColor','b');
% plot(MUE_C_3, '--ok', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
% title('MUE transmission rate','FontSize',12);%, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Transmission rate (b/s/Hz)','FontSize',12);%, 'FontWeight','bold');
% xlim([1 16]);
% ylim([0 4]);
legend({'qos', 'IL','CL'},'Interpreter','latex','FontSize',12);
% legend({'qos','IL+$\mathcal{X}_2$','CL+$\mathcal{X}_2$', 'CL+$\mathcal{X}_1$', '$CL+\mathcal{X}_3$'},'Interpreter','latex','FontSize',12);
% legend({'$\mathbf{\tilde{q}}$','IL+$\mathbf{X}_1$','CL+$\mathbf{X}_2$', 'X_3', 'X_4'},'FontSize',14, 'FontWeight','bold','Interpreter','latex');
% legend({'$\mathbf{\tilde{q}}$','IL+$\mathbf{X}_2$','CL+$\mathbf{X}_3$'},'FontSize',14, 'FontWeight','bold','Interpreter','latex');
% legend({'qos','proposed RF','proximity RF'},'FontSize',12);%, 'FontWeight','bold','Interpreter','latex');
%%
figure;
hold on;
grid on;
box on;
plot( ones(1,10)*0.50, '--k', 'LineWidth',1);
for i=1:10
    vec = C_FUE_Mat{i};
    vec_ref = C_FUE_Mat_1{i};
%     vec_ilq = C_FUE_Mat_2{i};
%     vec_4 = C_FUE_Mat_3{i};
    for j=1:size(vec,2)
%         plot(i,vec(j), 'sr', 'LineWidth',1.2,'MarkerSize',10, 'MarkerEdgeColor','r');
        plot(i,vec_ref(j), 'sb', 'LineWidth',1.2,'MarkerSize',10, 'MarkerEdgeColor','b');
%         plot(i,vec_ilq(j), '*g', 'LineWidth',1,'MarkerSize',10);
%         plot(i,vec_4(j), '*k', 'LineWidth',1,'MarkerSize',10);
    end
end
% plot(min_FUE, '--r', 'LineWidth',1.2,'MarkerSize',10);
plot(min_FUE_1, '--b', 'LineWidth',1.2,'MarkerSize',10);
% plot(min_FUE_2, '--g', 'LineWidth',1,'MarkerSize',10);
% plot(min_FUE_3, '--k', 'LineWidth',1,'MarkerSize',10);
% title('FUEs capacity','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Transmission rate(b/s/Hz)','FontSize',12);%, 'FontWeight','bold');
% xlim([2 16]);
% ylim([0 3.5]);
% legend({'qos','proposed RF','proximity RF'},'FontSize',12);%, 'FontWeight','bold');
%%
figure;

subplot(3,1,1);
hold on;
grid on;
box on;
plot(1:10, ones(1,10)*0.50, '--k', 'LineWidth',1);
errorbar(1:10, mean_FUE, max_FUE-min_FUE, '--or', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
%  xlim([1 16]);
% ylim([0 4.0]);
legend({'qos','IL'}, 'FontSize',12);
% legend({'qos', 'old pathloss'},'Interpreter','latex','FontSize',12);
% legend({'qos','IL+$\mathcal{X}_2$'},'FontSize',10, 'Interpreter','latex');


subplot(3,1,2);
hold on;
grid on;
box on;
plot(1:10, ones(1,10)*0.50, '--k', 'LineWidth',1);
errorbar(1:10, mean_FUE_1, max_FUE_1-min_FUE_1, '--ob', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
%  xlim([1 16]);
% ylim([0 4.0]);
% legend({'qos','proximity RF'},'FontSize',12);
legend({'qos', 'CL'},'Interpreter','latex','FontSize',12);
% legend({'qos','CL+$\mathcal{X}_2$'},'FontSize',10, 'Interpreter','latex');

subplot(3,1,3);
hold on;
grid on;
box on;
plot(1:10, ones(1,10)*.50, '--k', 'LineWidth',1);
errorbar(1:10, mean_FUE_2, max_FUE_2-min_FUE_2, '--og', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','g', 'MarkerEdgeColor','b');
xlim([1 10]);
% ylim([0 4.0]);
% legend('X_3');
legend({'qos', '$\rho$'},'Interpreter','latex','FontSize',12);
% legend({'qos','CL+$\mathcal{X}_1$'},'FontSize',10, 'Interpreter','latex');
%%
subplot(4,1,4);
hold on;
grid on;
box on;
plot(1:10, ones(1,10)*0.50, '--k', 'LineWidth',1);
errorbar(1:10, mean_FUE_3, max_FUE_3-min_FUE_3, '--ok', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
% xlim([2 15]);
% ylim([0 4.0]);
% legend({'qos','CL+$\mathcal{X}_3$'},'FontSize',10, 'Interpreter','latex');

supertitle('','FontSize',14, 'FontWeight','bold');

% xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
% ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
%%
figure;
hold on;
grid on;
box on;
% plot( ones(1,16)*2.0, '--k', 'LineWidth',1 );
plot(sum_FUE, '--or', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
plot(sum_FUE_1, '--ob', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
plot(sum_FUE_2, '--og', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','g', 'MarkerEdgeColor','b');
% plot(sum_FUE_3, '--ok', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
% title('SUM capacity of FUEs','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Sum transmission rate (b/s/Hz)','FontSize',12);%, 'FontWeight','bold');
% xlim([1 16]);
 ylim([0 15]);
legend({'IL','CL', '$\rho$'},'Interpreter','latex','FontSize',12);
% legend({'\alpha_1','\alpha_2', '\alpha_3'},'FontSize',14, 'FontWeight','bold');
% legend({'X_1','X_2', 'X_3', 'X_4'},'FontSize',14, 'FontWeight','bold');
% legend({'IL+$\mathbf{X}_2$','CL+$\mathbf{X}_3$'},'FontSize',14, 'FontWeight','bold','Interpreter','latex');
% zoomPlot to highlight a portion of the major plot 
% [p,z] = zoomPlot(x,y,[5 50],[0.33 0.35 0.3 0.55],[1 3]); 
% hold on 
% plot(f1) 
% legend hide

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
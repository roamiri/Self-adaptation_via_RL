%%
clear;
clc;
%%
T1 = 'oct4/ex/pro_ex_%d_%d.mat';
[MUE_C ,min_FUE ,sum_FUE ,mean_FUE ,max_FUE ,failed_FUE ,diff_FUE, P_sum_FUE] = performance(T1);
T2 = 'oct4/p9/pro_ex_9_%d_%d.mat';
[MUE_C(9) ,min_FUE ,sum_FUE(9) ,mean_FUE ,max_FUE ,failed_FUE(9) ,diff_FUE(9), P_sum_FUE(9)] = performance_iter(T2,9);
T3 = 'oct4/p10/pro_ex_10_1_%d.mat';
[MUE_C(10) ,sum_FUE(10) ,failed_FUE(10) ,diff_FUE(10), P_sum_FUE(10)] = performance_iter10(T3);
%%
% figure;
hold on;
grid on;
box on;
% plot(ones(1,10)*4.0, '--k', 'LineWidth',1);
plot(MUE_C, '--sk', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k');%, 'MarkerEdgeColor','b');

xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('MUE transmission rate (b/s/Hz)','FontSize',12);%, 'FontWeight','bold');
xlim([1 10]);
ylim([3 12]);
%%
% figure;
hold on;
grid on;
box on;
% plot( ones(1,16)*2.0, '--k', 'LineWidth',1 );
plot(sum_FUE, '--sk', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','k');%, 'MarkerEdgeColor','b');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Sum transmission rate (b/s/Hz)','FontSize',12);%, 'FontWeight','bold');
xlim([1 10]);
ylim([0 20]);
%%
% figure;
hold on;
grid on;
box on;
% plot( ones(1,16)*2.0, '--k', 'LineWidth',1 );
plot(failed_FUE./100., '--sk', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','k');%, 'MarkerEdgeColor','b');
% title('SUM capacity of FUEs','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Probability of FUEs with $\gamma_k<\Gamma_k$ ','Interpreter','latex','FontSize',12);%, 'FontWeight','bold');
xlim([1 10]);
ylim([0 1]);
%%
% figure;
hold on;
grid on;
box on;
% plot( ones(1,16)*2.0, '--k', 'LineWidth',1 );
plot(P_sum_FUE, '--sk', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','k');%, 'MarkerEdgeColor','b');
% title('SUM capacity of FUEs','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Sum power (mWatt)','FontSize',12);%, 'FontWeight','bold');
xlim([1 10]);
ylim([0 350]);
%%
function [MUE_C ,sum_FUE ,failed_FUE ,diff_FUE, P_sum_FUE] = performance_iter10(T)
    dum_final=0;
    dum=0;
    for j=1:1200
        s = sprintf(T,j);
        filename = strcat(s);
        if exist(s)
            load(filename);
            if dum<final.rsum
                dum=final.rsum;
                dum_final = final;
            end
        end
    end
    MUE_C = dum_final.r0; 
    sum_FUE = dum_final.rsum;
    C_FUE_Mat = dum_final.r;
    failed_FUE = 100*(sum(dum_final.r<0.5)/10);
    diff_FUE = 0;
    P_FUE_Mat_W = 10.^(dum_final.p/(10));
    P_sum_FUE = sum(P_FUE_Mat_W);
    fprintf('FBS 10 Done!\n');
end
%%
function [MUE_C ,min_FUE ,sum_FUE ,mean_FUE ,max_FUE ,failed_FUE ,diff_FUE, P_sum_FUE] = performance_iter(T,k)
   MUE_C = [];    
    min_FUE = [];
    sum_FUE = [];
    mean_FUE = [];
    max_FUE = [];
    failed_FUE = [];
    diff_FUE = []; % difference of rate of failed FUEs from QoS
    P_sum_FUE = [];
    P_FUE_Mat_W = cell(1,10);
    C_FUE_Mat = cell(1,40);
    allCombinations = cell(1,100);
    for i=1:100
        dum = 0.;
        dum_final = [];
%         fprintf('Comb num = %d\t', i);
        for j=1:200
            s = sprintf(T,i,j);
            filename = strcat(s);
            if exist(s)
                load(filename);
                if dum<final.rsum
                    dum=final.rsum;
                    dum_final = final;
                end
            end
        end
        allCombinations{1,i} = dum_final;
    end
    Cnt = 0;
    maxmue = 0.;
    maxfue = 0.;
    mue_C = 0.;
    minfue = 0.;
    sumfue = 0.;
    c_fue_vec = zeros(1,k);
    p_fue_vec = zeros(1,k);
    lowCnt = 0;
    failedFUE = 0;
    diffFUE = 0;
    for i=1:100
        final = allCombinations{1,i};
        if size(final)
            mue_C = mue_C + final.r0;
            sumfue = sumfue + final.rsum;
            c_fue_vec = c_fue_vec + final.r;
            failedFUE = failedFUE + sum(final.r<0.5);
            if sum((final.r<0.5)) > 0
                diffFUE = diffFUE + sum((final.r<0.5).*(0.5-final.r))./sum((final.r<0.5));
            end
            p_fue_vec = p_fue_vec + final.p;
            Cnt = Cnt+1;
        end
    end
    fprintf('FBS num = 9    Total Cnt = %d\n',Cnt);
    MUE_C = [MUE_C mue_C/Cnt]; 
    sum_FUE = [sum_FUE sumfue/Cnt];
    C_FUE_Mat = c_fue_vec./Cnt;
%         mean_FUE = [mean_FUE mean(C_FUE_Mat{i})];
%         max_FUE = [max_FUE max(C_FUE_Mat{i})];
%         min_FUE = [min_FUE min(C_FUE_Mat{i})];
    failed_FUE = [failed_FUE (failedFUE/(9*Cnt))*100];
    diff_FUE = [diff_FUE diffFUE./(Cnt)];
    P_FUE_Mat_W = 10.^(p_fue_vec./(10*Cnt));
    P_sum_FUE = [P_sum_FUE sum(P_FUE_Mat_W)];
end
%%
function [MUE_C ,min_FUE ,sum_FUE ,mean_FUE ,max_FUE ,failed_FUE ,diff_FUE, P_sum_FUE] = performance(T)
    MUE_C = [];    
    min_FUE = [];
    sum_FUE = [];
    mean_FUE = [];
    max_FUE = [];
    failed_FUE = [];
    diff_FUE = []; % difference of rate of failed FUEs from QoS
    P_sum_FUE = [];
    P_FUE_Mat_W = cell(1,10);
    C_FUE_Mat = cell(1,40);
    for i=1:10
        fprintf('FBS num = %d\t', i);
        maxmue = 0.;
        maxfue = 0.;
        mue_C = 0.;
        minfue = 0.;
        sumfue = 0.;
        c_fue_vec = zeros(1,i);
        p_fue_vec = zeros(1,i);
        Cnt = 0;
        lowCnt = 0;
        failedFUE = 0;
        diffFUE = 0;
        for j=1:100
    %         s = sprintf('Jun14/learn_rate/pro_IL_77_%d_%d.mat',i,j);
            s = sprintf(T,i,j);
    %         s = sprintf('Aug16/IL/pro_IL_77_%d_%d.mat',i,j);
            filename = strcat(s);
            if exist(s)
                load(filename);
    %                 mue_C  = QFinal.mue.C;
    %                 cc = sum(C(40000:size(C,2)))/(-40000+size(C,2)+1);
                    mue_C = mue_C + final.r0;
                    sumfue = sumfue + final.rsum;
                    c_fue_vec = c_fue_vec + final.r;
                    failedFUE = failedFUE + sum(final.r<0.5);
                    if sum((final.r<0.5)) > 0
                        diffFUE = diffFUE + sum((final.r<0.5).*(0.5-final.r))./sum((final.r<0.5));
                    end
                    p_fue_vec = p_fue_vec + final.p;
                    Cnt = Cnt+1;
            end
        end
        fprintf('Total Cnt = %d\n',Cnt);

        MUE_C = [MUE_C mue_C/Cnt]; 
        sum_FUE = [sum_FUE sumfue/Cnt];
        C_FUE_Mat{i} = c_fue_vec./Cnt;
%         mean_FUE = [mean_FUE mean(C_FUE_Mat{i})];
%         max_FUE = [max_FUE max(C_FUE_Mat{i})];
%         min_FUE = [min_FUE min(C_FUE_Mat{i})];
        failed_FUE = [failed_FUE (failedFUE/(i*Cnt))*100];
        diff_FUE = [diff_FUE diffFUE./(Cnt)];
        P_FUE_Mat_W{i} = 10.^(p_fue_vec./(10*Cnt));
        P_sum_FUE = [P_sum_FUE sum(P_FUE_Mat_W{i})];
    end
end
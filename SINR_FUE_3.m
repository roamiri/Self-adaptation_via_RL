function SINR = SINR_FUE_3(G, L, fbsNum, p_ar, MBS_P, sigma2)
    SINR = zeros(1,fbsNum);
    mbsP = 10^((MBS_P-30)/10);
    sigma = 10^((sigma2-30)/10);
    P_interf = 0.0;
    pAgent = zeros(1,fbsNum);
    for i=1:fbsNum
        pAgent(i) = 10.^((p_ar(i)-30)/10);
    end
    
    for i=1:fbsNum
        for j=1:fbsNum
            if i ~= j
                P_interf = P_interf + pAgent(j)*(G(j,i)/L(j,i));
            end
        end
        SINR(i) = (pAgent(i)*(G(i,i)/L(i,i)))/(mbsP*(G(fbsNum+1,i)/L(fbsNum+1,i))+P_interf+sigma);
    end
end
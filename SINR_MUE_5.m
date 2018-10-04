function SINR = SINR_MUE_5(G, L, fbsNum, p_ar, MBSP, sigma2)
    P_interf = 0.0;
    sigma = 10^((sigma2-30)/10);
    mbsP = 10^((MBSP-30)/10);
    for i=1:fbsNum
        pAgent = 10^((p_ar(i)-30)/10);
        P_interf = P_interf + pAgent*G(i,fbsNum+1)/L(i,fbsNum+1);
    end
    SINR = (mbsP*G(fbsNum+1,fbsNum+1)/L(fbsNum+1,fbsNum+1))/(P_interf+sigma);
end
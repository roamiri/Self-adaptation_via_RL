function [G, L] = measure_channel_3GPP(FBS,MBS,MUE,NumRealization)
    fbsNum = size(FBS,2);
    G = zeros(fbsNum+1, fbsNum+1);
    L = zeros(fbsNum+1, fbsNum+1);
    f=2.0; % Frequency in GHz
    for i=1:fbsNum
        xAgent = FBS{i}.X;
        yAgent = FBS{i}.Y;
        for j=1:fbsNum
            d = sqrt((xAgent-FBS{j}.FUEX)^2+(yAgent-FBS{j}.FUEY)^2);
            if i==j
                PL0 = 38.46+20*log10(d)+0.7+18.3+5.0;
                L(i,j) = 10^((PL0)/10);
            else
                PL0 = max(15.3+37.6*log10(d), 38.46+20*log10(d))+0.7+18.3+5.0+20.0;
                L(i,j) = 10^((PL0)/10);
            end        
        end
        d = sqrt((xAgent-MUE.X)^2+(yAgent-MUE.Y)^2);
        PL0 = max(15.3+37.6*log10(d), 38.46+20*log10(d))+0.7+18.3+5.0+20.0;
        L(i,fbsNum+1) = 10.^((PL0)/10);
        
        d = sqrt((MBS.X-FBS{i}.FUEX)^2+(MBS.Y-FBS{i}.FUEY)^2);
        PL_BS = 15.3+37.6*log10(d)+20;
        L(fbsNum+1,i) = 10^((PL_BS)/10);
    end
    d = sqrt((MBS.X-MUE.X).^2+(MBS.Y-MUE.Y).^2);
    PL_BS = 15.3+37.6*log10(d);
    L(fbsNum+1,fbsNum+1) = 10.^((PL_BS)/10);
    
    Hij = abs((1/sqrt(2)) * (randn(fbsNum+1, fbsNum+1, NumRealization)+1i*randn(fbsNum+1, fbsNum+1, NumRealization)));
    hij = Hij.^2;
    for i=1:fbsNum+1
        for j=1:fbsNum+1
            G(i,j)=(sum(hij(i,j,:))/NumRealization);
        end
    end
end
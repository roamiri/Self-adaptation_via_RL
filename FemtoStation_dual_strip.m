classdef FemtoStation_dual_strip
   properties
      X
      Y
      P
      P_index
      dBS
      dMUE
      dFUE
      FUEX
      FUEY
      M  % distance with MUE
      B  % distance with BS
%       dM1 = 15; dM2 = 30; dM3 = 45; 
      dM1 = 17.5; dM2 = 22.5; dM3 = 45; 
%       dB1 = 50; dB2 = 150; dB3 = 400;
      dB1 = 50; dB2 = 360; dB3 = 400;
      state = zeros(1,3)
      powerProfile = []
      C_FUE
      C_profile = []
      Q
      index = -1
      s_index = -1
      s_new
      stable = 0
      dr = [0]  %discounted_reward
      alpha
      Error
      apt_strip = 0
      SINR_min
      SINR_max
   end
   methods
      function obj = FemtoStation_dual_strip(xPos, yPos, BS, MUE, dFUE, aptStrip)
        obj.X = xPos;
        obj.Y = yPos;
        obj.dBS = sqrt((xPos-BS.X)^2 + (yPos-BS.Y)^2);
        obj.dMUE = nearest_MUE(xPos, yPos, MUE);% sqrt((xPos-MUE.X)^2 + (yPos-MUE.Y)^2); %distance to nearest MUE
        n = floor(2*rand);
        obj.dFUE = dFUE;
        obj.FUEX = xPos+dFUE*(-1)^n;
        obj.FUEY = yPos+dFUE*(-1)^n;
        obj.apt_strip = aptStrip;
      end
      
      function obj = setPower(obj,power)
%           obj.P = 10^((power-30)/10);
            obj.P = power;
%             obj.powerProfile = [obj.powerProfile power];
      end
      
      function obj = setQTable(obj, Q_init)
          obj.Q = Q_init;
      end
      function obj = setLearningRate(obj, a_init)
          obj.alpha = a_init;
      end
      function obj = setCapacity(obj,c)
        obj.C_FUE = c;
%         obj.C_profile = [obj.C_profile c];
      end
      function obj = getDistanceStatus(obj)
          if(obj.dMUE <= obj.dM1 )
              obj.state(2) = 0;
          elseif(obj.dMUE <= obj.dM2 )
              obj.state(2) = 1;
          elseif(obj.dMUE <= obj.dM3 )
              obj.state(2) = 2;
          else
              obj.state(2) = 3;
          end
          
          if(obj.dBS <= obj.dB1 )
              obj.state(3) = 0;
          elseif(obj.dBS <= obj.dB2 )
              obj.state(3) = 1;
          elseif(obj.dBS <= obj.dB3 )
              obj.state(3) = 2;
          else
              obj.state(3) = 3;
          end
          obj.index = 4*obj.state(2)+obj.state(3)+1;
          obj.s_index = 4*obj.state(2)+obj.state(3)+1;
      end
   end
end
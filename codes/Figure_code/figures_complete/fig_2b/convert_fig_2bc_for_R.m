DKL_1 = @(x) DJS(x,1,:);
DKL_2 = @(x) DJS(x,2,:);
DKL_3 = @(x) DJS(x,3,:);

for i=1:49
   DKL_m_1(i) = mean(DKL_1(i));
   DKL_m_2(i) = mean(DKL_2(i));
   DKL_m_3(i) = mean(DKL_3(i));

   DKL_e_1(i) = std(DKL_1(i))/sqrt(80);
   DKL_e_2(i) = std(DKL_2(i))/sqrt(80);
   DKL_e_3(i) = std(DKL_3(i))/sqrt(80);
end

N_vals = (4:2:100)';

fig_2b_R = [[N_vals;N_vals;N_vals],...
            [0.05*ones(49,1);0.1*ones(49,1);0.25*ones(49,1)],...
            [DKL_m_1';DKL_m_2';DKL_m_3'],...
            [DKL_e_1';DKL_e_2';DKL_e_3']];
        
        
% DKL_1 = @(x) DJS(x,1,:);
% DKL_2 = @(x) DJS(x,2,:);
% DKL_3 = @(x) DJS(x,3,:);
% 
% for i=1:49
%    DKL_m_1(i) = mean(DKL_1(i));
%    DKL_m_2(i) = mean(DKL_2(i));
%    DKL_m_3(i) = mean(DKL_3(i));
% 
%    DKL_e_1(i) = std(DKL_1(i))/sqrt(80);
%    DKL_e_2(i) = std(DKL_2(i))/sqrt(80);
%    DKL_e_3(i) = std(DKL_3(i))/sqrt(80);
% end
% 
% N_vals = (4:2:100)';
% 
% fig_2c_R = [[N_vals;N_vals;N_vals],...
%             [0.1*ones(49,1);0.2*ones(49,1);0.3*ones(49,1)],...
%             [DKL_m_1';DKL_m_2';DKL_m_3'],...
%             [DKL_e_1';DKL_e_2';DKL_e_3']];
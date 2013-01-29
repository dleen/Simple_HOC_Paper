Pt =[0.0666
    0.0777
    0.0769
    0.0721
    0.0659
    0.0603
    0.0548
    0.0493
    0.0447
    0.0404
    0.0365
    0.0331
    0.0299
    0.0269
    0.0246
    0.0223
    0.0203
    0.0183
    0.0167
    0.0151
    0.0138
    0.0125
    0.0114
    0.0103
    0.0094
    0.0086
    0.0078
    0.0070
    0.0063
    0.0058
    0.0053
    0.0048
    0.0044
    0.0040
    0.0036
    0.0033
    0.0030
    0.0027
    0.0024
    0.0022
    0.0020
    0.0018
    0.0016
    0.0015
    0.0013
    0.0012
    0.0011
    0.0010
    0.0009
    0.0008
    0.0007
    0.0006
    0.0006
    0.0005
    0.0004
    0.0004
    0.0004
    0.0003
    0.0003
    0.0002
    0.0002
    0.0002
    0.0002
    0.0001
    0.0001
    0.0001
    0.0001
    0.0001
    0.0001
    0.0001
    0.0001
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
    0.0000
         0
         0
         0
         0
         0
         0
         0
         0];



options = optimset('GradObj','off','LargeScale','on',...
'Display','off',...
'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-3,'TolX',1e-3);
% We set GradObj = on as we supply the gradient.
% We set LargeScale = on as that is the algorithm that uses the
% user supplied gradient.
% Display just shows some accuracy output.
% The final options are tolerances etc.

param_list_init = [-60.4,0.25];
% The minimization function
[param_list,fval,output] =...
fminunc(@(x)lin_func...
(x,Pt,10,6.295),param_list_init,options);

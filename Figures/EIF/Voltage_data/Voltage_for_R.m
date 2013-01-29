clear all

volt = importdata('./EIF_voltage_4.dat',' ');

volt(volt>100) = 100;

v_range = 31000:50999;

R_volt = volt(v_range,:);

dt = 5/100;

t = [0;cumsum(dt*ones(length(v_range)-1,1))];

%plot(t,R_volt)

volt_lab = ones(size(v_range'));

volt_R_data = [[R_volt(:,1);R_volt(:,2);R_volt(:,3);R_volt(:,4)],...
 [t;t;t;t],...
 [volt_lab;2*volt_lab;3*volt_lab;4*volt_lab]
];

bins = [];
for i=1:50
   bins = [bins;[1;2]]; 
end

bin_times_1 = [0; cumsum(10*ones(99,1))];
bin_times_2 = 9.95+bin_times_1;

bin_times = [bin_times_1,bin_times_2,bins];
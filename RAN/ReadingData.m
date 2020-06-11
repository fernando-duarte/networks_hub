num = xlsread('Data.xlsx','DataSample');
c=transpose(num(:,11));
w=transpose(num(:,10));
p=transpose(num(:,12)-num(:,10));

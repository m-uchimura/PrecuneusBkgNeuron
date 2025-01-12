function normV=calc_normvec(data)
datapoint=[[0 0 data(1)];
[0 1 data(2)];
[1 0 data(3)];
[1 1 data(4)];];
[coeff,score]=pca(datapoint);
basis = coeff(:,1:2);
normV = coeff(:,3);
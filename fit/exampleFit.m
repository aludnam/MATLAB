%this is in separate m file example_exp.m
% function y = example_exp(x,p);
% y = exp(p(1)*(x-p(2))+p(3);

x = 1:0.1:5;
y = example_exp(x, [1 5 1]); %creates exp function
y_exp = y + 0.05*randn(size(x)); %simulates noisy data

[pbest,perror,nchi]=nonlinft('example_exp' ,x,y_exp,ones(size(x)),[2 3 0.5],[1 1 1]); %fits the data wit exponential...
figure
plot(x,y,'xb',x,y_exp,'or',x,example_exp(x,pbest),'g'); %original data without noise - blue, noisy data - red, fit - green...
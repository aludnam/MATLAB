function [dist1, dist2, x_mu_rear, y_mu_rear] = rearange(x_mu,y_mu,x_known, y_known, niter)
% distacnce of the first point to the center point
distmat1 = [x_mu(1,:); y_mu(1,:)] - [x_known(1,:); y_known(1,:)];
d1 = distmat(distmat1);

% distacnce of the second point to the center point
distmat2 = [x_mu(2,:); y_mu(2,:)] - [x_known(1,:); y_known(1,:)];
d2 = distmat(distmat2);

%points (indexed as 1) is cleser to the center:

indx1c = find(d1<=d2);
indx2c = find(d2<d1);

x_mu_rear = zeros(2,niter);
y_mu_rear = zeros(2,niter);

% these stays were they were
x_mu_rear(1, indx1c) = x_mu(1,indx1c);
x_mu_rear(2, indx1c) = x_mu(2,indx1c);
y_mu_rear(1, indx1c) = y_mu(1,indx1c);
y_mu_rear(2, indx1c) = y_mu(2,indx1c);

% these are changed
x_mu_rear(1, indx2c) = x_mu(2,indx2c);
x_mu_rear(2, indx2c) = x_mu(1,indx2c);
y_mu_rear(1, indx2c) = y_mu(2,indx2c);
y_mu_rear(2, indx2c) = y_mu(1,indx2c);

dist1 = distmat([x_mu_rear(1,:); y_mu_rear(1,:)] - [x_known(1,:); y_known(1,:)]);
dist2 = distmat([x_mu_rear(2,:); y_mu_rear(2,:)] - [x_known(2,:); y_known(2,:)]);



function InpolySpline
% Synthetic spline
x_d_p = [-4:4];
y_d_p = [0 1 0 1 0 1 0 1 0];
x_q_p = [-4:.1:4];
y_q_p = spline(x_d_p,y_d_p,x_q_p);
x_d_m = [4:-1:-4];
y_d_m = [0 -1 -2 -1 0 -1 -2 -1 0];
x_q_m = [4:-0.1:-4];
y_q_m = spline(x_d_m,y_d_m,x_q_m);

% Random data
xrand=8.*(rand(1,1000)-0.5);
yrand=4.*(rand(1,1000)-0.5);

% Generate edge
xe=[x_q_p,x_q_m];
ye=[y_q_p,y_q_m];
id=inpolygon(xrand,yrand,xe,ye);

% Plot figure
figure(10); clf(10)
plot(x_q_p,y_q_p); hold on
plot(x_q_m,y_q_m); hold on
plot(xrand,yrand,'ob'); hold on
plot(xrand(id),yrand(id),'or')

end
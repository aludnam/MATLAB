function plotd (d,xlimvec)
subplot(2,1,1)
plot([0:length(d)-1],d);
makelegend(1, 'Ddivergence', 'iteration', 'Ddiv')
xlim(xlimvec);

subplot(2,1,2)
dd = diff(d);
dd(2:end+1)=dd;
dd(1)=0;
plot([0:length(d)-1],dd,'r')
makelegend(1, 'diff(Ddivergence)', 'iteration', 'diff(Ddiv)')
xlim(xlimvec);
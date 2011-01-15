n = 100;
h = waitbar(0,'Waiting ...');
for i=1:n
  waitbar(i/n);
  % here perform some stuff
end
close(h)

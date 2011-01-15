cmap1=zeros(256,3);
da=2*pi/256;
off=-0.5;
for i=1:256
    if i<40
      mult=(0.05+double(i)/40);
      if mult<0
          mult=0;
      end
      if mult>1
          mult=1;
      end
    else
      mult=1;
    end
    red=0.54*(1+cos(i*da+pi*2/3+off))*mult;  % red
    if red<0
        red=0;
    end
    if red>1
        red=1;
    end
    cmap1(i,1)=red;
    
    green=0.50*(1+cos(i*da+pi*4/3+off))*mult;  % green
    if green>1
        green=1;
    end
    if green<0
        green=0;
    end
    cmap1(i,2)=green;
    
    blue=0.5*(1+cos(i*da+off))*mult;   % blue
    if blue<0
        blue=0;
    end
    if blue>1
        blue=1;
    end
    cmap1(i,3)=blue;
end
rand1=zeros(256,256);
for i=1:256
    for j=1:256
        rand1(i,j)=double(j)/256*256*(11)/11;
    end
end
% image(rand1);
% axis equal
% colormap(cmap1);

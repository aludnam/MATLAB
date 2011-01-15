function [beadsAt,bestslice]=ExtractCoordinates(myres,s) 
%Error: Undefined function or method 'ExtractCoordiantes' for input arguments of type 'dip_image'

[mymax,maxpos]=max(myres,[],3);
mymax
fprintf('Please select bead coordinates by left clicking. Finish with right click\n');
beadsAt=dipgetcoords(100);
fprintf('%d beads selected\n',size(beadsAt,1)-1);
beadsAt(end,:)=[];

for b=1:size(beadsAt,1)
    mypos=beadsAt(b,:);
    mybead=extract(squeeze(mymax),[s s],mypos);

    mymask=mybead >= max(mybead);
    mybeadslice=extract(squeeze(maxpos),[s s],mypos);
    bestslice(b)= round(mean(mybeadslice,mymask));
    
%    Version 1
%     [mm,mpx]=max(mybead,[],1); mpx=squeeze(mpx);
%     [mm,mppx]=max(squeeze(mm),[],1);
%     bestslice(b,1) = double(mpx(mppx));
%     [mm,mpy]=max(mybead,[],2);mpy=squeeze(mpy);
%     [mm,mppy]=max(squeeze(mm),[],1);
%     bestslice(b,2) = double(mpy(mppy));


end


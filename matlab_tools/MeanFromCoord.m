function [sumBds,mybds,sumBds2,mybds2]=MeanFromCoord(myres,beadsAt,s,myres2)

mybds=newim([s s size(beadsAt,1)]);
mybds2=newim([s s size(beadsAt,1)]);

for b=1:size(beadsAt,1)
    mypos=beadsAt(b,:);
    mybead=extract(squeeze(myres),[s s],mypos);
    mybds(:,:,b-1)=DampEdge(mybead,0.4,2,1,2);
    mybead2=extract(squeeze(myres2),[s s],mypos);
    mybds2(:,:,b-1)=DampEdge(mybead2,0.4,2,1,2);
end

ref=squeeze(mybds(:,:,3));
for b=1:size(beadsAt,1)
    si=findshift(squeeze(mybds(:,:,b-1)),ref,'iter');
    mybds(:,:,b-1)=shift(squeeze(mybds(:,:,b-1)),-si);
    mybds2(:,:,b-1)=shift(squeeze(mybds2(:,:,b-1)),-si);
end

sumBds=squeeze(sum(mybds,[],3));
sumBds2=squeeze(sum(mybds2,[],3));
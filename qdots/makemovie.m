function h = makemovie(matrix, cmap)
% h = makemovie(matrix, cmap)
h=figure;
% Resize the figure window to size of movie required.
% 3) Record the size of the plot window:
winsize = get(h,'Position');
% 4) Adjust size of this window to include the whole figure window (if you require the axes, title and axis labels
%      in the movie):
winsize(1:2) = [0 0];
% 5) Set the number of frames:
numframes=size(matrix,3);
% 6) Create the MATLAB movie matrix:
A=moviein(numframes,h,winsize);
% 7) Fix the features of the plot window (ensures each frame of the movie
% is the same size):
set(h,'NextPlot','replacechildren')
% 8) Within a loop, plot each picture and save to MATLAB movie matrix:
for i=1:numframes
   imagesc(matrix(:,:,i)); % plot command
   set (gca, 'DataAspectRatio',[1 1 1]);
   if exist ('cmap', 'var')
       colormap(cmap);
   end   
   % add axis label, legends, titles, etc. in here
   A(:,i)=getframe(h,winsize);
end
% This procedure creates a movie stored in a special format that is only readable in MATLAB. The first thing you will want to do is to play the movie:
movie(h,A,30,3,winsize) 

    
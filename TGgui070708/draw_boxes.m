function draw_boxes(n_boxes,boxes_xy,rbox)
%n_boxes

for i=1:n_boxes
   x0_box=boxes_xy(i,1)-rbox-0.5;
   y0_box=boxes_xy(i,2)-rbox-0.5;
   x1_box=boxes_xy(i,1)+rbox+0.5;
   y1_box=boxes_xy(i,2)+rbox+0.5;
   if (x0_box<1) 
       x0_box=1;
   end
   if (y0_box<1) 
       y0_box=1;
   end
   if (x1_box>2048) 
       x1_box=2048;
   end
   if (y1_box>2048) 
       y1_box=2048;
   end
   if (boxes_xy(i,3)==1)
      line([x0_box x1_box x1_box x0_box x0_box],[y0_box y0_box y1_box y1_box y0_box],'Color','g','HitTest','off');
   end

   if (boxes_xy(i,3)==-2)
      line([x0_box x1_box x1_box x0_box x0_box],[y0_box y0_box y1_box y1_box y0_box],'Color','r','HitTest','off');
   end

   if (boxes_xy(i,3)==-1)
     line([x0_box x1_box x1_box x0_box x0_box],[y0_box y0_box y1_box y1_box y0_box],'Color','y','HitTest','off');
   end

   if (boxes_xy(i,3)==-3)  % overlapping
     line([x0_box x1_box x1_box x0_box x0_box],[y0_box y0_box y1_box y1_box y0_box],'Color','b','HitTest','off');
   end

   %textlabel=int2str(i);
   %text(x0_box,y0_box,textlabel,'Color',[1 0 0]);
end
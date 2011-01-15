function zvalues = interactivedendrogramvalues(fig)
% zvalues = interactivedendrogramvalues(fig)
% from http://www.mathworks.com/help/techdoc/ref/datacursormode.html#bsawkea-7
dcm_obj = datacursormode(fig);
stop = 0;
inp = input('How many nodes?\n');
for ii=1:inp
    
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
    
    disp('Click on a node, then press Return.')
    figure(fig)
    pause                            % Wait while the user does this.    
    c_info = getCursorInfo(dcm_obj);
    zvalues(ii)=c_info.Position(2);        
    fprintf('Zvalue:%g\n', zvalues(ii))
end


% while ~stop
%     set(dcm_obj,'DisplayStyle','datatip',...
%         'SnapToDataVertex','off','Enable','on')
%     
%     disp('Click on a node, then press Return.')
%     pause                            % Wait while the user does this.
%     
%     c_info = getCursorInfo(dcm_obj);
%     zvalues(ii)=c_info.Position(2);
%     ii=ii+1;
%     
%     inp = input('More? Yes = [anything but n]/ No = [n]\n','s');
%     
%     if strcmp(inp,'n')
%         stop=1;
%     end
% end

function winit_pix=init_w_meandata1res(dpixc,peval,res)
% winit_pix=init_w_meandata1res(dpixc,peval,res)
% first iteration initialize with mean(dpixc), then with results form
% previous comuopation (res.w)
if isempty(res);
    fprintf('Initializing W with: mean(dpixc,3)\n');
    a=mean(dpixc,3);
    winit_pix=init_w(peval.initw1_method,peval,a);    
    
else
    fprintf('Initializing W with: res.w\n');
    winit_pix=init_w(peval.initw_method,peval,res.w);    
end
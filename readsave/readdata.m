% readdata
mfprintf(peval.fid, 'Data read from:\n%s\n',datasource)
load (datasource);

[peval.nx, peval.ny, peval.nt]=size(dpixc);

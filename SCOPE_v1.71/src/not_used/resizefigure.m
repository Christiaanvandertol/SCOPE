function resizefigure(spfig,nx,ny,xo,yo,xi,yi, xend, yend)

if (nargin<9)
    xend = .97;
    yend = .97;
end

dx = (xend - xo - (nx-1) * xi)/nx;
dy = (yend - yo - (ny-1) * yi)/ny;
NoPlots = nx*ny;

for iy = ny:-1:1
    y = yo + (ny-iy) * (dy+yi);
    for ix = 1:nx
       	plotno = ix + (iy-1)*nx;
	    x    = xo + (ix-1) * (dx+xi);
   	    set(spfig(plotno),'Position',[x y dx dy])
    end % for ix
end % for iy
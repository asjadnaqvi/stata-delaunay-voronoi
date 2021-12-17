

clear all
set trace off
set scheme white_tableau, perm


cap cd "D:/Programs/Dropbox/Dropbox/STATA - DELAUNAY"
cap cd "C:\Users\asjad\Dropbox\STATA - DELAUNAY"


net d delaunay, from("C:\Users\asjad\Dropbox\STATA - DELAUNAY\installation")
net install delaunay, replace

help delaunay



set obs 2000

gen x = runiform(0,100)
gen y = runiform(0,100)




// create a donut
drop if sqrt((x-50)^2 + (y-50)^2) > 50
*drop if sqrt((x-50)^2 + (y-50)^2) < 20
local obs = _N 
di `obs'

gen id = _n
order id


delaunay x y, id(id) triangles hull voronoi

// all information is automatically stored in matrices
mat dir





**************************
// Delaunay triangles   //
**************************


	twoway ///
		(area tri_y tri_x, cmissing(n) nodropbase fc(gs14%10) lw(0.04) lc()) ///
		(line hull_y hull_x) ///
		(scatter y x , msize(0.2) mc(black) ) ///
		, ///
		legend(off) aspect(1) xsize(1) ysize(1)





**************************
// Voronoi tesselations //
**************************


	twoway ///
		(pcspike vor_y1 vor_x1 vor_y2 vor_x2, lw(0.08) lc()) ///
		(scatter y x , msize(0.1) mc(black) ) ///
				, ///
		legend(off)	aspect(1) xsize(1) ysize(1)

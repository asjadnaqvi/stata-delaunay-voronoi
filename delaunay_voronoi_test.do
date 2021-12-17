

clear all
set trace off

cd "D:/Programs/Dropbox/Dropbox/STATA - DELAUNAY"

net d delaunay, from("D:\Programs\Dropbox\Dropbox\STATA - DELAUNAY\installation")
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

mat drop _all

timer on 1
delaunay x y, id(id) triangles hull voronoi
timer off 1

timer on 2
voronoi
timer off 2


// all information is automatically stored in matrices
mat dir





**************************
// Delaunay triangles   //
**************************

timer on 3
	twoway ///
		(area tri_y tri_x, cmissing(n) nodropbase fc(gs14%10) lw(0.04) lc()) ///
		(line hull_y hull_x) ///
		(scatter y x , msize(0.2) mc(black) ) ///
		, ///
		legend(off) aspect(1) xsize(1) ysize(1)
timer off 3		

timer list


	twoway ///
		(area tri_y tri_x, cmissing(n) nodropbase fc(gs14%15) lw(0.04) lc()) ///
		(line hull_y hull_x, lw(thin)) ///
		(scatter y x , msize(0.15) mc(black) ) ///
		, ///
		title("Delaunay triangulation and convex hull") ///
		subtitle("Delaunay based on the S-hull algorithm implementation by Mapbox", size(2.2)) ///
		note("Total observation: `obs' points.. Start to add data to Stata in `r(t1)' seconds. Graph drawn in `r(t3)' seconds." "A project of the Stata Guide on Medium (https://medium.com/the-stata-guide).", size(1.5)) ///
		legend(off) aspect(1) xsize(1) ysize(1)
		
		graph export delaunay1_triangles.png, replace wid(2000)


**************************
// Voronoi tesselations //
**************************

timer on 4
	twoway ///
		(pcspike vor_y1 vor_x1 vor_y2 vor_x2, lw(0.08) lc()) ///
		(scatter y x , msize(0.1) mc(black) ) ///
				, ///
		legend(off)	aspect(1) xsize(1) ysize(1)
timer off 4		

timer list

	twoway ///
		(pcspike vor_y1 vor_x1 vor_y2 vor_x2, lw(0.08) lc()) ///
		(scatter y x , msize(0.15) mc(black) ) ///
				, ///
		title("Voronoi tesselations and box clipping") ///
		subtitle("Voronoi based on d3js implementation. Box clipping based on the Cohen–Sutherland algorithm", size(2.2)) ///
		note("Total observation: `obs' points. Start to add data to Stata in `r(t2)' seconds. Graph drawn in `r(t4)' seconds.", size(1.5)) ///
		legend(off)	aspect(1) xsize(1) ysize(1)
		
		graph export delaunay2_voronoi.png, replace wid(2000)
		
	//  "A project of the Stata Guide (https://medium.com/the-stata-guide)."	
*! Voronoi tesselations v1.2 04.09.2022: Export voronoi as files. Reshaping removed to optimize speed.
*! by Asjad Naqvi (asjadnaqvi@gmail.com, @AsjadNaqvi)
*
* Ver 1.11 05.03.2022: x and y flipped. faster export to Stata using views
* Ver 1.10 01.03.2022: polygons exported as shapes. rays fixed. other bug fixes.
* Ver 1.02 27.12.2021: rescale integration
* Ver 1.01 20.12.2021: minor fixes
* Ver 1.00 22.11.2021: first release


**************************
**************************
**					    **
**       Voronoi        **
**						**
**     Asjad Naqvi		**
**                      ** 
**    Last updated:     **
**     04 Sep 2022      **
**						** 
**    First release:    **
**     22 Nov 2021      **
**                      **
**************************
**************************



cap program drop voronoi

program define voronoi, eclass sortpreserve
	version 15
	
	syntax varlist(min = 2 max = 2 numeric), ///
		[ lines 			] ///   // export lines
		[ polygons 			] ///  	// export polygons
		[ OFFset(real 0.05) ]		// bounding box of the voronoi. default set at +5% * (max - min)
	
	
	mata: vorlines 	  = .
	mata: vorpolygons = .
	
	mata: voronoi_core(triangles, points, points2, halfedges, hull, `offset', vorlines, vorpolygons)
	
	if "`lines'" != "" {
			
		qui {
		preserve		
			clear
			getmata (vline_x1 vline_y1 vline_x2 vline_y2) = vorlines, force double replace		
			
			gen serial = _n
			
			lab var serial "Serial"
			lab var vline_x1 "Voronoi line: x1"
			lab var vline_y1 "Voronoi line: y1"
			lab var vline_x2 "Voronoi line: x2"
			lab var vline_y2 "Voronoi line: y2"	
			
			compress
			save _vorlines.dta, replace	
			
		restore
		}
		
	}
	
	
	if "`polygons'" != "" {
			vorpoly `varlist'
	}	
	
	
end	

************************
// 	   voronoi_core	  //		
************************


cap mata: mata drop voronoi_core()

mata // voronoi_core
function voronoi_core(triangles, points, points2, halfedges, hull, offs, vorl, vorp)
{
	// real matrix vorl, vorp
	
	coords = initialize(points2)

	triangles = select(triangles, (triangles[.,1] :< .)) // added 17.12.2021
	tri3 = colshape(triangles',3)'  // reshape triangles
	
	
	xmin = .
	xmax = .
	ymin = .
	ymax = .
	
	bounds(points2,xmin,xmax,ymin,ymax, offs) 
	
	// collect the voronoi centers in vorcenter
	num2 = rows(triangles) / 3  // drop the missing rows


	vorcenter = J(num2,2,.)
	for (i=1; i <= num2; i++) {
		vorcenter[i,.]  = circumcenter2(points2,triangles,i)
	}

	// triangle pairs from which the center points are extracted (point0, point1)
	voredges = J(rows(triangles),2,.)  // voronoi edge pairs indexed to vorcenter
	forEachVoronoiEdge(triangles,halfedges,voredges)
	voredges = select(voredges, (voredges[.,1] :< .)) // drop the missing rows
	
	// coordinates of the interior points (vorcenters of triangles in voredges)
	// each point of a triangle is associated with at most two points
	
	point0 = J(rows(triangles),2,.)
	point1 = J(rows(triangles),2,.)

			for (i=1; i <= rows(triangles); i++) {
				if (i < halfedges[i]) {
					point0[i,.] = vorcenter[triangleOfEdge(i),.]
					point1[i,.] = vorcenter[triangleOfEdge(halfedges[i]),.]
				}
			}

			
	point0 = select(point0, (point0[.,1] :< .)) // drop the missing rows
	point1 = select(point1, (point1[.,1] :< .)) 		


	// exterior cell rays
	
	hlen = hull[rows(hull), 1]
	
	p0 = hlen * 4
	p1 = hlen * 4
	
	x0 = coords[2 * hlen - 1, 1]
	x1 = coords[2 * hlen - 1, 1]
	
	y0 = coords[2 * hlen, 1]
	y1 = coords[2 * hlen, 1]
	
	vectors = J(rows(triangles) * 2, 1, 0)
	
		for (i=1; i<=rows(hull); i++) {			
			
			hme = hull[i,1]
			
			p0 = p1
			x0 = x1
			y0 = y1
			
			p1 = hme * 4
			
			x1 = coords[2 * hme - 1, 1]		
			y1 = coords[2 * hme    , 1]
			
			vectors[p0 + 2, 1] = y0 - y1
			vectors[p1 - 0, 1] = y0 - y1
			
			vectors[p0 + 3, 1] = x1 - x0
			vectors[p1 + 1, 1] = x1 - x0			
		}	
		
	// add the boundary edges
	pointh0 = J(rows(hull),2,.)
	pointh1 = J(rows(hull),2,.)
			
	h0 = hull[rows(hull),1]
	h1 = hull[rows(hull),1]

		for (i=1; i <= rows(hull); i++) {
			
			h0 = h1
			h1 = hull[i,1]
			t = findhulltri(h0, h1, tri3)
			
			vorc  = vorcenter[t,.] 
			
			v = h0 * 4
			p = project(vorc[1,1], vorc[1,2], vectors[v + 2, 1], vectors[v + 3, 1], xmin, xmax, ymin, ymax)
			
			if (p[1,1] != 0) {
			pointh0[i,.] = vorcenter[t,.]
			pointh1[i,.] = p
			}
		}
		

	// append with point type
	point0 =  point0 \ pointh0
	point1 =  (point1, J(rows(point1),1,1)) \ (pointh1, J(rows(pointh1),1,2))

	
	// clip the rays

	vorl = J(rows(point0),4,.)
		
		for (i=1; i <= rows(point0); i++) {			
			vorl[i,.] = clipline(point0[i,1], point0[i,2], point1[i,1], point1[i,2], xmin, xmax, ymin, ymax)	
		}

	
		for (i=1; i <= rows(vorl); i++) {	
			
			if (vorl[i,1]==vorl[i,3] & vorl[i,2]==vorl[i,4]) {
					
				vorl[i,1] = .
				vorl[i,2] = .
				vorl[i,3] = .
				vorl[i,4] = .
			}			
		}	

		
	xminr = .
	xmaxr = .
	yminr = .
	ymaxr = .
	
	bounds(points, xminr, xmaxr, yminr, ymaxr, offs) 
	
	
	// rescale back from [0,1] to the original domain
	vorl[.,1] = rescale2(vorl[.,1], xmin, xmax, xminr, xmaxr)
	vorl[.,2] = rescale2(vorl[.,2], ymin, ymax, yminr, ymaxr)
	vorl[.,3] = rescale2(vorl[.,3], xmin, xmax, xminr, xmaxr)
	vorl[.,4] = rescale2(vorl[.,4], ymin, ymax, yminr, ymaxr)			
		
	
	
	////////////////////////////////
	/// polygon procedures here  ///
	////////////////////////////////
	

	// interior points
	point0c = vorl[1::rows(voredges),1..2]
	point1c = vorl[1::rows(voredges),3..4]

	
	vorpoly = J(rows(voredges) * 4, 4, .)
	

	for (i=1; i<=rows(voredges); i++) {
	
		vorpoly[i * 4 - 3, 1] = CommonVals(tri3[.,voredges[i,1]], tri3[.,voredges[i,2]])[1,1]
		vorpoly[i * 4 - 2, 1] = CommonVals(tri3[.,voredges[i,1]], tri3[.,voredges[i,2]])[1,1]
		vorpoly[i * 4 - 1, 1] = CommonVals(tri3[.,voredges[i,1]], tri3[.,voredges[i,2]])[2,1]		
		vorpoly[i * 4    , 1] = CommonVals(tri3[.,voredges[i,1]], tri3[.,voredges[i,2]])[2,1]

		vorpoly[i * 4 - 3, 2..3] = point0c[i,.]  // line start point A 
		vorpoly[i * 4 - 2, 2..3] = point1c[i,.]  // line end   point A
		
		vorpoly[i * 4 - 1, 2..3] = point0c[i,.]  // line start point B
		vorpoly[i * 4    , 2..3] = point1c[i,.]  // line start point B		
		
	}
	

	vorpoly = uniqrows(vorpoly)
	vorpoly= sort(vorpoly, (1,2))
	
	// mark the hull points
	for (i = 1; i <= rows(vorpoly); i++ ) { 
		for (j = 1; j <= rows(hull); j++ ) { 
			if (vorpoly[i,1]==hull[j,1]) vorpoly[i,4]=1
		}
	}
	
	
// ray points

	vorpoly2 = J(rows(hull) * 4, 3, .)
	
	for (i=1; i<=rows(hull); i++) {
		
		i1 = i
		
		if (i1 > 1) {
			i2 = i1 - 1
		}
		else {
			i2 = rows(hull)
		}
		
		vorpoly2[i * 4 - 3, 1] = hull[i1, 1]
		vorpoly2[i * 4 - 2, 1] = hull[i1, 1]
		vorpoly2[i * 4 - 1, 1] = hull[i2, 1]
		vorpoly2[i * 4    , 1] = hull[i2, 1]

		
		vorpoly2[i * 4 - 3, 2..3] = pointh0[i,.]
		vorpoly2[i * 4 - 2, 2..3] = pointh1[i,.]
		
		vorpoly2[i * 4 - 1, 2..3] = pointh0[i,.]
		vorpoly2[i * 4    , 2..3] = pointh1[i,.]		
		
	}
	
	vorpoly2 = select(vorpoly2, (vorpoly2[.,2] :< .)) // drop the missing rows	
	
	// rescale them back
	vorpoly2[.,2] = rescale2(vorpoly2[.,2], xmin, xmax, xminr, xmaxr)
	vorpoly2[.,3] = rescale2(vorpoly2[.,3], ymin, ymax, yminr, ymaxr)
	
	vorpoly2 = vorpoly2, J(rows(vorpoly2),1,1)
	

	vorp = vorpoly \ vorpoly2
	vorp = uniqrows(vorp)
	vorp = vorp, J(rows(vorp),1,.)
	vorp = select(vorp, (vorp[.,2] :< .)) // drop the missing rows	
}

end





////////////////////////////////////
///   voronoi subroutines here   ///
////////////////////////////////////

********************
// 	   bounds	  //		
********************

cap mata: mata drop bounds()

mata // bounds
function bounds(points,xmin,xmax,ymin,ymax, offs)
{

	// 5% displacment based on (max - min)
	displacex = abs((max(points[.,1]) - min(points[.,1])) * offs)
	displacey = abs((max(points[.,2]) - min(points[.,2])) * offs)

	xmin 	  = min(points[.,1]) - displacex
	xmax 	  = max(points[.,1]) + displacex

	ymin 	  = min(points[.,2]) - displacey
	ymax 	  = max(points[.,2]) + displacey

	st_numscalar("xmin", xmin)
	st_numscalar("xmax", xmax)

	st_numscalar("ymin", ymin)
	st_numscalar("ymax", ymax)
	
}
end


************************
// 	 triangleOfEdge   //     // index
************************

cap mata: mata drop triangleOfEdge()
mata: // triangleOfEdge
function triangleOfEdge(x)
	{
		return(floor((x - 1) / 3) + 1)
	}
end


************************
// 	 edgesOfTriangle  //    
************************

cap mata: mata drop edgesOfTriangle()
mata: // edgesOfTriangle
function edgesOfTriangle(t)
	{
		return (3 * t - 2, 3 * t - 1, 3 * t)
	}
end

***************************
// 	 pointsOfTriangle    // 
***************************

cap mata: mata drop pointsOfTriangle()
mata:  // pointsOfTriangle
function pointsOfTriangle(triangles,t)
	{
		return (triangles[edgesOfTriangle(t)[1,1],1],triangles[edgesOfTriangle(t)[1,2],1],triangles[edgesOfTriangle(t)[1,3],1] )
	}
end



************************
// 	  circumcenter2	  //	 // another circumcenter. it just returns a different structure
************************    

cap mata: mata drop circumcenter2()

mata: // circumcenter2
function circumcenter2(points,triangles,i)
{
	real matrix myset
	real scalar t1, t2, t3
	
	myset = points[pointsOfTriangle(triangles,i),.]
		
	x1 = myset[1,1]
	y1 = myset[1,2]
	x2 = myset[2,1]
	y2 = myset[2,2]
	x3 = myset[3,1]
	y3 = myset[3,2]
	
	
	dx = x2 - x1
	dy = y2 - y1
	ex = x3 - x1
	ey = y3 - y1
	ab = (dx * ey - dy * ex) * 2
	
	
	// if (abs(ab) < 1e-9 ) {  // indexing issue
		// degenerate case 
	//	a = 1e9
	//	r = triangles[1] * 2
	//	a = a * sign((points[r - 1,1] - x1) * ey - (points[r,1] - y1) * ex)
	//	myx = (x1 + x3) / 2 - a * ey
	//	myy = (y1 + y3) / 2 + a * ex
	// }
	// else {	
	    d = 1 / ab
        bl = dx * dx + dy * dy
        cl = ex * ex + ey * ey
        myx = x1 + (ey * bl - dy * cl) * d
        myy = y1 + (dx * cl - ex * bl) * d
	// }
	
	return (myx, myy)
}
end



**************************
// 	 Voronoi edges      // the triangle center pairs that need to be connected
**************************


cap mata: mata drop forEachVoronoiEdge()

mata: // forEachVoronoiEdge
function forEachVoronoiEdge(triangles,halfedges, real matrix xx) 
	{
		for (i=1; i <= rows(triangles); i++) {
			if (i < halfedges[i]) {
				
				xx[i,.] = triangleOfEdge(i), triangleOfEdge(halfedges[i])

			}
		}
			
	}
end

************************
// 	  _project		  //
************************    

cap mata: mata drop project()

mata: // project
function project(x0, y0, vx, vy, xmin, xmax, ymin, ymax)
{
	real scalar t, c
	t = .   // just a very large number
	
	if (vy < 0) { // top
		
		if (y0 <= ymin) return(0,0)
		
		
		c = (ymin - y0) / vy
		if (c < t) {
			myy = ymin
			myx = x0 + (t = c) * vx
		}
    } 

	else { // bottom
	
		if (y0 >= ymax) return(0,0)
		
		
		c = (ymax - y0) / vy
		if (c < t) {
			myy = ymax
			
			myx = x0 + (t = c) * vx
		}
	}
    
	if (vx > 0) { // right
		
		if (x0 >= xmax) return(0,0)
		
		c = (xmax - x0) / vx
		if (c < t) {
			myx = xmax
			myy = y0 + (t = c) * vy
		}
    } 
	
	else { // left
		if (x0 <= xmin) return(0,0)
		
		c = (xmin - x0) / vx
		if (c < t) {
			myx = xmin
			myy = y0 + (t = c) * vy
		}
    }
	
	return (myx, myy)
	
}
end



************************
// 	  find hull tri	  //	
************************

cap mata: mata drop findhulltri()

mata: //  findhulltri
function findhulltri(x,y,tri3)
{	
	for (i=1; i <= cols(tri3); i++) {	
		if (max(x:==tri3[.,i])==1 & max(y:==tri3[.,i])==1)	return(i)
	}
}
end



*********************
// 	initialize	   //
*********************

cap mata: mata drop initialize()

mata: // initialize
real vector initialize(real matrix data)
{
	real scalar num, i
	real vector coords
	
	num  = rows(data)
	
	coords = J(num*2,1,.)
	
		for (i=1; i<=num; i++) {	
			coords[2 * i - 1, 1] = data[i,1]
			coords[2 * i    , 1] = data[i,2]
		}
		
	return(coords)
}
end

********************
// 	 rescale2	  // rescale a vector column, based on rescaled min/max, and original minmax	
********************

cap mata: mata drop rescale2()

mata: // rescale2
real vector rescale2(points, mymin, mymax, a, b)
{
	newpoints = (b - a) * (points :- mymin) :/ (mymax - mymin) :+ a
	
	return(newpoints)
}
end
			

*********************
// 	   clipline    //
*********************




cap mata: mata drop clipline()
mata:  // clipline
	function clipline(x1, y1, x2, y2, minX, maxX, minY, maxY)  
	{
		real scalar code1, code2, accept, x, y

		// Defining region codes 
		LEFT   = 1  // 0001 
		RIGHT  = 2  // 0010 
		BOTTOM = 4  // 0100 
		TOP    = 8  // 1000 
		
		
		code1 = computeCode(x1, y1, minX, maxX, minY, maxY)
		code2 = computeCode(x2, y2, minX, maxX, minY, maxY)
		
		accept = 0
		
		t = 0 // counter
		
		zz = 0
		while (zz == 0)  {
		
			// If both endpoints lie within rectangle 
			if (code1 == 0 & code2 == 0) {
				accept = 1
				break
			}
			// If both endpoints are outside rectangle 
			else if (code1==code2 & code1!=0 & code2!=0) {
				accept = 0
				
				x1 = .
				y1 = .
				x2 = .
				y2 = .
				
				break
			}
			
			// Some segment lies within the rectangle 			
			else {	
				x = 1
				y = 1
				
				if (code1 != 0) {
					codeout = code1
				}
				else {
					codeout = code2
				}
				
				// point is above the clip rectangle 
				if (codeout >= TOP) {   					// fix this line. what is top?
					x = x1 + ((x2 - x1) * (maxY - y1) / (y2 - y1))
					y = maxY
				}
				// point is below the clip rectangle
				else if (codeout >= BOTTOM) {
					x = x1 + ((x2 - x1) * (minY - y1) / (y2 - y1))
					y = minY
				}
				// point is to the right of the clip rectangle 
				else if (codeout >= RIGHT) {
					y = y1 + ((y2 - y1) * (maxX - x1) / (x2 - x1))
					x = maxX
				}
				// point is to the left of the clip rectangle 
				else if (codeout >= LEFT) {
					y = y1 + ((y2 - y1) * (minX - x1) / (x2 - x1))
					x = minX
				}	
				
			if (codeout == code1) { 
                x1 = x 
                y1 = y 
                code1 = computeCode(x1, y1, minX, maxX, minY, maxY) 
			}
			else {
                x2 = x 
                y2 = y 
                code2 = computeCode(x2, y2, minX, maxX, minY, maxY) 
			}
		
		t = t + 1
		if (t > 100) {
		//	printf("bad egg\n")
				x1 = .
				y1 = .
				x2 = .
				y2 = .
			break
			}
		
		}

	}
	
			
	return(x1,y1,x2,y2)
}
end




*********************
// 	 computeCode   //
*********************


cap mata: mata drop computeCode()
mata // computeCode
	function computeCode(xx, yy, minX, maxX, minY, maxY)  // triangleCenter,
	{
		real scalar code
		code = 0   // 1 = left, 2 = right, 4 = bottom, 8 = top (defined in binary)
			if (xx < minX) code = code + 1
			if (xx > maxX) code = code + 2		
			if (yy < minY) code = code + 4
			if (yy > maxY) code = code + 8	
		return(code)
	}
end



**** END OF CLIPLINE ****

*******************
// 	 CommonVals  //    // comapre two vectors and find their common values
*******************

** form this statalist: https://www.statalist.org/forums/forum/general-stata-discussion/mata/177560-common-elements-in-a-list

cap mata: mata drop CommonVals()

mata: // CommonVals
real colvector CommonVals(real colvector Z, real colvector Y) {
    
	real vector mylist
	real scalar i
	
	mylist = J(0,1,.)
    for (i=1; i<=rows(Z); i++) {
		if (any(Y:==Z[i])) {
			mylist=mylist \ Z[i]
		}
	}
    return(mylist)
}
end



**************************
// convert to polygons  //
**************************

cap program drop vorpoly
program define vorpoly, sortpreserve 
	version 15

	syntax varlist(min = 2 max = 2 numeric)
	
	gettoken y x : varlist
	
	
	////////////////////////
	// save the raw data  //
	////////////////////////
	
	qui {
		preserve
			tempfile _raw
			keep _ID `varlist'
			drop if _ID==.
			drop if `x'==.   
			ren _ID vpoly_id
			save `_raw', replace
		restore
	}
	
	
	///////////////////////////
	// clean up the polygons //
	///////////////////////////
	
	qui {
		preserve	
			clear
		
			getmata (vpoly_id vpoly_x vpoly_y hull angle) = vorpolygons, force double replace

			drop angle
			bysort vpoly_id: gen order = _n
			
			// add the corner here
			cap drop dist*

			bysort vpoly_id: egen double dist_L = min(vpoly_x - xmin)    if hull==1    // left-most   point on x-axis
			bysort vpoly_id: egen double dist_R = min(xmax    - vpoly_x) if hull==1    // right most  point on x-axis
			bysort vpoly_id: egen double dist_B = min(vpoly_y - ymin)	 if hull==1    // bottom-most point on y-axis
			bysort vpoly_id: egen double dist_T = min(ymax    - vpoly_y) if hull==1    // top-most    point on y-axis

			cap drop edge*

			egen edge_L = group(dist_L)
			egen edge_R = group(dist_R)
			egen edge_B = group(dist_B)
			egen edge_T = group(dist_T)


			recode edge* (2/. = .)	
			
			cap drop corner*
			gen corner=.

			replace corner = 1 if edge_B ==1 & edge_R ==1   // bottom right
			replace corner = 2 if edge_B ==1 & edge_L ==1   // bottom left
			replace corner = 3 if edge_T ==1 & edge_R ==1   // top right
			replace corner = 4 if edge_T ==1 & edge_L ==1   // top left

			tab corner
			
			expand 2 if corner!=. & order==1, gen(markme)

			replace vpoly_x = xmax if corner==1 & markme==1
			replace vpoly_y = ymin if corner==1 & markme==1

			replace vpoly_x = xmin if corner==2 & markme==1
			replace vpoly_y = ymin if corner==2 & markme==1

			replace vpoly_x = xmax if corner==3 & markme==1
			replace vpoly_y = ymax if corner==3 & markme==1

			replace vpoly_x = xmin if corner==4 & markme==1
			replace vpoly_y = ymax if corner==4 & markme==1

			summ order, meanonly
			replace order = r(max) + 1 if markme==1

			sort vpoly_id order
			drop order edge* dist* corner markme hull

			merge m:1 vpoly_id using `_raw'

			gen angle = atan2(vpoly_y - `y', vpoly_x - `x')  // fix the order of drawing the points
			gsort vpoly_id -angle

			by vpoly_id: gen order = _n		
			
			egen tag = tag(vpoly_id)
			expand 2 if tag==1, gen(temp)
			replace order=0 if temp==1
			replace vpoly_x = . if temp==1
			replace vpoly_y = . if temp==1
			
			drop tag temp `varlist' _m angle

			sort vpoly_id order
			ren order shape_order
			ren vpoly_id _ID
			ren vpoly_x _X
			ren vpoly_y _Y

			*mkmat _ID _X _Y shape_order, mat(vorpoly)
			compress
			save _vorpoly.dta, replace	
			
		restore
	}
	
end


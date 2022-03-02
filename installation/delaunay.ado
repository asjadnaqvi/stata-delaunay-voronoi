*! Delaunay triangulation
*! by Asjad Naqvi (asjadnaqvi@gmail.com, @AsjadNaqvi)
*!
*! Ver 1.10 01.03.2022: missing triangles fixed. rays fixed. 
*  Ver 1.02 27.12.2021: rescale added. if made optional
*  Ver 1.01 20.12.2021: if/in & notes by wbuchanan
*  Ver 1.00 22.11.2021: first release


**************************
**************************
**                      **
**       Delaunay       **
**						**
**     Asjad Naqvi      **
**                      ** 
**    Last updated:     **
**     01 Mar 2022      **
**			            ** 
**    First release:    **
**     22 Nov 2021      **
**			            ** 
**************************
**************************


cap prog drop delaunay

prog def delaunay, eclass sortpreserve

	version 15
	
	syntax varlist(min = 2 max = 2 numeric) [if] [in],	 ///   
		[ REScale TRIangles Hull VORonoi(namelist min=1 max=2) ] [OFFset(real 0.05)] 
	

	di "Delaunay: Initializing"
	
	gettoken x y : varlist
	
		
	// generate an internal _id variable
		
		cap confirm variable _id
		
		if !_rc {
			drop if _id==.
			drop _id
		}

		gen _id = _n
		lab var _id "observation id"
	
	marksample touse, strok
	
	

	mata: points   = select(st_data(., ("`x'", "`y'")), st_data(., "`touse'"))
	
	if "`rescale'" != "" {	
		mata: points2 = points  // create a copy
		mata: points2[.,1] = rescale(points[.,1], 0, 1)
		mata: points2[.,2] = rescale(points[.,2], 0, 1)
	}
	else {
		mata: points2 = points
	}
	
	mata: eps          = 1e-20
	mata: edgestack    = J(512, 1, .) 
	mata: coords       = initialize(points2)
	mata: num          = rows(points2)   // updated
	mata: maxTriangles = max(((2 * num) - 5, 0)) 
	
	// core arrays
	mata: triangles = J(maxTriangles * 3, 1, .)   
	mata: halfedges = J(maxTriangles * 3, 1, .)
	mata: hull      = J(num, 1, .)  
		
	// arrays for tracking 
	mata: hashSize  = ceil(sqrt(num))
	mata: hullPrev  = J(num, 1, .)  
	mata: hullNext  = J(num, 1, .)  
	mata: hullTri   = J(num, 1, .)  
	mata: hullHash  = J(hashSize, 1, -1) 
	mata: ids   	= J(num, 1, .)
	mata: dists     = J(num, 1, .)	
	
	// temporary arrays for sorting points
	mata: ids   = J(num, 1, .)
	mata: dists = J(num, 1, .)
	
	
	di "Starting core routines"

	// run the core routine 

	mata: _delaunay_core(coords, ids, dists, triangles, halfedges, hull, hullNext,  ///   
				  hullPrev, hullTri, hullHash, hashSize, eps, edgestack)

	**************************
	***  export geometry   ***
	**************************

	di "Exporting geometry"
	
	
	// delaunay triangles
	if "`triangles'" != ""  add_triangles
	
	// convex hull
	if "`hull'" 	 != ""  add_hull
	
	// voronoi lines and/or polygons
	
	// if "`voronoi'" 	 == "" di as err `" You need to specify {ul:lines} or {ul:polygon} or both as options."'
	
	local vor_valid = `" "lines", "polygons" "'
	capt assert inlist( "`voronoi'", `vor_valid')

	/*
		di "`_rc'"
		
		if _rc {
		
		di as err `" Valid options for voronoi() are: `vor_valid' or both. "'
		exit
		}
	*/
	
	
	if "`voronoi'" 	!= ""  {
		if "`offset'" != "" {
			voronoi `varlist', `voronoi' offset(`offset')
		}
		else {
			voronoi `varlist', `voronoi'
		}
	}
	
	di "Done!"
		
	
end



************************
// 	 _delaunay_core	  //  main routine
************************

cap mata: mata drop _delaunay_core()

mata: // _delaunay_core
void _delaunay_core(coords,ids,dists,triangles,halfedges,hull,hullNext,hullPrev,hullTri,hullHash,hashSize,eps,edgestack)
{
	num = rows(coords) / 2
	
	// populate an array of point indices; calculate input data bbox
	
	minX =  1e16
	minY =  1e16
	maxX = -1e16   // arbitrary large numbers representing infinity
	maxY = -1e16

	for (i=1; i<= num; i++) {	
		
		x = coords[2 * i - 1, 1] 
		y = coords[2 * i    , 1]
				
		if (x < minX) minX = x
		if (y < minY) minY = y
		if (x > maxX) maxX = x
		if (y > maxY) maxY = y
		
		ids[i,1] = i
	}

	cx = (minX + maxX) / 2
	cy = (minY + maxY) / 2	
		
	minDist = .
	i0 = 0
	i1 = 0
	i2 = 0
	
	// pick a seed point close to the center (first observation)

	for (i=1; i<=num; i++) {	
		dval = dist(cx, cy, coords[2 * i - 1], coords[2 * i])
	
		if (dval < minDist) {
			i0 = i
			minDist = dval
		}
	}

	i0x = coords[2 * i0 - 1, 1]
	i0y = coords[2 * i0    , 1]
	minDist = .


	// find the point closest to the seed (second observation)
	for (i=1; i<=num; i++) {
		if (i == i0) continue
		
		dval = dist(i0x, i0y, coords[2 * i - 1], coords[2 * i])

		if ((dval < minDist) & (dval > 0)) {
			i1 = i
			minDist = dval
		}
	}
	    
	i1x = coords[2 * i1 - 1, 1]
	i1y = coords[2 * i1    , 1]
	minRadius = .

	// find the third point which forms the smallest circumcircle with the first two (third observation)
	
	for (i=1; i<=num; i++) {
		if (i == i0 | i == i1) continue
		
		r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i - 1, 1], coords[2 * i, 1])
		
		if (r < minRadius) {
			i2 = i
			minRadius = r
		}
	}
		
	i2x = coords[2 * i2 - 1, 1]
	i2y = coords[2 * i2    , 1]		

	if (minRadius == .) {  // a very large number
		
		// calculate the distances ot the points with the first point	
		for (i=1; i<=num; i++) {
			dists[i,1] = (coords[2 * i - 1, 1] - coords[1, 1]) | (coords[2 * i, 1] - coords[2, 1])  
		}
		
		// and order them
		_quicksort(ids, dists, 1, num - 1)  
		
		hull =  J(num, 1, .)
		j = 1
		d0 = -1e20  // minus infinity

		for (i=1; i <= num; i++) {   
			id = ids[i,1]

			if (dists[id,1] > d0) {
				hull[j, 1] = id
				j = j + 1
				d0 = dists[id,1]
			}
		}
	
		// initiatte the hull, triangle, and halfedges matrices 
		hull = hull[1..j,1]  
		triangles =  .  	
		halfedges =  .	
	}


	// swap the order of the seed points for counter-clockwise orientation

	if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
		i = i1
		x = i1x
		y = i1y
		
		i1  = i2
		i1x = i2x
		i1y = i2y
		
		i2  = i
		i2x = x
		i2y = y		
	}

	c0x = .   // blanks to be evaluated in circumcenter below
	c0y = .
		
	circumcenter(i0x, i0y, i1x, i1y, i2x, i2y, c0x, c0y)

		
	for (i=1; i <= num; i++) {   // changed from <= to <
		dists[i,1] = dist(coords[2 * i - 1,1], coords[2 * i,1], c0x, c0y)
	}
	
	// sort the points by distance from the seed triangle circumcenter
	_quicksort(ids, dists, 1, num-1)  // num changed to num - 1

	
	// set up the seed triangle as the starting hull
	hullStart = i0
	hullSize  = 3		
		
	hullNext[i0] = i1
	hullPrev[i2] = i1	
	
	hullNext[i1] = i2 
	hullPrev[i0] = i2		
	
	hullNext[i2] = i0 
	hullPrev[i1] = i0
	

	// initialize the hull triangles
	hullTri[i0] = 1
	hullTri[i1] = 2
	hullTri[i2] = 3

		
	// manual index fix to match the py output
	hullHash[hashKey(i0x,i0y,cx,cy,hashSize) + 1] = i0   
	hullHash[hashKey(i1x,i1y,cx,cy,hashSize) + 1] = i1
	hullHash[hashKey(i2x,i2y,cx,cy,hashSize) + 1] = i2
		
	trianglesLen = 1
	

	// dump is just a dummy to prevent the function from returning a value on the screen
	dump = addTriangle(triangles,trianglesLen,halfedges, i0, i1, i2, -1, -1, -1)   

	xp=.
	yp=.

	
	///////////////////////////////////
	///       main loop below       ///
	///////////////////////////////////
	
	for (k=1; k <= rows(ids); k++) {   // evalaute the sorted points
		
		i = ids[k,1]
		x = coords[2 * i - 1, 1]
		y = coords[2 * i	, 1]
			
		// skip near-duplicate points 
		if ((k > 1) & (abs(x - xp) <= eps) & (abs(y - yp) <= eps)) continue
		
		xp = x
		yp = y
		
		
		// skip seed triangle points
		if ((i == i0) | (i == i1) | (i == i2)) continue
		
		
		
		// find a visible edge on the convex hull using edge hash
		start = 1
		key = hashKey(x, y, c0x, c0y, hashSize)

		for (j=1; j <= hashSize; j++) {			
			start = hullHash[(mod((key + j - 1), hashSize) + 1),1] 				

			if ((start != -1) & (start != hullNext[max((start,1)), 1])) break  // added max to prevent the 0 argument
		}
		

		start = hullPrev[start,1]
		emp = start                		
		
		www = 0 // while		
		
		
		while (www == 0) {  // here is where points are being dropped
			q = hullNext[emp,1]
			if (orient(x, y, coords[2 * emp - 1], coords[2 * emp], coords[2 * q - 1], coords[2 * q])) break		
			emp = q
			
			if (emp == start) {
				emp = -1
				break
			}
		}


		if (emp == -1) continue	
		
		
		// add the first triangle from the point
		
		t = addTriangle(triangles, trianglesLen, halfedges, emp, i, hullNext[emp,1], -1, -1, hullTri[emp,1])
		
		// recursively flip triangles from the point until they satisfy the Delaunay condition
		hullTri[i  , 1] = legalize(t + 2, coords, halfedges, edgestack, triangles, hullStart, hullTri, hullPrev)			
		hullTri[emp, 1] = t 
		hullSize = hullSize + 1	
		
		// walk forward through the hull, adding more triangles and flipping recursively
		num = hullNext[emp, 1]
	
		www = 0 // for the infinite while
		while (www == 0) {
			
			q  = hullNext[num, 1]	

			if (orient(x, y, coords[2 * num - 1,1], coords[2 * num,1], coords[2 * q - 1,1], coords[2 * q,1]) == 0) break					
			t = addTriangle(triangles,trianglesLen, halfedges,num, i, q, hullTri[i,1], -1, hullTri[num,1])	
					
			hullTri[i,1] = legalize(t + 2, coords,halfedges,edgestack,triangles,hullStart,hullTri,hullPrev)
			hullNext[num,1] = num // mark as removed		
			
			hullSize = hullSize - 1
			num = q		
		}
	
		// walk backward from the other side, adding more triangles and flipping
		if (emp == start) {
			www = 0 // while
			while (www == 0) {
				q = hullPrev[emp,1]
				
				if ((orient(x, y, coords[2 * q - 1], coords[2 * q], coords[2 * emp - 1], coords[2 * emp])) == 0) break
				t = addTriangle(triangles, trianglesLen, halfedges, q, i, emp, -1, hullTri[emp,1], hullTri[q,1])
				
				dump = legalize(t + 2, coords, halfedges, edgestack, triangles, hullStart, hullTri, hullPrev)
				
				hullTri[q,1] = t
				hullNext[emp,1] = emp // mark as removed
				hullSize = hullSize - 1
				emp = q	
				
			}
		}		
		
		
		// update the hull indices
		hullStart 			= emp
		hullPrev[i,1] 		= emp
		hullNext[emp,1] 	= i 
		hullPrev[num,1] 	= i
		hullNext[i,1]   	= num

		// save the two new edges in the hash table
		hullHash[hashKey(x,y,cx,cy,hashSize) + 1,1] = i				
		hullHash[hashKey(coords[2 * emp - 1,1], coords[2 * emp,1],cx,cy,hashSize) + 1,1] = emp		

	}
	
	
	hull = J(hullSize,1,.)
	emp = hullStart
        
	for (i=1; i<=hullSize; i++) {
		hull[i,1] = emp
		emp = hullNext[emp,1]
	}		
}
end


****************************************************
****************  SUBROUTINES   ********************
****************************************************

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
// 	   _hashKey	  //		
********************

cap mata: mata drop hashKey()

mata // _hashKey
real scalar hashKey(x,y,cx,cy,hashSize)
{
	return(mod(floor(pseudoAngle(x - cx, y - cy) * hashSize),hashSize))
}
end



************************
// 	   pseudoAngle	  //
************************

cap mata: mata drop pseudoAngle()

mata: // pseudoAngle
function pseudoAngle(real scalar dx, real scalar dy)
{

	real scalar p
	p = dx / (abs(dx) + abs(dy))

	if (dy > 0) {
        return((3 - p) / 4)
	}
    else {
        return((1 + p) / 4)
	}
}
end



************************
// 	   _legalize	  //		
************************

cap mata: mata drop legalize()

mata // _legalize
function legalize(real scalar a, real vector coords, halfedges, edgestack, triangles, hullStart, hullTri, hullPrev)
{
	i = 1
	ar = 0
	
	www = 0		
	while (www==0) {
		
		b = halfedges[a,1]
	
            /*
              if the pair of triangles doesn't satisfy the Delaunay condition
              (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
              then do the same check/flip recursively for the new pair of triangles
             
                        pl                    pl
                       /||\                  /  \
                    al/ || \bl            al/    \a
                     /  ||  \              /      \
                    /  a||b  \    flip    /___ar___\
                  p0\   ||   /p1   =>   p0\---bl---/p1
                     \  ||  /              \      /
                    ar\ || /br             b\    /br
                       \||/                  \  /
                        pr                    pr
             
            */
	
		a0 =  a - mod(a - 1, 3)  	// mod(a,3)
		ar = a0 + mod(a + 1, 3)  	//+ mod((a + 2),3)
	
		if (b == -1) { //  convex hull edge
			if (i == 1) break
			i = i - 1
			a = edgestack[i,1]  // define in the other file. empty array
			continue
		}

		b0 = b  - mod(b - 1,3)
		al = a0 + mod(a, 3)
		bl = b0 + mod((b + 1),3)
			
		
		p0 = triangles[ar,1]
		pr = triangles[a ,1]
		pl = triangles[al,1]
		p1 = triangles[bl,1]
		
		
		illegal = inCircle(coords[2 * p0 - 1, 1], coords[2 * p0, 1], coords[2 * pr - 1], coords[2 * pr, 1], coords[2 * pl - 1, 1], coords[2 * pl, 1], coords[2 * p1 - 1, 1], coords[2 * p1, 1])
		
		
		
		if (illegal) {
		
			triangles[a,1] = p1
			triangles[b,1] = p0
			
			hbl = halfedges[bl, 1]		
			
			// edge swapped on the other side of the hull (rare); fix the halfedge reference

			if (hbl == -1) {
				eme = hullStart  
				
				yyy = 0
				while (yyy==0) {
					if (hullTri[eme, 1] == bl) {
						hullTri[eme, 1] = a
						break
					}
					
					eme = hullPrev[eme, 1]
					if (eme == hullStart) break
				}
			}
			
			_link(halfedges,a, hbl)
			_link(halfedges,b, halfedges[ar, 1])
			_link(halfedges,ar, bl)
		
			brme = b0 + mod(b,3) 
		
			if (i < rows(edgestack)) {
				edgestack[i,1] = brme
				i = i + 1
			}
		}	
		else {
			if (i == 1) break
			i = i - 1		
			a = edgestack[i,1]
		}
	}
	return(ar) 
}
end



********************
// 	   _link	  //		
********************

cap mata: mata drop _link()

mata: // _link
void _link(real vector halfedges, real scalar a, real scalar b)
{
	halfedges[a,1] = b
	if (b != -1) {
		halfedges[b,1] = a
	}
}
end


************************
// 	   _addTriangle	  //		
************************

cap mata: mata drop addTriangle()

mata: // addTriangle
function addTriangle(triangles, trianglesLen, halfedges, i0, i1, i2, a, b, c)
{
        real scalar t
		t = trianglesLen

        triangles[t    ] = i0
        triangles[t + 1] = i1
        triangles[t + 2] = i2

		_link(halfedges, t    , a)
        _link(halfedges, t + 1, b)
        _link(halfedges, t + 2, c)

        trianglesLen = trianglesLen + 3
        return(t)
}
end


****************
// 	   dist	  //
****************

cap mata: mata drop dist()

mata:  // dist
real vector dist(ax, ay, bx, by)
{
	real scalar dx, dy
	dx = ax - bx
	dy = ay - by

	return(dx*dx + dy*dy)
}
end


************************
// 	   orientIfSure	  //
************************

cap mata: mata drop orientIfSure()

mata: // orientIfSure
function orientIfSure(px, py, rx, ry, qx, qy)
{
	real scalar left, right
	
	left  = (ry - py) * (qx - px)
    right = (rx - px) * (qy - py)

    
	// J. Shewchuk error bound check
	if (abs(left - right) >= (3.3306690738754716e-16 * abs(left + right))) {
        return(left - right)
	}
    else {
        return(0)
	}
}
end


********************
// 	   orient	  //		
********************

cap mata: mata drop orient()

mata: // orient
function orient(px, py, rx, ry, qx, qy)
{
	real scalar x
	x = (orientIfSure(px, py, rx, ry, qx, qy) < 0  | orientIfSure(rx, ry, qx, qy, px, py) < 0  |  orientIfSure(qx, qy, px, py, rx, ry)  < 0 )
	
	return(x)  // a true false statement
}
end



********************
// 	   inCircle	  //		
********************

cap mata: mata drop inCircle()

mata: // inCircle
function inCircle(ax, ay, bx, by, cx, cy, px, py)
{
    real scalar dx, dy, ex, ey, fx, fy, ap, bp, cp, x
	
	dx = ax - px
    dy = ay - py
    ex = bx - px
    ey = by - py
    fx = cx - px
    fy = cy - py

    ap = (dx * dx) + (dy * dy)
    bp = (ex * ex) + (ey * ey)
    cp = (fx * fx) + (fy * fy)

	x = (dx * (ey * cp - bp * fy)) - (dy * (ex * cp - bp * fx)) + (ap * (ex * fy - ey * fx))

	return(x < 0) 

}
end


************************
// 	   circumradius	  //		
************************

cap mata: mata drop circumradius()

mata: // circumradius
function circumradius(ax, ay, bx, by, cx, cy)
{
    real scalar dx, dy, ex, ey, bl, cl, d, x, y
	
	dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = (dx * dx) + (dy * dy)
    cl = (ex * ex) + (ey * ey)
    d = 0.5 / ((dx * ey) - (dy * ex))

    x = ((ey * bl) - (dy * cl)) * d
    y = ((dx * cl) - (ex * bl)) * d

    return(x*x + y*y)
}
end


************************
// 	   circumcenter	  //		
************************

cap mata: mata drop circumcenter()

mata: // circumcenter
function circumcenter(ax, ay, bx, by, cx, cy, real scalar retx, real scalar rety)
{
	real scalar dx, dy, ex, ey, bl, cl, d
	
	dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = (dx * dx) + (dy * dy)
    cl = (ex * ex) + (ey * ey)
    d = 0.5 / ((dx * ey) - (dy * ex))

    retx = ax + ((ey * bl) - (dy * cl)) * d  // return x
    rety = ay + ((dx * cl) - (ex * bl)) * d  // return y
	return
}
end


************************
// 	   quicksort   	  //   
************************

cap mata: mata drop _quicksort()  // function updated

mata: // quicksort
void _quicksort(ids, dists, left, right)
{
    
	real scalar i, temp, tempDist, j, medianme, xxx, yyy, zzz
	
//	if ((right - left) <= 20) {  // taken out completely. was causing points to be dropped or not sorted correctly
		for (i = left + 1; i <= right + 1; i++) {
            
			temp = ids[i,1]
            tempDist = dists[temp,1]
            j = i - 1
            
			while ((j >= left) & (dists[ids[max((j,1)),1],1] > tempDist)) {  // if j hits zero, the second condition cannot be checked. so adding the max condition here
                ids[j + 1,1] = ids[j,1]
                j = j - 1
			}
			
            ids[j + 1,1] = temp
		}
}
end



************************
// 	   swapme     	  //  
************************

cap mata: mata drop swapme()

mata: // swapme
function swapme(arr,i,j)
{
	real scalar tmp1
	
	tmp1     = arr[i,1]
    arr[i,1] = arr[j,1]
    arr[j,1] = tmp1
}
end


********************
// 	 rescale	  // rescale a vector column	
********************

cap mata: mata drop rescale()

mata: // rescale
real vector rescale(points, a, b)
{
	newpoints = (b - a) * (points :- colmin(points)) :/ (colmax(points) - colmin(points)) :+ a
	return(newpoints)
}
end



********************************
***    END OF SUBROUTINES    ***
********************************



*****************************************
***  routines for adding data back    ***
*****************************************


************************
// 	   return index	  //   for exporting and merging with original data
************************

cap mata: mata drop returnindex()

mata: // returnindex
function returnindex(triangles4)
{

real vector myindex

myindex = J(rows(triangles4),1,.)

	for (i = 1; i <= rows(triangles4); i++ ) { 
		//myindex[i,1] = i
		myindex[i,1] = floor((i - 1) / 4) + 1 
	}

return(myindex)
}
end




************************
// 	   triangles  	  //  
************************

cap mata: mata drop fixtriangles()

mata: // fix triangles
function fixtriangles(inputtri,points)
{
	// duplicate the first observation

	exporttri = inputtri, J(rows(inputtri),2,.)

		for (i = 1; i <= rows(exporttri); i++ ) { 
			if (exporttri[i,2] != .) {
				exporttri[i,3..4] = points[exporttri[i,2],.]
			}
		}
	
		
	return(exporttri)	
}
end


cap program drop add_triangles
program define add_triangles
	version 15


	// these intermediate names will eventually be removed
	mata: triangles = select(triangles, (triangles[.,1] :< .)) // drop the missing rows
	mata: triangles2 = colshape(triangles',3) // // convert triplets to quadruplets. need to add a black row for drawing
	mata: triangles3 = triangles2, J(rows(triangles2),1,.) 
	mata: triangles4 = colshape(triangles3,1)
	mata: triindex   = returnindex(triangles4)
	mata: triangles5 = triindex,triangles4	
	mata: mytriangles = fixtriangles(triangles5,points)
	
	
	mata: st_matrix("triangles", mytriangles)
	mat colnames triangles = "tri_num" "tri_id" "tri_x" "tri_y"

	cap drop tri* // make sure the variables are clear
	
	qui svmat triangles, n(col)
		lab var tri_num "Triangle: number"
		lab var tri_id  "Triangle: _id"
		lab var tri_x   "Triangle: x"
		lab var tri_y   "Triangle: y"

	// drop the junk
	mata: mata drop triangles2 triangles3 triangles4 triangles5 
	cap drop tri_num  // dud
	
end



********************
// 	   hull  	  //  
********************



cap mata: mata drop fixhull()

mata: // generate hull
function fixhull(hull,points)
{
	// duplicate the first observation
	hull2 = hull \ . 
	hull2[rows(hull2),1] = hull[1,1]
	hull3 = J(rows(hull2),1,.), hull2, J(rows(hull2),2,.)


		for (i = 1; i <= rows(hull2); i++ ) { 
			hull3[i,1] = i
			hull3[i,3..4] = points[hull3[i,2],.]	
		}	
	return(hull3)	
}
end

**********************
// export the hull  //
**********************

cap program drop add_hull
program define add_hull
	version 15


cap drop hull* // make sure the variables are clear

	mata: myhull = fixhull(hull,points)
	mata st_matrix("hull", myhull)
	mat colnames hull = "hull_num" "hull_id" "hull_x" "hull_y"
	svmat hull, n(col)

		lab var hull_num "Hull: number"
		lab var hull_id  "Hull: _id"
		lab var hull_x   "Hull: x-coord"
		lab var hull_y   "Hull: y-coord"
		
		cap drop hull_num  // don't really need it
		
end




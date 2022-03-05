![StataMin](https://img.shields.io/badge/stata-2015-blue)


# Stata delaunay voronoi

Implementation of the [S-Hull Delaunay](http://www.s-hull.org/) algorithm in Stata. The Convex Hull is generated as a residual from the sweeping search algorithm, and the Voronoi tessellations are recovered as a dual to the Delaunay. 

*Note:* 
-   This is a beta release and still needs to be improved. It has been uploaded here for testing purposes only. Feedback is always welcome!
-   The package now requires **gtools** for fast reshaping (ssc install gtools, replace)


## Install the package:

Please install and replace the package you want to use it. Small tweaks are continuously being made. The following command can be used to install it directly in Stata:

```applescript
net install delaunay, from("https://raw.githubusercontent.com/asjadnaqvi/stata-delaunay-voronoi/main/installation/") replace force
```

The force option ensures that the package files are replaced even if Stata thinks they are the same. 


The syntax is:

```applescript
delaunay y x [if] [in], [rescale triangles hull voronoi(lines polygons) offset(value)]
```

where `y` and `x` are coordinates or data points. If `y` and `x` do not have similar value ranges, then rescale normalizes to the same interval calculates the triangles and rescales them back. The next set of options export the `triangles`, `hull`, and `voronoi` back to Stata.

See the help files for details:

```applescript
help delaunay
```

The summary of options is as follows:


| Option | Description |
| --- |--- |
| rescale | If the `y` and `x` coordinates do not have similar value ranges, then rescale normalizes to the same interval calculates the triangles and rescales them back. |
| triangles | exports back the Delaunay triangles as shapes. |
| hull | exports back the hull as line coordinates. |
| voronoi() | exports back the Voronoi as *lines* or *polygons* or both. |
| offset()  | Overwrites the clipping box for the Voronoi. Default is 5% over the (max - min) range. |


*   All options are optional. If none are specified, then nothing is added back to the data. But all information is stored as Stata matrices (see `mat dir`). The command automatically adds an identifier variable `_id` that can be used to trace back triangles and Voronoi tessellations back to the original observation.


*   Export the values back to Stata will increase the number of observations.


## Test the package:

This section shows code to test the package. In order to replicate the figures as they are, you need to set the following installed:

```applescript
ssc install gtools, replace

ssc install colrspace, replace	
ssc install palettes , replace	

ssc install schemepack, replace  
set scheme white_tableau

set seed 7543223
```


Generate some random data or use your own coordinates:

```applescript
set obs 1000  // have tested up to 10k observations

gen x = runiform(0,100)
gen y = runiform(0,100)


// create a donut (if you want)
drop if sqrt((x-50)^2 + (y-50)^2) > 50
drop if sqrt((x-50)^2 + (y-50)^2) < 10
```


Export everything back to Stata:

```applescript
delaunay x y, triangles hull voronoi(lines poly) offset(0.1)
```

See the stored matrices by typing `mat dir`.


## Convex hull

```applescript
	twoway ///
		(line hull_y hull_x, lw(0.3) lc(orange)) ///
		(scatter y x , msize(0.4) mc(black) ) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off) aspect(1) xsize(1) ysize(1) 
```

<img src="/figures/delaunay_hull.png" height="600">


## Delaunay triangles

```applescript
	twoway ///
		(area tri_y tri_x, cmissing(n) nodropbase fc(gs14%50) lw(0.1) lc(orange)) ///
		(scatter y x , msize(0.4) mc(black) ) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off) aspect(1) xsize(1) ysize(1) 
```


<img src="/figures/delaunay_triangles.png" height="600">

## Voronoi tessellations

### Lines

```applescript
	twoway ///
		(pcspike vline_y1 vline_x1 vline_y2 vline_x2, lw(0.1) lc(orange)) ///
		(scatter y x , msize(0.4) mc(black) ) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off)	aspect(1) xsize(1) ysize(1) 
```

<img src="/figures/delaunay_voronoi_lines.png" height="600">

### Polygons

```applescript
	twoway ///
		(area vpoly_y vpoly_x, cmissing(n) nodropbase fc(gs14%50) lw(0.1) lc(orange)) ///
		(scatter y x , msize(0.4) mc(black)) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off) aspect(1) xsize(1) ysize(1) 
```

<img src="/figures/delaunay_voronoi_polygons.png" height="600">


### Color fill polygons

Let's assign 6 groups to polygons

```applescript
cap drop tag temp
cap drop group		
		
egen tag = tag(vpoly_id)
gen temp = runiformint(1,6) if tag==1

bysort vpoly_id: egen group = max(temp)
cap drop tag temp
```

These could be defined by unique attributes like school types, or shop types, or location types. We can now color specific groups separately:

```applescript
local vpoly		
		
levelsof group, local(lvls)
local items = `r(r)'

foreach x of local lvls {	
	colorpalette HSV intense, n(`items') nograph
	local vpoly `vpoly' (area vpoly_y vpoly_x if group==`x' , cmissing(n) nodropbase fi(100) fc("`r(p`x')'%90") lw(0.03) lc(black)) ||
	
}
		
		
	twoway ///
		`vpoly' ///
		(scatter y x , msize(0.15) mc(black) ) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off)	aspect(1) xsize(1) ysize(1) 	
```

<img src="/figures/delaunay_voronoi_polygons_colors.png" height="600">


## Versions


### 1.10 (01 Mar 2022)

Major update with several bug fixes and added features.

*   Triangles being missed in triangulation fixed
*   Voronoi rays fixed
*   Voronoi polygon feature added
*   Offset feature added
*   Several code improvements



### 1.02 (29 Dec 2021)

Added the **rescale** option and removed the id requirement. The program now generates its own _id variable. From the help file: 

> Delaunay triangles, and subsequently Voronoi tessellations, are not agnostic about the scale of the x and y-axis. They were designed to deal with physical geometry and therefore expect x and y values to be on a similar scale. If we are working with data where one variable is several times the magnitude of the other, then the command will correctly execute the triangles but they will be stretched in one direction. The *rescale* option normalizes both the x and y variables on a common range, calculates the triangles and rescales them back to provide reasonable looking triangles.

Let's test the following example where the yaxis is scaled up:

```
clear


set obs 20

set seed 103

gen x = runiform(1, 5)
gen y = runiform(1, 5)


foreach i in 1 5 10 {

cap drop y2	
	gen y2 = y * `i'


// without rescaling

delaunay x y2, tri hull vor 

	twoway (scatter y2 x, msize(small)) ///
		(line hull_y hull_x, lw(thin)) ///
		(area tri_y tri_x, cmissing(n) nodropbase fc(gs14%10) lw(0.04) lc()), ///
		legend(off) aspect(1) ///
		title("y = `i'x (no rescaling)")


	twoway (scatter y2 x, msize(small)) ///
		(pcspike vor_y1 vor_x1 vor_y2 vor_x2, lw(0.2) lc()) ///
		, legend(off) aspect(1) xlabel(#10) ylabel(#10) ///
		title("y = `i'x (no rescaling)")


// with rescaling

delaunay x y2, tri hull vor rescale

	twoway (scatter y2 x, msize(small)) ///
		(line hull_y hull_x, lw(thin)) ///
		(area tri_y tri_x, cmissing(n) nodropbase fc(gs14%10) lw(0.04) lc()), ///
		legend(off) aspect(1) ///
		title("y = `i'x (with rescaling)")


	twoway (scatter y2 x, msize(small)) ///
		(pcspike vor_y1 vor_x1 vor_y2 vor_x2, lw(0.2) lc()) ///
		, legend(off) aspect(1) xlabel(#10) ylabel(#10) ///
		title("y = `i'x (with rescaling)")
}
```

We can compare the triangles:

<img src="/figures/triangles_norescale_1.png" height="180"><img src="/figures/triangles_norescale_5.png" height="180"><img src="/figures/triangles_norescale_10.png" height="180">

<img src="/figures/triangles_rescale_1.png" height="180"><img src="/figures/triangles_rescale_5.png" height="180"><img src="/figures/triangles_rescale_10.png" height="180">


and the Voronoi tessellations:

<img src="/figures/voronoi_norescale_1.png" height="180"><img src="/figures/voronoi_norescale_5.png" height="180"><img src="/figures/voronoi_norescale_10.png" height="180">

<img src="/figures/voronoi_rescale_1.png" height="180"><img src="/figures/voronoi_rescale_5.png" height="180"><img src="/figures/voronoi_rescale_10.png" height="180">


### 1.01 (19 Dec 2021)

The [if] [in] options were added to the program based on wbuchanan's suggestions. Other code improvements. Some redundancy removed from programs.

### 1.00 (5 Dec 2021)

First release.


## The to do list

1.   Add checks and error returns.
2.   Get rid of Mata junk. 
3.   Optimize export back to Stata from Mata. Currently it is very slow.
4.   Add e-class locals.



![StataMin](https://img.shields.io/badge/stata-2015-blue)


# Stata delaunay voronoi

Implementation of the S-Hull Delaunay triangulation algorithm in Stata. The algorithm also generates the convex hull. The Voronoi are recovered as a dual to the Delaunay.

*Note:* This is a beta release and still needs to be improved. It has been uploaded here for testing purposes only. Feedback is always welcome!


## Install the package:

Please install and replace the package everytime  you want to use it. Small tweaks are continuously being made. The following command can be used to install it directly in Stata:

```applescript
net install delaunay, from("https://raw.githubusercontent.com/asjadnaqvi/stata-delaunay-voronoi/main/installation/") force
```

The force option ensures that the files are replaced even if Stata thinks they are the same. For older versions see the `versions` folder for roll-back options.


The syntax is:

```applescript
delaunay x y [if] [in], [rescale] [triangles] [hull] [voronoi]
```

where `x` and `y` are coordinates. If the x and y coordinates do not have similar value ranges, then rescale normalizes to the same interval calculates the triangles and rescales them back. The last three options export the `triangles`, `hull`, and `voronoi` back to Stata for plotting.

See the help files for details:

```applescript
help delaunay
```

## Test the package:

You can use the above do file for test the code. A sample code is provided here:

Generate some random data or use your own coordinates:

```applescript
set obs 2000  // have tested up to 10k observations

gen x = runiform(0,100)
gen y = runiform(0,100)


// create a donut (if you want)
drop if sqrt((x-50)^2 + (y-50)^2) > 50
```


Export everything back to Stata:

```applescript
delaunay x y, triangles hull voronoi
```

Plot the triangles and the hull:

```applescript
	twoway ///
		(area tri_y tri_x, cmissing(n) nodropbase fc(gs14%10) lw(0.04) lc()) ///
		(line hull_y hull_x) ///
		(scatter y x , msize(0.2) mc(black) ) ///
		, ///
		legend(off) aspect(1) xsize(1) ysize(1)
```


<img src="/figures/delaunay1_triangles.png" height="500">

Plot the Voronoi tessellations:

```applescript
	twoway ///
		(pcspike vor_y1 vor_x1 vor_y2 vor_x2, lw(0.08) lc()) ///
		(scatter y x , msize(0.1) mc(black) ) ///
				, ///
		legend(off)	aspect(1) xsize(1) ysize(1)
```

<img src="/figures/delaunay2_voronoi.png" height="500">

## Versions

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

<img src="/figures/triangles_norescale_1.png" height="220"><img src="/figures/triangles_rescale_1.png" height="220">
<img src="/figures/triangles_norescale_5.png" height="220"><img src="/figures/triangles_rescale_5.png" height="220">
<img src="/figures/triangles_norescale_10.png" height="220"><img src="/figures/triangles_rescale_10.png" height="220">
<img src="/figures/triangles_norescale_100.png" height="220"><img src="/figures/triangles_rescale_100.png" height="220">

and the Voronoi tessellations:

<img src="/figures/voronoi_norescale_1.png" height="220"><img src="/figures/voronoi_rescale_1.png" height="220">
<img src="/figures/voronoi_norescale_5.png" height="220"><img src="/figures/voronoi_rescale_5.png" height="220">
<img src="/figures/voronoi_norescale_10.png" height="220"><img src="/figures/voronoi_rescale_10.png" height="220">
<img src="/figures/voronoi_norescale_100.png" height="220"><img src="/figures/voronoi_rescale_100.png" height="220">

### 1.01 (19 Dec 2021)

The [if] [in] options were added to the program based on wbuchanan's suggestions. Other code improvements. Some redundancy removed from programs.

### 1.00 (5 Dec 2021)

First release.

## Known issues

1.   For some point combinations, the last point is being skipped from triangles.  
2.   For some Voronoi lines on the edges, the infinite rays are not being calculated.


## In the pipeline

1.   The above errors.  
2.   Add [if] [in] options. DONE. 
3.   Add an option to check and correct indices.  
4.   Get rid of Mata junk 
5.   Separate the Mata calculations from export back to Stata. 
6.   Call the *clipline* command from within the program, and add box options.  
7.   Convert Voronoi lines to shapes for more interesting visualizations.  
8.   Add e-class locals.



![StataMin](https://img.shields.io/badge/stata-2015-blue) ![issues](https://img.shields.io/github/issues/asjadnaqvi/stata-delaunay-voronoi) ![license](https://img.shields.io/github/license/asjadnaqvi/stata-delaunay-voronoi) ![Stars](https://img.shields.io/github/stars/asjadnaqvi/stata-delaunay-voronoi) ![version](https://img.shields.io/github/v/release/asjadnaqvi/stata-delaunay-voronoi) ![release](https://img.shields.io/github/release-date/asjadnaqvi/stata-delaunay-voronoi)


# DELAUNAY

Implementation of the [S-Hull Delaunay](http://www.s-hull.org/) triangulation algorithm in Stata. The Convex Hull is generated as a residual from the triangles. Voronoi tessellations are recovered as a dual to the triangles. 

The code is based on Mapbox's [Delaunator](https://github.com/mapbox/delaunator) and D3's [Delaunay/Voronoi](https://github.com/d3/d3-delaunay) implementations.


## Install the package

Please install and replace the package you want to use it. Small tweaks are continuously being made. The following command can be used to install it directly in Stata:

```applescript
net install delaunay, from("https://raw.githubusercontent.com/asjadnaqvi/stata-delaunay-voronoi/main/installation/") replace force
```


The syntax is:

```applescript
delaunay y x [if] [in], [rescale triangles hull voronoi(lines polygons) offset(value) replace addbox]
```

where `y` and `x` are coordinates or data points. If `y` and `x` have completely different scales, then use `rescale` to normalize the values to generate reasonable-looking triangles and tesselations. Note that `triangles` and `voronoi()` export calculated data as files, and `hull` is added back to the original dataset as variables.

See the help files for details:

```applescript
help delaunay
```

The summary of options is as follows:


| Option | Description |
| --- |--- |
| rescale | If the `y` and `x` coordinates do not have similar value ranges, then rescale normalizes to the same interval calculates the triangles and rescales them back. |
| triangles | exports the Delaunay triangles in the *_triangles.dta* file. |
| hull | exports back the hull as point coordinates. |
| voronoi(lines polygons) | exports the Voronoi as *lines* or *polygons* or both as *_vorlines.dta* or *_vorpoly.dta* files. |
| offset()  | Overwrites the clipping box for the Voronoi. Default is 5% over the (max - min) range. |
| replace  | Replace the exported files/variables. |
| addbox  | An experimental option to add a bounding box to the triangles.  |


**NOTE:** Running the command automatically adds an identifier variable `_ID` that can be used to trace back triangles and Voronoi tessellations to the original observations.



## Test the package

This section shows code to test the package. In order to replicate the figures as they are, you need to set the following installed:

```applescript
ssc install colrspace, replace	
ssc install palettes , replace	

ssc install schemepack, replace  
set scheme white_tableau
```


Generate some random data or use your own coordinates:

```applescript
set obs 1000 
set seed 1337

gen x = runiform(0,100)
gen y = runiform(0,100)


// create a donut (if you want)
drop if sqrt((x-50)^2 + (y-50)^2) > 50
drop if sqrt((x-50)^2 + (y-50)^2) < 20
```

Higher observations result in a higher processing time. For example, 12k data points roughly take a minute for generating the triangles and tessellations.

Run the command:

```applescript
delaunay y x, tri hull vor(lines polygons) replace
```



### Convex hull

```applescript
	twoway ///
		(line hull_Y hull_X, lw(0.3) lc(orange)) ///
		(scatter y x , msize(0.4) mc(black) ) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off) aspect(1) xsize(1) ysize(1) 
```

<img src="/figures/delaunay_hull.png" height="600">


### Delaunay triangles

```applescript
use _triangles.dta, clear

	twoway ///
		(area _Y _X, cmissing(n) nodropbase fc(gs14%50) lw(0.1) lc(orange)) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off) aspect(1) xsize(1) ysize(1) 
```


<img src="/figures/delaunay_triangles.png" height="600">

### Voronoi tessellations

**Lines**

```applescript
use _vorlines.dta, clear

	twoway ///
		(pcspike vline_y1 vline_x1 vline_y2 vline_x2, lw(0.1) lc(orange)) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off)	aspect(1) xsize(1) ysize(1) 
```

<img src="/figures/delaunay_voronoi_lines.png" height="600">

**Polygons**

```applescript
use _vorpoly.dta, clear
	twoway ///
		(area _Y _X, cmissing(n) nodropbase fc(gs14%50) lw(0.1) lc(orange)) ///
		, ///
		xlabel(#10, nogrid) ylabel(#10, nogrid) ///
		legend(off) aspect(1) xsize(1) ysize(1) 
```

<img src="/figures/delaunay_voronoi_polygons.png" height="600">


### Color fill tessellations

Let's assign 6 groups to polygons

```applescript
clear
set obs 1000 
set seed 1337    

gen x = runiform(0,100)
gen y = runiform(0,100)


drop if sqrt((x-50)^2 + (y-50)^2) > 50
drop if sqrt((x-50)^2 + (y-50)^2) < 20

delaunay y x, vor(polygons) replace

gen group = runiformint(1,6)
```

These could be defined by unique attributes like school types, or shop types, or location types. We can now color specific groups separately:

```applescript
spmap group using _vorpoly, id(_ID) fcolor(Pastel2) legend(off)
```

<img src="/figures/delaunay_voronoi_polygons_colors.png" height="600">


## Getting creative with Delaunay

In this example, we use coordinates of the following picture of Dali:

<img src="/figures/dali.jpg" height="600">

Import and set up the file:

```applescript
import delim using "https://github.com/asjadnaqvi/stata-delaunay-voronoi/blob/main/data/dali.csv?raw=true"
cap drop v1

// flip the yaxis. error in export
qui summ y, meanonly
replace y = r(max) - y + 1
drop if value==1
```

Generate a sample of points:

```applescript
cap drop sample
gen sample = runiform() > 0.80	// Lower values = large sample = more processing time. Use carefully. 
keep if sample==1
count							// see how many values you are processing
```

Run the script:

```applescript
delaunay y x, tri addbox vor(lines) offset(0) replace
```


```applescript
use _vorlines, clear

	twoway ///
		(pcspike vline_y1 vline_x1 vline_y2 vline_x2, lw(0.06) lc(black)) ///
		, ///
		xlabel(, nogrid) ylabel(, nogrid) ///
		xscale(off) yscale(off)	 ///
		xtitle("") ytitle("") ///
		legend(off)	xsize(3.13) ysize(4)
```

<img src="/figures/dali_voronoi.png" height="600">

```applescript
use _triangles, clear

	twoway ///
		(line _Y _X, cmissing(n) nodropbase fi(100) fc(none) lw(0.03) lc(black)) ///
		, ///
		xlabel(, nogrid) ylabel(, nogrid) ///
		xscale(off) yscale(off)	 ///
		xtitle("") ytitle("") ///
		legend(off)	xsize(3.13) ysize(4)
```

<img src="/figures/dali_triangles.png" height="600">

### R script for converting images to coordinates

If you want to play with other pictures, you can use the following R script to process any image:

```R
library(imager)
library(dplyr)
library(scales)

# Convert to grayscale
load.image("dali.jpg") %>% grayscale() -> x

# Filter image. A higher threshold of "black" = more data points (use carefully)
x %>%
  threshold("30%") %>% 
  as.cimg() %>% 
  as.data.frame() -> df

write.csv(df,"dali.csv", row.names = TRUE)
```

Share your images if you use this script!

## Versions

### 1.2 (04 Sep 2022)

Major update (the package is no longer beta):

*   Voronoi and triangles are now exported as files. Hull is still added directly to the original file.
*   An `_ID` variable is added to the original file, which can be used to trace values in the exported files.
*   the `_vorpoly.dta` is a shapefile that can be used with the `spmap` command.
*   Major code optimizations. All the requirements to use reshape taken out.

### 1.11 beta (05 Mar 2022)

Minor update:

*   The input options x and y have been flipped to conform with Stata conventions.
*   Exporting the geometry back to Stata is now faster.
*   Some error messages added. Package check for `gtools` added.

### 1.10 beta (01 Mar 2022)

Major update with several bug fixes and added features:

*   Triangles being missed in triangulation fixed
*   Voronoi rays fixed
*   Voronoi polygon feature added
*   Offset feature added
*   Several code improvements



### 1.02 beta (29 Dec 2021)

Added the **rescale** option and removed the id requirement. The program now generates its own _id variable. From the help file: 

> Delaunay triangles, and subsequently Voronoi tessellations, are not agnostic about the scale of the x and y-axis. They were designed to deal with physical geometry and therefore expect x and y values to be on a similar scale. If we are working with data where one variable is several times the magnitude of the other, then the command will correctly execute the triangles but they will be stretched in one direction. The *rescale* option normalizes the x and y variables on a common range, does the calculations, and rescales them back for exporting.

*NOTE: The code below uses the old syntax and will not work with v1.2. It still needs to be updated, so it's current purpose is mostly illustrative.*

Let's test the following example where the yaxis is scaled up:

```
clear


set obs 20

set seed 103

gen x = runiform(1, 5)
gen y = runiform(1, 5)


foreach i in 1 5 10 {

cap drop y2	
	gen y2 = y * `i'  // scale up the y-axis


// without rescaling

delaunay y2 x, tri hull vor 

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

delaunay y2 x, tri hull vor rescale

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


### 1.01 beta (19 Dec 2021)

The [if] [in] options were added to the program based on wbuchanan's suggestions. Other code improvements. Some redundancy removed from programs.

### 1.00 beta (5 Dec 2021)

First release.


## The to do list

1.   Add checks and error returns.
2.   Get rid of Mata junk. 
3.   Optimize export back to Stata from Mata. Currently it is very slow.
4.   Add e-class locals.



{smcl}
{* 05March2022}{...}
{hi:help delaunay}{...}
{right:{browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":delaunay v1.11 (GitHub)}}

{hline}

{title:DELAUNAY (beta release)}

{p 4 4 2}
{cmd:delaunay} implements the {browse "http://www.s-hull.org":S-Hull algorithm} to generate Delaunay triangles. The Convex hull is a derived as a residual. The Voronoi tessellations are generated as a dual to the triangles. 

{p 4 4 2}
The Delaunay and Voronoi core routines are based on {browse "https://github.com/mapbox/delaunator":Mapbox} and {browse "https://github.com/d3/d3-delaunay":d3} implementations respectively. 

{p 4 4 2}
The S-hull is sweep-hull algorithm that iterates across radially sorted points and flips triangles as it moves along. It is exceptionally fast as compared to {browse "https://github.com/mapbox/delaunator#performance":other implementations}.


{marker syntax}{title:Syntax}
{p 8 15 2}

{cmd:delaunay} {it:y x} [if] [in] [, {cmdab:res:cale} {cmdab:tri:angles} {cmdab:h:ull} {cmdab:vor:onoi}({it:lines} {it:{ul:poly}gons}) {ul:off}set({it:value})]

{p 4 4 2}
The options inside square brackets [] are optional and are mostly used for exporting the geometries back to Stata. These are described as follows:

{synoptset 36 tabbed}{...}
{synopthdr}
{synoptline}

{p2coldent : {opt delaunay y x}}The command requires an {it:(y,x)} coordinate pair stored as two numeric variables. Please make sure that the order of the variables is correct.{p_end}

{p2coldent : {opt res:cale}} Delaunay triangles and Voronoi tessellations are not agnostic about the scale of the axes. They are designed to work with physical geometry and therefore expect x and y axes to be on a similar scale. If we are working with data where one variable is several times the magnitude of the second, then the command will correctly execute the triangles but they might be stretched in one direction. The {it:rescale} option normalizes both the x and y variables on the [0,1] range, calculates the triangles and rescales them back to provide reasonable looking triangles. {p_end}

{p2coldent : {opt tri:angles}}Export the triangles back to Stata as three new variables {it:tri_id, tri_x, tri_y}. The first variable, {it:tri_id} is the the points {it:id}, and the remaining two are the coordinates. The data structure is fastest option to draw shapes in Stata. These can be plotted using the {it:twoway area} option.{p_end}

{p2coldent : {opt h:ull}}Exports the convex hull back to Stata as three new variables {it:hull_id, hull_x, hull_y}. The first variable, {it:hull_id} corresponds to the point {it:_id}. The remaining two variables are the coordinates. They starting point is repeated in order to draw a fully enclosed hull. The hull can be plotted using the {it:twoway line} option.{p_end}

{p2coldent : {opt vor:onoi}({it:lines} {it:{opt poly:gons}})}Exports the Voronoi tesselations back to Stata as lines or polygons or both. The {it:lines} option generates four new variables {it:(vline_x1,vline_y1,vline_x2,vline_y2)} which are the 
start and end coordinates of a line. The {it:{opt poly:gons}} option generates three variables {it:vpoly_id,vpoly_x,vpoly_y}, which are the {it:_id} of the point, and the {it:x} and {it:y} coordinates of the shape.

{p2coldent : {opt off:set}(value)}Voronoi tesselations, have a bounding box defined as {it:xmax + xdisp}, {it:xmin - xdisp}, {it:ymax + ydisp}, {it:ymin - ydisp}, where {it:xdisp = (xmax - xmin) * offset} and {it:ydisp = (ymax - ymin) * offset}.  
The default offset is 0.05 or 5% displacement. This can be overwritten using the {opt off:set}(value) option. Values in the range of 5-50% work well. {p_end}


{synoptline}
{p2colreset}{...}

{p 4 4 2}
The program will generate a unique identifier {it:_id}. The triangles, the convex hull, and the voronoi geometry are also stored as Stata matrices ({cmd:mat dir}). 

{p 4 4 2}
If the {it:triangle}, {it:hull} and/or {it:voronoi} options are specified, then the coordinates of the vertices or edges are returned together with their respective ids which can be matched to the {it:_id} variable. 

{p 4 4 2}
For code examples, please see the {browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":{ul:Delaunay GitHub}} repository.

{hline}

{title:Dependencies}

- The generation of Voronoi polygons requires {browse "https://gtools.readthedocs.io/en/latest/":gtools} for fast reshaping. Make sure you have it installed:
{p 4 4 2}
{stata ssc install gtools, replace}


{title:To do list}

- Add more error checks.
- Optimize Mata to Stata export options.
- Remove redundancies in code.


{title:Version history}

- {bf:1.11}: Order of x and y flipped. Faster export to Stata. Error messages added.
- {bf:1.10}: Delaunay point skipping fixed. Voronoi rays fixed. Voronoi polygons and offset added.
- {bf:1.02}: Rescale function added. id option removed. The program generates its own _id variable.
- {bf:1.01}: Minor code cleanups. [if] [in] conditions (thanks to wbuchanan!). 
- {bf:1.00}: First release

{hline}

{title:Package details}

Version      : {bf:delaunay} v1.11
This release : 05 Mar 2022
First release: 05 Dec 2021
Repository   : {browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":GitHub}
Keywords     : Stata, delaunay, voronoi, convex hull, s-hull algorithm
License      : {browse "https://opensource.org/licenses/MIT":MIT}

Author       : {browse "https://github.com/asjadnaqvi":Asjad Naqvi} (WU Wien/IIASA, Austria)
E-mail       : asjadnaqvi@gmail.com
Twitter      : {browse "https://twitter.com/AsjadNaqvi":@AsjadNaqvi}


{title:Want to cite this package?}

{p 4 8 2}Naqvi, A. (2022). Delaunay v1.11: A Stata package for Delaunay triangles, Convex Hulls, and Voronoi tessellations. {browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":https://github.com/asjadnaqvi/stata-delaunay-voronoi}.


{title:Acknowledgments}

{p 4 4 2}
The package greatly benefitted from discussions with William Buchanan, Eric Melse, and the broad Stata community.


{title:References}

{p 4 8 2}Sinclair, D.A. (2016). S-hull: a fast radial sweep-hull routine for Delaunay triangulation. {it:arXiv} {browse "https://arxiv.org/abs/1604.01428v1":1604.01428}.

{p 4 8 2}Mapbox (2021). {browse "https://github.com/mapbox/delaunator":Delaunator library}.

{p 4 8 2}D3 (2021). {browse "https://github.com/d3/d3-delaunay":d3-delaunay}.

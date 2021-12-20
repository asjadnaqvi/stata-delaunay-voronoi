{smcl}
{* 13Dec2021}{...}
{hi:help delaunay}{...}
{right:{browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":delaunay v1.01 (GitHub)}}

{hline}

{title:DELAUNAY (beta release)}

{p 4 4 2}
{cmd:delaunay} implements the S-Hull algorithm to generate Delaunay triangles. The Voronoi triangles are generated as a dual
to the triangulation. The script is derived from Mapbox and d3js implementations. The S-hull implementation is incredibly
fast. It is now the industry benchmark for visualizations. Exporting and rendering the triangles on the screen takes time
especially if there are more than 5,000 triangles. The Voronoi algorithm is derived from the script by Mike Bostock, the main
person behind the D3js visualization toolkit. Line clipping is based on the {cmd:clipline} command. 



{marker syntax}{title:Syntax}
{p 8 15 2}
{cmd:delaunay} {it:x y} [if] [in], id({it:id}) [{cmdab:tri:angles}] [{cmdab:h:ull}] [{cmdab:vor:onoi}]

{synoptset 36 tabbed}{...}
{synopthdr}
{synoptline}

{p2coldent : {opt delauay x y}}The command requires an {it:(x,y)} coordinate pair stored as two numeric variables. Please make the order sequence of the variables is correct.{p_end}

{p2coldent : {opt id(id)}}This variable contains the sequence order of the coordinate pairs. This should be sorted and without gaps. (Still need to add options to check and parse this correctly in Mata){p_end}

{p2coldent : {opt tri:angles}}Export the triangles back to Stata as three new variables {it:tri_id, tri_x, tri_y}. The first variable, {it:tri_id} is the the points {it:id}, and the remaining two are the coordinates. The data structure is set up to generate shapes and (AFAIK), it is the fastest option in Stata. These can be plotted using the {it:twoway area} option.{p_end}

{p2coldent : {opt h:ull}}Exports the convex hull back to Stata as four new variables {it:hull_num, hull_id, hull_x, hull_y}. The first variable {it:hull_num} is an ordered serial number which can be used for sorting. The second variable, {it:hull_id} is the points {it:id}, and the remaining two are the coordinates. They starting point is repeated to allow us to generate a line that fully encloses the hull. The hull can be plotted using the {it:twoway line} option.{p_end}

{p2coldent : {opt vor:onoi}}Exports the Voronoi lines back to Stata as four new variables {it:vor_x1,vor_y1,vor_x2,vor_y2}. Voronoi grids are extracted as lines that are then clipped on the minmax values of the x,y space using the {cmd clipline} command. A 10% buffer is added to expand the box. The data structure is set up to generate a line. The Voronoi lines can be plotted using the {it:twoway pcspike} option.

{synoptline}
{p2colreset}{...}

{p 4 4 2}
If no options are specified, then nothing is returned. The triangles, the convex hull, and the voronoi lines are stored as Stata matrices (mat dir), that can be extracted.


{title:Known issues}

{p 4 4 2}
1. For some point combinations, the last point is being skipped from triangles.
2. For some voronoi lines on the edge, the infinite rays are not being calculated.

{title:To be fixed}

{p 4 4 2}
1. The above errors.
2. Add if options. DONE.
3. Add an option to check and correct indices.
4. Get rid of Mata junk
5. Separate the Mata calculations from export back to Stata.
6. Call the {cmd clipline} command from within the program, and add box options.
7. Convert Voronoi lines to shapes for more interesting visualizations.
8. Add e-class locals.

{hline}


Keywords: Stata, delaunay, voronoi, convex hull, s-hull algorithm
Version: {bf:delaunay} version 1.01
This version: 20 Dec 2021
First release: 05 Dec 2021
License: {browse "https://opensource.org/licenses/MIT":MIT}

Author: {browse "https://github.com/asjadnaqvi":Asjad Naqvi} (WU Wien and IIASA, Austria)
E-mail: asjadnaqvi@gmail.com
Twitter: {browse "https://twitter.com/AsjadNaqvi":@AsjadNaqvi}


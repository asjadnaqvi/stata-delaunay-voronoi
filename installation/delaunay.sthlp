{smcl}
{* 04September2022}{...}
{hi:help delaunay}{...}
{right:{browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":delaunay v1.2 (GitHub)}}

{hline}

{title:DELAUNAY}

{p 4 4 2}
{cmd:delaunay} implements the {browse "http://www.s-hull.org":S-Hull algorithm} (Sinclair 2016) for generating the Delaunay triangles. The Convex hull is a derived as a residual. The Voronoi tessellations are generated as a dual to the triangles. 

{p 4 4 2}
The Delaunay and Voronoi routines are based on {browse "https://github.com/mapbox/delaunator":Mapbox} and {browse "https://github.com/d3/d3-delaunay":D3} implementations respectively. 
S-hull, or sweeping-hull, is a divide-and-conquer algorithm that iterates over radially-sorted points, relative to an intial seed point, and generates and flips triangles as it moves along. 
It is exceptionally fast as compared to {browse "https://github.com/mapbox/delaunator#performance":other implementations}.


{marker syntax}{title:Syntax}

{p 8 15 2}
{cmd:delaunay} {it:y x} [if] [in] [, {cmdab:res:cale} {cmdab:tri:angles} {cmd:hull} {cmdab:vor:onoi}({it:lines} {it:polygons}) {cmdab:off:set}({it:value}) {cmd:replace} {cmd:addbox}]



{p 4 4 2}
The options are described as follows:

{synoptset 36 tabbed}{...}
{synopthdr}
{synoptline}

{p2coldent : {opt delaunay y x}}The command requires an {it:(y,x)} coordinate pair stored as two numeric variables.
Please make sure that the order of the variables is correct.
Running the command will automatically add the {it:_ID} variable to the dataset.{p_end}

{p2coldent : {opt res:cale}} Delaunay triangles and Voronoi tessellations are not agnostic about the scale of the axes. They are designed for physical geometry where the relation of the {it:x} and {it:y} scales matter.
If we are working with data where one variable is several times the magnitude of the other, then the command will correctly execute the triangles but they might be stretched in one direction.
The {cmd:rescale} option normalizes both the {it:x} and {it:y} variables on the [0,1] range, 
calculates the triangles and rescales them back to provide reasonable looking triangles.{p_end}

{p2coldent : {opt tri:angles}}Export the triangles in the {it:_triangles.dta} files. Each triangle is identified by a {it:_group} variable. Each coordinate has an {it:_ID} variable that allows it to be matched to the original data.{p_end}

{p2coldent : {opt hull}}Exports the convex hull back to Stata as three new variables {it:hull_ID, hull_X, hull_Y}. {it:hull_ID} corresponds to the variable {it:_ID}. The remaining two variables are the coordinates.
The starting point is repeated at the end, in order to draw a fully enclosed hull in Stata.{p_end}

{p2coldent : {opt vor:onoi}({it:lines} {it:polygons})}Exports the Voronoi tesselations back to Stata as lines or polygons or both. 
The {it:lines} option generates a {it:_vorlines.dta} file, which contains four variables {it:vline_x1}, {it:vline_y1}, {it:vline_x2}, {it:vline_y2} which are the end coordinates of each line. 
The {it:polygons} option generates a {it:_vorpoly.dta} file, which contains three variables {it:_ID}, {it:_X}, {it:_Y}. This is essentially a shapefile can be used with the main data using the {stata help spmap:spmap} command.{p_end}

{p2coldent : {opt off:set}(value)}Voronoi tesselations, have a bounding box defined as {it:xmax + xdisp}, {it:xmin - xdisp}, {it:ymax + ydisp}, {it:ymin - ydisp}, where {it:xdisp = (xmax - xmin) * offset} and {it:ydisp = (ymax - ymin) * offset}.  
The default offset is 0.05 or +5% displacement.{p_end}

{p2coldent : {opt replace}}This option replaces the existing files and hull variables. If {cmd:replace} is not specified, and the files already exist, the command will throw an error.
This ensures that files are not accidentally overwritten.{p_end}

{p2coldent : {opt addbox}}This is an experimental command that is used for adding a bounding box for triangles.
This is still being tested and it is highly recommended to use this only with triangles. It might result errors with Voronoi exports.{p_end}

{synoptline}
{p2colreset}{...}

{p 4 4 2}
NOTE: The program will generate a unique identifier, {it:_ID}, in the main file. If this variable already exists, the {cmd:replace} option will overwrite it. 


{p 4 4 2}
For code examples, please see the {browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":{ul:GitHub}} page.
If you find and error, please open an {browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi/issues":Issue}. 


{title:Version history}

- {bf:1.2} : Stable release. Triangles and Voronoi are exported as files. The command is now much faster.
- {bf:1.11}: Order of x and y flipped. Faster export to Stata. Error messages added.
- {bf:1.10}: Delaunay point skipping fixed. Voronoi rays fixed. Voronoi polygons and offset added.
- {bf:1.02}: Rescale function added. id option removed. The program generates its own _id variable.
- {bf:1.01}: Minor code cleanups. [if] [in] conditions. 
- {bf:1.00}: First release


{title:Package details}

Version      : {bf:delaunay} v1.2
This release : 04 Sep 2022
First release: 05 Dec 2021
Repository   : {browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":GitHub}
Keywords     : Stata, delaunay, voronoi, convex hull, s-hull, algorithm,
License      : {browse "https://opensource.org/licenses/MIT":MIT}

Author       : {browse "https://github.com/asjadnaqvi":Asjad Naqvi}
E-mail       : asjadnaqvi@gmail.com
Twitter      : {browse "https://twitter.com/AsjadNaqvi":@AsjadNaqvi}


{title:Want to cite this package?}

{p 4 8 2}Naqvi, A. (2022). Delaunay v1.2: A Stata package for Delaunay triangles, Convex hulls, and Voronoi tessellations. {browse "https://github.com/asjadnaqvi/stata-delaunay-voronoi":https://github.com/asjadnaqvi/stata-delaunay-voronoi}.


{title:Acknowledgments}

{p 4 4 2}
The package greatly benefited from discussions with Eric Melse and William Buchanan.


{title:References}

{p 4 8 2}Sinclair, D.A. (2016). S-hull: a fast radial sweep-hull routine for Delaunay triangulation. {it:arXiv} {browse "https://arxiv.org/abs/1604.01428v1":1604.01428}.

{p 4 8 2}Mapbox (2022). {browse "https://github.com/mapbox/delaunator":Delaunator library}.

{p 4 8 2}D3 (2022). {browse "https://github.com/d3/d3-delaunay":d3-delaunay}.

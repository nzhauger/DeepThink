MOOSalogScraper takes in an alog file and grabs the NAV X, Y and Z data (and the times associated) 
from the uSimMarine MOOSapp.
If we want to take more data like heading, we can make this a function with the argument being the variable we wanted
and then output the array of data.
that would be more scalable than writing each individual variable we want and writing elseifs for each.

As of 7FEB25 here's some issues I immediately see:
- I preset the axis but it can probably be constructed better
- The ship is made of polygons, but I didn't make it 3D. I wonder if we can take polygon points from the Mechanical
  Engineer's model and make a polygon that way? It's not hard, just time consuming to make a patch object manually out of points.
-Nick

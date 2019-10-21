isolines documentation
**********************


**Street Networks Isolines and Isochrones**

Create Isolines and Isochrones for street networks.

**isolines** is a Python library for creating street networks isolines (equal-distance)
and isochrones (equal-time) polygons with just one line of code.

**isolines** allows you to create isolines:
from OSM data **without any data input** or **use your own data**:  your own network and source locations

Main Features
^^^^^^^^^^^^^
- Generation of isolines/isochroones for areas, segments, and point locations
- Augmentation of the network for precise isoline/isochroone delineation
- Based on a **concave hull** algorithm by default (for walking networks, otherwise convex hull is the default
  and concave hull is optional)
- Source location can be either an address string to be geocoded using OSM Nominatim or a shapely geometry
- Accepts various graph inputs:  the ``GpdIsoliner`` class excepts edges GeoDataFrames and the
  ``NxIsoliner`` class accepts NetworkX graphs as input
- Downloads a street network from OpenStreetMap (OSM) when no graph input is provided
  (the ``OsmIsoliner`` class)

Other Features
^^^^^^^^^^^^^^
  - Easy built-in visualization
  - Dynamically tweaking isolines/isochroones delineation by re-setting the default
    knn parameter that controls the concave hull heuristic
  - Dynamically adding isolines/isochroones to an existing instance
  - Extraction of the augmented output graph
  - Generation of nodes and edges GeoDataFrames

Examples
^^^^^^^^

 - Create isolines/isochroones for complex geometries (polygons and linestrings) as well as simple point geometries::

    import isolines as il
    iso = il.OSMIsoliner(
    'Prospect Park, Brooklyn, NYC, USA',
     metric = 'time',
     values=[3, 6, 9]
     )

    iso.get_isolines()


+---+-------------------------------+-----+
|   | geometry                      | time|
+---+-------------------------------+-----+
| 0 | POLYGON ((-73.9703 40.6600,...| 3   |
+---+-------------------------------+-----+
| 1 | POLYGON ((-73.9718 40.6589,...| 6   |
+---+-------------------------------+-----+
| 2 | POLYGON ((-73.9728 40.6572,...| 9   |
+---+-------------------------------+-----+

 plot::

    iso.plot_isolines()



.. image:: /../figs/prospect_park.png
   :align: center
   :alt:



- The isolines/isochroones boundaries are not confined to the existing network nodes. The network is amended to include
  new source and target nodes based on the input geometry and distances/times specified, so that large edges are cut in
  the respective location of new source/ target nodes to yield a more realistic isolines/isochroones boundary::

    iso = il.OSMIsoliner(Point(179370.985,664422.488), values=[500, 250],crs = 2039)
    iso.plot_isolines(plot_nodes=True, cmap = 'GnBu')


.. image:: /../figs/habima1.png
   :align: center
   :alt:

- The output isolines/isochroones are based on a **concave hull** (rather than convex hull) heuristic by default
  resulting in a more realistically shaped and accurate isolines/isochroones::

    iso = il.OSMIsoliner(Point(179370.985,664422.488), values=[500, 250],crs = 2039, knn = 20)
    iso.plot_isolines(plot_nodes=True, cmap = 'GnBu')

.. image:: /../figs/habima2.png
   :align: center
   :alt:

- isolines lets you either download a network from **OpenStreetMap (OSM)**
  or use an edges **geopandas GeoDataFrame**
  or use a **NetworkX graph**::

    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import LineString
    df = gpd.GeoDataFrame.from_file('../data/tl_2019_36047_edges/tl_2019_36047_edges.shp')
    # pre-process: add edges in the opposite direction
    df2 = df.copy()
    df2['TNIDF'] = df['TNIDT'].copy()
    df2['TNIDT'] = df['TNIDF'].copy()
    # reverse the line geometry coordinate sequence
    df2['geometry'] = df['geometry'].apply(lambda x: LineString(x.coords[::-1]))
    df = pd.concat([df, df2]).reset_index(drop = True)
    self = il.GDFIsoliner('Prospect Park, Brooklyn, NYC, USA',
                                edges = df,
                                network_type = 'walk',
                                metric = 'time',
                                values=[3, 6, 8, 16, 20],
                                edge_idcol = 'TLID',
                                fromcol = 'TNIDF',
                                tocol = 'TNIDT',
                                sample= 400,
                                knn = 50
                                     )

    self.plot_isolines()


.. image:: /../figs/tiger.png
   :align: center
   :alt:
- Add isolines to an existing instance::

    isochrones_from_poly = il.OSMIsoliner('Prospect Park, Brooklyn,NYC, USA',
                                     network_type = 'walk',
                                     values=[500, 1500, 2500, 3500],
                                     ft= True,
                                     sample = 600)
    isochrones_from_poly.plot_isolines()

.. image:: /../figs//prospect_ft1.png
   :align: center
   :alt:

add isolines::

  isochrones_from_poly.add_isolines([1000, 2000, 3000, 4000])
  isochrones_from_poly.plot_isolines()


.. image:: /../figs//prospect_ft2.png
   :align: center
   :alt:

- Extract nodes and edges GeoDataFrames, NetworkX graphs::

    nodes = isochrones.get_nodes()
    edges = isochrones.get_edges()
    G = isochrones.get_graph()












.. toctree::
   :maxdepth: 2

   examples
   api
   tutorial

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



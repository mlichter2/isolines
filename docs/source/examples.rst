.. _examples:

********
EXAMPLES
********


Create isolines/isochroones for complex geometries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you only want to output a GeoDataFrame of isolines/isochrones 
use the isolines function::

    import isolines as il
    df = il.isolines('Prospect Park, Brooklyn, NYC, USA', metric = 'time',values=[3, 6, 9])
    df


+---+-------------------------------+-----+
|   | geometry                      | time|
+---+-------------------------------+-----+
| 0 | POLYGON ((-73.9715 40.6492,...| 3   |
+---+-------------------------------+-----+
| 1 | POLYGON ((-73.9708 40.6470,...| 6   |
+---+-------------------------------+-----+
| 2 | POLYGON ((-73.9704 40.644,...| 9   |
+---+-------------------------------+-----+

However, if you which to explore your result with built-in visualization and be able to 
amend them and perform EDA use one of the following classes: GpdIsolimes, NXIsoliner, OsmIsolines::

    iso = il.OsmIsolines('Prospect Park, Brooklyn, NYC, USA', metric = 'time', values=np.arange(2.5, 22.5, 2.5), unit = 'ft', sample = 600)
    iso.plot_isolines(figsize = (10, 10))




.. image:: /../figs/prospect_park.png
   :align: center
   :alt:


Augmentation of the network for precise isoline/isochroone delineation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The isolines/isochroones boundaries are not confined to the existing network nodes. The network is amended to  include
 new source and target nodes based on the input geometry and distances/times specified, so that large edges are cut in
 the respective location of new source/ target nodes to yield a more realistic isolines/isochroones boundary.
 the built-in visualization can can help explore and refine these boundaries::

    iso = il.OsmIsolines('Bozeman High School, Bozeman, Montana, USA',values=[250, 500],sample = 200)
    iso.plot_isolines(plot_nodes=True, plot_source_nodes=True,figsize = (10,10))

.. image:: /../figs/bozeman1.png
   :align: center
   :alt:

Isoline boundary adjustment
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The output isolines/isochroones are based on a **concave hull** (rather than convex hull) heuristic by default
resulting in a more realistically shaped and accurate isolines/isochroones. These can be tweaked to get even more
refined delineation if needed::

    iso.change_isolines(knn = 15)
    iso.plot_isolines(plot_nodes=True, plot_source_nodes=True,figsize = (10,10))

.. image:: /../figs/bozeman2.png
   :align: center
   :alt:


Use different data sources
^^^^^^^^^^^^^^^^^^^^^^^^^^
isolines lets you either download a network from **OpenStreetMap (OSM)**
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
    tiger = il.GpdIsolines('Prospect Park, Brooklyn, NYC, USA',
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

    tiger.plot_isolines()


.. image:: /../figs/tiger.png
   :align: center
   :alt:

Add isolines to an existing instance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    isochrones = il.OsmIsolines('Prospect Park, Brooklyn,NYC, USA',
                                 network_type = 'walk',
                                 values=[500, 1500, 2500, 3500],
                                 unit= 'ft',
                                 sample = 600)
    isochrones.plot_isolines()


.. image:: /../figs/habima1.png
   :align: center
   :alt:

add isolines::

    isochrones.add_isolines([1000, 2000, 3000, 4000])
    isochrones.plot_isolines()


.. image:: /../figs/habima2.png
   :align: center
   :alt:

Generate nodes and edges GeoDataFrames, extract NetworkX graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    nodes = isochrones.get_nodes()
    edges = isochrones.get_edges()
    G = isochrones.get_graph()



.. toctree::
   :maxdepth: 2

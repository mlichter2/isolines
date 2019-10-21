from isolines.utils import *
import osmnx as ox
import networkx as nx
import numpy as np
from shapely.geometry import *
from shapely import wkt
import geopandas as gpd
import pandas as pd
from geopy.geocoders import Nominatim
import rtree
import time
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from dask import delayed
import warnings


class OsmIsolines(object):
    """


   :param location: an address string, or a shapely geometry
    the source location

   :param network_type: string
    {'walk', 'drive'}
    the type of mobility, default is ``'walk'``

   :param metric: string
    {'distance', 'time'}
    determines if the output will be isolines based on distances or isochrones based on time, default is ``'distance'``

   :param values: list or numpy array
    a list or array of time values in minutes or distance values in meters or feet (if unit = True), default is [500]

   :param speed: number
    if the ``metric`` parameter is set to ``'time'`` this will be the speed of walking im meters per
    minute or feet per minute (if ``unit =='ft'``), if network_type is set to ``'drive'``. this
    will be kilometers per hour or mile per hour (if ``unit =='ft'``)

   :param address_type: string
    {'all', 'point'}
    for a geocoded address this determines if to get back a point address or multipoint,
    polygon, multipolygon, linestring, multilinestring when relevant

   :param sample: number
    if a geocoded location geometry is polygon, multipolygon, linestring, multilinestring it is
    sampled in evenly spaced intervals (the boundary of polygons or the line itself) and sample points from it are
    added as nodes in the network, the default sample interval length is 100 m. or 330 ft.

   :param crs: int
    the crs of the input location and other inputs if relevant (input edges dataframe or
    networkX graph) - default is 4326

   :param convex: boolean
    whether the isoline/isochrone will be constructed using a convex or concave hull
    algorithms - default is ``False`` if ``network_type == 'walk'`` and ``True`` if ``network type == 'drive'``

   :param knn: int
    if the isoline/isochrone boundary is concave hull - the hull algorithm is based on
    k-nearest-neighbors (knn) heuristic. If no knn is provided, knn will be determined based on the number of points.
    smaller knn - will produce a more and elaborate isoline/isochrone

   :param unit: str
    {'all', 'point'}
    whether measurement are in meters or feet, default is meters

   :param validate_geometry: boolean
    if ``True`` the isolines topological validity and counter clock wise (ccw) coordinate sequence direction if
    ensured. default is ``True``

   :param par: boolean
    whether data processing is parallelized, default is ``True``

   :param smooth: boolean
    whether to smooth the isoline boundary, default is ``False``

   :param sigma: number
    standard deviation for Gaussian kernel if ``smooth ==True``, default is 2

   :param num_points_factor: int
    the number of points determines the density of vertices for smoothing if ``smooth ==True``, default is 5

   :return: None


    """

    def __init__(self,
                 location,
                 **kwargs):
        #: input location (str or shapely geometry)
        self.location = location
        #: network type ('walk' or 'drive')
        self.network_type = 'walk'
        #: 'distance' for isolines, 'time for isochrones
        self.metric = 'distance'
        #: an array of time/distance values
        self.values = [500]
        #: speed of walking or driving if the metric is 'time'. m. per min. or ft. per min if network type is walk
        #: and the units are m. or ft. respectfuly, : km. per hour. or miles per hour if network type is drive and the
        #: units are m. or ft. respectfuly
        self.speed = 84
        #: spacing interval of polygon/polyline boundary
        self.sample = 100
        #: crs
        self.crs = 4326
        #: whether the boundary delineation should be convex or concave
        self.convex = None
        #: the number of nearest neighbors used in the concave hull algorithm
        self.knn = None
        #: should the location be a point or any other type of geometry
        self.address_type = 'all'
        #: the units of measure
        self.unit = 'm'
        #: whether to make sure the output isolines are topologically valid
        self.validate_geometry = True
        #: the edge/graph id column if relevant
        self.edge_idcol = 'id'
        #: the edge from column if relevant
        self.fromcol = 'u'
        #: the edge to column if relevant
        self.tocol = 'v'
        #: the edge geometry column if relevant
        self.edge_geometrycol = 'geometry'
        #: the edge weight column if relevant
        self.weightcol = 'length'
        #: oneway attribute of the edges
        self.onewaycol = 'oneway'
        #: use dask parallelization
        self.par = True
        #: whether to smooth the output isolines boundary
        self.smooth = False
        #: standard deviation for Gaussian kernel if smooth ==True
        self.sigma = 2
        #: int the number of points determine the density of vertices for smoothing if smooth ==True
        self.num_points_factor = 5
        self.args = kwargs
        #: a dictionary of target nodes for every value
        self.iso_nodes = self.iso_nodes = {}

        t = time.time()

        self._set_params(self.args)

        self._set_distances()

        t = time.time()

        if isinstance(self.location, (str)):
            self._geocode_address()

        t = time.time()
        if isinstance(self.location, (Point,
                                      MultiPoint,
                                      Polygon,
                                      MultiPolygon,
                                      LineString,
                                      MultiLineString)):

            self._create_graph_from_geometry()

        else:
            raise TypeError('The Location provided is not valid')

        t = time.time()
        self.source_points = self._create_source_points()

        t = time.time()
        self._add_source_nodes()

        t = time.time()
        self.isos = gpd.GeoDataFrame(crs={'init': 'epsg:' + str(self.utm)})
        self._create_isolines(self.knn)

    def _create_isoline(self, idx, knn):
        t = time.time()
        nodes_coords = self._extract_target_nodes_coords(idx, self.distances, self.values)

        t = time.time()
        if len(nodes_coords) == 1:
            iso = self._create_iso(nodes_coords[0], knn=knn)
        else:
            iso = MultiPolygon([self._create_iso(i, knn=knn) for i in nodes_coords])

        return idx, iso

    def _change_isoline(self, idx, knn):
        nodes_coords = self.iso_nodes[self.values[idx]]

        # take care of disconnected components
        wcc = list(nx.weakly_connected_components(self.GP))

        if len(wcc) > 1:
            polygon = self.geom.buffer(self.distance + 100)
            if isinstance(polygon, MultiPolygon):

                filter = [[Point(x).intersects(p) for x in nodes_coords] for p in polygon]
                nodes_coords = [nodes_coords[f] for f in filter]

                # nodes_coords =[item for sublist in nodes_coords  for item in sublist]
                nodes_coords = [i for i in nodes_coords if i.shape[0] > 0]

                if len(nodes_coords) == 1:
                    iso = self._create_iso(nodes_coords[0], knn=knn)
                else:
                    iso = MultiPolygon([self._create_iso(i, knn=knn) for i in nodes_coords])
            else:
                iso = self._create_iso(nodes_coords, knn=knn)

        else:
            iso = self._create_iso(nodes_coords, knn=knn)

        return idx, iso

    def _create_isolines(self, knn):

        isos = []
        if self.par == True:
            for idx, d in enumerate(self.values):
                iso = delayed(self._create_isoline)(idx, knn)
                isos.append(iso)

            results = delayed()(isos)
            isos = results.compute()
            self.results = results
        else:
            for idx, d in enumerate(self.values):
                iso = self._create_isoline(idx, knn)
                isos.append(iso)

        self.isos = gpd.GeoDataFrame({
            'geometry': [i[1] for i in isos],
            self.metric: [self.values[i[0]] for i in isos]
        }, crs={'init': 'epsg:' + str(self.utm)})

        self.isos['geometry'] = self.isos['geometry'].apply(
            lambda x: make_ccw(make_valid_polygon(x)))
        self.isos = self.isos.sort_values(self.metric).reset_index(drop=True)

        if self.smooth:
            self.isos['geometry'] = self.isos['geometry'].apply(
                lambda x: gaussian_smooth_geom(x, self.sigma, self.num_points_factor))
        try:

            self.isos['geometry'] = self.isos.apply(lambda x: x['geometry'].intersection(
                self.isos.iloc[x.name + 1]['geometry']) if x.name < len(self.isos) - 1 else x['geometry'], axis=1)
        except:
            pass

    def change_isolines(self, **kwargs):
        """

        :param knn:  int
         change the k nearest neighbors controlling the concave hull boundary
        :param smooth: boolean
         whether to smooth the boundary or not
        :param sigma: number
         standard deviation for Gaussian kernel
        :param num_points_factor: int
         the number of points determine the density of vertices  - resolution
        :return: geopandas GeoDataFrame
        """

        for name in kwargs:
            setattr(self, name, kwargs[name])

        isos = []
        if self.par == True:
            for idx, d in enumerate(self.values):
                iso = delayed(self._change_isoline)(idx, self.knn)
                isos.append(iso)

            results = delayed()(isos)
            isos = results.compute()
        else:
            for idx, d in enumerate(self.values):
                iso = self._change_isoline(idx, self.knn)
                isos.append(iso)

        self.isos = gpd.GeoDataFrame({
            'geometry': [i[1] for i in isos],
            self.metric: [self.values[i[0]] for i in isos]
        }, crs={'init': 'epsg:' + str(self.utm)})

        self.isos['geometry'] = self.isos['geometry'].apply(
            lambda x: make_ccw(make_valid_polygon(x)))

        self.isos = self.isos.sort_values(self.metric).reset_index(drop=True)

        if self.smooth:
            self.isos['geometry'] = self.isos['geometry'].apply(
                lambda x: gaussian_smooth_geom(x, self.sigma, self.num_points_factor))
        try:

            self.isos['geometry'] = self.isos.apply(lambda x: x['geometry'].intersection(
                self.isos.iloc[x.name + 1]['geometry']) if x.name < len(self.isos) - 1 else x['geometry'], axis=1)
        except:
            pass

    def add_isolines(self, values):

        """
        Adds isolines/isochrones to an existing instance

        :param values: list/array
         values to add
        :return: None

        """

        values = np.array(values)

        if self.metric == 'time':
            if self.network_type == 'walk':

                if self.unit == 'ft':
                    if 'speed' not in self.args.keys():
                        self.speed = 276

                    distances = values * self.speed / 3.28084
                else:

                    distances = values * self.speed

            else:
                if self.unit == 'ft':

                    distances = (values / 60) * self.speed * 1609.344

                else:

                    distances = (values / 60) * self.speed * 1000

        else:
            if self.unit == 'ft':
                distances = values * 0.3048
            else:
                distances = values
        max_val = distances.max()

        if max_val > self.distance:
            values = values[distances <= self.distance]
            distances = distances[distances <= self.distance]

            if distances.shape[0] > 0:

                warnings.warn("Some of the values are larger than the maximum value set on the Isolines (" + str(
                    round(self.unit_dist, 2)) + " " + str(
                    self.unit) + "). New Isolines will be calculated for the maximum value set on the instance. "
                                 "Alternatively, re-initiate the Isolines with a larger maximum value.",
                              Warning)

            else:

                warnings.warn(
                    "All of the values are larger than the maximum value set on the Isolines instance. No new "
                    "Isolines will be calculated. To calclulate larger isolines, re-initiate the Isolines with a "
                    "larger maximum value.",
                    Warning)

        isos = []
        for idx in range(0, values.shape[0]):
            nodes_coords = self._extract_target_nodes_coords(idx, distances, values)
            if len(nodes_coords) == 1:
                iso = self._create_iso(nodes_coords[0], knn=self.knn)
            else:
                iso = MultiPolygon([self._create_iso(i, knn=self.knn) for i in nodes_coords])

            isos.append(iso)

        self.isos = pd.concat([self.isos,
                               gpd.GeoDataFrame({
                                   'geometry': isos,
                                   self.metric: values
                               })])
        self.isos = self.isos.sort_values(self.metric).reset_index(drop=True)

        self.isos['geometry'] = self.isos.apply(lambda x: x['geometry'].intersection(
            self.isos.iloc[x.name + 1]['geometry']) if x.name < len(self.isos) - 1 else x['geometry'], axis=1)

        if self.validate_geometry == True:
            self.isos['geometry'] = self.isos['geometry'].apply(
                lambda x: make_ccw(make_valid_polygon(x)))

        self.distances = np.unique(np.append(self.distances, distances))

    def get_isolines(self):
        """
         Returns a geopandas GeoDataFrame of the outpot isolines/isochrones
          :return: geopandas GeoDataFrame
        """
        result = self.isos.to_crs({'init': 'epsg:' + str(self.crs)}).copy()
        if self.validate_geometry == True:
            result['geometry'] = result['geometry'].apply(
                lambda x: make_valid_polygon(x))
        return result

    def get_nodes(self):
        """
        Generates a geopandas GeoDataFrame of the network nodes
         :return: geopandas GeoDataFrame
        """
        edges = self.get_edges()
        ids = edges.groupby(self.fromcol)[self.fromcol].head(1).tolist() + edges.groupby(self.tocol)[
            self.tocol].head(1).tolist()
        nodes_coords = np.array(
            edges.groupby(self.fromcol)['geometry'].head(1).apply(lambda x: x.coords[0]).tolist() +
            edges.groupby(self.tocol)['geometry'].head(1).apply(lambda x: x.coords[-1]).tolist())

        nodes = pd.DataFrame({self.edge_idcol: ids,
                              'lon': nodes_coords[:, 0],
                              'lat': nodes_coords[:, 1]}).drop_duplicates()
        nodes = gpd.GeoDataFrame(nodes,
                                 geometry=[Point(i) for i in zip(nodes['lon'], nodes['lat'])])
        nodes.crs = edges.crs

        return nodes

    def get_edges(self):
        """
        Generates an edges GeoDataFrame
         :return: geopandas GeoDataFrame
        """
        self.fromcol = 'source'
        self.tocol = 'target'
        edges = nx.to_pandas_edgelist(self.GP).copy()
        edges = gpd.GeoDataFrame(edges, geometry=[LineString(x) for x in edges["geometry"]],
                                 crs={'init': 'epsg:' + str(self.utm)})
        edges = edges.to_crs({'init': 'epsg:' + str(self.crs)})
        return edges

    def get_graph(self):
        """
        Returns the network graph used to generate the isolines/isochroones
         :return: NetworkX MultiDiGraph
        """
        edges = self.get_edges().to_crs({'init': 'epsg:' + str(self.crs)}).copy()
        G = self.GP.copy()
        idx = set()
        group = pd.DataFrame(
            {'count': edges.groupby([self.fromcol, self.tocol]).size()})
        for c in group['count'].unique():
            df = edges.iloc[~edges.index.isin(idx)].groupby(
                [self.fromcol, self.tocol]).head(1)

            idx.update(df.index)
            edges_names = list(zip(df[self.fromcol],
                                   df[self.tocol],
                                   np.full(len(df), c - 1)
                                   ))

            geom_dict = dict(zip(edges_names, df['geometry']))

            nx.set_edge_attributes(G, name='geometry', values=geom_dict)

        return G

    def get_sample_points(self):
        """
        Returns the sample points of the input location
         :return: geopandas GeoDataFrame
        """
        return self.get_projected_source_nodes().to_crs({'init': 'epsg:' + str(self.crs)})

    def get_location(self):
        """
         Returns the input location geometry
          :return: shapely geometry
         """
        return self.location

    def get_projected_isolines(self):
        """
         Returns a geopandas GeoDataFrame of the outpot isolines/isochrones projected to UTM
          :return: geopandas GeoDataFrame
         """
        return self.isos

    def get_projected_nodes(self):
        """
        Generates a geopandas GeoDataFrame of the network nodes projected to UTM
         :return: geopandas GeoDataFrame
        """
        edges = self.get_projected_edges()
        ids = edges.groupby(self.fromcol)[self.fromcol].head(1).tolist() + edges.groupby(self.tocol)[
            self.tocol].head(1).tolist()
        nodes_coords = np.array(
            edges.groupby(self.fromcol)['geometry'].head(1).apply(lambda x: x.coords[0]).tolist() +
            edges.groupby(self.tocol)['geometry'].head(1).apply(lambda x: x.coords[-1]).tolist())

        nodes = pd.DataFrame({self.edge_idcol: ids,
                              'x': nodes_coords[:, 0],
                              'y': nodes_coords[:, 1]}).drop_duplicates()
        nodes = gpd.GeoDataFrame(nodes,
                                 geometry=[Point(i) for i in zip(nodes['x'], nodes['y'])])
        nodes.crs = edges.crs

        return nodes

    def get_projected_edges(self):
        """
        Generates a geopandas GeoDataFrame of the network edges projected to UTM
         :return: geopandas GeoDataFrame
        """
        self.fromcol = 'source'
        self.tocol = 'target'
        edges = nx.to_pandas_edgelist(self.GP).copy()
        edges = gpd.GeoDataFrame(edges, geometry=edges["geometry"], crs={'init': 'epsg:' + str(self.utm)})
        return edges

    def get_projected_graph(self):
        """
        Returns the network graph used to generate the isolines/isochroones projected to UTM
         :return: NetworkX MultiDiGraph
        """
        return self.GP

    def get_projected_source_nodes(self):
        """
        Returns the added source nodes based on the input location projected to UTM
         :return: geopandas GeoSeries
        """
        return gpd.GeoSeries([Point(i) for i in self.source_points], crs={'init': 'epsg:' + str(self.utm)})

    def get_projected_location(self):
        """
         Returns the input location geometry projected to UTM
          :return: shapely geometry
         """
        return self.geom

    def plot_isolines(self,
                      cmap='plasma_r',
                      plot_nodes=False,
                      plot_source_nodes=True,
                      plot_location=True,
                      plot_samples_points=False,
                      figsize=(15, 15)):
        """

        :param cmap: matplotlib cmap
        :param plot_nodes: boolean
         whether to plot the network nodes
        :param plot_source_nodes: boolean
         whether to plot the added new source network nodes
        :param plot_location: boolean
         whether to plot the input location geometry
        :param plot_samples_points: boolean
         whether to plot the sample location points
        :param figsize: tuple
         figure size

        """

        if cmap is None:
            cmap = 'plasma_r'  # 'YlOrRd'
        cmap = plt.get_cmap(cmap)

        legend_elements = []
        fig, ax = plt.subplots(figsize=figsize)
        if self.metric == 'time':
            txt = ' min.'
        else:
            if self.unit == 'ft':
                txt = ' ft.'
            else:
                txt = ' m.'
        isolines = self.get_isolines()
        isolines.sort_values(self.metric, ascending=False).plot(alpha=0.5,
                                                                ax=ax,
                                                                column=self.metric,
                                                                legend=True,
                                                                categorical=True,
                                                                cmap=cmap)
        d = len(isolines) - 1
        if d <= 0:
            d = 1
        legend_elements += [Patch(facecolor=cmap(ii / d), alpha=0.5,
                                  label=str(i) + txt) for ii, i in enumerate(isolines[self.metric])]

        if plot_location:
            gpd.GeoSeries([self.location]).plot(ax=ax, color='#2bbed8')
            if isinstance(self.location, (Point, MultiPoint)):
                legend_elements += [Line2D([0], [0], color='w', marker='o', markerfacecolor='#2bbed8', markersize=5,
                                           label='Input Location')]
            else:

                legend_elements += [Patch(facecolor='#2bbed8', label='Input Location')]

        if plot_nodes:
            nodes = self.get_nodes()
            nodes.plot(color='#7bfc96', ax=ax)
            legend_elements += [Line2D([0], [0], color='w', marker='o', markerfacecolor='#7bfc96', markersize=5,
                                       label='Original Nodes')]

            new_nodes = pd.concat([gpd.GeoDataFrame([(Point(i), j) for i in self.iso_nodes[j]]) for ji, j in
                                   enumerate(self.values)])
            new_nodes = gpd.GeoDataFrame({'geometry': new_nodes[0],
                                          self.metric: new_nodes[1]}, crs={'init': 'epsg:' + str(self.utm)})

            cmap2 = plt.get_cmap('rainbow')

            new_nodes = new_nodes.to_crs({'init': 'epsg:' + str(self.crs)})
            new_nodes_added = new_nodes[new_nodes.disjoint(MultiPoint(nodes['geometry'].tolist()))]
            new_nodes = new_nodes[new_nodes.intersects(MultiPoint(nodes['geometry'].tolist()))]

            new_nodes = new_nodes.sort_values(self.metric, ascending=False)
            new_nodes_added = new_nodes_added.sort_values(self.metric, ascending=False)

            d = len(new_nodes[self.metric].unique()) - 1
            if d <= 0:
                d = 1

            for ii, i in enumerate(new_nodes[self.metric].unique()):
                df = new_nodes[new_nodes[self.metric] == i]
                if len(df) > 0:
                    c = cmap2(ii / d)
                    df.plot(ax=ax,

                            marker='o',
                            color=c,
                            )

            if len(new_nodes[self.metric].unique()) > 1:

                c_list = [j / d for j in range(0, d + 1)]
                c_list = list(reversed(c_list))
            else:
                c_list = [0]

            for ii, i in enumerate(reversed(new_nodes[self.metric].unique())):
                c = cmap2(c_list[ii])

                legend_elements += [
                    Line2D([0], [0], color='w', marker='o', markerfacecolor=c, markersize=5, alpha=0.5,

                           label='Orig. target Nodes ' + str(i) + txt)]

            d = len(new_nodes_added[self.metric].unique()) - 1
            if d <= 0:
                d = 1

            for ii, i in enumerate(new_nodes_added[self.metric].unique()):
                df = new_nodes_added[new_nodes_added[self.metric] == i]
                if len(df) > 0:
                    c = cmap2(ii / d)
                    df.plot(ax=ax,
                            markersize=150,
                            marker='o',
                            color='none',
                            edgecolor=c,
                            linewidth=2,
                            alpha=0.3)

            if len(new_nodes_added[self.metric].unique()) > 1:

                c_list = [j / d for j in range(0, d + 1)]
                c_list = list(reversed(c_list))
            else:
                c_list = [0]

            for ii, i in enumerate(reversed(new_nodes_added[self.metric].unique())):
                c = cmap2(c_list[ii])

                legend_elements += [Line2D([0], [0], markersize=10,
                                           marker='o',
                                           color=c,
                                           markerfacecolor='w',
                                           alpha=1,
                                           label='Added Target Nodes ' + str(i) + txt)]

        if plot_source_nodes:
            nodes_dict = nx.get_node_attributes(self.GP, 'coords')
            df = gpd.GeoSeries([Point(nodes_dict[i]) for i in nodes_dict.keys() if i in self.new_nodes],
                               crs={'init': 'epsg:' + str(self.utm)})
            df.to_crs({'init': 'epsg:' + str(self.crs)}).plot(ax=ax, markersize=100, color='#ff7272')
            legend_elements += [Line2D([0], [0], markersize=5, color='w',
                                       markerfacecolor='#ff7272', marker='o',
                                       label='Added Source Nodes')]

        if plot_samples_points:
            sample_points = self.get_sample_points()
            if len(sample_points) > 0:
                sample_points.plot(ax=ax, markersize=100, color='#706072')
                legend_elements += [Line2D([0], [0], markersize=5, color='w',
                                           markerfacecolor='#706072', marker='o',
                                           label='Sampled Points')]

        self.get_edges().plot(edgecolor='grey', ax=ax, alpha=0.3)
        legend_elements += [Line2D([0], [0], color='grey', label='Edges')]

        ax.legend(handles=legend_elements, )

        plt.show()

    def _check_params(self):
        if isinstance(self.location,
                      (str,
                       Point,
                       MultiPoint,
                       LineString,
                       MultiLineString,
                       Polygon,
                       MultiPolygon)) == False:
            raise TypeError("""location must be one of :str, Point, MultiPoint, LineString, MultiLineString, Polygon,
                       MultiPolygon but a {} was provided.
                       """.format(type(self.location)))

        if self.network_type not in {'walk', 'drive'}:
            raise ValueError(""" network_type must be 'walk' or 'drive',  but {} was provided.
                       """.format(self.network_type))

        if self.metric not in {'distance', 'time'}:
            raise ValueError("""metric  must be 'distance' or 'time',  but {} was provided.
                       """.format(self.metric))

        if isinstance(self.values,
                      (np.ndarray,
                       list,
                       set,
                       tuple)) == False:
            raise TypeError("""values must be an iterable but a {} was provided.
                       """.format(type(self.values)))

        # self.values = np.array(self.values)

        if isinstance(self.speed,
                      (float, int, np.integer, np.float)) == False:
            raise TypeError("""speed must be a number but a {} was provided.
                       """.format(type(self.speed)))

        if isinstance(self.sample,
                      (float, int, np.integer, np.float)) == False:
            raise TypeError("""sample must be a number but a {} was provided.
                       """.format(type(self.sample)))

        if isinstance(self.crs,
                      (int, np.integer)) == False:
            raise TypeError("""crs must be an integer but a {} was provided.
                       """.format(type(self.crs)))
        if self.convex is None:
            pass
        else:
            if isinstance(self.convex,
                          (bool)) == False and self.knn is not None:
                raise TypeError("""convex must be a boolean or None but a {} was provided.
                           """.format(type(self.convex)))

        if self.convex is None:
            pass
        else:
            if isinstance(self.knn,
                          (int, np.integer)) == False and self.knn is not None:
                raise TypeError("""knn must be an integer or None but a {} was provided.
                           """.format(type(self.knn)))

        if self.address_type not in {'all', 'point'}:
            raise ValueError("""address_type  must be 'all' or 'point',  but {} was provided.
                       """.format(self.address_type))

        if self.unit not in {'ft', 'm'}:
            raise ValueError("""unit  must be 'ft' or 'm',  but {} was provided.
                                   """.format(self.unit))

        if isinstance(self.validate_geometry,
                      (bool)) == False:
            raise TypeError("""validate_geometry must be a boolean but a {} was provided.
                       """.format(type(self.validate_geometry)))

        if isinstance(self.edge_idcol,
                      (str)) == False:
            raise TypeError("""edge_idcol must be a string but a {} was provided.
                       """.format(type(self.edge_idcol)))

        if isinstance(self.fromcol,
                      (str)) == False:
            raise TypeError("""fromcol must be a string but a {} was provided.
                       """.format(type(self.fromcol)))

        if isinstance(self.tocol,
                      (str)) == False:
            raise TypeError("""tocol must be a string but a {} was provided.
                       """.format(type(self.tocol)))

        if isinstance(self.edge_geometrycol,
                      (str)) == False:
            raise TypeError("""edge_geometrycol must be a string but a {} was provided.
                       """.format(type(self.edge_geometrycol)))

        if isinstance(self.weightcol,
                      (str)) == False:
            raise TypeError("""weightcol must be a string but a {} was provided.
                       """.format(type(self.weightcol)))

        if isinstance(self.onewaycol,
                      (str)) == False:
            raise TypeError("""onewaycol  must be a string but a {} was provided.
                       """.format(type(self.onewaycol)))

    def _set_params(self, kwargs):
        for name in kwargs:
            setattr(self, name, kwargs[name])

        self.values = np.array(self.values)

        self._check_params()

        self.iso_nodes = {}
        self.args = kwargs

    def _set_distances(self):
        if self.metric == 'time':
            if self.network_type == 'walk':

                if self.unit == 'ft':
                    if 'speed' not in self.args.keys():
                        self.speed = 276
                    self.distance = self.values.max() * self.speed / 3.28084
                    self.distances = self.values * self.speed / 3.28084
                else:
                    self.distance = self.values.max() * self.speed
                    self.distances = self.values * self.speed

            else:
                if self.unit == 'ft':
                    self.distance = (self.values.max() / 60) * \
                                    self.speed * 1609.344
                    self.distances = (self.values / 60) * self.speed * 1609.344

                else:
                    self.distance = (self.values.max() / 60) * \
                                    self.speed * 1000
                    self.distances = (self.values / 60) * self.speed * 1000

        else:

            if self.unit == 'ft':
                self.distance = self.values.max() * 0.3048
                self.distances = self.values * 0.3048
            else:
                self.distance = self.values.max()
                self.distances = self.values

        if self.unit == 'ft':
            self.unit_dist = self.distance * 3.28084

        else:
            self.unit_dist = self.distance

    def _geocode_address(self):

        geolocator = Nominatim(user_agent="isolines")

        if self.address_type == 'all':
            geocode = geolocator.geocode(self.location, geometry='wkt')
            if geocode is None:
                raise ValueError("Could not locate " + self.location)
            else:
                geocode = wkt.loads(geocode.raw['geotext'])
        else:
            geocode = geolocator.geocode(self.location)
            if geocode is None:
                raise ValueError("Could not locate " + self.location)
            else:
                geocode = Point([geocode.longitude, geocode.latitude])

        self.address = self.location
        self.location = geocode

    def _load_graph_from_geometry(self):

        self.geom, self.utm = convert_to_utm(
            self.location, crs=self.crs)

        polygon = self.geom.buffer(self.distance + 100)
        polygon = project_geom(
            polygon, self.utm, 4326)

        start = time.time()

        print('Downloading graph from OSM using osmnx, maximum distance is {} {}.'.format(
            round(self.unit_dist, 2), self.unit))

        if isinstance(polygon, MultiPolygon):
            graphs = [ox.graph_from_polygon(p, network_type=self.network_type) for p in polygon]
            for g, _ in enumerate(graphs):
                if g == 0:
                    G = graphs[g]
                else:
                    G = nx.compose(G, graphs[g])
        else:

            G = ox.graph_from_polygon(polygon, network_type=self.network_type)

        print('Finished downloading graph from OSM using osmnx, time elapsed: {} seconds.'.format(
            time.time() - start))

        nodes, edges = ox.save_load.graph_to_gdfs(G)

        edges = edges.to_crs({'init': 'epsg:' + str(self.utm)})
        edges['length'] = edges['geometry'].length
        self.edge_idcol = 'osmid'

        return edges

    def _multidigraph_from_edges(self, edges):
        if self.weightcol not in edges.columns and self.weightcol == 'length':
            edges['length'] = edges['geometry'].length
        elif self.weightcol not in edges.columns:
            raise AttributeError("The weight column {} is not in the edges dataframe columns. ".format(self.weightcol))
        self.GP = nx.MultiDiGraph()
        self.GP.add_weighted_edges_from(zip(edges[self.fromcol],
                                            edges[self.tocol],
                                            edges[self.weightcol]))
        edges['geometry'] = edges['geometry'].apply(lambda x: np.array(x.coords))

        edges['coordsf'] = edges['geometry'].apply(lambda x: x[0])
        edges['coordst'] = edges['geometry'].apply(lambda x: x[-1])
        from_dict = dict(zip(edges[self.fromcol], edges['coordsf']))
        to_dict = dict(zip(edges[self.tocol], edges['coordst']))

        idx = set()
        group = pd.DataFrame(
            {'count': edges.groupby([self.fromcol, self.tocol]).size()})
        for c in group['count'].unique():
            df = edges.iloc[~edges.index.isin(idx)].groupby(
                [self.fromcol, self.tocol]).head(1)

            idx.update(df.index)
            edges_names = list(zip(df[self.fromcol],
                                   df[self.tocol],
                                   np.full(len(df), c - 1)
                                   ))

            # length_dict = dict(zip(edges_names, df[self.weightcol]))
            id_dict = dict(zip(edges_names, df[self.edge_idcol]))

            geom_dict = dict(zip(edges_names, df['geometry']))

            nx.set_edge_attributes(self.GP, name='geometry', values=geom_dict)
            nx.set_edge_attributes(
                self.GP, name=self.edge_idcol, values=id_dict)
            if self.onewaycol in edges.columns:
                oneway_dict = dict(zip(edges_names, df[self.onewaycol]))
                nx.set_edge_attributes(
                    self.GP, name=self.onewaycol, values=oneway_dict)

        nx.set_node_attributes(self.GP, name='coords', values=from_dict)
        nx.set_node_attributes(self.GP, name='coords', values=to_dict)

    def _create_graph_from_geometry(self):
        edges = self._load_graph_from_geometry()
        self._multidigraph_from_edges(edges)

    def _create_source_points(self):

        source_points = []

        if self.unit == 'ft' and 'sample' in self.args.keys():
            self.sample = self.sample * 0.3048

        if isinstance(self.geom, (Point)):
            source_points = [self.geom.coords[0]]

        elif isinstance(self.geom, (MultiPoint)):
            source_points = [i.coords[0] for i in self.geom]

        elif isinstance(self.geom, (Polygon, MultiPolygon)):
            if self.sample == 0:
                raise ValueError('Sample has to be larger than 0')

            else:
                if isinstance(self.geom, (Polygon)):
                    ext = self.geom.exterior

                    source_points = [ext.interpolate(d).coords[0] for d in range(
                        0, int(ext.length) + 1, int(self.sample))]
                else:
                    for g in self.geom:
                        ext = g.exterior
                        source_points += [ext.interpolate(d).coords[0] for d in
                                          range(0, int(ext.length) + 1, self.sample)]

        elif isinstance(self.geom, (LineString, MultiLineString)):
            if self.sample == 0:
                raise ValueError('Sample has to be larger than 0')


            else:
                if isinstance(self.geom, (LineString)):

                    source_points = [self.geom.interpolate(d).coords[0] for d in range(0,
                                                                                       int(
                                                                                           self.geom.length) + 1,
                                                                                       self.sample)]
                else:
                    for g in self.geom:
                        source_points += [g.interpolate(d).coords[0] for d in range(
                            0, int(g.length) + 1, int(self.sample))]

        return source_points

    def _insert_new_source(self, tp, rt_idx, edges_dict):
        edges_to_remove = []
        d = np.inf
        closest = list(rt_idx.nearest(Point(tp).buffer(10).bounds, 1))
        for idx in closest:
            s = edges_dict[idx]
            dd = LineString(self.GP.edges[s]['geometry']).distance(Point(tp))
            if dd < d:
                idxi = idx
                d = dd

        source_edge = edges_dict[idxi]

        source_line_coords = self.GP.edges[source_edge]['geometry']

        source_line = LineString(source_line_coords)

        dist_from = source_line.project(Point(tp))
        dist_to = source_line.length - dist_from
        sp = source_line.interpolate(dist_from)
        source_node_coords = sp.coords[0]

        source_node = random_string(15)
        nodes_dict = {source_node: source_node_coords}
        self.new_nodes.append(source_node)
        start_node = source_edge[0]
        end_node = source_edge[1]
        num_node = source_edge[0]
        oneway = False
        if self.onewaycol in self.GP.edges[edges_dict[0]].keys():
            oneway = self.GP.edges[source_edge][self.onewaycol]

        if self.onewaycol in self.GP.edges[edges_dict[0]].keys() and self.network_type == 'drive' \
                and oneway == True:

            lines_dict = {
                (start_node, source_node, 0):
                    cut_segment_at_point(source_line_coords, source_node_coords),

                (source_node, end_node, 0):
                    cut_segment_at_point(source_line_coords[::-1], source_node_coords)[::-1],

            }

            id_dict = {(start_node, source_node, 0): random_string(15),
                       (source_node, end_node, 0): random_string(15),
                       }

            oneway_dict = {(start_node, source_node, 0): oneway,
                           (source_node, end_node, 0): oneway,
                           }

            self.GP.add_weighted_edges_from([(start_node, source_node, dist_from),
                                             (source_node, end_node, dist_to)])
            # remove original edges
            edges_to_remove.append((start_node, end_node))

        else:
            ## manage circular edges:

            if start_node == end_node:

                lines_dict = {
                    (start_node, source_node, 0):
                        cut_segment_at_point(source_line_coords, source_node_coords),
                    (source_node, start_node, 0):
                        cut_segment_at_point(source_line_coords, source_node_coords)[::-1],
                    (source_node, end_node, 1):
                        cut_segment_at_point(source_line_coords[::-1], source_node_coords)[::-1],
                    (end_node, source_node, 1):
                        cut_segment_at_point(source_line_coords[::-1], source_node_coords),
                }

                id_dict = {(start_node, source_node, 0): random_string(15),
                           (source_node, start_node, 0): random_string(15),
                           (source_node, end_node, 1): random_string(15),
                           (end_node, source_node, 1): random_string(15),
                           }

                # add edges
                self.GP.add_weighted_edges_from([(start_node, source_node, dist_from),
                                                 (source_node, start_node, dist_from),
                                                 (source_node, end_node, dist_to),
                                                 (end_node, source_node, dist_to)])
                # remove original edges

                edges_to_remove.append((start_node, end_node, num_node))
                edges_to_remove.append((end_node, start_node, num_node))


            else:

                lines_dict = {
                    (start_node, source_node, 0):
                        cut_segment_at_point(source_line_coords, source_node_coords),
                    (source_node, start_node, 0):
                        cut_segment_at_point(source_line_coords, source_node_coords)[::-1],
                    (source_node, end_node, 0):
                        cut_segment_at_point(source_line_coords[::-1], source_node_coords)[::-1],
                    (end_node, source_node, 0):
                        cut_segment_at_point(source_line_coords[::-1], source_node_coords),
                }

                id_dict = {(start_node, source_node, 0): random_string(15),
                           (source_node, start_node, 0): random_string(15),
                           (source_node, end_node, 0): random_string(15),
                           (end_node, source_node, 0): random_string(15),
                           }

                oneway_dict = {(start_node, source_node, 0): oneway,
                               (source_node, start_node, 0): oneway,
                               (source_node, end_node, 0): oneway,
                               (end_node, source_node, 0): oneway,
                               }

                # add edges
                self.GP.add_weighted_edges_from([(start_node, source_node, dist_from),
                                                 (source_node, start_node, dist_from),
                                                 (source_node, end_node, dist_to),
                                                 (end_node, source_node, dist_to)])
                # remove original edges

                edges_to_remove.append((start_node, end_node, num_node))
                edges_to_remove.append((end_node, start_node, num_node))

        nx.set_edge_attributes(
            self.GP, name='geometry', values=lines_dict)
        nx.set_edge_attributes(
            self.GP, name=self.edge_idcol, values=id_dict)

        nx.set_edge_attributes(
            self.GP, name=self.onewaycol, values=oneway_dict)
        nx.set_node_attributes(self.GP, name='coords', values=nodes_dict)

        return edges_to_remove

    def _add_source_nodes(self):

        self.new_nodes = []
        edges_dict = dict(zip(range(len(self.GP.edges)), self.GP.edges))
        rt_idx = rtree.index.Index()
        [rt_idx.insert(i[0], (i[1])) for i in \
         [(i, coords_bounds(self.GP.edges[edges_dict[i]]['geometry'])) for i in edges_dict.keys()]]
        edges_to_remove = []

        for tpi, tp in enumerate(self.source_points):
            # result = delayed(self._insert_new_source)( tp, rt_idx, edges_dict)
            # edges_to_remove.append(result)
            result = self._insert_new_source(tp, rt_idx, edges_dict)
            edges_to_remove += result

        # results = delayed()(edges_to_remove)
        # edges_to_remove = results.compute()

        self.GP.remove_edges_from(edges_to_remove)

    def _compute_new_nodes(self, tpi, distance):
        nodes_coords = np.empty([0, 2])
        pred, dist = nx.dijkstra_predecessor_and_distance(
            self.GP, tpi, weight='weight')
        lower_weights = set([i for i in dist if dist[i] <= distance])
        lower_weights_pred = [pred[i] for i in dist if dist[i] > distance]
        lower_weights_pred = set([item for sublist in lower_weights_pred for item in sublist]) \
            .intersection(lower_weights)

        for node in lower_weights_pred:
            neighs = dict(self.GP[node])
            for i in neighs.keys():
                if i not in lower_weights:
                    ii = dict(neighs[i])
                    for j in ii.keys():
                        vertices = ii[j]['geometry']

                        new_line = LineString(vertices)

                        end_point_dist = distance - dist[node]

                        new_end_points = new_line.interpolate(end_point_dist).coords[0]

                        nodes_coords = np.append(nodes_coords, new_end_points)

        coords_dict = nx.get_node_attributes(self.GP, 'coords')
        filtered_nodes_coords = np.array(
            [coords_dict[i] for i in lower_weights])

        nodes_coords = np.append(nodes_coords,
                                 filtered_nodes_coords)
        return nodes_coords

    def _concat_nodes(self, nodes_coords):
        nodes_coords = np.hstack(nodes_coords)
        nodes_coords = nodes_coords.reshape(
            (int(nodes_coords.shape[0] / 2), 2))
        nodes_coords = np.unique(nodes_coords, axis=0)

        return nodes_coords

    def _extract_target_nodes_coords(self, idx, distances_arr, values_arr):
        distance = distances_arr[idx]

        # get the graph weakly connected components
        wcc = list(nx.weakly_connected_components(self.GP))

        nodes_coords = []

        if self.par == True:
            list_ = []
            for cc in wcc:
                cc_list = []
                tpis = set(self.new_nodes).intersection(cc)
                if len(tpis) > 0:
                    for tpi in tpis:
                        result = delayed(self._compute_new_nodes)(tpi, distance)
                        cc_list.append(result)

                    concat = delayed(self._concat_nodes)(cc_list)
                    results = concat.compute()
                    self.results2 = concat
                    list_.append(results)
                    nodes_coords.append(results)
            self.iso_nodes[values_arr[idx]] = np.vstack(list_)


        else:
            list_ = []
            for cc in wcc:
                cc_list = []
                tpis = set(self.new_nodes).intersection(cc)
                if len(tpis) > 0:
                    for tpi in tpis:
                        result = self._compute_new_nodes(tpi, distance)
                        cc_list.append(result)

                    results = self._concat_nodes(cc_list)
                    list_.append(results)
                    nodes_coords.append(results)
            self.iso_nodes[values_arr[idx]] = np.vstack(list_)

        return nodes_coords

    def _create_iso(self, coords, knn):

        if self.convex is None:
            if self.network_type == 'walk':
                self.convex = False
            else:
                self.convex = True
        if self.convex:
            return MultiPoint(coords).convex_hull
        else:

            if knn is None:
                if coords.shape[0] < 15:
                    knn = 10
                elif coords.shape[0] < 20:
                    knn = int(coords.shape[0] / 2)
                elif coords.shape[0] < 40:
                    knn = int(coords.shape[0] / 3)
                else:
                    knn = int(coords.shape[0] / 5)

                if knn < 3:
                    knn = 3

            return ConcaveHull(coords, knn).concave_hull


class GpdIsolines(OsmIsolines):
    """

   :param location: an address string, or a shapely geometry
    the source location

   :param edges: geopandas GeoDataFrame
    geopandas GeoDataFrame with the network's edges

   :param edge_idcol: string
    id column name of a edges geodataframe (default is 'id')

   :param fromcol: string
    a 'from' node column name of an edges geodataframe (default is 'u')

   :param tocol: string
    a 'to' node column name of an edges geodataframe (default is 'v')

   :param edge_geometrycol: string
    geometry column name of the edges geodataframe  (default is 'geometry')

   :param weightcol: string
    a weight column name of the edges geodataframe (default is 'length')

   :param onewaycol: string
    a oneway column name (specifies if a street is oneway - relevant only if network_type is
    'drive') of the edges geodataframe (default is 'oneway')

   :param network_type: string
    {'walk', 'drive'}
    the type of mobility, default is ``'walk'``

   :param metric: string
    {'distance', 'time'}
    determines if the output will be isolines based on distances or isochrones based on time, default is ``'distance'``

   :param values: list or numpy array
    a list or array of time values in minutes or distance values in meters or feet (if unit = True), default is [500]

   :param speed: number
    if the ``metric`` parameter is set to ``'time'`` this will be the speed of walking im meters per
    minute or feet per minute (if ``unit =='ft'``), if network_type is set to ``'drive'``. this
    will be kilometers per hour or mile per hour (if ``unit =='ft'``)

   :param address_type: string
    {'all', 'point'}
    for a geocoded address this determines if to get back a point address or multipoint,
    polygon, multipolygon, linestring, multilinestring when relevant

   :param sample: number
    if a geocoded location geometry is polygon, multipolygon, linestring, multilinestring it is
    sampled in evenly spaced intervals (the boundary of polygons or the line itself) and sample points from it are
    added as nodes in the network, the default sample interval length is 100 m. or 330 ft.

   :param crs: int
    the crs of the input location and other inputs if relevant (input edges dataframe or
    networkX graph) - default is 4326

   :param convex: boolean
    whether the isoline/isochrone will be constructed using a convex or concave hull
    algorithms - default is ``False`` if ``network_type == 'walk'`` and ``True`` if ``network type == 'drive'``

   :param knn: int
    if the isoline/isochrone boundary is concave hull - the hull algorithm is based on
    k-nearest-neighbors (knn) heuristic. If no knn is provided, knn will be determined based on the number of points.
    smaller knn - will produce a more and elaborate isoline/isochrone

   :param unit: str
    {'all', 'point'}
    whether measurement are in meters or feet, default is meters

   :param validate_geometry: boolean
    if ``True`` the isolines topological validity and counter clock wise (ccw) coordinate sequence direction if
    ensured. default is ``True``

   :param par: boolean
    whether data processing is parallelized, default is ``True``

   :param smooth: boolean
    whether to smooth the isoline boundary, default is ``False``

   :param sigma: number
    standard deviation for Gaussian kernel if ``smooth ==True``, default is 2

   :param num_points_factor: int
    the number of points determines the density of vertices for smoothing if ``smooth ==True``, default is 5



   :return: None



        """

    def __init__(self,
                 location,
                 **kwargs):

        super().__init__(location,
                         **kwargs)

    def _set_params(self, kwargs):
        if 'edges' not in kwargs.keys():
            raise ValueError("You must supply an edges GeoDatatFrame")

        for name in kwargs:
            setattr(self, name, kwargs[name])

        self.values = np.array(self.values)

        self._check_params()

        self.iso_nodes = {}
        self.args = kwargs

    def _load_graph_from_geometry(self):

        self.geom, self.utm = convert_to_utm(
            self.location, crs=self.crs)
        if self.edge_geometrycol != 'geometry':
            self.edges['geometry'] = self.edges[self.edge_geometrycol]

        edges = self.edges.to_crs({'init': 'epsg:' + str(self.utm)})
        del self.edges
        edges = edges.reset_index(drop=True)

        edges = edges.iloc[list(edges.sindex.intersection(self.geom.buffer(self.distance + 100).bounds))]

        edges['length'] = edges['geometry'].length
        if self.edge_idcol not in edges.columns:
            edges = edges.reset_index(drop=True)
            edges[self.edge_idcol] = edges.index

        return edges


class NxIsolines(GpdIsolines):
    """

   :param location: an address string, or a shapely geometry
    the source location

   :param graph: NetworkX graph

   :param edge_idcol: string
    id attribute name of the graph (default is 'id')

   :param fromcol: string
    a 'from' attribute name of the graph (default is 'u')

   :param tocol: string
    a 'to' node attribute name of an the graph (default is 'v')

   :param edge_geometrycol: string
    geometry attribute name of the graph  (default is 'geometry')

   :param weightcol: string
    a weight attribute name of the graph (default is 'length')

   :param onewaycol: string
    a oneway attribute name (specifies if a street is oneway - relevant only if network_type is
    'drive') of the edges geodataframe (default is 'oneway')

   :param network_type: string
    {'walk', 'drive'}
    the type of mobility, default is ``'walk'``

   :param metric: string
    {'distance', 'time'}
    determines if the output will be isolines based on distances or isochrones based on time, default is ``'distance'``

   :param values: list or numpy array
    a list or array of time values in minutes or distance values in meters or feet (if unit = True), default is [500]

   :param speed: number
    if the ``metric`` parameter is set to ``'time'`` this will be the speed of walking im meters per
    minute or feet per minute (if ``unit =='ft'``), if network_type is set to ``'drive'``. this
    will be kilometers per hour or mile per hour (if ``unit =='ft'``)

   :param address_type: string
    {'all', 'point'}
    for a geocoded address this determines if to get back a point address or multipoint,
    polygon, multipolygon, linestring, multilinestring when relevant

   :param sample: number
    if a geocoded location geometry is polygon, multipolygon, linestring, multilinestring it is
    sampled in evenly spaced intervals (the boundary of polygons or the line itself) and sample points from it are
    added as nodes in the network, the default sample interval length is 100 m. or 330 ft.

   :param crs: int
    the crs of the input location and other inputs if relevant (input edges dataframe or
    networkX graph) - default is 4326

   :param convex: boolean
    whether the isoline/isochrone will be constructed using a convex or concave hull
    algorithms - default is ``False`` if ``network_type == 'walk'`` and ``True`` if ``network type == 'drive'``

   :param knn: int
    if the isoline/isochrone boundary is concave hull - the hull algorithm is based on
    k-nearest-neighbors (knn) heuristic. If no knn is provided, knn will be determined based on the number of points.
    smaller knn - will produce a more and elaborate isoline/isochrone

   :param unit: str
    {'all', 'point'}
    whether measurement are in meters or feet, default is meters

   :param validate_geometry: boolean
    if ``True`` the isolines topological validity and counter clock wise (ccw) coordinate sequence direction if
    ensured. default is ``True``

   :param par: boolean
    whether data processing is parallelized, default is ``True``

   :param smooth: boolean
    whether to smooth the isoline boundary, default is ``False``

   :param sigma: number
    standard deviation for Gaussian kernel if ``smooth ==True``, default is 2

   :param num_points_factor: int
    the number of points determines the density of vertices for smoothing if ``smooth ==True``, default is 5



   :return: None
    """

    def __init__(self,
                 location,
                 **kwargs):

        super().__init__(location,
                         **kwargs)

    def _set_params(self, kwargs):
        if 'graph' not in kwargs.keys():
            raise ValueError("You must supply a NetworkX graph")

        for name in kwargs:
            setattr(self, name, kwargs[name])
        self.G = self.graph

        del self.graph

        self.values = np.array(self.values)

        self._check_params()

        self.iso_nodes = {}
        self.args = kwargs

        self.fromcol = 'source'
        self.tocol = 'target'
        self.edges = nx.to_pandas_edgelist(self.G)
        if isinstance(self.edges.loc[0, self.edge_geometrycol], np.ndarray):
            self.edges = gpd.GeoDataFrame(
                self.edges, geometry=[LineString(i) for i in self.edges[self.edge_geometrycol]])
        elif isinstance(self.edges.loc[0, self.edge_geometrycol], (LineString, MultiLineString)):

            self.edges = gpd.GeoDataFrame(
                self.edges, geometry=self.edges[self.edge_geometrycol])
        else:
            raise TypeError("""geometry must be either a 2d numpy array of line coordinate sequence or a shapely 
            LineString/MultiLineString but instead a {} was provided""".format(
                type(self.edges.loc[0, self.edge_geometrycol])))
        self.edges.crs = {'init': 'epsg:' + str(self.crs)}

        if self.edge_idcol not in self.edges.columns:
            self.edges = self.edges.reset_index(drop=True)
            self.edges[self.edge_idcol] = self.edges.index


def isolines(location, **kwargs):
    """

   :param location: an address string, or a shapely geometry
    the source location

   :param graph: NetworkX graph

   :param edges: geopandas GeoDataFrame

   :param edge_idcol: string
    id attribute/column name of the graph/edges GeoDataFrame (default is 'id')

   :param fromcol: string
    a 'from' attribute/column name of the graph/edges GeoDataFrame (default is 'u')

   :param tocol: string
    a 'to' nodeattribute/column name of an the graph/edges GeoDataFrame (default is 'v')

   :param edge_geometrycol: string
    geometry attribute/column name of the graph/edges GeoDataFrame (default is 'geometry')

   :param weightcol: string
    a weight attribute/column name of the graph/edges GeoDataFrame (default is 'length')

   :param onewaycol: string
    a oneway attribute/column name (specifies if a street is oneway - relevant only if network_type is
    'drive') of the graph/edges GeoDataFrame (default is 'oneway')

   :param network_type: string
    {'walk', 'drive'}
    the type of mobility, default is ``'walk'``

   :param metric: string
    {'distance', 'time'}
    determines if the output will be isolines based on distances or isochrones based on time, default is ``'distance'``

   :param values: list or numpy array
    a list or array of time values in minutes or distance values in meters or feet (if unit = True), default is [500]

   :param speed: number
    if the ``metric`` parameter is set to ``'time'`` this will be the speed of walking im meters per
    minute or feet per minute (if ``unit =='ft'``), if network_type is set to ``'drive'``. this
    will be kilometers per hour or mile per hour (if ``unit =='ft'``)

   :param address_type: string
    {'all', 'point'}
    for a geocoded address this determines if to get back a point address or multipoint,
    polygon, multipolygon, linestring, multilinestring when relevant

   :param sample: number
    if a geocoded location geometry is polygon, multipolygon, linestring, multilinestring it is
    sampled in evenly spaced intervals (the boundary of polygons or the line itself) and sample points from it are
    added as nodes in the network, the default sample interval length is 100 m. or 330 ft.

   :param crs: int
    the crs of the input location and other inputs if relevant (input edges dataframe or
    networkX graph) - default is 4326

   :param convex: boolean
    whether the isoline/isochrone will be constructed using a convex or concave hull
    algorithms - default is ``False`` if ``network_type == 'walk'`` and ``True`` if ``network type == 'drive'``

   :param knn: int
    if the isoline/isochrone boundary is concave hull - the hull algorithm is based on
    k-nearest-neighbors (knn) heuristic. If no knn is provided, knn will be determined based on the number of points.
    smaller knn - will produce a more and elaborate isoline/isochrone

   :param unit: str
    {'all', 'point'}
    whether measurement are in meters or feet, default is meters

   :param validate_geometry: boolean
    if ``True`` the isolines topological validity and counter clock wise (ccw) coordinate sequence direction if
    ensured. default is ``True``

   :param par: boolean
    whether data processing is parallelized, default is ``True``

   :param smooth: boolean
    whether to smooth the isoline boundary, default is ``False``

   :param sigma: number
    standard deviation for Gaussian kernel if ``smooth ==True``, default is 2

   :param num_points_factor: int
    the number of points determines the density of vertices for smoothing if ``smooth ==True``, default is 5

   :param edge_idcol: string
    id column of a edges geodataframe if an edges geodataframe is provided (default is
     'id')

   :return: geopandas GeoDataFrame of isolines/isochrones
    """

    if 'edges' in kwargs.keys():
        return GpdIsolines(location, **kwargs).get_isolines()
    elif 'graph' in kwargs.keys():
        return NxIsolines(location, **kwargs).get_isolines()
    else:
        return OsmIsolines(location, **kwargs).get_isolines()

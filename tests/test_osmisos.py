from shapely import wkt
from isolines.utils import *
import isolines as il
import math

def test_make_valid_polygon():
    poly = Point([0, 0]).buffer(50)
    poly = Polygon(poly.exterior.coords[:] + [[49, 2], [51, -2], [50, 0]])
    poly2 = Point([25, 25]).buffer(50)
    hole = Point([0, 0]).buffer(25)
    hole2 = Point([25, 25]).buffer(25)
    geom = MultiPolygon([Polygon(poly.exterior, [hole.exterior]),
                         Polygon(poly2.exterior, [hole2.exterior])])

    result = make_valid_polygon(geom)

    for i in result:
        assert i.is_valid


def test_location_input():

    try:
        il.OsmIsolines(4)
        assert False
    except TypeError:
        assert True

def test_crs_input():

    try:
        convert_to_utm(Point([0,0]), 00000)
        assert False
    except ValueError:
        assert True

def test_convert_to_utm():

    p, utm = convert_to_utm(Point([3, 0]), 4326)
    assert math.isclose(p.coords[0][0], 500000)
    assert math.isclose(p.coords[0][1], 0)
    assert  utm == 32631

    p, utm = convert_to_utm(Point([3, -0.00000000001]), 4326)
    assert math.isclose(p.coords[0][0], 500000)
    assert math.isclose(p.coords[0][1], 10000000)

    assert utm == 32731

    p, utm = convert_to_utm(Point([0,0]), 32636)
    assert utm == 32636

    p, utm = convert_to_utm(Point([0, 0]), 32736)
    assert utm == 32736

    p, utm = convert_to_utm(Point([0, 0]),2039)
    assert utm == 32636


def test_is_ccw():
    mpoly = 'MULTIPOLYGON (((667768.2818781252 3548514.477686874, 667334.8150165473 3548556.992710626, 667235.0777269548 3548617.734395519, 666834.8150165473 3548656.992710626, 666702.0430262301 3548737.853083225, 666493.7239733685 3548978.472464215, 666426.9510792443 3549136.252650722, 666275.1079055682 3549460.079263736, 666262.1208103356 3549499.01289328, 666077.3820849422 3549925.35078082, 666070.2926062428 3549996.434951082, 666097.9917770844 3550091.779910524, 666294.5489067507 3550769.423171774, 666403.8706166596 3550993.006219771, 666887.0932789206 3551296.862138772, 667214.1617881672 3551443.037086175, 667586.4814016479 3551490.367944148, 667680.7941483615 3551500.293769633, 667683.0227187943 3551499.875369023, 667728.7917326988 3551491.261606091, 667989.1036878376 3551377.988873653, 668086.4814016479 3551390.367944148, 668180.7941483615 3551400.293769633, 668183.0227187943 3551399.875369023, 668228.7917326988 3551391.261606091, 668541.0251886041 3551255.395634471, 668844.3821052716 3550928.867344437, 669353.567889452 3550617.198068652, 669604.9869401017 3549795.315142757, 669192.24521816 3549495.585266457, 668723.3477656682 3549049.894101015, 668579.6068137914 3548970.384763314, 668127.2354101954 3548642.754718971, 667768.2818781252 3548514.477686874), (666676.1030731541 3550261.073593734, 666707.8310840167 3550370.457938075, 666698.0582113585 3550372.995574201, 666656.7181028681 3550361.841724961, 666676.1030731541 3550261.073593734), (667870.0941545123 3550510.46245681, 667836.9665606197 3550461.691891549, 667636.8941544265 3550430.822491099, 667590.0527067006 3550515.865839677, 667584.7572832993 3550523.744795364, 667510.9833372888 3550546.970432545, 667343.3271147688 3550235.274872193, 667198.0582113585 3550272.995574201, 667156.7181028681 3550261.841724961, 667177.0830697499 3550155.979315843, 667229.6980328183 3549909.538396385, 667305.4026837416 3549775.85955404, 667424.5592640812 3549547.70189763, 667551.1897058861 3549332.605444953, 667624.0805913871 3549369.072791826, 667818.9040590749 3549449.000046127, 667962.2411685224 3549623.816262954, 667987.4892217622 3549729.719455983, 667959.2242101671 3549971.470488441, 668241.1116720552 3550060.143698678, 668298.6632920186 3550061.091587366, 668103.6206043581 3550363.217212968, 667906.3279906875 3550490.69571075, 667870.0941545123 3550510.46245681)), ((668368.2818781252 3549114.477686874, 667934.8150165473 3549156.992710626, 667835.0777269548 3549217.734395519, 667434.8150165473 3549256.992710626, 667302.0430262301 3549337.853083225, 667093.7239733685 3549578.472464215, 667026.9510792443 3549736.252650722, 666875.1079055682 3550060.079263736, 666862.1208103356 3550099.01289328, 666677.3820849422 3550525.35078082, 666670.2926062428 3550596.434951082, 666697.9917770844 3550691.779910524, 666894.5489067507 3551369.423171774, 667003.8706166596 3551593.006219771, 667487.0932789206 3551896.862138772, 667814.1617881672 3552043.037086175, 668186.4814016479 3552090.367944148, 668280.7941483615 3552100.293769633, 668283.0227187943 3552099.875369023, 668328.7917326988 3552091.261606091, 668589.1036878376 3551977.988873653, 668686.4814016479 3551990.367944148, 668780.7941483615 3552000.293769633, 668783.0227187943 3551999.875369023, 668828.7917326988 3551991.261606091, 669141.0251886041 3551855.395634471, 669444.3821052716 3551528.867344437, 669953.567889452 3551217.198068652, 670204.9869401017 3550395.315142757, 669792.24521816 3550095.585266457, 669323.3477656682 3549649.894101015, 669179.6068137914 3549570.384763314, 668727.2354101954 3549242.754718971, 668368.2818781252 3549114.477686874), (667276.1030731541 3550861.073593734, 667307.8310840167 3550970.457938075, 667298.0582113585 3550972.995574201, 667256.7181028681 3550961.841724961, 667276.1030731541 3550861.073593734), (668470.0941545123 3551110.46245681, 668436.9665606197 3551061.691891549, 668236.8941544265 3551030.822491099, 668190.0527067006 3551115.865839677, 668184.7572832993 3551123.744795364, 668110.9833372888 3551146.970432545, 667943.3271147688 3550835.274872193, 667798.0582113585 3550872.995574201, 667756.7181028681 3550861.841724961, 667777.0830697499 3550755.979315843, 667829.6980328183 3550509.538396385, 667905.4026837416 3550375.85955404, 668024.5592640812 3550147.70189763, 668151.1897058861 3549932.605444953, 668224.0805913871 3549969.072791826, 668418.9040590749 3550049.000046127, 668562.2411685224 3550223.816262954, 668587.4892217622 3550329.719455983, 668559.2242101671 3550571.470488441, 668841.1116720552 3550660.143698678, 668898.6632920186 3550661.091587366, 668703.6206043581 3550963.217212968, 668506.3279906875 3551090.69571075, 668470.0941545123 3551110.46245681)))'

    mpoly = wkt.loads(mpoly)
    mpoly = make_ccw(mpoly)

    for i in mpoly:
        assert i.exterior.is_ccw
        for n in i.interiors:
            assert n.is_ccw

def test_init_params():
    try:
        il.GpdIsolines('times square', values=[500])
        assert False
    except ValueError:
        assert True

    try:
        il.NxIsolines('times square', values=[500])
        assert False
    except ValueError:
        assert True

def test_geocode():
    try:
        iso = il.OsmIsolines(';oiruetoueo', values=[500])
        assert False
    except ValueError:
        assert True

    iso = il.OsmIsolines('times square', values=[500], address_type='point')

    assert isinstance(iso.location, Point) == True
    iso2 = il.OsmIsolines('times square', values=[500], address_type='all')

    assert isinstance(iso2.location, Point) == False

    iso3 = il.GpdIsolines('times square', edges=iso2.get_edges(), fromcol='source',
                                        tocol='target', values=[500], address_type='point')

    assert isinstance(iso3.location, Point) == True
    iso4 = il.GpdIsolines('times square', edges=iso3.get_edges(), fromcol='source',
                                        tocol='target', values=[500], address_type='all')

    assert isinstance(iso4.location, Point) == False

    iso5 = il.NxIsolines('times square', graph=iso4.get_graph(), fromcol='source',
                                       tocol='target', values=[500], address_type='point')

    assert isinstance(iso5.location, Point) == True
    iso6 = il.NxIsolines('times square', graph=iso5.get_graph(), fromcol='source',
                                       tocol='target', values=[500], address_type='all')

    assert isinstance(iso6.location, Point) == False

def test_distances():
    poly = 'POLYGON ((34.7789012 32.0723507, 34.7791154 32.0722353, 34.7793334 32.0721058, 34.7795058 32.0720034, 34.7798462 32.0718012, 34.7798952 32.0717779, 34.7799146 32.0717792, 34.7799342 32.0717872, 34.7799749 32.0718304, 34.7800152 32.0718938, 34.7800494 32.0719833, 34.7800951 32.0721526, 34.7801415 32.0721419, 34.7802657 32.0725926, 34.7802395 32.0725983, 34.7802115 32.0726044, 34.7803126 32.0729483, 34.7802872 32.0729538, 34.7800336 32.0730091, 34.7794776 32.0731262, 34.7793901 32.0731446, 34.7793928 32.0731335, 34.7794169 32.0730342, 34.7792506 32.0724367, 34.7792089 32.0724654, 34.7791549 32.0724809, 34.7791208 32.0724538, 34.7790698 32.0724201, 34.7790322 32.0724077, 34.7789862 32.0723998, 34.7789538 32.0724051, 34.7789236 32.0724102, 34.7789209 32.0724021, 34.7789012 32.0723507), (34.7794869 32.0722471, 34.7796115 32.0726659, 34.779777 32.0726305, 34.7796524 32.0722117, 34.7794869 32.0722471), (34.7796538 32.0728024, 34.7796934 32.0729576, 34.7798518 32.0729285, 34.7798122 32.0727734, 34.7796538 32.0728024), (34.7796569 32.0720271, 34.7796577 32.0720468, 34.7796659 32.0720653, 34.7796805 32.0720807, 34.7797 32.0720914, 34.7797226 32.0720964, 34.7797458 32.0720951, 34.7797674 32.0720877, 34.7797852 32.0720749, 34.7797973 32.0720581, 34.7798026 32.0720388, 34.7797999 32.0720176, 34.7797889 32.0719984, 34.7797709 32.0719834, 34.779748 32.0719746, 34.7797228 32.0719728, 34.7796985 32.0719784, 34.7796779 32.0719907, 34.7796633 32.0720081, 34.7796569 32.0720271))'

    poly = wkt.loads(poly)
    iso = il.OsmIsolines(poly, metric='time', values=[5, 10], sample=100, speed=10, network_type='drive')
    assert math.isclose(iso.distances[0], 10 / 12 * 1000)
    assert math.isclose(iso.distances[1], 10 / 6 * 1000)
    assert iso.utm == 32636

    # if added isolines are all above the network max distance

    iso.add_isolines([15])

    iso.add_isolines([7.5, 2.5, 20])
    # ensures isolines and index are sorted in ascending order
    assert all(iso.isos['time'][i] <= iso.isos['time'][i + 1] for i in range(len(iso.isos['time']) - 1))
    assert all(iso.isos.index[i] <= iso.isos.index[i + 1] for i in range(len(iso.isos.index) - 1))

    iso2 = il.GpdIsolines(poly,
                                        edges=iso.get_edges(),
                                        fromcol='source',
                                        tocol='target',
                                        metric='time',
                                        values=[5, 10],
                                        sample=100, speed=10,
                                        network_type='drive')
    assert math.isclose(iso2.distances[0], 10 / 12 * 1000)
    assert math.isclose(iso2.distances[1], 10 / 6 * 1000)
    assert iso2.utm == 32636

    bounds = poly.bounds
    utm = get_utm_zone(bounds[0], bounds[1])
    poly_utm = project_geom(poly, 4326, str(utm))

    iso3 = il.NxIsolines(poly_utm, graph=iso.GP, metric='time', values=[5, 10],
                                       sample=100, speed=10, network_type='drive', crs = iso.utm)

    assert math.isclose(iso3.distances[0], 10 / 12 * 1000)
    assert math.isclose(iso3.distances[1], 10 / 6 * 1000)
    assert iso3.utm == 32636

    iso4 = il.NxIsolines(poly_utm, graph=iso.GP, metric='time', values=[5, 10, 60],
                                       sample=100, speed=10, network_type='drive', crs=iso.utm, unit='ft')
    assert round(iso4.distances[0]) == round(10 * 1609.344 / (60 / 5))
    assert round(iso4.distances[1]) == round(10 * 1609.344 / (60 / 10))
    assert round(iso4.distances[2]) == round(10 * 1609.344)
    assert iso4.utm == 32636

    # distance - feet to meter

    iso = il.OsmIsolines(poly, metric='distance', values=[100, 200], unit='ft', network_type='walk')

    assert round(iso.distances[0], 2) == 30.48
    assert round(iso.distances[1], 2) == 60.96

    iso.add_isolines([150])

    assert round(iso.distances[1], 2) == 45.72

    iso = il.OsmIsolines(poly, metric='time', values=[15], unit='ft', network_type='walk')

    assert 2 * int(iso.distances[0] / 2) == 1260

    iso.add_isolines([10, 5])

    assert 2 * int(iso.distances[0] / 2) == 1260 / 3
    assert 2 * int(iso.distances[1] / 2) == 1260 / 1.5

    iso = il.OsmIsolines(poly, metric='time', values=[10], unit='ft', network_type='drive', speed=20)

    assert 10 * int(iso.distances[0] / 10) == 10 * int(5363.3 / 10)

    iso.add_isolines([2.5, 5, 7.5])

    assert round(10 * int(iso.distances[0] / 10)) == round(10 * int(5363.3 / 10) / 4)
    assert round(10 * int(iso.distances[1] / 10)) == round(10 * int(5363.3 / 10) / 2)
    assert round(10 * int(iso.distances[2] / 10)) == round(10 * int(5363.3 / 10) / 1.3333)

def test_data_types():
    isocrones = il.OsmIsolines('times square',
                                     network_type='walk',
                                     metric='time',
                                     values=[5],
                                     unit='ft',
                                     address_type='point')
    assert isinstance(isocrones.iso_nodes, dict)
    assert isinstance(isocrones.values, np.ndarray)


def test_short_distances():
    multipolygon = 'MULTIPOLYGON (((34.7889012 32.0823507, 34.7891154 32.0822353, 34.7893334 32.0821058, 34.7895058 32.0820034, 34.7898462 32.0818012, 34.7898952 32.0817779, 34.7899146 32.0817792, 34.7899342 32.0817872, 34.7899749 32.0818304, 34.7900152 32.0818938, 34.7900494 32.0819833, 34.79009509999999 32.0821526, 34.7901415 32.0821419, 34.7902657 32.0825926, 34.7902395 32.0825983, 34.7902115 32.0826044, 34.7903126 32.0829483, 34.79028719999999 32.0829538, 34.7900336 32.0830091, 34.7894776 32.0831262, 34.7893901 32.0831446, 34.78939279999999 32.0831335, 34.7894169 32.0830342, 34.7892506 32.0824367, 34.7892089 32.0824654, 34.7891549 32.0824809, 34.7891208 32.0824538, 34.7890698 32.0824201, 34.7890322 32.0824077, 34.7889862 32.0823998, 34.78895379999999 32.0824051, 34.7889236 32.0824102, 34.7889209 32.0824021, 34.7889012 32.0823507), (34.7894869 32.0822471, 34.7896115 32.08266589999999, 34.789777 32.0826305, 34.7896524 32.08221169999999, 34.7894869 32.0822471), (34.7896538 32.0828024, 34.7896934 32.0829576, 34.7898518 32.0829285, 34.7898122 32.0827734, 34.7896538 32.0828024), (34.7896569 32.0820271, 34.7896577 32.0820468, 34.7896659 32.0820653, 34.7896805 32.0820807, 34.7897 32.0820914, 34.7897226 32.0820964, 34.7897458 32.0820951, 34.7897674 32.0820877, 34.7897852 32.08207489999999, 34.7897973 32.0820581, 34.78980259999999 32.0820388, 34.7897999 32.0820176, 34.7897889 32.0819984, 34.7897709 32.0819834, 34.789748 32.0819746, 34.7897228 32.0819728, 34.7896985 32.0819784, 34.7896779 32.0819907, 34.7896633 32.0820081, 34.7896569 32.0820271)), ((34.7789012 32.0723507, 34.7791154 32.0722353, 34.7793334 32.0721058, 34.7795058 32.0720034, 34.7798462 32.0718012, 34.7798952 32.0717779, 34.7799146 32.0717792, 34.7799342 32.0717872, 34.7799749 32.0718304, 34.7800152 32.0718938, 34.7800494 32.0719833, 34.7800951 32.0721526, 34.7801415 32.0721419, 34.7802657 32.0725926, 34.7802395 32.0725983, 34.7802115 32.0726044, 34.7803126 32.0729483, 34.7802872 32.0729538, 34.7800336 32.0730091, 34.7794776 32.0731262, 34.7793901 32.0731446, 34.7793928 32.0731335, 34.7794169 32.0730342, 34.7792506 32.0724367, 34.7792089 32.0724654, 34.7791549 32.0724809, 34.7791208 32.0724538, 34.7790698 32.0724201, 34.7790322 32.0724077, 34.7789862 32.0723998, 34.7789538 32.0724051, 34.7789236 32.0724102, 34.7789209 32.0724021, 34.7789012 32.0723507), (34.7794869 32.0722471, 34.7796115 32.0726659, 34.779777 32.0726305, 34.7796524 32.0722117, 34.7794869 32.0722471), (34.7796538 32.0728024, 34.7796934 32.0729576, 34.7798518 32.0729285, 34.7798122 32.0727734, 34.7796538 32.0728024), (34.7796569 32.0720271, 34.7796577 32.0720468, 34.7796659 32.0720653, 34.7796805 32.0720807, 34.7797 32.0720914, 34.7797226 32.0720964, 34.7797458 32.0720951, 34.7797674 32.0720877, 34.7797852 32.0720749, 34.7797973 32.0720581, 34.7798026 32.0720388, 34.7797999 32.0720176, 34.7797889 32.0719984, 34.7797709 32.0719834, 34.779748 32.0719746, 34.7797228 32.0719728, 34.7796985 32.0719784, 34.7796779 32.0719907, 34.7796633 32.0720081, 34.7796569 32.0720271)))'
    multipolygon = wkt.loads(multipolygon)
    iso = il.OsmIsolines(multipolygon[0].centroid,
                         metric='distance',
                         values=[10, 50, 100],
                         unit='ft', network_type='walk',
                         address_type='point',
                         knn=2)
    iso.plot_isolines(plot_nodes=True, figsize=(10, 10))

def test_multigeometry():
    multipolygon = 'MULTIPOLYGON (((34.7889012 32.0823507, 34.7891154 32.0822353, 34.7893334 32.0821058, 34.7895058 32.0820034, 34.7898462 32.0818012, 34.7898952 32.0817779, 34.7899146 32.0817792, 34.7899342 32.0817872, 34.7899749 32.0818304, 34.7900152 32.0818938, 34.7900494 32.0819833, 34.79009509999999 32.0821526, 34.7901415 32.0821419, 34.7902657 32.0825926, 34.7902395 32.0825983, 34.7902115 32.0826044, 34.7903126 32.0829483, 34.79028719999999 32.0829538, 34.7900336 32.0830091, 34.7894776 32.0831262, 34.7893901 32.0831446, 34.78939279999999 32.0831335, 34.7894169 32.0830342, 34.7892506 32.0824367, 34.7892089 32.0824654, 34.7891549 32.0824809, 34.7891208 32.0824538, 34.7890698 32.0824201, 34.7890322 32.0824077, 34.7889862 32.0823998, 34.78895379999999 32.0824051, 34.7889236 32.0824102, 34.7889209 32.0824021, 34.7889012 32.0823507), (34.7894869 32.0822471, 34.7896115 32.08266589999999, 34.789777 32.0826305, 34.7896524 32.08221169999999, 34.7894869 32.0822471), (34.7896538 32.0828024, 34.7896934 32.0829576, 34.7898518 32.0829285, 34.7898122 32.0827734, 34.7896538 32.0828024), (34.7896569 32.0820271, 34.7896577 32.0820468, 34.7896659 32.0820653, 34.7896805 32.0820807, 34.7897 32.0820914, 34.7897226 32.0820964, 34.7897458 32.0820951, 34.7897674 32.0820877, 34.7897852 32.08207489999999, 34.7897973 32.0820581, 34.78980259999999 32.0820388, 34.7897999 32.0820176, 34.7897889 32.0819984, 34.7897709 32.0819834, 34.789748 32.0819746, 34.7897228 32.0819728, 34.7896985 32.0819784, 34.7896779 32.0819907, 34.7896633 32.0820081, 34.7896569 32.0820271)), ((34.7789012 32.0723507, 34.7791154 32.0722353, 34.7793334 32.0721058, 34.7795058 32.0720034, 34.7798462 32.0718012, 34.7798952 32.0717779, 34.7799146 32.0717792, 34.7799342 32.0717872, 34.7799749 32.0718304, 34.7800152 32.0718938, 34.7800494 32.0719833, 34.7800951 32.0721526, 34.7801415 32.0721419, 34.7802657 32.0725926, 34.7802395 32.0725983, 34.7802115 32.0726044, 34.7803126 32.0729483, 34.7802872 32.0729538, 34.7800336 32.0730091, 34.7794776 32.0731262, 34.7793901 32.0731446, 34.7793928 32.0731335, 34.7794169 32.0730342, 34.7792506 32.0724367, 34.7792089 32.0724654, 34.7791549 32.0724809, 34.7791208 32.0724538, 34.7790698 32.0724201, 34.7790322 32.0724077, 34.7789862 32.0723998, 34.7789538 32.0724051, 34.7789236 32.0724102, 34.7789209 32.0724021, 34.7789012 32.0723507), (34.7794869 32.0722471, 34.7796115 32.0726659, 34.779777 32.0726305, 34.7796524 32.0722117, 34.7794869 32.0722471), (34.7796538 32.0728024, 34.7796934 32.0729576, 34.7798518 32.0729285, 34.7798122 32.0727734, 34.7796538 32.0728024), (34.7796569 32.0720271, 34.7796577 32.0720468, 34.7796659 32.0720653, 34.7796805 32.0720807, 34.7797 32.0720914, 34.7797226 32.0720964, 34.7797458 32.0720951, 34.7797674 32.0720877, 34.7797852 32.0720749, 34.7797973 32.0720581, 34.7798026 32.0720388, 34.7797999 32.0720176, 34.7797889 32.0719984, 34.7797709 32.0719834, 34.779748 32.0719746, 34.7797228 32.0719728, 34.7796985 32.0719784, 34.7796779 32.0719907, 34.7796633 32.0720081, 34.7796569 32.0720271)))'
    multipolygon = wkt.loads(multipolygon)
    iso = il.OsmIsolines(multipolygon,
                         metric='distance',
                         values=[500, 1000, 1500],
                         unit='ft', network_type='walk',
                         )

def test_change_isolines_multi():
    multipolygon = 'MULTIPOLYGON (((34.7889012 32.0823507, 34.7891154 32.0822353, 34.7893334 32.0821058, 34.7895058 32.0820034, 34.7898462 32.0818012, 34.7898952 32.0817779, 34.7899146 32.0817792, 34.7899342 32.0817872, 34.7899749 32.0818304, 34.7900152 32.0818938, 34.7900494 32.0819833, 34.79009509999999 32.0821526, 34.7901415 32.0821419, 34.7902657 32.0825926, 34.7902395 32.0825983, 34.7902115 32.0826044, 34.7903126 32.0829483, 34.79028719999999 32.0829538, 34.7900336 32.0830091, 34.7894776 32.0831262, 34.7893901 32.0831446, 34.78939279999999 32.0831335, 34.7894169 32.0830342, 34.7892506 32.0824367, 34.7892089 32.0824654, 34.7891549 32.0824809, 34.7891208 32.0824538, 34.7890698 32.0824201, 34.7890322 32.0824077, 34.7889862 32.0823998, 34.78895379999999 32.0824051, 34.7889236 32.0824102, 34.7889209 32.0824021, 34.7889012 32.0823507), (34.7894869 32.0822471, 34.7896115 32.08266589999999, 34.789777 32.0826305, 34.7896524 32.08221169999999, 34.7894869 32.0822471), (34.7896538 32.0828024, 34.7896934 32.0829576, 34.7898518 32.0829285, 34.7898122 32.0827734, 34.7896538 32.0828024), (34.7896569 32.0820271, 34.7896577 32.0820468, 34.7896659 32.0820653, 34.7896805 32.0820807, 34.7897 32.0820914, 34.7897226 32.0820964, 34.7897458 32.0820951, 34.7897674 32.0820877, 34.7897852 32.08207489999999, 34.7897973 32.0820581, 34.78980259999999 32.0820388, 34.7897999 32.0820176, 34.7897889 32.0819984, 34.7897709 32.0819834, 34.789748 32.0819746, 34.7897228 32.0819728, 34.7896985 32.0819784, 34.7896779 32.0819907, 34.7896633 32.0820081, 34.7896569 32.0820271)), ((34.7789012 32.0723507, 34.7791154 32.0722353, 34.7793334 32.0721058, 34.7795058 32.0720034, 34.7798462 32.0718012, 34.7798952 32.0717779, 34.7799146 32.0717792, 34.7799342 32.0717872, 34.7799749 32.0718304, 34.7800152 32.0718938, 34.7800494 32.0719833, 34.7800951 32.0721526, 34.7801415 32.0721419, 34.7802657 32.0725926, 34.7802395 32.0725983, 34.7802115 32.0726044, 34.7803126 32.0729483, 34.7802872 32.0729538, 34.7800336 32.0730091, 34.7794776 32.0731262, 34.7793901 32.0731446, 34.7793928 32.0731335, 34.7794169 32.0730342, 34.7792506 32.0724367, 34.7792089 32.0724654, 34.7791549 32.0724809, 34.7791208 32.0724538, 34.7790698 32.0724201, 34.7790322 32.0724077, 34.7789862 32.0723998, 34.7789538 32.0724051, 34.7789236 32.0724102, 34.7789209 32.0724021, 34.7789012 32.0723507), (34.7794869 32.0722471, 34.7796115 32.0726659, 34.779777 32.0726305, 34.7796524 32.0722117, 34.7794869 32.0722471), (34.7796538 32.0728024, 34.7796934 32.0729576, 34.7798518 32.0729285, 34.7798122 32.0727734, 34.7796538 32.0728024), (34.7796569 32.0720271, 34.7796577 32.0720468, 34.7796659 32.0720653, 34.7796805 32.0720807, 34.7797 32.0720914, 34.7797226 32.0720964, 34.7797458 32.0720951, 34.7797674 32.0720877, 34.7797852 32.0720749, 34.7797973 32.0720581, 34.7798026 32.0720388, 34.7797999 32.0720176, 34.7797889 32.0719984, 34.7797709 32.0719834, 34.779748 32.0719746, 34.7797228 32.0719728, 34.7796985 32.0719784, 34.7796779 32.0719907, 34.7796633 32.0720081, 34.7796569 32.0720271)))'
    multipolygon = wkt.loads(multipolygon)
    iso = il.OsmIsolines(multipolygon,
                         metric='distance',
                         values=[500, 1000, 1500],
                         unit='ft', network_type='walk',
                         )

    iso.change_isolines(smooth=True, knn = 5)

def test_change_isolines():
    multipolygon = 'MULTIPOLYGON (((34.7889012 32.0823507, 34.7891154 32.0822353, 34.7893334 32.0821058, 34.7895058 32.0820034, 34.7898462 32.0818012, 34.7898952 32.0817779, 34.7899146 32.0817792, 34.7899342 32.0817872, 34.7899749 32.0818304, 34.7900152 32.0818938, 34.7900494 32.0819833, 34.79009509999999 32.0821526, 34.7901415 32.0821419, 34.7902657 32.0825926, 34.7902395 32.0825983, 34.7902115 32.0826044, 34.7903126 32.0829483, 34.79028719999999 32.0829538, 34.7900336 32.0830091, 34.7894776 32.0831262, 34.7893901 32.0831446, 34.78939279999999 32.0831335, 34.7894169 32.0830342, 34.7892506 32.0824367, 34.7892089 32.0824654, 34.7891549 32.0824809, 34.7891208 32.0824538, 34.7890698 32.0824201, 34.7890322 32.0824077, 34.7889862 32.0823998, 34.78895379999999 32.0824051, 34.7889236 32.0824102, 34.7889209 32.0824021, 34.7889012 32.0823507), (34.7894869 32.0822471, 34.7896115 32.08266589999999, 34.789777 32.0826305, 34.7896524 32.08221169999999, 34.7894869 32.0822471), (34.7896538 32.0828024, 34.7896934 32.0829576, 34.7898518 32.0829285, 34.7898122 32.0827734, 34.7896538 32.0828024), (34.7896569 32.0820271, 34.7896577 32.0820468, 34.7896659 32.0820653, 34.7896805 32.0820807, 34.7897 32.0820914, 34.7897226 32.0820964, 34.7897458 32.0820951, 34.7897674 32.0820877, 34.7897852 32.08207489999999, 34.7897973 32.0820581, 34.78980259999999 32.0820388, 34.7897999 32.0820176, 34.7897889 32.0819984, 34.7897709 32.0819834, 34.789748 32.0819746, 34.7897228 32.0819728, 34.7896985 32.0819784, 34.7896779 32.0819907, 34.7896633 32.0820081, 34.7896569 32.0820271)), ((34.7789012 32.0723507, 34.7791154 32.0722353, 34.7793334 32.0721058, 34.7795058 32.0720034, 34.7798462 32.0718012, 34.7798952 32.0717779, 34.7799146 32.0717792, 34.7799342 32.0717872, 34.7799749 32.0718304, 34.7800152 32.0718938, 34.7800494 32.0719833, 34.7800951 32.0721526, 34.7801415 32.0721419, 34.7802657 32.0725926, 34.7802395 32.0725983, 34.7802115 32.0726044, 34.7803126 32.0729483, 34.7802872 32.0729538, 34.7800336 32.0730091, 34.7794776 32.0731262, 34.7793901 32.0731446, 34.7793928 32.0731335, 34.7794169 32.0730342, 34.7792506 32.0724367, 34.7792089 32.0724654, 34.7791549 32.0724809, 34.7791208 32.0724538, 34.7790698 32.0724201, 34.7790322 32.0724077, 34.7789862 32.0723998, 34.7789538 32.0724051, 34.7789236 32.0724102, 34.7789209 32.0724021, 34.7789012 32.0723507), (34.7794869 32.0722471, 34.7796115 32.0726659, 34.779777 32.0726305, 34.7796524 32.0722117, 34.7794869 32.0722471), (34.7796538 32.0728024, 34.7796934 32.0729576, 34.7798518 32.0729285, 34.7798122 32.0727734, 34.7796538 32.0728024), (34.7796569 32.0720271, 34.7796577 32.0720468, 34.7796659 32.0720653, 34.7796805 32.0720807, 34.7797 32.0720914, 34.7797226 32.0720964, 34.7797458 32.0720951, 34.7797674 32.0720877, 34.7797852 32.0720749, 34.7797973 32.0720581, 34.7798026 32.0720388, 34.7797999 32.0720176, 34.7797889 32.0719984, 34.7797709 32.0719834, 34.779748 32.0719746, 34.7797228 32.0719728, 34.7796985 32.0719784, 34.7796779 32.0719907, 34.7796633 32.0720081, 34.7796569 32.0720271)))'
    multipolygon = wkt.loads(multipolygon)
    iso = il.OsmIsolines(multipolygon[0],
                         metric='distance',
                         values=[500, 1000, 1500],
                         unit='ft', network_type='walk',
                         )

    iso.change_isolines(smooth=True, knn = 5)



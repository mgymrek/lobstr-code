import os
import sys


def usage():
    print """
usage: python getYsearchLink.py  < genotypes.tab> 

takes in tab-delimited file with columns:

marker  genotype

and prints out ysearch link


file should include row containing header, followed by one row per marker
all markers without values are converted to -100, so they will be null on ysearch



"""

markerToQuery={
"DYS393":1,
"DYS390":2,
"DYS19/394":3,
"DYS19b":4,
"DYS391":5,
"DYS385a":6,
"DYS385b":7,
"DYS426":8,
"DYS388":9,
"DYS439":10,
"DYS389-1":11,
"DYS392":12,
"DYS389-2":13,
"DYS458":14,
"DYS459a":15,
"DYS459b":16,
"DYS455":17,
"DYS454":18,
"DYS447":19,
"DYS437":20,
"DYS448":21,
"DYS449":22,
"DYS464a":23,
"DYS464b":24,
"DYS464c":25,
"DYS464d":26,
"DYS464e":27,
"DYS464f":28,
"DYS464g":29,
"DYS460":30,
"GATA-H4":31,
"YCAIIa":32,
"YCAIIb":33,
"DYS456":34,
"DYS607":35,
"DYS576":36,
"DYS570":37,
"CDYa":38,
"CDYb":39,
"DYS442":40,
"DYS438":41,
"DYS531":54,
"DYS578":55,
"DYS395S1a":56,
"DYS395S1b":57,
"DYS590":58,
"DYS537":59,
"DYS641":60,
"DYS472":61,
"DYS406S1":62,
"DYS511":63,
"DYS425":42,
"DYS413a":64,
"DYS413b":65,
"DYS557":66,
"DYS594":67,
"DYS436":68,
"DYS490":69,
"DYS534":70,
"DYS450":71,
"DYS444":49,
"DYS481":72,
"DYS520":73,
"DYS446":51,
"DYS617":74,
"DYS568":75,
"DYS487":76,
"DYS572":77,
"DYS640":78,
"DYS492":79,
"DYS565":80,
"DYS461":43,
"DYS462":44,
"GATA-A10":45,
"DYS635":46,
"GAAT1B07":47,
"DYS441":48,
"DYS445":50,
"DYS452":52,
"DYS463":53,
"DYS434":81,
"DYS435":82,
"DYS485":83,
"DYS494":84,
"DYS495":85,
"DYS505":86,
"DYS522":87,
"DYS533":88,
"DYS549":89,
"DYS556":90,
"DYS575":91,
"DYS589":92,
"DYS636":93,
"DYS638":94,
"DYS643":95,
"DYS714":96,
"DYS716":97,
"DYS717":98,
"DYS726":99,
"DXY156-Y":100
}

def readMarkers(fname):
    markers = {}
    try:
        f = open(fname,"r")
    except:
        print "Error... couldn't open %s"%fname
        sys.exit(1)
    line = f.readline() # header
    line = f.readline()
    while(line != ""):
        marker,allele = line.split("\t")[0:2]
        marker = marker.strip('"')
        print marker, allele
        try:
            markers[markerToQuery[marker]] = int(allele)
        except:
            markers[markerToQuery[marker]] = -100
        line = f.readline()
    return markers


def getURL(markers):
    """ convert markers to a URL for ysearch """
    searchString = ""
    for i in range(1,101):
        searchString += "L%s=%s&"%(i,markers.get(i,-100))  
    url = "http://www.ysearch.org/search_results.asp?uid=&freeentry=true&"+searchString + "min_markers=8&mismatch_type=absolute&mismatches_max=0&mismatches_sliding_starting_marker=8"
    return url
    
def main():
    try:
        fname = sys.argv[1]
    except: 
        usage()
        sys.exit(1)
    if fname == "--help":
        usage()
        sys.exit(1)
    markers = readMarkers(fname)
    url = getURL(markers)
    print url

main()

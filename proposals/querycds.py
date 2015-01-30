# -*- coding: utf-8 -*-
"""
querycds, an object-oriented, pure Python query tool for
SIMBAD and Vizier by Brett Morris (UW). 

Notes
-----
Vizier queries based in part on QUERYVIZIER.pro (IDL), available at: 
http://idlastro.gsfc.nasa.gov/ftp/pro/sockets/queryvizier.pro

Vizier API reference: http://vizier.u-strasbg.fr/vizier/vizHelp/args.htx
"""

#from urllib import urlopen
import urllib
import urllib2
from urllib2 import urlopen
import time

def usemirror(simbad='cfa', vizier='cadc'):
    '''
    Set base URLs for queries
    Options: 'u-strasbg', 'cfa'
    '''
    global SIMBADbaseurl, vizierbaseurl
    if simbad.lower() == 'u-strasbg':
        SIMBADbaseurl = 'http://simbad.u-strasbg.fr/simbad/sim-id?Ident='
    elif simbad.lower() == 'cfa':
        SIMBADbaseurl = 'http://simbad.harvard.edu/simbad/sim-id?Ident='
    else:
        raise ValueError('SIMBAD mirror must be either "cfa" or "u-strasbg"')

    if vizier.lower() == 'u-strasbg':
        vizierbaseurl = 'http://vizier.u-strasbg.fr/viz-bin/'
    elif vizier.lower() == 'cfa':
        vizierbaseurl = 'http://vizier.cfa.harvard.edu/viz-bin/'
    elif vizier.lower() == 'cadc':
        vizierbaseurl = 'http://vizier.hia.nrc.ca/viz-bin/'
    else:
        raise ValueError('Vizier mirror must be either "cfa", "u-strasbg", or "cadc"')

def getmirrorspeeds(timeoutduration=2):
    '''
    Check how long it takes to access each mirror, return the name of the mirror
    with the fastest response.
    '''

    # Check SIMBAD Mirrors
    print 'Finding fastest SIMBAD mirror...'
    star = 'Rigel'
    baseurls = ['http://simbad.u-strasbg.fr/simbad/sim-id?Ident=', \
                'http://simbad.harvard.edu/simbad/sim-id?Ident=']
    mirrors = ['u-strasbg', 'cfa']
    pingtimes = []
    for mirror, baseurl in zip(mirrors, baseurls):
        print 'Pinging mirror: %s' % mirror
        try:
            starttime = time.time()
            txt = urlopen('%s%s&submit=SIMBAD+search' % (baseurl, star), timeout=timeoutduration)
            endtime = time.time()
            elapsedtime = endtime-starttime
        except urllib2.URLError:
            elapsedtime = timeoutduration
        pingtimes.append(elapsedtime)
    fastestSIMBADmirror = sorted(zip(pingtimes, mirrors))[0][1]
    print 'Fastest responding SIMBAD mirror: %s' % fastestSIMBADmirror

    # Check Vizier Mirrors
    print 'Finding fastest Vizier mirror...'
    baseurls = ['http://vizier.u-strasbg.fr/viz-bin/', \
                'http://vizier.cfa.harvard.edu/viz-bin/',\
                'http://vizier.hia.nrc.ca/viz-bin/']
    mirrors = ['u-strasbg', 'cfa', 'cadc']
    pingtimes = []
    catalog = 'SDSS-DR9'
    radiusMin = 1
    radiusMax = 2
    identifier = 'Vega'
    for mirror, baseurl in zip(mirrors, baseurls):
        print 'Pinging mirror: %s' % mirror
        try:
            starttime = time.time()
            txt = urlopen(baseurl + "asu-tsv/?" + \
                     urllib.urlencode({'-source':catalog, \
                     '-c':identifier, '-c.rm':"%s,%s" % (radiusMin,radiusMax),\
                     '-out.max':'unlimited'}), timeout=timeoutduration)
            endtime = time.time()
            elapsedtime = endtime-starttime
        except urllib2.URLError:
            elapsedtime = timeoutduration
        pingtimes.append(elapsedtime)
    fastestViziermirror = sorted(zip(pingtimes, mirrors))[0][1]
    print 'Fastest responding Vizier mirror: %s' % fastestViziermirror
    return {'simbad':fastestSIMBADmirror, 'vizier':fastestViziermirror}


## Find which mirror is fastest, set global variables (and each method's 
## defaults) based on the result
#fastestmirrors = getmirrorspeeds()
#usemirror(**fastestmirrors)
usemirror(simbad='u-strasbg', vizier='u-strasbg')

def txt2url(string):
    '''
    Take a string of text and replace some special characters
    with the equivalent URL encoding
    '''
    #return string.replace('+','%2b').replace(' ','%20')
    return urllib.quote_plus(string)

def sortSources(band, names, N_comps):
    '''
    Sort the lists of sources by their brightnesses in the `band` band. 
    
    Parameters
    ----------
    band : list of floats
        Magnitudes in one bandpass   
        
    names : list of strings
        Source identifiers for each magnitude in `band` 
        
    N_comps : int
        Number of comparison stars to sort into dictionary
    '''
    unsortedlists = zip(band, names)
    band, names = [list(eachlist) for eachlist in zip(*sorted(unsortedlists))]
    brightest_comps_sources = names[:N_comps]
    brightest_comps_mags = band[:N_comps]
    return brightest_comps_sources, brightest_comps_mags

def checksimbadforentry(identifier2MASS,SIMBADbaseurl=SIMBADbaseurl):
    '''
    Enter 2MASS identifier string (HHMMSSss+/-DDMMSSs), return the URL if there
    is a SIMBAD entry for the object, return None otherwise. 
    '''

    HTMLidentifier = '2MASS+J'+identifier2MASS.replace('+', '%2B').replace('-','%2D')
    URL = '%s%s&submit=SIMBAD+search' % (SIMBADbaseurl, HTMLidentifier)
    rawhtml = urlopen(URL).read()
    #print rawhtml
    insimbad = 'Identifier not found in the database' not in rawhtml


    if insimbad: 
        return URL
    else: 
        return None

def string2float(string):
    if '(~)' in string:
        return 0.0
    else: 
        return float(string)

class star:
    def __init__(self,identifier,SIMBADbaseurl=SIMBADbaseurl, vizierbaseurl=vizierbaseurl, manualdict=None):
        '''
        Initialize a SIMAD query object. Give an identifier which SIMBAD
        can resolve. Results may be misleading if the star is a member of a binary.

        Parameters
        ----------
        identifier : str
            The name of the star or planet to query for

        closest_SIMBADmirror : string
            Name the closest SIMBAD mirror to you, either 'cfa' or 'u-strasbg'
        '''
        if manualdict is None:
            self.dictionary = {}
            self.querySimbad(identifier,baseurl=SIMBADbaseurl)
            self.getRADec()
            self.getFluxes()
            self.getSpectralType()
        else:
            if '+' in identifier:
                idRA, idDec = identifier.split('+')
                idDec = '+'+idDec
            elif '-' in identifier:     
                idRA, idDec = identifier.split('-')
                idDec = '-'+idDec
            ralist = [idRA[0:2], idRA[2:4], idRA[4:9]]
            declist = [idDec[0:3], idDec[3:5], idDec[5:10]]
            RA_LIST = map(float, ralist)
            DEC_LIST = map(float, declist)
            RA_STRING = ' '.join(ralist)
            DEC_STRING = ' '.join(declist)

            self.dictionary = {'RA_LIST':RA_LIST, 'DEC_LIST':DEC_LIST, 'RA_STRING':RA_STRING, 'DEC_STRING':DEC_STRING, 'ID':identifier}
        #self.query2MASS_target()
        #self.query2MASS_comparisons(baseurl=vizierbaseurl)
        #self.querySDSS_comparisons()

    def __str__(self):
        longestkey = max([len(key) for key in self.dictionary.keys()]) + 1
        longestvalue = max([len(str(value)) for key, value in zip(self.dictionary.keys(), self.dictionary.values()) if 'comp' not in key]) + 1
        string = ''
        string += '{0:<{keyfill}s}{1:>{valuefill}s}\n'.format('KEYS', 'VALUES', keyfill=longestkey, valuefill=longestvalue)
        string += '{0:<{keyfill}s}{1:>{valuefill}s}\n'.format('----', '------', keyfill=longestkey, valuefill=longestvalue)

        for key in self.dictionary:
            if 'comp' in key:
                string += '{0:<{keyfill}s}{1:>{valuefill}s}\n'.format(key, '...', keyfill=longestkey, valuefill=longestvalue)                
            elif type(self.dictionary[key]) == str or type(self.dictionary[key]) == list:
                string += '{0:<{keyfill}s}{1:>{valuefill}s}\n'.format(key, self.dictionary[key], keyfill=longestkey, valuefill=longestvalue)
            elif type(self.dictionary[key]) == float:
                string += '{0:<{keyfill}s}{1:>{valuefill}f}\n'.format(key, self.dictionary[key], keyfill=longestkey, valuefill=longestvalue)
            elif type(self.dictionary[key]) == bool:
                string += '{0:<{keyfill}s}{1:>{valuefill}s}\n'.format(key, str(self.dictionary[key]), keyfill=longestkey, valuefill=longestvalue)
        return string

    def querySimbad(self,identifier,baseurl=SIMBADbaseurl,format='&output.format=ASCII'):
        '''
        Query Simbad for `identifier`, save the output data in a dictionary
        for quick lookups. Automatically remove the trailing ' b' to designate
        planets in `identifier`s.

        Parameters
        ----------
        identifier : str
            The name of the star or planet to query for

        baseurl : str, optional
            The first part of the URL address for the query. For the CfA mirror,
            use "http://simbad.harvard.edu/simbad/sim-id?Ident=" (default),
            for the Universite De Strasbourg, use
            "http://simbad.u-strasbg.fr/simbad/sim-id?Ident="

        format : str, optional
            Return results in ASCII format (default).
        '''
        # If identifier is a planet name (ending in a lower case letter)
        #  separated by a space from the rest of the name, then remove that
        #  last letter and rename the identifier for just the star.
        if identifier.split(' ')[-1].islower() and len(identifier.split(' ')[-1]) == 1:
            identifier_list = identifier.split(' ')[:-1]
            identifier = ' '.join(identifier_list).strip()

        # If identifier is a binary star, identify just the name of the binary
        #  system, not the particular member that's a planet host. Test case: WASP-77 A b
        if identifier.split(' ')[-1].isupper() and len(identifier.split(' ')[-1]) == 1:   
            identifier_list = identifier.split(' ')[:-1]
            identifier = ' '.join(identifier_list).strip()        
        print 'Querying Simbad for '+identifier+'...'

        searchURL = baseurl+txt2url(identifier)+format
        self.rawResult = urlopen(searchURL).read()

        ## Catch failed queries
        if self.rawResult.startswith('!!'):
            ## If failed query contained a '-', replace it with a space
            ##  and try again
            if '-' in identifier:
                corrected_identifier = identifier.replace('-',' ')
                print 'Retrying query with : "%s"...' % corrected_identifier
                self.rawResult = urlopen(baseurl+corrected_identifier+format).read()

            ## If query is still failing raise error, otherwise notify the
            ##  user of what identifier actually worked
            if self.rawResult.startswith('!!'):
                #print 'Query failed:\n'+self.rawResult
                #return
                raise ValueError('Identifier "%s" could not be resolved' % identifier)
            else:
                print ('Identifier "%s" resolved upon automatically correcting' + \
                      ' to "%s".')  % (identifier,corrected_identifier)

        ## If target has been successfully resolved, now continue...
        ## Store the successfully resolved identifier
        self.dictionary['ID'] = self.rawResult.splitlines()[2]

    def getRADec(self):
        '''
        Parse the SIMBAD results for the IRCS RA and Dec of the object
        '''
        for row in self.rawResult.splitlines():
            if row.startswith('Coordinates(ICRS') or row.startswith('nullCoordinates(ICRS'):
                ## Raw coordinate string comes after the colon and is separated
                ##  by a double white space
                raw_coord_str = row.split(':')[1].split('  ')

                ## If SIMBAD has coordinates:
                if raw_coord_str[0] != ' No Coord.':
                    ## Raw coordinate string has the RA and DEC separated by a double space
                    self.dictionary['RA_STRING'] = raw_coord_str[0].strip()
                    self.dictionary['RA_LIST'] = map(string2float, \
                                           self.dictionary['RA_STRING'].split(' '))
    
                    ## Raw coordinate string has the DEC with some extra stuff on the end,
                    ##  so split up the string by single spaces and take the first three
                    ##  values, which are (deg, min, sec)
                    self.dictionary['DEC_STRING'] =  ' '.join(
                                           raw_coord_str[1].split(' ')[:3]).strip()
                    self.dictionary['DEC_LIST'] = map(string2float, \
                          self.dictionary['DEC_STRING'].replace('+','').split(' '))

                ## If SIMBAD does not have coordinates:
                else: 
                    print "SIMBAD has no coordinates for %s" % self.dictionary['ID']

    def getFluxes(self):
        '''
        Parse the SIMBAD results for the fluxes of the object
        '''
        for row in self.rawResult.splitlines():
            if row.startswith('Flux ') and len(row) > 5:
                ## Second column has filter name, fourth column has flux
                filtername = row.split(' ')[1]
                filterflux = row.split(' ')[3]
                self.dictionary[filtername] = float(filterflux)


    def getSpectralType(self):
        '''
        Parse the SIMBAD results for the spectral type
        '''
        for row in self.rawResult.splitlines():
            if row.startswith('Spectral type:'):
                # Spectral type is in the third row divided by spaces
                self.dictionary['SPEC_TYPE'] = row.split(' ')[2]

    def query2MASS_target(self,catalog='2MASS-PSC',baseurl=vizierbaseurl,radiusMax=0.01):
        '''
        Query Vizier for comparison stars in 2MASS-PSC by target identifier 
        resolved with SIMBAD,  
        
        Parameters
        ----------
        catalog : str
            Default = '2MASS-PSC'
            
        radiusMax : float
            The maximal radial distance from the 2MASS-PSC RA/Dec, in arcminutes.
            
        Returns
        -------
        source : str
            The 2MASS source designation (i.e. "Name") of the target star
            
        mags : list of floats
            The corresponding [J,H,K] band fluxes for the target      
                
        '''
        identifier = self.dictionary['ID']
        
        print "Querying Vizier for "+str(identifier)+" (as target)..."
        radiusMin = 0.0
        searchURL = baseurl + "asu-tsv/?"  + \
                     urllib.urlencode({'-source':catalog, \
                     '-c':identifier, '-c.rm':"%s,%s" % (radiusMin,radiusMax),\
                     '-out.max':'unlimited'})
        rawViz = urlopen(searchURL).read()

        ## Look to see if data were retrieved 
        data_row = 0
        for i in range(len(rawViz.splitlines())):
            if rawViz.splitlines()[i].startswith('RAJ2000'):
                data_row = i+3
        ## If data were found...
        if data_row != 0:
            rowdata = rawViz.splitlines()[data_row].split('\t')
            RA = float(rowdata[0])
            Dec = float(rowdata[1])
            source = rowdata[2]  # 2MASS source designation
            J = float(rowdata[3]) # J magnitude
            H = float(rowdata[5]) # H magnitude
            K = float(rowdata[7]) # K magnitude

            self.dictionary['2MASS_identifier'] = source
            self.dictionary['2MASS_J'] = J
            self.dictionary['2MASS_H'] = H
            self.dictionary['2MASS_K'] = K
            self.dictionary['2MASS_RA'] = RA
            self.dictionary['2MASS_DEC'] = Dec
        else: 
            print "2MASS point source not resolved with Vizier for target %s" % identifier

    def query2MASS_comparisons(self,catalog='2MASS-PSC',baseurl=vizierbaseurl,radiusMin=0.01,radiusMax=6,Ncomps=5):
        '''
        Query Vizier for comparison stars in 2MASS-PSC by target identifier 
        resolved with SIMBAD. 
        
        Finds the J, H, and K band magnitudes and source identifiers for the 
        `Ncomps` brightest stars within the annulus `radiusMin` to `radiusMax` 
        (in arcminutes) centered on the the target star from the 2MASS Point
        Source Catalog. These targets are then added to the object dictionary
        star object's dictionary, in this format for each band (the "J"s can be
        replaced with H and K, as well):
        
        Comparison star magnitudes: starobject.dictionary['comp_J_mags']
        Comparison star source IDs: starobject.dictionary['comp_J_sources']
        
                
        Parameters
        ----------
        catalog : str
            Default = '2MASS-PSC'
            
        radiusMin : float
            The minimum radial distance from the target to search for a
            comparison star, in arminutes. 
            
        radiusMax : float
            The maximal radial distance from the target to search for a
            comparison star, in arcminutes.
            
        Ncomps : int
            The number of comparison stars to return data for
        
        '''
        identifier = self.dictionary['ID']
        
        print "Querying Vizier for "+str(identifier)+" comparisons in catalog: "+str(catalog)+"..."
        searchURL = baseurl + "asu-tsv/?" + \
                     urllib.urlencode({'-source':catalog, \
                     '-c':identifier, '-c.rm':"%s,%s" % (radiusMin,radiusMax),\
                     '-out.max':'unlimited'})
        rawViz = urlopen(searchURL).read()

        ## Look to see if data were retrieved 
        data_row = 0
        for i in range(len(rawViz.splitlines())):
            if rawViz.splitlines()[i].startswith('RAJ2000'):
                data_row = i+3
        ## If data were found...
        if data_row != 0:
            first_data_row = len(rawViz.splitlines())
            for i in range(len(rawViz.splitlines())):
                if rawViz.splitlines()[i].startswith('RAJ2000'):
                    first_data_row = i+3
                elif i > first_data_row and len(rawViz.splitlines()[i]):
                    last_data_row = i-1

            sources = []
            RAs = []
            Decs = []
            J = []
            H = []
            K = []
            for row in rawViz.splitlines()[first_data_row:last_data_row]:
                rowdata = row.split('\t')
                try: 
                    RAs.append(float(rowdata[0]))
                    Decs.append(float(rowdata[1]))
                    sources.append(rowdata[2].strip())  # 2MASS source designation
                    J.append(float(rowdata[3])) # J magnitude
                    H.append(float(rowdata[5])) # H magnitude
                    K.append(float(rowdata[7])) # K magnitude
                except ValueError:
                    pass

            bands = [J, H, K]
            bands_str = ['J','H','K']
            for band, band_str in zip(bands, bands_str):
                brightest_comps_sources, brightest_comps_mags = sortSources(band, sources, Ncomps)
                self.dictionary['comp_'+band_str] = {}
                self.dictionary['comp_'+band_str+'_mags'] = brightest_comps_mags
                self.dictionary['comp_'+band_str+'_sources'] = brightest_comps_sources
                for name, mag in zip(brightest_comps_sources, brightest_comps_mags):
                    self.dictionary['comp_'+band_str][name] = mag

    def querySDSS_target(self,catalog='SDSS-DR9',baseurl=vizierbaseurl,radiusMin=0.0,radiusMax=10./60):
        '''
        Query Vizier for comparison stars in SDSS-DR9 by target identifier 
        resolved with SIMBAD.
        
        Finds the Sloan u, g, r, i, z band magnitudes and source identifiers 
        for the `Ncomps` brightest stars within the annulus `radiusMin` to 
        `radiusMax` (in arcminutes) centered on the the target star from the 
        2MASS Point Source Catalog. These targets are then added to the object 
        dictionary star object's dictionary, in this format for each band (the 
        "J"s can be replaced with H and K, as well):
        
        Comparison star magnitudes: starobject.dictionary['comp_r_mags']
        Comparison star source IDs: starobject.dictionary['comp_r_sources']
        
        Parameters
        ----------
        catalog : str
            Default = 'SDSS-DR9'
            
        radiusMin : float
            The minimum radial distance from the target to search for a
            comparison star, in arminutes. 
            
        radiusMax : float
            The maximal radial distance from the target to search for a
            comparison star, in arcminutes.
            
        Ncomps : int
            The number of comparison stars to return data for
        '''
        #import ephem
        identifier = self.dictionary['ID']
        #catalog = 'V/139/sdss9'
        searchURL = baseurl + "asu-tsv/?" + \
                     urllib.urlencode({'-source':catalog, \
                     '-c':identifier, '-c.rm':"%s,%s" % (radiusMin,radiusMax),\
                     '-out.max':'unlimited', '-c.eq':'J2000'})
        rawViz = urlopen(searchURL).read()
        ## Look to see if data were retrieved 
        data_row = 0
        for i in range(len(rawViz.splitlines())):
            if rawViz.splitlines()[i].startswith('mode'):
                data_row = i+3
        ## If data were found...
        last_data_row = len(rawViz.splitlines()[i])-2
        if data_row != 0:
            self.dictionary['inSDSS'] = True
            first_data_row = len(rawViz.splitlines())
            for i in range(len(rawViz.splitlines())):
                if rawViz.splitlines()[i].startswith('mode'):
                    first_data_row = i+3
                elif i > first_data_row and len(rawViz.splitlines()[i]) > 0:
                    last_data_row = i#-0
            sources = []
            u = []
            g = []
            r = []
            i = []
            z = []
            for row in rawViz.splitlines()[first_data_row:last_data_row]:
                rowdata = row.split('\t')
                try: 
                    sources.append(rowdata[3].strip())  # 2MASS source designation
                    u.append(float(rowdata[10])) # u magnitude
                    g.append(float(rowdata[12])) # g magnitude
                    r.append(float(rowdata[14])) # r magnitude
                    i.append(float(rowdata[16])) # i magnitude
                    z.append(float(rowdata[18])) # z magnitude
                except ValueError:
                    pass
            bands = [u, g, r, i, z]
            bands_str = ['u','g','r','i','z']
            for band, band_str in zip(bands, bands_str):
                #brightest_comps_sources, brightest_comps_mags = sortSources(band, sources, Ncomps)
                if len(band) > 0:
                    self.dictionary['SDSS_'+band_str] = min(band)#brightest_comps_mags
                #self.dictionary['SDSS_'+band_str+'_mags'] = brightest_comps_mags
                #self.dictionary['SDSS_'+band_str+'_sources'] = brightest_comps_sources
                #for name, mag in zip(brightest_comps_sources, brightest_comps_mags):
                #    self.dictionary['comp_'+band_str][name] = mag
        else:
             self.dictionary['inSDSS'] = False

    def querySDSS_comparisons(self,catalog='SDSS-DR9',baseurl=vizierbaseurl,radiusMin=0.01,radiusMax=6,Ncomps=5):
        '''
        Query Vizier for comparison stars in SDSS-DR9 by target identifier 
        resolved with SIMBAD.
        
        Finds the Sloan u, g, r, i, z band magnitudes and source identifiers 
        for the `Ncomps` brightest stars within the annulus `radiusMin` to 
        `radiusMax` (in arcminutes) centered on the the target star from the 
        2MASS Point Source Catalog. These targets are then added to the object 
        dictionary star object's dictionary, in this format for each band (the 
        "J"s can be replaced with H and K, as well):
        
        Comparison star magnitudes: starobject.dictionary['comp_r_mags']
        Comparison star source IDs: starobject.dictionary['comp_r_sources']
        
        Parameters
        ----------
        catalog : str
            Default = 'SDSS-DR9'
            
        radiusMin : float
            The minimum radial distance from the target to search for a
            comparison star, in arminutes. 
            
        radiusMax : float
            The maximal radial distance from the target to search for a
            comparison star, in arcminutes.
            
        Ncomps : int
            The number of comparison stars to return data for
        '''
        identifier = self.dictionary['ID']
        print "Querying Vizier for "+str(identifier)+" comparisons in catalog: "+str(catalog)+"..."
        searchURL = baseurl + "asu-tsv/?" + \
                     urllib.urlencode({'-source':catalog, \
                     '-c':identifier, '-c.rm':"%s,%s" % (radiusMin,radiusMax),\
                     '-out.max':'unlimited', '-c.eq':'J2000'})
        rawViz = urlopen(searchURL).read()

        ## Look to see if data were retrieved 
        data_row = 0
        for i in range(len(rawViz.splitlines())):
            if rawViz.splitlines()[i].startswith('mode'):
                data_row = i+3
        ## If data were found...
        #last_data_row = rawViz.splitlines()
        if data_row != 0:
            first_data_row = len(rawViz.splitlines())
            for i in range(len(rawViz.splitlines())):
                if rawViz.splitlines()[i].startswith('mode'):
                    first_data_row = i+3
                elif i > first_data_row and len(rawViz.splitlines()[i]) > 0:
                    last_data_row = i-1

            sources = []
            u = []
            g = []
            r = []
            i = []
            z = []
            for row in rawViz.splitlines()[first_data_row:last_data_row]:
                rowdata = row.split('\t')
                try: 
                    sources.append(rowdata[3].strip())  # 2MASS source designation
                    u.append(float(rowdata[10])) # u magnitude
                    g.append(float(rowdata[12])) # g magnitude
                    r.append(float(rowdata[14])) # r magnitude
                    i.append(float(rowdata[16])) # i magnitude
                    z.append(float(rowdata[18])) # z magnitude
                except ValueError:
                    pass
            bands = [u, g, r, i, z]
            bands_str = ['u','g','r','i','z']
            for band, band_str in zip(bands, bands_str):
                brightest_comps_sources, brightest_comps_mags = sortSources(band, sources, Ncomps)
                self.dictionary['comp_'+band_str] = {}
                self.dictionary['comp_'+band_str+'_mags'] = brightest_comps_mags
                self.dictionary['comp_'+band_str+'_sources'] = brightest_comps_sources
                for name, mag in zip(brightest_comps_sources, brightest_comps_mags):
                    self.dictionary['comp_'+band_str][name] = mag

    def comps_report(self):
        '''
        Return sorted lists of the brightness of comparison stars in each
        queried bandpass, meant for printing to screen for review.
        '''
        allkeys = self.dictionary.keys()
        bands = [key.split('comp_')[1].split('_mags')[0] for key in allkeys if key.startswith('comp_') and key.endswith('_mags')]
        
        printstring = ""
        for band in bands:
            title = "\n"+band+"-band: Source ID, mags"
            printstring += title
            printstring += '\n'+(len(title)-1)*'-'+'\n'
            for name, mag in zip(self.dictionary['comp_'+band+'_sources'], self.dictionary['comp_'+band+'_mags']):
                printstring += name+'    '+str(mag)+'\n'
        return printstring
        
    def SDSScomps_nearest(self, band, threshold, sledgehammer=False):
        '''
        Return the comparisons in `band` within the nearest `threshold` magnitudes
        
        Parameters
        ----------
        band : str
            The name of the band to query
            
        threshold : float
            Number of magnitudes to return comparisons from

        sledgehammer : bool (optional)
            If True will return an empty dictionary if no SDSS data are available, 
            if False will raise an error if no SDSS data are available.
        '''
        targetkey = 'SDSS_%s' % band
        compkey = 'comp_%s' % band

        if sledgehammer:
            if self.dictionary.get(targetkey, None) is None or self.dictionary.get(compkey, None) is None:
                return {}
        else:
            if self.dictionary.get(targetkey, None) is None: 
                raise ValueError('Target mags from SDSS not retrieved. Use querycds.querySDSS_target() first.')
            elif self.dictionary.get(compkey, None) is None: 
                raise ValueError('Comparison star mags from SDSS not retrieved. Use querycds.querySDSS_comparisons() first.')
        nearestdict = {}
        for key in self.dictionary[compkey]:
            if abs(self.dictionary[compkey][key] - self.dictionary[targetkey]) <= threshold:
                nearestdict[key] = self.dictionary[compkey][key]
        return nearestdict
        
        
        
        
        

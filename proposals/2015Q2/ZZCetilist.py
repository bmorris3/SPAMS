import numpy as np
from matplotlib import pyplot as plt
import querycds
import ephem
import datetime

# Set up observatory parameters
observatory_minHorizon = '25:00:00'
apo = ephem.Observer()
apo.lat =  '32.0:46.0:49.0' 
apo.long = '-105.0:49.0:13.0' 
apo.elevation = 2788.0 # m
apo.temp = 10.0  ## Celsius
apo.horizon = observatory_minHorizon

def aboveHorizon(star, observatory, time, observatory_minHorizon=observatory_minHorizon):
    '''
    Take querycds star object, PyEphem observatory object, and datetime object. 
    Return true if the star is above the horizon at time `time`.
    '''
    RA = star.dictionary['RA_STRING'].strip().replace(' ',':')
    Dec = star.dictionary['DEC_STRING'].strip().replace(' ',':')
    observatory.date = time.isoformat(' ')
    star = ephem.FixedBody()
    star._ra = ephem.hours(RA)
    star._dec = ephem.degrees(Dec)
    star.compute(observatory)
    altitude = float(repr(star.alt))/(2*np.pi) * 360    ## Convert altitudes to degrees
    if altitude > float(ephem.degrees(observatory_minHorizon))*(180/np.pi): 
        return True
    else: 
        return False
    
def nighthours(observatory, time, twilight=-12, observatory_minHorizon=observatory_minHorizon):
    '''
    For time in datetime object `time` at observatory` observatory, 
    return a list of datetime objects for every hour during that night.
    '''
    # Temporarily set horizon to some twilight horizon (civil, nautical, astronomical)
    observatory.horizon = str(twilight)+':00:00'
    observatory.date = time.isoformat(' ')
    sunrise_hh, sunrise_mm, sunrise_ss = map(int, str(observatory.previous_rising(ephem.Sun())).split(' ')[1].split(':'))
    sunset_hh, sunset_mm, sunset_ss = map(int, str(observatory.next_setting(ephem.Sun())).split(' ')[1].split(':'))
    
    sunrise = sunrise_hh + sunrise_mm/60.0
    sunset = sunset_hh + sunset_mm/60.0
    Nhours = int(sunrise - sunset)
    # Set horizon back to old value
    observatory.horizon = observatory_minHorizon
    nighttimes = [datetime.datetime(time.year, time.month, time.day, i, sunset_mm) for i in range(Nhours+1)]
    return nighttimes

def getconstellation(starobj, observatory):
    RA = starobj.dictionary['RA_STRING'].strip().replace(' ',':')
    Dec = starobj.dictionary['DEC_STRING'].strip().replace(' ',':')
    star = ephem.FixedBody()
    star._ra = ephem.hours(RA)
    star._dec = ephem.degrees(Dec)
    star.compute(observatory)
    return ephem.constellation(star)

def hrsabovehoriz(star, observatory, date):
    '''
    Compute the number of hours a star spends above the horizon
    '''
    return sum([aboveHorizon(star, observatory, eachhour) for eachhour in nighthours(observatory, date)]) - 1

def bestmonths(star, observatory, threshold):
    '''
    Find the available months when the target is visible for more than `threshold` hours a night
    '''
    startdate = datetime.datetime(2015,1,1,7)
    Ndays = 365
    observablemonths = []
    Nnightsobservable = 0
    for i in range(Ndays):
        testdate = startdate + datetime.timedelta(days=i)
        Nhoursabovehorizon = hrsabovehoriz(s, apo, testdate)
        if Nhoursabovehorizon > threshold:
            observablemonths.append(testdate.strftime('%B'))
            Nnightsobservable += 1
    observablemonths = list(set(observablemonths))
    return observablemonths, Nnightsobservable

#stars = ['WD 2148+286','LFT 0487','GJ 3435',\
         #'WD 0046+051','GJ 3121','G191-B2B',\
         #'WD 1004+665']
filename= 'newZZC.txt' #import list of stars
#reading data into file
starlist= np.genfromtxt(filename,dtype=str, delimiter=' ')
#stars=starlist[:,0]
stars=starlist
Numberofgoodcomparisonstars=[]
goodstars=[]
renamethese=[]
radiusMax = 1.5
threshold=8
for star in stars:
    try:
        #get rid of weirdness in names
        star=star.replace('\xe2\x88\x92','-')
        #search for the star
        s = querycds.star(star, manualdict=True)
        s.querySDSS_target()
        s.querySDSS_comparisons(radiusMax=radiusMax, Ncomps=250)
        #Find the number of nearby comparisons
        Nnearbycomparisons = len(s.SDSScomps_nearest('g', 2, sledgehammer=True))
        #const, constellation = getconstellation(s, apo)
        #identify the months when the target is visible
        goodmonths, Nnights = bestmonths(s, apo, threshold)
        #add the number of comparisons for this star to the list
        Numberofgoodcomparisonstars.append(Nnearbycomparisons)
        if (('January' in goodmonths) or ('February' in goodmonths) or ('March' in goodmonths)) and Nnearbycomparisons >0:
            goodstars.append(star)
        #print 'Best months to observe:', ', '.join(goodmonths)
        #print 'Number of nights observable:', Nnights
        #print 'Number of available comparison stars:', Nnearbycomparisons
        #print 'Target in constellation:', constellation, '\n'
    except (ValueError,KeyError):
        renamethese.append(star)



with open('ZZcetifinal.txt', 'w') as outputfile:
    outputfile.write('\n'.join(goodstars))





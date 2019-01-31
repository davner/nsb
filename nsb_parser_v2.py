#!/usr/bin/env python3

#!python
#
#Usage:Matches Mid-Observation Times with Best Auxillary Sensor Times (SQMs and Boltwood)
#Assumes sensor observations have frequency of less than one-half the exposure time of images.
#E.G. if the exposure time is 10 seconds- the frequency of sensor observations should be
#separated by less than 5 seconds.

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta
import glob, pdb
import argparse


filefilter = '*.fts'

def output_filename(files, ftype):
    """
    ftype: '_nsb_sqm_tel_', '_nsb_sqm_zenith_', '_nsb_bolt'
    """

    hdu = fits.open(files[0])
    header = hdu[0].header
    date_obs = header['DATE-OBS']
    t_delta = TimeDelta(25200, format = 'sec')
    local_date_obs = Time(date_obs) - t_delta
    mod_date_obs = local_date_obs.value.split('T')[0].replace('-','')
    Filter = header['FILTER'] if 'sqm' in ftype else ''
    fname = mod_date_obs + ftype + Filter + '.csv'
    hdu.close()

    return fname


def Create_NSB_SQM_Tel(DataArray, files):
    fname = output_filename(files, '_nsb_sqm_tel_')
    with open(fname, 'w') as f:
        header = ('filename,Date-Obs,sqm_ut,az,elv,counts,reading\n')
        f.write(header)
        for i in range(DataArray.shape[0]):
            row = (DataArray[i,0]+','+
                   DataArray[i,2]+','+
                   DataArray[i,17]+','+
                   DataArray[i,4]+','+
                   DataArray[i,5]+','+
                   DataArray[i,18]+','+
                   DataArray[i,19]+'\n')
            f.write(row)
    f.close()
    return()

def Create_NSB_SQM_Zenith(DataArray, files):
    fname = output_filename(files, '_nsb_sqm_zenith_')
    with open(fname, 'w') as f:
        header = ('filename,Date-Obs,sqm_ut,counts,reading\n')
        f.write(header)
        for i in range(DataArray.shape[0]):
            row = (DataArray[i,0]+','+
                   DataArray[i,2]+','+
                   DataArray[i,20]+','+
                   DataArray[i,21]+','+
                   DataArray[i,22]+'\n')
            f.write(row)
    f.close()
    return()

#For clarity, making a function to populate the SQM data
#instead of having this block repeated twice for the two time cases in MapTimesSQM
def PopulateSQMData(Files_Times_Data,sqmtime,splitrow,j,SQMUnit):



    if ('Tele' in SQMUnit):
        Files_Times_Data[j,17] = str(sqmtime)     #Telescope SQM time
        Files_Times_Data[j,18] = splitrow[3]      #Telescope SQM counts
        Files_Times_Data[j,19] = splitrow[5]      #Telescope SQM reading


    if ('Zen' in SQMUnit):
        Files_Times_Data[j,20] = str(sqmtime)     #Zenith SQM time
        Files_Times_Data[j,21] = splitrow[3]      #Zenith SQM counts
        Files_Times_Data[j,22] = splitrow[5]      #Zenith SQM reading



    return(Files_Times_Data)



def Create_NSB_Boltwood(DataArray, files):
    fname = output_filename(files, '_nsb_bolt')
    with open(fname, 'w') as f:
        header = ('filename,Date-Obs,bolt_ut,sky_t,amb_t,wind,hum,dew,cloud\n')
        f.write(header)
        for i in range(DataArray.shape[0]):
            row = (DataArray[i,0]+','+
                   DataArray[i,2]+','+
                   DataArray[i,6]+','+
                   DataArray[i,10]+','+
                   DataArray[i,11]+','+
                   DataArray[i,12]+','+
                   DataArray[i,15]+','+
                   DataArray[i,16]+','+
                   DataArray[i,7]+'\n')
            f.write(row)
    f.close()
    return()


#For clarity, making a function to populate the boltwood data
#instead of having this block repeated twice for the two time cases in MapTimesBoltwood
def PopulateBoltwoodData(Files_Times_Data,btime,splitrow,j):


    Files_Times_Data[j,6] = str(btime)     #boltwood time
    Files_Times_Data[j,7] = splitrow[5]    #cloud flag "C"
    Files_Times_Data[j,8] = splitrow[6]    #wind flag "W"
    Files_Times_Data[j,9] = splitrow[7]    #Rain flag "R"
    Files_Times_Data[j,10] = splitrow[10]   #sky temperature
    Files_Times_Data[j,11] = splitrow[11]   #ambient temperature
    Files_Times_Data[j,12] = splitrow[12]   #wind speed
    Files_Times_Data[j,13] = splitrow[13]  #wet condition "w"
    Files_Times_Data[j,14] = splitrow[14]  #rain condition "r"
    Files_Times_Data[j,15] = splitrow[15]  #humidity
    Files_Times_Data[j,16] = splitrow[16]  #dew point


    return(Files_Times_Data)



#My goodness..this one really makes me mad
#hh:mm:60.00 exists! Time methods expect 0..59 for seconds....
def checkfor60sec(timepart):
    checkthestring = timepart.split(':')
    if (float(checkthestring[2]) >= 60.00):
       fixedtime = checkthestring[0]+':'+checkthestring[1]+':00.00'
    else:
       fixedtime = timepart #no change

    return(fixedtime)

def SortDataFilesChronologically(filefilter):
    FitfileNames = [] #list of file names
    JD_Obs = []  #list of times from Fits Header in JD format
    Date_Obs = []  #list of Date-Obs times from Fits Header in isot str format
    Az = [] #azimuth of telescope observation from Fits Header
    Alt = [] #altitude of telescope observation from Fits Header
    Exposure = [] #Exposure time of image. Necessary to computer mid-exposure time
    currentfile = fits.PrimaryHDU()  #variable to hold header
    for file in glob.glob(filefilter):  #create a list of files from directory
        FitfileNames.append(file)
    for file in range(len(FitfileNames)):  #Open each file header and create list of JD numeric times
        currentfile = fits.open(FitfileNames[file])
        header = currentfile[0].header
        JD_Obs.append(header['JD-Obs'])
        Date_Obs.append(header['Date-Obs'])
        Exposure.append(header['EXPOSURE'])
        Az.append(header['AZIMUTH'])
        Alt.append(header['ALTITUDE'])

    Boltwood_info = np.zeros([len(FitfileNames),11]) #stores boltwood data once correct time is found- to be populated later
                                                     #12 elements because it is the chosen boltwood time plus 11 other
                                                     #boltwood measurements
    SQM_info = np.zeros([len(FitfileNames),6])       #stores SQM data once correct time is found- to be populated later
                                                     #6 elements because it is the chosen SQM time, counts and sky
                                                     #brightness. These 3 measurements for each device (Telescope and Zenith)
    #creating array of FileNames,JD-Obs, Date-Obs, Azimuth, Altitude, Boltwood_Chosen time, Boltwood data[0], ... Boltwood data[10], SQMTel time, SQMTel counts, SQMTel reading, SQMzenith time, SQMzenith counts, SQMzenith reading,
    Files_Times_Data = np.column_stack((FitfileNames,JD_Obs,Date_Obs,Exposure,Az,Alt,Boltwood_info,SQM_info))
    #sort array based on column with JD-Obs time information
    Files_Times_Data = Files_Times_Data[Files_Times_Data[:,1].astype(np.float).argsort()]




    return(Files_Times_Data, FitfileNames)

def MapTimesSQM(Files_Times_Data,SQMUnit):

    StartTime = Time(Files_Times_Data[0,2],format='isot')
    EndTime =  Time(Files_Times_Data[(Files_Times_Data.shape[0])-1,2],format='isot')

    with open(SQMUnit, "r") as f:
      i = 0
      imagedex = 0
      LastSQMTime = 0 #This is a flag that when tripped will not compare any more image times with SQM times
      rows = f.readlines()
      sqm_datetime = prior_sqm_datetime = Time((Files_Times_Data[0,2]),format='isot')
      while (LastSQMTime == 0):
         #There are occassionally some returns with no strings
         #This protects .stripping the '\n' and leaving []
         if(rows[i] != '\n'):
            rows[i] = rows[i].strip()
            splitrow = rows[i].split(';')
            #If it is not a blank row, then check to see if the row is fully written with the correct number of elements
            if((len(splitrow) == 6) and ('#' not in splitrow[0])):
                 #splitrow[1] = checkfor60sec(splitrow[1]) #splitrow[1] is the time part, check for silly 60sec issue
                 sqm_datetime = Time((splitrow[0]),format='isot') # Nicely SQM times are in UT like the data file Date-Obs
                 DataFileObs = Time(Files_Times_Data[imagedex,2],format='isot')+TimeDelta(float(Files_Times_Data[0,3])/2,format='sec')
                 if (sqm_datetime > StartTime) and (LastSQMTime == 0):
                         if DataFileObs < sqm_datetime:
                            futuretime_difference = sqm_datetime - DataFileObs
                            pasttime_difference = DataFileObs - prior_sqm_datetime
                            if futuretime_difference <= pasttime_difference:
                               PopulateSQMData(Files_Times_Data,sqm_datetime,splitrow,imagedex,SQMUnit)
                            else:
                                PopulateSQMData(Files_Times_Data,prior_sqm_datetime,splitrow,imagedex,SQMUnit)
                            imagedex = imagedex + 1
                            if (sqm_datetime > EndTime): #when this condition is met, it is the very last SQM time that occurred
                                                              #after the last data frame was taken. Thus the flag is tripped and subsequently
                                                              #no more evaluations will happen
                                LastSQMTime = 1
            prior_sqm_datetime = sqm_datetime   #current time becomes previous time for the next loop
            i = i + 1
         else:
            i= i + 1 #if nothing else happens, i must be incremented
         if i == (len(rows)):  #error trapping
            print('Something bad happened here. The SQM times were never between the start and end times of the data files.\nMake certain you are using the correct boltwood file.')
            break
         if (len(splitrow) > 6):
            print('Please check the file at line ',i+1,' , there appears to be extra information there...')
            break
    return(Files_Times_Data)



def MapTimesBoltwood(Files_Times_Data,Boltwoodfile):

    StartTime = Time(Files_Times_Data[0,2],format='isot')
    EndTime =  Time(Files_Times_Data[(Files_Times_Data.shape[0])-1,2],format='isot')
    UToffset = 0.2916666666 #fraction of a day offset from UT (7/24)
    dotstring = '.'

    with open(Boltwoodfile, "r") as f:
      rows = f.readlines()
      prior_boltwood_datetime = Time((Files_Times_Data[0,2]),format='isot')
      i = 0
      imagedex = 0
      LastBoltTime = 0 #This is a flag that when tripped will not compare any more image times with boltwood times
      while (LastBoltTime == 0):
         #There are occassionally some returns with no strings
         #This protects .stripping the '\n' and leaving []
         if(rows[i] != '\n'):
            rows[i] = rows[i].strip()
            splitrow = rows[i].split()
            splitrow[1] = checkfor60sec(splitrow[1]) #splitrow[1] is the time part, check for silly 60sec issue
            boltwood_datetime = Time((splitrow[0] +'T'+ splitrow[1]),format='isot') #puts time in isot format
            boltwood_datetime = boltwood_datetime + TimeDelta(UToffset, format='jd') #this TimeDelta coverts boltwood time to UT
            #DataFileObs is the mid-exposure time. TimeDelta here adds 1/2 the exposure time to the Date-Obs time. [0,3] is exposure time.
            DataFileObs = Time(Files_Times_Data[imagedex,2],format='isot')+TimeDelta(float(Files_Times_Data[0,3])/2,format='sec')
            #There exist rows with less than 4 elements- but I check splitrow[3] explicitly...
            #so there must be at least 4 elements of the row to proceed any further
            if(len(splitrow) > 4):
                 if (boltwood_datetime > StartTime) and (LastBoltTime == 0) and (splitrow[3] == '~D'):
                         if DataFileObs < boltwood_datetime:
                            futuretime_difference = boltwood_datetime - DataFileObs
                            pasttime_difference = DataFileObs - prior_boltwood_datetime
                            print(dotstring,end="")
                            if futuretime_difference <= pasttime_difference:
                               PopulateBoltwoodData(Files_Times_Data,boltwood_datetime,splitrow,imagedex)
                            else:
                                PopulateBoltwoodData(Files_Times_Data,prior_boltwood_datetime,splitrow,imagedex)
                            imagedex = imagedex + 1
                            if (boltwood_datetime > EndTime): #when this condition is met, it is the very last boltwood time that occurred
                                                              #after the last data frame was taken. Thus the flag is tripped and subsequently
                                                              #no more evaluations will happen
                                LastBoltTime = 1
            if(prior_boltwood_datetime == StartTime):
                print('Searching Boltwood file forward in time for the time that matches the first image Date-Obs.\n')
            prior_boltwood_datetime = boltwood_datetime   #current time becomes previous time for the next loop
            i = i + 1
         else:
            i= i + 1 #if nothing else happens, i must be incremented
         if i == (len(rows)):  #error trapping
            print('ERROR: \n'+
                  'Something bad happened here. The boltwood times were never between the start and end times of the data files.\nMake certain you are using the correct boltwood file.')
            break
    return(Files_Times_Data)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-ls', nargs=1, type=str, help='Choose particular files to operate on. e.g. *.v.fts or *.fits')
    parser.add_argument('-bolt','--processbolt', help ='Use this option to search the directory and process the Boltwood file found there.', action='store_true')
    parser.add_argument('-sqm','--processsqm', help = 'Use this option to search for and process all SQM files in this directory.',action='store_true' )

    args = parser.parse_args()

    print('\nOrganizing Data Files Chronologically...\n')
    filefilter = args.ls[0]
    SortedDataFiles, filenames = SortDataFilesChronologically(filefilter)
    print('Files sorted.\n')

    if args.processbolt:
        Boltwoodfile = glob.glob('Boltwood*.txt')
        print('Mapping Boltwood Times....\n')
        DataArray = MapTimesBoltwood(SortedDataFiles,Boltwoodfile[0])
        print('\nBoltwood Times Mapped to Observation Times.\n')
        print('All done mapping times!\n')
        print('Writing CSV file...\n')
        Create_NSB_Boltwood(DataArray, filenames)


    if args.processsqm:
        SQMFiles = glob.glob('*SQM.dat')
        for datafile in SQMFiles:
            DataArray = MapTimesSQM(SortedDataFiles,datafile)
            if ('Tele' in datafile):
                print('Mapping Telescope SQM Times...\n')
            if ('Zen' in datafile):
                print('Mapping Zenith SQM Times...\n')
        print('All done mapping times!\n')
        print('Writing CSV files...\n')
        Create_NSB_SQM_Tel(DataArray, filenames)
        Create_NSB_SQM_Zenith(DataArray, filenames)


    print('FINISHED!')



#This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()

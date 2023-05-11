"""The madrigal_amisr module is used to automate updating Madrigal with
SRI formated hdf5 files.  Written by Bill Rideout; specification by
Bill Rideout and Todd Valentic.  Meant to be used in realtime by a data
transport application.

This module will automatically create a new Madrigal experiment if needed.
Also allows for switching back to previously existing experiments.  The start
method is used either to start a new experiment, or to switch to an old one.
The update method is used to add data.

Modified 2007-10-25 to also add batch updating.
Modified 2008-01-11 to match latest fitter output.
Modified 2008-01-22 to allow for plot directory creation and plots with same name
Modified 2018-04-18 to use Madrigal 3
Modified 2018-05-11 ??? brideout
Modified 2020-07-31 add chi-squared and SNR - brideout
Modified 2021-10-28 remove unused code, refactor, removed extraneous
        dependencies - Ashton S. Reimer
Modified 2021-11-03 removed multiprocessing, optimized and simplified codebase
        - Ashton S. Reimer
Modified 2021-11-18 added multiprocessing on a per file basis, instead of doing
        multiple files in parallel - Ashton S. Reimer
Modified 2022-08-14 Changed experiments for experiments0 as the folder to hold
        the data - Pablo M. Reyes
Modified 2023-02-15 change madFilename to include information about the radar
        mode, the integration time and the type of process,
        e.g. _lp_1min, or _vvels_5min - Pablo M. Reyes
Modified 2023-02-16 bring experimentsDirNum outside uploadExperiment() so that
        one can specify the folder experiments0 like the place to save madrigal
        files - Pablo M. Reyes 
Modified 2023-03-20 Creating definition kindat2fname(kindat) in order to help
        with filename formatting. - P. M. Reyes
        # e.g. _lp_fit_5min
        # e.g. _bc_nenotr_5min
        # e.g. _5min-lp_vvels_5min
        # e.g. _1min-lp_vvels_5min
        Adding file_version to createNewExperimentFromIni and uploadExperiment
        in order to add a file_version.
"""
#.......10........20........30........40........50........60........70........80

import os
import sys
import datetime
import traceback
import configparser
import logging
import glob
import shutil
import struct
import distutils.dir_util
import multiprocessing as mp
from multiprocessing import shared_memory
import dill # needed because we can't pickle
import tempfile

import tables
import numpy as np

import madrigal.metadata
import madrigal.cedar
import madrigal.admin

# global flag
WRITEHEADER = 1
POOLSIZE = 8

it2min = {
        1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7, 8:8, 9:9, 10:10, 11:11,
        12:12, 13:13, 14:14, 15:15, 16:16, 17:17, 18:18, 19:19, 20:20,
        25 : 30,
        30 : 45,
        35 : 60,
        40 : 90,
        45 : 120,
        60 : 45/60,
        64 : 30/60,
        68 : 15/60,
        72 : 10/60,
        76 : 5/60,
        80 : 4/60,
}

def kindat2fname(kindat):
    pc = kindat//10000
    pt = (kindat - 10000 * pc) // 100
    it = (kindat - 10000 * pc - 100 * pt)

    if pc == 100:
        pc_desc = "nenotr"
    elif pc == 200:
        pc_desc = "fit"
    elif pc == 300:
        pc_desc = "vvels"
    else:
        raise Exception(f'process type pc={pc} not available. '\
                f'kindat = {kindat}. '\
                f'valid values: 100,200,300')

    if pc in [100,200, 300]:
        if pt == 1:
            pt_desc = "lp" # long pulse
        elif pt == 2:
            pt_desc = "ac" # alternating code
        elif pt == 3:
            pt_desc = "bc" # binary code
        else:
            raise Exception(f'pulse type pt={pt} not available. '\
                f'kindat = {kindat}. '\
                    f'valid values: 1,2,3')

    intg_min = it2min[it]
    if intg_min < 1:
        it_desc = f"{int(intg_min * 60):02d}sec"
    elif intg_min >= 1:
        it_desc = f"{intg_min:02d}min"


    return f"_{pt_desc}_{pc_desc}_{it_desc}" 
    # e.g. _lp_fit_05min
    # e.g. _bc_nenotr_05min
    # e.g. _lp_vvels_05min
    # e.g. _lp_vvels_05min

def update_typetab(kindat,ckindat,typetab_file='/opt/madrigal/madrigal3/metadata/typeTab.txt'):
    # read the existing kindats from the typeTab
    with open(typetab_file,'r') as f:
        lines = f.readlines()

    existing_kindats = list()
    existing_ckindats = list()

    for line in lines:
        temp_k, temp_ck = line.strip('\n').split(',')
        existing_kindats.append(int(temp_k))
        existing_ckindats.append(temp_ck)

    # check to see if input kindat is in typeTab

    if not kindat in existing_kindats:
        print("Adding kindat %s to typeTab.txt")
        existing_kindats.append(kindat)
        existing_ckindats.append(ckindat)
        # sort them by kindat
        sorted_kindats = sorted(existing_kindats)
        sorted_ckindats = [y for _, y in sorted(zip(existing_kindats,
                           existing_ckindats), key=lambda pair: pair[0])]
        with open(typetab_file,'w') as f:
            for kindat,ckindat in zip(sorted_kindats,sorted_ckindats):
                f.write('%s,%s\n' % (str(kindat),ckindat))

    return


def update_cachedfiles(inst,kindat,cachedfiles_file='/opt/madrigal/madrigal3/cachedFiles.ini'):
    # if the kindat is for resolved velocities, we need to specify the array splitting
    mapping = {300: "{'array':'cgm_lat'}",  # latitude binning resolved velocities
              }

    proc_config = int(kindat) // 10000

    # if the kindat isn't a derived data product kindat we're done
    if not proc_config in list(mapping.keys()):
        return

    # check to make sure the kindat is in the cachedFiles.ini
    configfile = configparser.ConfigParser(interpolation=None)
    configfile.read(cachedfiles_file)
    existing_kindats = list(set([x.split('_')[0] for x in configfile.options(str(inst))]))

    # if not, add it:
    if not kindat in existing_kindats:
        configfile.set(str(inst),"%s_params" % kindat,value='')
        configfile.set(str(inst),"%s_format" % kindat,value=mapping[proc_config])

    with open(cachedfiles_file,'w') as f:
        configfile.write(f)

    return

def uploadMadrigalFile(madAdminObj,expPath,fullMadFilename,fileDesc,
        category,fileAnalyst,fileAnalystEmail):
# addMadrigalFile adds a new file to an experiment using metadata read from madFilename.
# Inputs:
#    expDir - full path to experiment directory (as returned by createMadriogalExperiment)
#    madFilename - full path to the complete Madrigal file.  Basename will be maintained.
#    permission - 0 (public) or 1 (private).
#    fileDesc - file description
#    category - 1=default, 2=variant, 3=history, or 4=realtime. Default is 1 (default file)
#    kindat - if not None (the default), use this kindat instead of what is found in the file.
#    notify - if True (the default), send a message to all registered users.  If False, do not.
#    fileAnalyst - full name of file Analyst.  Default is ''
#    fileAnalystEmail - email of file Analyst.  Default is ''
#    createCachedText - if True, add cached text file in overview/<basename>.txt.gz.  If False,
#        no cached file.
#    createCachedNetCDF4 - if True, add cached netCDF4 file in overview/<basename>.nc.  If False,
#        no cached file.
#    updateToMad3 - if False (the default), error raised if madFilename non-Hdf5 file. If True, try to
#        convert madFilename to Madrigal with .hdf5 extension before loading.
#    acceptOldSummary - if True, accept an old summary file. Used mainly for upgrading to Madrigal 3.
#        Default is False.
    madAdminObj.addMadrigalFile(expDir = expPath,
            madFilename = fullMadFilename,
            permission = 0,
            fileDesc = fileDesc,
            category = category,
            kindat = None,
            notify = True,
            fileAnalyst = fileAnalyst,
            fileAnalystEmail = fileAnalystEmail,
            createCachedText = True,
            createCachedNetCDF4 = True,
            updateToMad3 = False,
            acceptOldSummary = False)


def uploadMadrigalExp(madAdminObj,fullMadFilename,expTitle,
        fileDesc,category,optChar,experimentsDirNum,PI,PIEmail,
        fileAnalyst,fileAnalystEmail):
    """
    uploadMadrigalExp uploads an already created Madrigal file to the
    corresponding Madrigal local path.
    """
# createMadrigalExperiment creates a new experiment on Madrigal using metadata read from madFilename.
# Inputs:
#    madFilename - full path to the complete Madrigal file.  Basename will be maintained.
#    expTitle - experiment title
#    permission - 0 (public) or 1 (private) or -1 (ignore).
#    fileDesc - file description
#    instCode - instrument code.  If default (None), instrument code is taken from file,
#          but error    is thrown if more than one kinst found.
#    category - 1=default, 2=variant, 3=history, or 4=realtime. Default is 1 (default file)
#    optChar - optional character to be added to experiment directory if no dirName
#              given.  If dirName argument given, this argument ignored.  optChar
#              is used if the default directory name DDmmmYY is used for
#              more than one experiment created for a given instrument on a given day.
#              For example, if --optChar=h for a MLH experiment on September 12, 2005,
#              then the experiment directory created would be experiments/2005/mlh/12sep05h.
#    dirName - directory name to use for experiment.  If None (the default), the directory
#              name will be the default name DDmmmYY[optChar].  Cannot contain "/"
#    kindat - if not None (the default), use this kindat instead of what is found in the file.
#    experimentsDirNum - the number to be appended to the experiments directory, if experiments
#              directory being used is of the form experiments[0-9]* instead of just
#              experiments.  For example, if experimentsDirNum is 7, then the experiment
#              would be created in MADROOT/experiments7 instead of MADROOT/experiments.
#    PI- full name of principal investigator.  The default is ''
#    PIEmail - email of principal investigator.  The default is ''
#    fileAnalyst -full name of file analyst.  The default is ''
#    fileAnalystEmail - email of file analyst,.  The default is ''
#    createCachedText - if True, add cached text file in overview/<basename>.txt.gz.  If False,
#        no cached file.
#    createCachedNetCDF4 - if True, add cached netCDF4 file in overview/<basename>.nc.  If False,
#        no cached file.
#    notify - if True (the default), send a message to all registered users.  If False, do not.
#    updateToMad3 - if False (the default), error raised if madFilename non-Hdf5 file. If True, try to
#        convert madFilename to Madrigal with .hdf5 extension before loading.
    expPath = madAdminObj.createMadrigalExperiment(
            madFilename = fullMadFilename,
            expTitle = expTitle,
            permission = 0,
            fileDesc = fileDesc,
            instCode = None,
            category = category,
            optChar = optChar,
            dirName = None,
            kindat = None,
            experimentsDirNum = experimentsDirNum,
            PI = PI,
            PIEmail = PIEmail,
            fileAnalyst = fileAnalyst,
            fileAnalystEmail = fileAnalystEmail,
            createCachedText = True,
            createCachedNetCDF4 = True,
            notify = True,
            updateToMad3 = False)

    return expPath

def get_unique_fname(expPath, plotBasenames = []):
    """
    Routine to get a unique plotX.html file name.
    """
    plotFiles = glob.glob(os.path.join(expPath, 'plot*.html'))

    for plotFile in plotFiles:
        plotBasenames.append(os.path.basename(plotFile))

    plotNum = 0
    while True:
        plotName = 'plot%.3i.html' % (plotNum)

        if plotName not in plotBasenames:
            break

        plotNum += 1

    return plotName, plotBasenames

def createMad3File(args):
    """createMad3File is the method creates a single Madrigal 3 in parallel.

    args:   hdf5Type, hdf5Filename, instrument, kindat, 
            fullMadFilename, thisLowerRange, thisUpperRange, writeHeader, iniData,
            fileSection, principleInvestigator, expPurpose, expMode,
            cycleTime, correlativeExp, sciRemarks
    """
    hdf5Type, hdf5Filename, instrument, kindat, fullMadFilename, thisLowerRange, \
    thisUpperRange, writeHeader, iniData, fileSection, principleInvestigator, expPurpose, \
    expMode, cycleTime, correlativeExp, sciRemarks = args

    print("Creating Madrigal hdf5 file...")

    with tempfile.NamedTemporaryFile(delete=True,dir='/tmp') as tf:
        tempfile_name = tf.name + ".hdf5"

    fileHandler = hdf5Handler(hdf5Type)
    fileHandler.createMadrigalFile(hdf5Filename,instrument,kindat,None,
            tempfile_name,thisLowerRange,thisUpperRange)

    # header
    if writeHeader:
        try: kindatDesc = iniData.get(fileSection, 'extend_ckindat')
        except: kindatDesc = None
        try: analyst = iniData.get(fileSection, 'analyst')
        except: analyst = None
        try: comments = iniData.get(fileSection, 'comments')
        except: comments = None
        try: history = iniData.get(fileSection, 'history')
        except: history = None

        if not hdf5Filename is None:
            if not comments is None:
                comments = comments + '\n' + 'SOURCE_FILE %s' % (
                        os.path.basename(hdf5Filename))
            else:
                comments = 'SOURCE_FILE %s' % (
                        os.path.basename(hdf5Filename))

        print("Creating Madrigal Catalog Header...")
        now = datetime.datetime.now()
        catHeadObj = madrigal.cedar.CatalogHeaderCreator(tempfile_name)
        # principleInvestigator - Names of responsible Principal Investigator(s)
        #                         or others knowledgeable about the experiment.
        # expPurpose - Brief description of the experiment purpose
        # expMode - Further elaboration of meaning of MODEXP; e.g. antenna
        #           patterns and pulse sequences.
        # cycleTime - Minutes for one full measurement cycle
        # correlativeExp - Correlative experiments (experiments with related data)
        # sciRemarks - scientific remarks
        # instRemarks - instrument remarks

        catHeadObj.createCatalog(principleInvestigator=principleInvestigator,
                expPurpose=expPurpose,
                expMode=expMode,
                cycleTime=cycleTime,
                correlativeExp=correlativeExp,
                sciRemarks=sciRemarks)
        # kindatDesc - description of how this data was analyzed (the kind of data)
        # analyst - name of person who analyzed this data
        # comments - additional comments about data (describe any instrument-specific parameters)
        # history - a description of the history of the processing of this file
        catHeadObj.createHeader(kindatDesc=kindatDesc,
                analyst=analyst,
                comments=comments,
                history=history)
        catHeadObj.write()

        shutil.copyfile(tempfile_name, fullMadFilename)
        os.remove(tempfile_name)
        print("This took %s seconds" % (str(datetime.datetime.now()-now)))


def parseExpId(expId):
    """parseExpId parses an experiment id in the form YYYYMMDD.<inst_code>.<number>, and
    returns a tuple of (datetime, YYYYMMSS string, instId, optional char associated with number,
    and full path to the Madrigal experiment.

    Inputs:  expId - experiment id string in the form YYYYMMDD.<inst_code>.<number>, where
    the date represents the first day of the experiment, <inst_code> is the instrument
    code, and the trailing <number> is between 0 and 26

    Returns: a tuple with 5 items: 1. datetime represented by YYYYMMDD, 2. YYYYMMSS string
    itself, 3) the inst id, 4) the optional char associated with the number (0='', 1='a', ...
    26='z'), and 5) the string representing the full path to the Madrigal experiment in form
    $MADROOT/experiments/YYYY/<3_letter_inst>/DDmmmYYY<char>.

    Raises ValueError if expId not in proper format, instrument code not found,
    or trailing number < 0 or > 26.
    """
    madDBObj = madrigal.metadata.MadrigalDB()
    madInstObj = madrigal.metadata.MadrigalInstrument(madDBObj)

    try:
        year = int(expId[0:4])
        month = int(expId[4:6])
        day = int(expId[6:8])
    except:
        traceback.print_exc()
        raise ValueError('expId not in form YYYYMMDD.<inst_code>.<number>: <%s>' % (str(expId)))

    if year < 1900:
        raise ValueError('expId <%s> has too early year %i' % (str(expId), year))

    try:
        thisDate = datetime.datetime(year, month, day)
    except:
        traceback.print_exc()
        raise ValueError('expId not in form YYYYMMDD.<inst_code>.<number>: <%s>' % (str(expId)))

    try:
        items = expId.split('.')
        instCode = int(items[1])
        num = int(items[2])
    except:
        traceback.print_exc()
        raise ValueError('expId not in form YYYYMMDD.<inst_code>.<number>: <%s>' % (str(expId)))

    if num == 0:
        extChar = ''
    else:
        if num < 26:
            extChar = chr(96 + num)
        else:
            extChar = chr(96 + int(num/26)) + chr(96 + num%26)


    # get 3 letter instrument mnemonic
    mnem = madInstObj.getInstrumentMnemonic(instCode)
    if mnem == None:
        raise ValueError('unknown instrument code in expId: <%i>' % (instCode))


    dirName = os.path.join(madDBObj.getMadroot(),'experiments0','%04i' % year,mnem,'%s%s' % (thisDate.strftime('%d%b%y').lower(), extChar))

    return((thisDate, items[0], instCode, extChar, dirName))


class BatchExperiment:
    """BatchExperiment is a class to create and update AMISR Madrigal experiments

    """
    # defines length of line in Cedar catalog/header file

    __CEDAR_LEN__ = 80

    def uploadExperiment(self,iniFile,plotsdir='plots',file_version=1,
            removeTmpFiles=True, experimentsDirNum=None):
        """
        file_version : the 3 digit version attached to each file
        removeTmpFiles : When True, remove tmp file after uploading to Madrigal
        """

        # create needed Madrigal objects
        self.madDBObj = madrigal.metadata.MadrigalDB()
        madExpObj = madrigal.metadata.MadrigalExperiment(self.madDBObj)
        self.madInstObj = madrigal.metadata.MadrigalInstrument(self.madDBObj)  

        # read ini file
        self.__iniData__ = configparser.ConfigParser()
        self.__iniData__.read(iniFile)

        # DEFAULT
        ExperimentName = self.__iniData__.get('DEFAULT','ExperimentName')

        # get reqiured experiment info
        expTitle = self.__iniData__.get('Experiment', 'title')
        self.instrument = int(self.__iniData__.get('Experiment', 'instrument'))
        logFile = self.__iniData__.get('Experiment', 'logFile')
        expId = self.__iniData__.get('Experiment', 'expId')
        OutPath = self.__iniData__.get('Experiment','OutPath')

        # PI, Analyst
        PI = self.__iniData__.get('Experiment', 'pi')
        PIEmail = self.__iniData__.get('Experiment', 'PIEmail')
        fileAnalyst = self.__iniData__.get('Experiment', 'fileAnalyst')
        fileAnalystEmail = self.__iniData__.get('Experiment', 'fileAnalystEmail')
        # parse the expId
        try:
            items = expId.split('.')
            date = int(items[0])
            num = int(items[1])
            expId = items[0] + '.' + str(self.instrument) + '.' + items[1]
        except:
            traceback.print_exc()
            raise ValueError('expId not in form YYYYMMDD.<inst_code>.<number>: <%s>' % (str(expId)))

        # find the number of files being created
        numFiles = 0
        while True:
            try:
                self.__iniData__.get('File%i' % (numFiles + 1), 'hdf5Filename')
                numFiles += 1
            except configparser.NoSectionError:
                break

        # get optional character, if any
        optChar = parseExpId(expId)[3]

        # next find the time range in the data
        #firstTime = None
        #lastTime = None

        #for fileNum in range(numFiles):
        #    self.fileSection = 'File%i' % (fileNum + 1)
        #    hdf5Filename = self.__iniData__.get(self.fileSection, 'hdf5Filename')
        #    hdf5Type = self.__iniData__.get(self.fileSection, 'type')
        #    fileHandler = hdf5Handler(hdf5Type)
        #    startTime, endTime = fileHandler.getStartEndTimes(hdf5Filename)

        #    print(startTime, endTime, hdf5Filename)

        #    if firstTime == None:
        #        firstTime = startTime
        #    elif firstTime > startTime:
        #        firstTime = startTime

        #    if lastTime == None:
        #        lastTime = endTime
        #    elif lastTime < endTime:
        #        lastTime = endTime


        # create new madrigal file name template
        instMnemonic = self.madInstObj.getInstrumentMnemonic(self.instrument).lower()
        #madFilenameTemplate = '%s%02i%02i%02i.' % (instMnemonic, firstTime.year % 100,
        #                                           firstTime.month, firstTime.day)
        #madFilenameTemplate = '%s%04i%02i%02i' % (instMnemonic, firstTime.year,
        #                                           firstTime.month, firstTime.day)
        madFilenameTemplate = '%s%s' % (instMnemonic, ExperimentName)

        madAdminObj = madrigal.admin.MadrigalDBAdmin()
        for fileNum in range(numFiles):
            self.fileSection = 'File%i' % (fileNum + 1)
            kindat = int(self.__iniData__.get(self.fileSection, 'kindat'))

            if type(file_version) == int:
                ver2use = file_version
            elif type(file_version) == type(None):
                # here should be the last version created
                allversions = sorted(glob.glob(os.path.join(OutPath,
                    madFilenameTemplate + kindat2fname(kindat) + ".???.h5")))
                if len(allversions)>0:
                    ver2use = int(allversions[-1].split('.')[-2]) + 1
                else:
                    ver2use = 1
            else:
                raise Exception("file_version needs to be int or None.")

            madFilename = madFilenameTemplate + kindat2fname(kindat)\
                + f'.{ver2use:03d}.h5'
            fullMadFilename = os.path.join(OutPath,madFilename)

            print(f"working on file: {fullMadFilename}")
            if not os.path.exists(fullMadFilename):
                raise Exception(f"file {fullMadFilename} does not exist.")

            hdf5Type = self.__iniData__.get(self.fileSection, 'type')
            status = self.__iniData__.get(self.fileSection, 'status')
            category = int(self.__iniData__.get(self.fileSection, 'category'))
            fileDesc=status

            shutil.copyfile(fullMadFilename, os.path.join('/tmp',madFilename))
            tmpfullMadFilename=os.path.join('/tmp',madFilename)

            if fileNum==0:# create the experiment
                try:

                    expPath = uploadMadrigalExp(madAdminObj,tmpfullMadFilename,
                            expTitle,fileDesc,category,optChar,experimentsDirNum,
                            PI,PIEmail,fileAnalyst,fileAnalystEmail)
                except IOError:
                    x,y,z = sys.exc_info()
                    print(y)

                    expPath = str(y).split()[1]

                    info = input('Okay to Remove ' + expPath + '? type Yes: ')

                    if info=='Yes':
                        distutils.dir_util.remove_tree(expPath+'/',verbose=1)
                        expPath = uploadMadrigalExp(madAdminObj,tmpfullMadFilename,
                            expTitle,fileDesc,category,optChar,experimentsDirNum,
                            PI,PIEmail,fileAnalyst,fileAnalystEmail)
                    else:
                        raise IOError(y)
            else:
                uploadMadrigalFile(madAdminObj,expPath,tmpfullMadFilename,fileDesc,
                        category,fileAnalyst,fileAnalystEmail)

            # see if links to images are desired
            numLinks = 0
            image_dict = dict()
            while True:
                try:
                    imageTitle = self.__iniData__.get(self.fileSection,
                            'imageTitle%i' % (numLinks + 1))
                    image = self.__iniData__.get(self.fileSection,
                            'image%i' % (numLinks + 1))
                    #if 'Geometry Plot' in imageTitle:
                    if imageTitle not in image_dict.keys():
                        image_dict.update({imageTitle:[]})
                    image_dict[imageTitle].append(image)
                    logging.info(f'Including image {image}')
                    numLinks += 1
                except (configparser.NoSectionError, configparser.NoOptionError):
                    break
            if len(image_dict.keys())>0:
                self.createPlotLinks(madFilename,image_dict, expPath)

            if removeTmpFiles:
                logging.info(f'{tmpfullMadFilename} has been uploaded. Removing tmp file ...')
                os.remove(tmpfullMadFilename)

            self.expPath = expPath

        # create link to ISR AMISR database
        self.createLinkBack2amisr(self.instrument,ExperimentName,expPath)
        #return(expPath)



    def createNewExperimentFromIni(self, iniFile, skip_existing=False,
               skip_doc_plots=False, file_version = 1):
        """createNewExperimentFromIni will try to create a single new Madrigal experiment using
        information parsed from an input ini file.

        This method will also allow importing summary plots created outside this script into Madrigal.
        It will also allow the option of the automatic creation of individual record plots.

        file_version : the 3 digit version attached to each file

        Example ini file:

        [Experiment]

        title:              World Day
        instrument:         61
        logFile:            $BASE/amisr/pokerflat/face1/20071011-120000/log.txt
        expId:              20071011.61.000
        pi:                 Craig Heiselman
        modexp:             This is a one line experiment title
        cmodexp:            In this section you can write a multi-line description
                            of your experiment.  An ini file recognizes multiline
                            descriptions as long as every continuation line begins
                            with whitespace (as in this example).

        [File1]
        hdf5Filename:       $BASE/amisr/pokerflat/face1/20071011-120000/20071011.004_lp_3min.h5
        kindat:             5950
        type:               standard
        createRecPlots:     True
        imageTitle1:        Electron density - long pulse - beam 1
        image1:             /tmp/elecDensBeam1_lp.jpg
        imageTitle2:        Electron density - long pulse - beam 2
        image2:             /tmp/elecDensBeam2_lp.jpg
        ckindat:            In this section you can write a multi-line description
                            of how this particular file was created.  An ini file
                            recognizes multiline descriptions as long as every continuation line begins
                            with whitespace (as in this example)

        [File2]
        hdf5Filename:       $BASE/amisr/pokerflat/face1/20071011-120000/20071011.004_ac_3min.h5
        kindat:             5951
        type:               standard
        createRecPlots:     False
        imageTitle1:        Electron density - alternating code - beam 1
        image1:             /tmp/elecDensBeam1_ac.jpg
        imageTitle2:        Electron density - alternating code - beam 2
        image2:             /tmp/elecDensBeam2_ac.jpg

        [File3]
        hdf5Filename:       $BASE/amisr/pokerflat/face1/20071011-120000/20071011.004_lp_3min-vvels.h5
        kindat:             5953
        type:               velocity


        The input ini file is made up of an [Experiment] section and one or more [File*] section.
        No time values are used here. since the experiment times will be determined from the data itself.
        The names title, instrument, logFile, and expId are all required in the Experiment section.
        The expId name must be in the form YYYYMMDD.<inst_code>.<number>, where
        the date represents the first day of the experiment, <inst_code> is the instrument
        code, and the trailing <number> is between 0 and 26.

        In addition to the required fields in the Experiment section, there are also some optional
        fields designed to add experiment level information to the files' catalog record:

            pi - the experiment principle investigator
            modexp - a short experiment title
            cmodexp - a full description of the experiment.  These fields describe the experiment
                as a whole, and so will be the same for each file.  The field cmodexp will
                typically be multiple lines.  It is legal to have multiline values in an ini
                file as long as each new line begins with some white space.  

        For each file to be added to the experiment, a section called File* is required, and the
        numbering must be in increasing order starting with 1.  The names hdf5Filename, kindat, and type
        are required.  It is highly recommended that every file in an experiment have a unique
        kindat, because the kindat description is how the user determines the differences between
        files.  Madrigal administrators can always add additional kindats by editing the
        $MADROOT/metadata.typeTab.txt file (kindat descriptions must not contain commas).
        type deterimines the type of hdfs being loaded.  Presently supported types are standard,
        velocity, and uncorrected_ne_only.

        In addition to the required fields in the File* section, there are also some optional fields:
            createRecPlots -If set to True, the createRecPlots name will allow the creation of
                individual record plots. If not given or False, these plots will not be created.
                If type != standard, createRecPlots is ignored, since only works with standard data.

            imageTitle%i and image%i - must be given as a pair with matching numbers.  Allows images
                relating to that experiment to be imported into Madrigal.  The user will see a link to
                the imported image with text set by imageTitle.

            ckindat - a description of how this particular file was processed.  Will be included in the
                header record prepended to the file.  The field ckindat will typically be multiple lines.
                It is legal to have multiline values in an ini file as long as each new line begins
                with some white space.

            lowerRange - sets a lower range cutoff in km for uncorrected_ne_only files
            upperRange - sets a upper range cutoff in km for uncorrected_ne_only files
        """

        # create needed Madrigal objects
        self.madDBObj = madrigal.metadata.MadrigalDB()
        madExpObj = madrigal.metadata.MadrigalExperiment(self.madDBObj)
        self.madInstObj = madrigal.metadata.MadrigalInstrument(self.madDBObj)        

        # read ini file
        self.__iniData__ = configparser.ConfigParser()
        self.__iniData__.read(iniFile)

        # DEFAULT
        ExperimentName = self.__iniData__.get('DEFAULT','ExperimentName')

        # get reqiured experiment info
        expTitle = self.__iniData__.get('Experiment', 'title')
        self.instrument = int(self.__iniData__.get('Experiment', 'instrument'))
        logFile = self.__iniData__.get('Experiment', 'logFile')
        expId = self.__iniData__.get('Experiment', 'expId')
        OutPath = self.__iniData__.get('Experiment','OutPath')

        # PI, Analyst
        PI = self.__iniData__.get('Experiment', 'pi')
        PIEmail = self.__iniData__.get('Experiment', 'PIEmail')
        fileAnalyst = self.__iniData__.get('Experiment', 'fileAnalyst')
        fileAnalystEmail = self.__iniData__.get('Experiment', 'fileAnalystEmail')
        # parse the expId
        try:
            items = expId.split('.')
            date = int(items[0])
            num = int(items[1])
            expId = items[0] + '.' + str(self.instrument) + '.' + items[1]
        except:
            traceback.print_exc()
            raise ValueError('expId not in form YYYYMMDD.<inst_code>.<number>: <%s>' % (str(expId)))

        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(message)s',
                            filename=logFile,
                            filemode='w')

        logging.info('Creating exp using ini file %s with instrument %i and title <%s> and expId <%s>' % \
                     (iniFile,  self.instrument, expTitle, expId))

        # find the number of files being created
        numFiles = 0
        while True:
            try:
                self.__iniData__.get('File%i' % (numFiles + 1), 'hdf5Filename')
                numFiles += 1
            except configparser.NoSectionError:
                break

        if numFiles == 0:
            raise IOError('No File* section specified in ini file')

        # next find the time range in the data
        # firstTime = None
        # lastTime = None
        # for fileNum in range(numFiles):
        #     self.fileSection = 'File%i' % (fileNum + 1)
        #     hdf5Filename = self.__iniData__.get(self.fileSection, 'hdf5Filename')
        #     hdf5Type = self.__iniData__.get(self.fileSection, 'type')
        #     fileHandler = hdf5Handler(hdf5Type)
        #     startTime, endTime = fileHandler.getStartEndTimes(hdf5Filename)

        #     print(startTime, endTime, hdf5Filename)

        #     if firstTime == None:
        #         firstTime = startTime
        #     elif firstTime > startTime:
        #         firstTime = startTime

        #     if lastTime == None:
        #         lastTime = endTime
        #     elif lastTime < endTime:
        #         lastTime = endTime

        # create new madrigal file name template
        instMnemonic = self.madInstObj.getInstrumentMnemonic(self.instrument).lower()
        #madFilenameTemplate = '%s%02i%02i%02i.' % (instMnemonic, firstTime.year % 100,
        #                                           firstTime.month, firstTime.day)
        #madFilenameTemplate = '%s%04i%02i%02i' % (instMnemonic, firstTime.year,
        #                                           firstTime.month, firstTime.day)
        madFilenameTemplate = '%s%s' % (instMnemonic, ExperimentName)

        # header
        try: principleInvestigator = self.__iniData__.get('Experiment', 'pi')
        except: principleInvestigator = None
        try: expPurpose = self.__iniData__.get('Experiment', 'modexp')
        except: expPurpose = None
        try: expMode = self.__iniData__.get('Experiment', 'cmodexp')
        except: expMode = None
        try: cycleTime = self.__iniData__.get('Experiment', 'cycletime')
        except: cycleTime = None
        try: correlativeExp = self.__iniData__.get('Experiment', 'correxp')
        except: correlativeExp = None
        try: sciRemarks = self.__iniData__.get('Experiment', 'remarks')
        except: sciRemarks = None

        args = [] # each will be a tuple of (hdf5Type,hdf5Filename, self.instrument, kindat, 
                     #    fullMadFilename, thisLowerRange, thisUpperRange, writeHeader, iniData,
                     #  fileSection, principleInvestigator, expPurpose, expMode,
                     # cycleTime, correlativeExp, sciRemarks)

        # loop through all the files, and add to madrigal
        plots_links_dict = {}
        for fileNum in range(numFiles):
            self.fileSection = 'File%i' % (fileNum + 1)
            kindat = int(self.__iniData__.get(self.fileSection, 'kindat'))
            if type(file_version) == int:
                madFilename = madFilenameTemplate + kindat2fname(kindat)\
                    + f'.{file_version:03d}.h5'
                fullMadFilename = os.path.join(OutPath,madFilename)
            elif type(file_version) == type(None):
                for fcount in range(1,1000):
                    madFilename = madFilenameTemplate + kindat2fname(kindat)\
                                + f'.{fcount:03d}.h5'
                    fullMadFilename = os.path.join(OutPath,madFilename)
                    if not os.path.exists(fullMadFilename):
                        break
            else:
                raise Exception("file_version needs to be int or None.")

            print(f"working on file: {fullMadFilename}")

            hdf5Filename = self.__iniData__.get(self.fileSection, 'hdf5Filename')
            hdf5Type = self.__iniData__.get(self.fileSection, 'type')
            ckindat = self.__iniData__.get(self.fileSection, 'ckindat')

            # update the madrigal/metadata/typeTab.txt and madrigal/cachedFiles.ini files
            update_typetab(kindat, ckindat)
            update_cachedfiles(self.instrument,kindat)

            try:
                thisLowerRange = float(self.__iniData__.get(
                    self.fileSection, 'lowerRange'))
            except:
                thisLowerRange = None
            try:
                thisUpperRange = float(self.__iniData__.get(
                    self.fileSection, 'upperRange'))
            except:
                thisUpperRange = None

            args.append((hdf5Type, hdf5Filename, self.instrument, kindat,
                fullMadFilename, thisLowerRange, thisUpperRange, WRITEHEADER,
                self.__iniData__, self.fileSection, principleInvestigator,
                expPurpose, expMode, cycleTime, correlativeExp, sciRemarks))

            #logging.info(f'adding file from {str(firstTime)} to {str(lastTime)}'
            logging.info(f'adding experiment {ExperimentName} '
                    f' using {hdf5Filename}, kindat {kindat:d}, type {hdf5Type}'
                    f', lowerRange={str(thisLowerRange)}'
                    f', upperRange={str(thisUpperRange)}')

            # link plots associated with the file
            numLinks = 0
            image_dict = dict()
            while True:
                try:
                    imageTitle = self.__iniData__.get(self.fileSection,
                        'imageTitle%i' % (numLinks + 1))
                    image = self.__iniData__.get(self.fileSection,
                        'image%i' % (numLinks + 1))
                    if imageTitle not in image_dict.keys():
                        image_dict.update({imageTitle:[]})
                    image_dict[imageTitle].append(image)
                    numLinks += 1
                except (configparser.NoSectionError, configparser.NoOptionError):
                    break
            if len(image_dict.keys())>0:
                plots_links_dict.update({madFilename:[image_dict, OutPath]})

        # all file info has been gathered - kick off multiprocessing
        print('Processing %i files' % (len(args)))
        for i,arg in enumerate(args):
            print("...working on file %d: %s" % (i+1,arg[1]))
            if skip_existing and os.path.exists(arg[4]):
                print('Skipping file %s - already exists\n' % arg[4])
                continue
            else:
                createMad3File(arg)

        if not skip_doc_plots:
            print('Creating plot links...')
            for madfname, [image_dict, OutPath] in plots_links_dict.items():
                self.createPlotLinks(madfname,image_dict, OutPath)
            self.createLinkBack2amisr(self.instrument,ExperimentName,OutPath)
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # def createNewExperimentFromIni

    def createLinkBack2amisr(self,instrument,Experiment,expPath):
        """
        Create a link to the amisr.com page of the current Experiment

        Inputs:
            instrument -  int, instrument code e.g. 61 for PFISR
            Experiment - YYYYMMDD.XXX unique AMISR experiment name
            expPath - path to experiment directory
        """
        # get a unique filename
        plotName, plotBasenames = get_unique_fname(expPath)
        baseurl = "https://data.amisr.com/database/"
        amisr_isr_database = f"{baseurl}{instrument}/experiment/{Experiment}/3/"
        output_file = os.path.join(expPath, plotName)
        logging.info(f"Saving file {output_file}")
        with open(output_file, 'w') as f:
            f.write(f"""<html><head>
              <TITLE>AMISR ISR Database</TITLE>
              <meta http-equiv="refresh" content="0; url='{amisr_isr_database}'" />
             </head>
             <body>
               <p>You will be redirected to data.amisr.com soon!</p>
             </body>
            </html>
            """)

    def createPlotLinks(self,hdf5Filename,image_dict, expPath):
        """createPlotLinks is a method to create an html file associated
        with each hdf5 file and links to groups of plots.

        Inputs:
            hdf5Filename - the file that the plots will be associated with
            image_dict - a dictionary with plot types as keys and plot images
                as items
            expPath - path to experiment directory 
        """


        # get a unique filename
        plotName, plotBasenames = get_unique_fname(expPath)

        if not os.path.exists(os.path.join(expPath, 'plots')):
            os.mkdir(os.path.join(expPath, 'plots'))


        f_links_html = f"""<html><head>
        <TITLE>Plots associated with : {hdf5Filename}</TITLE>
        """
        f_links_html += """
        <style>
            div.gallery {
              margin: 2px;
              border: 1px solid #ccc;
              width: 400px;
              float: left;
            }

            div.gallery:hover {
              border: 1px solid #777;
            }

            div.gallery img {
              width: 45%;
              height: auto;
            }

            div.desc {
              width: 55%;
              float: right;
              font-size : 14px;
              text-align: center;
            }
            div.rowDiv {
              overflow: hidden; /* add this to contain floated children */
            }
            h4 {
               margin-top: 2px ;
               margin-bottom: 0 ;
            }
        </style>
        </head> <body>
        """
        f_links_html += f'<h4>Plots associated with : {hdf5Filename}</h4>\n'
        for imageTitle,imagefiles in image_dict.items():
            if "Geometry Plot" in imageTitle:
                assert len(image_dict) >1, "Need other plots, not only Geometry."
                assert len(imagefiles) == 1, f"So far, {imageTitle} expects 1 fig"
                logging.info("Creating link for imageTitle.")
                self.createLink(hdf5Filename,imageTitle,imagefiles[0],expPath,
                                  plotBasenames = [plotName])
                continue
            f_links_html += f'    <div class="rowDiv">\n'
            f_links_html += f'      <h4>{imageTitle.split(" ",1)[1]}</h4>\n'
            for imageFile in imagefiles:
                imgFileBase = os.path.basename(imageFile)
                # change base with madrigal fname
                uniquename = imgFileBase.split(imageTitle.split(" ",1)[0],1)[-1]
                imgFileBase = hdf5Filename.rsplit('.',1)[0] + " " + uniquename
                outName = os.path.join(expPath, 'plots',imgFileBase)

                if os.path.exists(outName):
                    imageFileTrailing = imgFileBase.rsplit('.',1)
                    imgNum=1

                    while True:
                        outName = '%s-%i.%s' % (imageFileTrailing[0],imgNum,
                                                imageFileTrailing[1])
                        if not os.path.exists(os.path.join(expPath, 'plots',
                                                                outName)):
                            break
                        imgNum+=1

                    outName = os.path.join(expPath,'plots',outName)

                shutil.copyfile(imageFile, outName)
                try:
                    os.chmod(outName,0o664)
                except:
                    pass
                image_link = os.path.basename(outName)
                f_links_html +=  '   <div class="gallery">\n'
                f_links_html += f'    <a target="_blank" href="plots/{image_link}">\n'
                f_links_html += f'      <img src="plots/{image_link}">\n'
                f_links_html +=  '      </a>\n'
                f_links_html += f'      <div class="desc">{image_link}</div>\n'
                f_links_html += '   </div>\n'
            f_links_html += '    </div>\n\n'

        f_links_html +="""
        </body></html>
        """

        with open(os.path.join(expPath, plotName), 'w') as f:
            f.write(f_links_html)

#.......10........20........30........40........50........60........70........80
    def createLink(self,hdf5Filename,imageTitle,imageFile,expPath,
            plotBasenames = []):
        """createLink is a method to create a new html file in expPath with a
           link to an image file.

        Inputs:
            hdf5Filename - the file that the plots will be associated with
            imageTitle - title of plot (will show up in Madrigal)
            imageFile - external image to be copied and displayed in Madrigal
            expPath - path to experiment directory
            plotBasenames - list of plot names to avoid
        """

        templateHtml = """<html><head>
        <TITLE>%s</TITLE>
        </head> <body><img src="plots/%s"></body></html>
        """


        # get a unique filename
        plotName, plotBasenames  = get_unique_fname(expPath, plotBasenames)

        if not os.path.exists(os.path.join(expPath, 'plots')):
            os.mkdir(os.path.join(expPath, 'plots'))

        imgFileBase = os.path.basename(imageFile)
        # change base with madrigal fname
        uniquename = imgFileBase.split(imageTitle.split(" ",1)[0],1)[-1]
        if uniquename[0] not in ['-','_']:
            uniquename = "_" + uniquename
        imgFileBase = hdf5Filename.rsplit('.',1)[0] + uniquename
        outName = os.path.join(expPath, 'plots',imgFileBase)

        if os.path.exists(outName):
            imageFileTrailing = imgFileBase.rsplit('.',1)
            imgNum=1

            while True:
                outName = '%s-%i.%s' % (imageFileTrailing[0],imgNum,
                                                imageFileTrailing[1])
                if not os.path.exists(os.path.join(expPath, 'plots',outName)):
                    break
                imgNum+=1

            outName = os.path.join(expPath,'plots',outName)

        shutil.copyfile(imageFile, outName)
        try:
            os.chmod(outName,0o664)
        except:
            pass

        output_file = os.path.join(expPath, plotName)
        logging.info(f"Saving file {output_file}")
        plot_title = hdf5Filename.rsplit('.',1)[0] + " " + imageTitle.split(" ",1)[1]
        with open(output_file, 'w') as f:
            f.write(templateHtml % (plot_title, os.path.basename(outName)))


class analyzeHdf5:
    """analyzeHdf5 is a class to analyze a SRI-formated hdf5 file containing standard ISR parameters
    """
    def __init__(self,hdf5File):
        """__init__ gets summary information about hdf5File containing uncorrected electron density.
        """
        self.__startTime = None # will be set to the earliest datetime
        self.__endTime = None # will be set to the latest datetime

        # read in all required data
        with tables.open_file(hdf5File) as hdfObj:
            start_time = hdfObj.root.Time.UnixTime[0,0]
            end_time = hdfObj.root.Time.UnixTime[-1,1]

        self.__startTime = datetime.datetime.utcfromtimestamp(start_time)
        self.__endTime = datetime.datetime.utcfromtimestamp(end_time)

        if self.__startTime > self.__endTime:
            raise Exception("Problem with start and end times: %s, %s, %s" % (self.__startTime,self.__endTime,hdf5File))


    def getStartEndTimes(self):
        return (self.__startTime, self.__endTime)


class hdf5ToMadrigal:
    """hdf5ToMadrigal is a class to turn a standard SRI-formated hdf5 file into a Madrigal file
    """

    def __init__(self,hdf5File,kinst,kindat,cedarObj,madrigalFile):
        """__init__ will write or update a Madrigal file using data in hdf5File

        Inputs:
            hdf5File - full path to hdf5 file with ISR data in SRI format
            kinst - instrument code (integer)
            kindat - data file kindat (integer)
            cedarObj - existing madrigal.cedar.MadrigalCedarFile to append data to.
                       If None, new madrigal.cedar.MadrigalCedarFile
                       created using madrigalFile.
            madrigalFile - name of Madrigal file to create or append to.

        Sets attributes self.numRecs, self.numTimes, self.numBeams
        """
        print("doing the fitted data stuff")
        now = datetime.datetime.now()

        # hard-coded indices defined by the format of the hdf file
        o_index = 0
        e_index = -1
        fractIndex = 0
        tempIndex = 1
        colIndex = 2
        velIndex = 3

        # parameter boundaries
        minTemp = 100.0 # MJN changed 06/05/2008
        maxTemp = 10000.0
        maxTempErr = 32765.0 # limit on Cedar format
        minTempErr = 1.0 # limit on Cedar format
        minNe = 1.0E9
        maxNe = 1.0E13
        maxNeErr = 3.2E13
        minNeErr = 1.0
        maxVo = 32765.0 # limit on Cedar format
        maxVoErr = 32765.0 # limit on Cedar format
        minVoErr = 0.01
        maxFract = 1.0
        minFract = 0.0

        # read in all required data
        with tables.open_file(hdf5File) as hdfObj:
            # beam codes
            beamCodes = hdfObj.root.BeamCodes
            beamCodeArray = beamCodes.read()
            self.numBeams = beamCodeArray.shape[0]

            # geomag
            plat = hdfObj.root.Geomag.MagneticLatitude
            platArray = plat.read()
            plong = hdfObj.root.Geomag.MagneticLongitude
            plongArray = plong.read()

            # ranges
            ranges = hdfObj.root.FittedParams.Range
            rangeArray = ranges.read()
            numRanges = rangeArray.shape[1]
            RangeTime=0
            if rangeArray.ndim==3:
                RangeTime=1

            # electron density (ne)
            ne = hdfObj.root.FittedParams.Ne
            neArray = ne.read()
            self.numTimes = neArray.shape[0]

            # error in electron density
            dne = hdfObj.root.FittedParams.dNe
            dneArray = dne.read()

            # ion info
            fits = hdfObj.root.FittedParams.Fits
            fitsArray = fits.read()

            # ion error info
            errors = hdfObj.root.FittedParams.Errors
            errArray = errors.read()

            # time info
            days = hdfObj.root.Time.Day
            dayArray = days.read()
            months = hdfObj.root.Time.Month
            monthArray = months.read()
            years = hdfObj.root.Time.Year
            yearArray = years.read()
            dtimes = hdfObj.root.Time.dtime
            dtimeArray = dtimes.read()

            # number of tx, rx
            numTxAeu = hdfObj.root.ProcessingParams.AeuTx
            numTxAeuArray = numTxAeu.read()
            numRxAeu = hdfObj.root.ProcessingParams.AeuRx
            numRxAeuArray = numRxAeu.read()

            # power info
            txPower = hdfObj.root.ProcessingParams.TxPower
            txPowerArray = txPower.read()

            # baud length
            baudLength = hdfObj.root.ProcessingParams.BaudLength.read()

            # pulse length
            pulseLength = hdfObj.root.ProcessingParams.PulseLength.read()

            baudCount = int(pulseLength/baudLength)
            if baudCount <= 0:
                baudCount = 1

            # tx freq
            txFreq = hdfObj.root.ProcessingParams.TxFrequency.read()

            # rx freq
            rxFreq = hdfObj.root.ProcessingParams.RxFrequency.read()

            # chisq
            chisq = hdfObj.root.FittedParams.FitInfo.chi2
            chisqArray = chisq.read()
            
            # snr
            snr = hdfObj.root.NeFromPower.SNR
            snrArray = snr.read()

        # create cedarObj if needed
        if cedarObj == None:
            cedarObj = madrigal.cedar.MadrigalCedarFile(madrigalFile,True,arraySplitParms=['beamid'])

        # create all data records 
        # loop first through num records, then through num beams
        print("Building list of data...")
        inputs = list()
        for recIndex in range(self.numTimes):
            # get start and end times for this record
            startYear = int(yearArray[recIndex][0])
            endYear = int(yearArray[recIndex][1])
            startMonth = int(monthArray[recIndex][0])
            endMonth = int(monthArray[recIndex][1])
            startDay = int(dayArray[recIndex][0])
            endDay = int(dayArray[recIndex][1])
            startDtime = dtimeArray[recIndex][0]
            endDtime = dtimeArray[recIndex][1]
            startHour = int(startDtime)
            endHour = int(endDtime)
            startMin = int(startDtime*60.0 - startHour*60.0)
            endMin = int(endDtime*60.0 - endHour*60.0)
            startSec = int(startDtime*3600.0) % 60
            endSec = int(endDtime*3600.0) % 60
            startTime = datetime.datetime(startYear, startMonth, startDay, startHour, startMin, startSec)
            endTime = datetime.datetime(endYear, endMonth, endDay, endHour, endMin, endSec)

            for beamIndex in range(self.numBeams):
                beamId = beamCodeArray[beamIndex][0]
                az = beamCodeArray[beamIndex][1]
                el = beamCodeArray[beamIndex][2]
                txpower = txPowerArray[recIndex]/1000.0
                numtxaeu = numTxAeuArray[recIndex]
                numrxaeu = numRxAeuArray[recIndex]
                if RangeTime:
                    rangeValue = rangeArray[recIndex][0]
                else:
                    rangeValue = rangeArray[0]

                neValue = neArray[recIndex][beamIndex]
                dneValue = dneArray[recIndex][beamIndex]
                fitsValue = fitsArray[recIndex][beamIndex]
                errValue = errArray[recIndex][beamIndex]
                chisqValue = chisqArray[recIndex][beamIndex]
                snrValue = snrArray[recIndex][beamIndex]
                platValue = platArray[beamIndex]
                plongValue = plongArray[beamIndex]

                limits = (minTemp,maxTemp,maxTempErr,minTempErr,minNe,maxNe,maxNeErr,minNeErr,
                          maxVo,maxVoErr,minVoErr,maxFract,minFract)

                indices = (o_index,e_index,fractIndex,tempIndex,colIndex,velIndex)

                inputs.append((kinst,kindat,startTime,endTime,az,el,beamId,txpower,numtxaeu,
                               numrxaeu,baudCount,pulseLength,txFreq,rxFreq,rangeValue,neValue,
                               dneValue,fitsValue,errValue,chisqValue,snrValue,platValue,
                               plongValue,limits,indices)
                             )


        print("There are %d jobs to execute..." % (len(inputs)))

        # do one to get the size
        dataRec = set_fitted_data_rec(inputs[0])
        encoded = dill.dumps(dataRec)
        datarecsize = len(encoded)
        dtype = '|S%d' % (datarecsize)
        output = np.zeros((len(inputs),), dtype=dtype)

        shm = shared_memory.SharedMemory(create=True, size=output.nbytes)
        shm_name = shm.name
        output_array = np.ndarray((len(inputs),),dtype=dtype,buffer=shm.buf)

        worker_pool = mp.Pool(POOLSIZE)
        jobs = [worker_pool.apply_async(mp_wrapper,args=(set_fitted_data_rec,i,args,shm_name,len(inputs),dtype)) for i,args in enumerate(inputs)]

        results = [p.get() for p in jobs]

        for i in range(len(inputs)):
            if i % 100 == 0:
                print("Dumping data: %d/%d" % (i + 1,len(inputs)))
            data_record = dill.loads(output_array[i])
            cedarObj.append(data_record)
            # dump records every 100
            if i % 50 == 0: 
                cedarObj.dump()

        # dump remaining records
        cedarObj.dump()
        shm.close()
        shm.unlink()
        print("This took %s seconds" % (str(datetime.datetime.now()-now)))
        print("Closing out the file...")
        now = datetime.datetime.now()
        cedarObj.close()
        self.numRecs = self.numTimes * self.numBeams
        print("This took %s seconds" % (str(datetime.datetime.now()-now)))


class hdf5VelocityToMadrigal:
    """hdf5VelocityToMadrigal is a class to turn a SRI-formated hdf5 file with vector velocities
    into a Madrigal file
    """

    def __init__(self,hdf5File,kinst,kindat,cedarObj,madrigalFile):
        """__init__ will write or update a Madrigal file using data in hdf5File containing vector velocities
        Inputs:
            hdf5File - full path to hdf5 file with vector velocity data in SRI format
            kinst - instrument code (integer)
            kindat - data file kindat (integer)
            cedarObj - existing madrigal.cedar.MadrigalCedarFile to append data to.
                       If None, new madrigal.cedar.MadrigalCedarFile
                       created using madrigalFile.
            madrigalFile - name of Madrigal file to create or append to.

        Sets attribute self.numRecs
        """
        print("doing the fitted data stuff")
        now = datetime.datetime.now()

        # hard-coded indices defined by the format of the hdf file
        north_index = 0
        east_index = 1
        parallel_index = 2

        maxV = 32765.0 # limit on Cedar format
        maxVErr = 32765.0 # limit on Cedar format
        minVErr = 1.0
        maxE = 32765.0*1e-5 # limit on Cedar format
        maxEErr = 32765.0*1e-5 # limit on Cedar format
        minEErr = 0.01*1e-5
        maxDirErr = 360.0

        # read in all required data
        with tables.open_file(hdf5File) as hdfObj:
            # time info
            unixTimes = hdfObj.root.Time.UnixTime
            timeArray = unixTimes.read()
            self.numRecs = timeArray.shape[0]

            # altitude info
            minAlt = hdfObj.root.ProcessingParams.MinAlt.read()
            maxAlt = hdfObj.root.ProcessingParams.MaxAlt.read()

            # nmeas - number of measurments
            nmeas = hdfObj.root.VectorVels.Nmeas
            nmeasArray = nmeas.read()
            params = ['nsmpta']

            try:
                byGeo = hdfObj.root.ProcessingParams.GeographicBinning.read()
            except:
                byGeo = 0

            # integration time
            integrationSecs = hdfObj.root.ProcessingParams.IntegrationTime.read()

            # magnetic latitude or latitude bins
            if byGeo!=0:
                if byGeo==2:
                    platArray = hdfObj.root.VectorVels.Latitude.read()
                    params.extend(['gdlat'])
                else:
                    platArray = hdfObj.root.VectorVels.MagneticLatitude.read()      
                    params.extend(['cgm_lat'])                  
                params.extend(['vi2','dvi2','vi1', 'dvi1', 'vi3', 'dvi3'])

            else:            
                try:
                    platArray = hdfObj.root.VectorVels.Plat.read()
                except:
                    platArray = hdfObj.root.VectorVels.MagneticLatitude.read()            
                params.extend(['cgm_lat','vipn','dvipn','vipe', 'dvipe', 'vi6', 'dvi6'])
            numPlat = platArray.shape[0]

            # velocity vectors and errors
            try:
                vestArray = hdfObj.root.VectorVels.Vest.read()
                dvestArray = hdfObj.root.VectorVels.errVest.read()
            except:
                vestArray = None
                dvestArray = None         
                xxxxxx               

            dvestArray=np.maximum(dvestArray,minVErr)
            # velocity magnitude, dir and errors
            try:
                vmagArray = hdfObj.root.VectorVels.Vmag.read()
                dvmagArray = hdfObj.root.VectorVels.errVmag.read()
                vdirArray = hdfObj.root.VectorVels.Vdir.read()
                dvdirArray = hdfObj.root.VectorVels.errVdir.read()
                if byGeo:
                    params.extend(['magvel','dmagvel','gnangle','dgnangle'])
                else:
                    params.extend(['magvel','dmagvel','nangle','dnangle'])            
            except:
                vmagArray = None
                dvmagArray = None
                vdirArray = None
                dvdirArray = None
                xxxxxx

            # electric field vectors and errors
            try:
                eestArray = hdfObj.root.VectorVels.Eest.read()
                deestArray = hdfObj.root.VectorVels.errEest.read()

                if byGeo:
                    params.extend(['en','den','ee', 'dee','eu','deu'])
                else:
                    params.extend(['epn','depn','epe', 'depe'])

            except:
                eestArray = None
                deestArray = None
                xxxxxx

            # electric field magnitude, dir and errors

            try:
                emagArray = hdfObj.root.VectorVels.Emag.read()
                demagArray = hdfObj.root.VectorVels.errEmag.read()
                edirArray = hdfObj.root.VectorVels.Edir.read()
                dedirArray = hdfObj.root.VectorVels.errEdir.read()
                if byGeo:
                    params.extend(['magef','dmagef','geangle','dgeangle'])
                else:
                    params.extend(['magef','dmagef','eangle','deangle'])  

            except:
                emagArray = None
                demagArray = None
                edirArray = None
                dedirArray = None
                xxxxxx

        # create cedarObj if needed
        if cedarObj == None:
            cedarObj = madrigal.cedar.MadrigalCedarFile(madrigalFile, True)

        # create all data records 
        for recIndex in range(self.numRecs):
            print("Doing record %d of %d" % (recIndex+1,self.numRecs))
            # get start and end times for this record
            startTime = datetime.datetime.utcfromtimestamp(timeArray[recIndex][0])
            endTime = datetime.datetime.utcfromtimestamp(timeArray[recIndex][1])
            dataRec = madrigal.cedar.MadrigalDataRecord(kinst,
                                                        kindat,
                                                        startTime.year,
                                                        startTime.month,
                                                        startTime.day,
                                                        startTime.hour,
                                                        startTime.minute,
                                                        startTime.second,
                                                        startTime.microsecond/10000,
                                                        endTime.year,
                                                        endTime.month,
                                                        endTime.day,
                                                        endTime.hour,
                                                        endTime.minute,
                                                        endTime.second,
                                                        endTime.microsecond/10000,
                                                        ('altb', 'alte', 'inttms'),
                                                        params,
                                                        numPlat,
                                                        ind2DList=['cgm_lat'])

            # set 1d values
            dataRec.set1D('altb', minAlt/1000.0) # m -> km
            dataRec.set1D('alte', maxAlt/1000.0) # m -> km
            dataRec.set1D('inttms', integrationSecs)

            # set 2d values
            for platIndex in range(numPlat):
                for parmIndex in range(len(params)):
                    # number of samples
                    if params[parmIndex] in ('nsmpta'): 
                        try:
                            if np.isnan(nmeasArray[recIndex][platIndex]):
                                raise ValueError()
                            dataRec.set2D('nsmpta', platIndex,nmeasArray[recIndex][platIndex])                
                        except:
                            dataRec.set2D('nsmpta', platIndex, 'missing')

                    # lat / plat
                    elif params[parmIndex] in ('gdlat','cgm_lat'):  
                        try:
                            if np.isnan(platArray[platIndex][0]) or np.isnan(platArray[platIndex][1]):
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex,(platArray[platIndex][0]+platArray[platIndex][1])/2.0)
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')                    

                    # vipn or vn
                    elif params[parmIndex] in ('vipn','vi2'):
                        try:
                            if np.isnan(vestArray[recIndex][platIndex][north_index]) \
                                or np.isnan(dvestArray[recIndex][platIndex][north_index]) \
                                or np.abs(dvestArray[recIndex][platIndex][north_index])<minVErr \
                                or np.abs(dvestArray[recIndex][platIndex][north_index])>maxVErr \
                                or np.abs(vestArray[recIndex][platIndex][north_index])>maxV:
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, vestArray[recIndex][platIndex][north_index])

                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # dvipn or dvn
                    elif params[parmIndex] in ('dvipn','dvi2'):
                        try:
                            if np.isnan(vestArray[recIndex][platIndex][north_index]) \
                                or np.isnan(dvestArray[recIndex][platIndex][north_index]) \
                                or np.abs(dvestArray[recIndex][platIndex][north_index])<minVErr \
                                or np.abs(dvestArray[recIndex][platIndex][north_index])>maxVErr \
                                or np.abs(vestArray[recIndex][platIndex][north_index])>maxV:
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, dvestArray[recIndex][platIndex][north_index])

                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # vipe or ve

                    elif params[parmIndex] in ('vipe','vi1'):                    
                        try:
                            if np.isnan(vestArray[recIndex][platIndex][east_index]) \
                                or np.isnan(dvestArray[recIndex][platIndex][east_index]) \
                                or np.abs(dvestArray[recIndex][platIndex][east_index])<minVErr \
                                or np.abs(dvestArray[recIndex][platIndex][east_index])>maxVErr \
                                or np.abs(vestArray[recIndex][platIndex][east_index])>maxV:                        
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, vestArray[recIndex][platIndex][east_index])

                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # dvipe or dve
                    elif params[parmIndex] in ('dvipe','dvi1'):                    
                        try:
                            if np.isnan(vestArray[recIndex][platIndex][east_index]) \
                                or np.isnan(dvestArray[recIndex][platIndex][east_index]) \
                                or np.abs(dvestArray[recIndex][platIndex][east_index])<minVErr \
                                or np.abs(dvestArray[recIndex][platIndex][east_index])>maxVErr \
                                or np.abs(vestArray[recIndex][platIndex][east_index])>maxV:                        
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, dvestArray[recIndex][platIndex][east_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')                

                    # vi6 or vi3
                    elif params[parmIndex] in ('vi6','vi3'):                       
                        try:
                            if np.isnan(vestArray[recIndex][platIndex][parallel_index]) \
                                or np.isnan(dvestArray[recIndex][platIndex][parallel_index]) \
                                or np.abs(dvestArray[recIndex][platIndex][parallel_index])<minVErr \
                                or np.abs(dvestArray[recIndex][platIndex][parallel_index])>maxVErr \
                                or np.abs(vestArray[recIndex][platIndex][parallel_index])>maxV:
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, vestArray[recIndex][platIndex][parallel_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # dvi6 or dvi3
                    elif params[parmIndex] in ('dvi6','dvi3'):                         
                        try:
                            if np.isnan(vestArray[recIndex][platIndex][parallel_index]) \
                                or np.isnan(dvestArray[recIndex][platIndex][parallel_index]) \
                                or np.abs(dvestArray[recIndex][platIndex][parallel_index])<minVErr \
                                or np.abs(dvestArray[recIndex][platIndex][parallel_index])>maxVErr \
                                or np.abs(vestArray[recIndex][platIndex][parallel_index])>maxV:
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, dvestArray[recIndex][platIndex][parallel_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # epn or en
                    elif params[parmIndex] in ('epn','en'):
                        try:
                            if np.isnan(eestArray[recIndex][platIndex][north_index]) \
                                or np.isnan(deestArray[recIndex][platIndex][north_index]) \
                                or np.abs(deestArray[recIndex][platIndex][north_index])<minEErr \
                                or np.abs(deestArray[recIndex][platIndex][north_index])>maxEErr \
                                or np.abs(eestArray[recIndex][platIndex][north_index])>maxE:
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, eestArray[recIndex][platIndex][north_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # depn or den
                    elif params[parmIndex] in ('depn','den'):
                        try:
                            if np.isnan(eestArray[recIndex][platIndex][north_index]) \
                                or np.isnan(deestArray[recIndex][platIndex][north_index]) \
                                or np.abs(deestArray[recIndex][platIndex][north_index])<minEErr \
                                or np.abs(deestArray[recIndex][platIndex][north_index])>maxEErr \
                                or np.abs(eestArray[recIndex][platIndex][north_index])>maxE: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, deestArray[recIndex][platIndex][north_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')
                    # epe or ee
                    elif params[parmIndex] in ('epe','ee'):                    
                        try:
                            if np.isnan(eestArray[recIndex][platIndex][east_index]) \
                                or np.isnan(deestArray[recIndex][platIndex][east_index]) \
                                or np.abs(deestArray[recIndex][platIndex][east_index])<minEErr \
                                or np.abs(deestArray[recIndex][platIndex][east_index])>maxEErr \
                                or np.abs(eestArray[recIndex][platIndex][east_index])>maxE:                         
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, eestArray[recIndex][platIndex][east_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # depe or dee
                    elif params[parmIndex] in ('depe','dee'):                    
                        try:
                            if np.isnan(eestArray[recIndex][platIndex][east_index]) \
                                or np.isnan(deestArray[recIndex][platIndex][east_index]) \
                                or np.abs(deestArray[recIndex][platIndex][east_index])<minEErr \
                                or np.abs(deestArray[recIndex][platIndex][east_index])>maxEErr \
                                or np.abs(eestArray[recIndex][platIndex][east_index])>maxE: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, deestArray[recIndex][platIndex][east_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')    

                    # eu
                    elif params[parmIndex] in ('eu'):
                        try:
                            if np.isnan(eestArray[recIndex][platIndex][parallel_index]) \
                                or np.isnan(deestArray[recIndex][platIndex][parallel_index]) \
                                or np.abs(deestArray[recIndex][platIndex][parallel_index])<minEErr \
                                or np.abs(deestArray[recIndex][platIndex][parallel_index])>maxEErr \
                                or np.abs(eestArray[recIndex][platIndex][parallel_index])>maxE: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, eestArray[recIndex][platIndex][parallel_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')          

                    # deu
                    elif params[parmIndex] in ('deu'):
                        try:
                            if np.isnan(eestArray[recIndex][platIndex][parallel_index]) \
                                or np.isnan(deestArray[recIndex][platIndex][parallel_index]) \
                                or np.abs(deestArray[recIndex][platIndex][parallel_index])<minEErr \
                                or np.abs(deestArray[recIndex][platIndex][parallel_index])>maxEErr \
                                or np.abs(eestArray[recIndex][platIndex][parallel_index])>maxE: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, deestArray[recIndex][platIndex][parallel_index])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')                                          

                    # magvel
                    elif params[parmIndex] in ('magvel'):
                        try:
                            if np.isnan(vmagArray[recIndex][platIndex]) \
                                or np.isnan(dvmagArray[recIndex][platIndex]) \
                                or np.abs(dvmagArray[recIndex][platIndex])<minVErr \
                                or np.abs(dvmagArray[recIndex][platIndex])>maxVErr \
                                or np.abs(vmagArray[recIndex][platIndex])>maxV: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, vmagArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                    # dmagvel
                    elif params[parmIndex] in ('dmagvel'):
                        try:
                            if np.isnan(vmagArray[recIndex][platIndex]) \
                                or np.isnan(dvmagArray[recIndex][platIndex]) \
                                or np.abs(dvmagArray[recIndex][platIndex])<minVErr \
                                or np.abs(dvmagArray[recIndex][platIndex])>maxVErr \
                                or np.abs(vmagArray[recIndex][platIndex])>maxV:                         
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, dvmagArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')                            
                    # magef
                    elif params[parmIndex] in ('magef'):
                        try:
                            if np.isnan(emagArray[recIndex][platIndex]) \
                                or np.isnan(demagArray[recIndex][platIndex]) \
                                or np.abs(demagArray[recIndex][platIndex])<minEErr \
                                or np.abs(demagArray[recIndex][platIndex])>maxEErr \
                                or np.abs(emagArray[recIndex][platIndex])>maxE: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, emagArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')
                    # dmagef
                    elif params[parmIndex] in ('dmagef'):
                        try: 
                            if np.isnan(emagArray[recIndex][platIndex]) \
                                or np.isnan(demagArray[recIndex][platIndex]) \
                                or np.abs(demagArray[recIndex][platIndex])<minEErr \
                                or np.abs(demagArray[recIndex][platIndex])>maxEErr \
                                or np.abs(emagArray[recIndex][platIndex])>maxE: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, demagArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')
                   # nangle or gnangle
                    elif params[parmIndex] in ('nangle','gnangle'):
                        try: 
                            if np.isnan(vdirArray[recIndex][platIndex]) or np.isnan(dvdirArray[recIndex][platIndex]) or dvdirArray[recIndex][platIndex]>maxDirErr or dvdirArray[recIndex][platIndex]<1e-2: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, vdirArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')

                   # dnangle or dgnangle
                    elif params[parmIndex] in ('dnangle','dgnangle'):
                        try:
                            if np.isnan(vdirArray[recIndex][platIndex]) or np.isnan(dvdirArray[recIndex][platIndex]) or dvdirArray[recIndex][platIndex]>maxDirErr or dvdirArray[recIndex][platIndex]<1e-2: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, dvdirArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')
                   # eangle or geangle
                    elif params[parmIndex] in ('eangle','geangle'):
                        try:
                            if np.isnan(edirArray[recIndex][platIndex]) or np.isnan(dedirArray[recIndex][platIndex]) or dedirArray[recIndex][platIndex]>maxDirErr or dedirArray[recIndex][platIndex]<1e-2: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, edirArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')                            
                   # deangle or dgeangle
                    elif params[parmIndex] in ('deangle','dgeangle'):
                        try:
                            if np.isnan(edirArray[recIndex][platIndex]) or np.isnan(dedirArray[recIndex][platIndex]) or dedirArray[recIndex][platIndex]>maxDirErr or dedirArray[recIndex][platIndex]<1e-2: 
                                raise ValueError()
                            dataRec.set2D(params[parmIndex], platIndex, dedirArray[recIndex][platIndex])
                        except:
                            dataRec.set2D(params[parmIndex], platIndex, 'missing')
                    else:
                        raise ValueError('Unhandled parameter %s' % (params[parmIndex]))

            # append new data record
            cedarObj.append(dataRec)
            # dump records every 100
            if recIndex % 100 == 0: 
                cedarObj.dump()

        # dump remaining records
        cedarObj.dump()
        print("This took %s seconds" % (str(datetime.datetime.now()-now)))
        print("Closing out the file...")
        now = datetime.datetime.now()
        cedarObj.close()
        print("This took %s seconds" % (str(datetime.datetime.now()-now)))


class hdf5UncorrectedToMadrigal:
    """hdf5ToMadrigal is a class to turn a SRI-formated hdf5 file with uncorrected electron
    density into a Madrigal file
    """
    def __init__(self,hdf5File,kinst,kindat,cedarObj,madrigalFile,lowerRange=None,upperRange=None):

        """__init__ will write or update a Madrigal file with uncorrected electron
        density using data in hdf5File

        Inputs:
            hdf5File - full path to hdf5 file with ISR data in SRI format
            kinst - instrument code (integer)
            kindat - data file kindat (integer)
            cedarObj - existing madrigal.cedar.MadrigalCedarFile to append data to.
                       If None, new madrigal.cedar.MadrigalCedarFile
                       created using madrigalFile.
            madrigalFile - name of Madrigal file to create or append to.
            lowerRange - lower range cutoff.  If None (the default), no lower range cutoff.
            upperRange - upper range cutoff.  If None (the default), no upper range cutoff.

        Sets attributes self.numRecs, self.numTimes, self.numBeams

        Raises IOError if no ranges fit in range limits
        """
        print("doing the uncorrected stuff")
        now = datetime.datetime.now()

        # hard-coded indices defined by the format of the hdf file
        o_index = 0
        e_index = 2
        fractIndex = 0
        tempIndex = 1
        colIndex = 2
        velIndex = 3

        # read in all required data
        with tables.open_file(hdf5File) as hdfObj:
            # beam codes
            beamCodes = hdfObj.root.BeamCodes
            beamCodeArray = beamCodes.read()
            self.numBeams = beamCodeArray.shape[0]

            # ranges (in meters)
            ranges = hdfObj.root.NeFromPower.Range
            rangeArray = ranges.read()
            if rangeArray.ndim==2:
                RangeTime=0
                numRanges = rangeArray.shape[1]
                lowerRangeIndex = None
                if lowerRange == None:
                    lowerRangeIndex = 0

                upperRangeIndex = None
                if upperRange == None:
                    upperRangeIndex = numRanges

                for i in range(numRanges):
                    if lowerRangeIndex == None and rangeArray[0][i]/1000.0 >=  lowerRange:
                        lowerRangeIndex = i
                    if upperRangeIndex == None and rangeArray[0][-1 - i]/1000.0 <  upperRange:
                        upperRangeIndex = numRanges - i

                if lowerRangeIndex == None:
                    # no ranges accepted
                    raise IOError('No valid ranges found between limits %s and %s' % (str(lowerRange), str(upperRange)))

                if upperRangeIndex == None:
                    upperRangeIndex = numRanges

                numUsedRanges = upperRangeIndex-lowerRangeIndex
                if numUsedRanges <= 0:
                    # no ranges accepted
                    raise IOError('No valid ranges found between limits %s and %s' % (str(lowerRange), str(upperRange)))
            else:

                RangeTime=1

            # uncorrected electron density (pop)
            pop = hdfObj.root.NeFromPower.Ne_NoTr
            popArray = pop.read()
            self.numTimes = popArray.shape[0]

            # error in electron density
            dpop = hdfObj.root.NeFromPower.dNeFrac
            dpopArray = dpop.read()

            # time info
            days = hdfObj.root.Time.Day
            dayArray = days.read()
            months = hdfObj.root.Time.Month
            monthArray = months.read()
            years = hdfObj.root.Time.Year
            yearArray = years.read()
            dtimes = hdfObj.root.Time.dtime
            dtimeArray = dtimes.read()

            # number of tx, rx
            numTxAeu = hdfObj.root.ProcessingParams.AeuTx
            numTxAeuArray = numTxAeu.read()
            numRxAeu = hdfObj.root.ProcessingParams.AeuRx
            numRxAeuArray = numRxAeu.read()

            # power info
            txPower = hdfObj.root.ProcessingParams.TxPower
            txPowerArray = txPower.read()

            # baud length
            baudLength = hdfObj.root.ProcessingParams.BaudLength.read()

            # pulse length
            pulseLength = hdfObj.root.ProcessingParams.PulseLength.read()
            baudCount = int(pulseLength/baudLength)
            if baudCount <= 0:
                baudCount = 1

            # tx freq
            txFreq = hdfObj.root.ProcessingParams.TxFrequency.read()
            # rx freq
            rxFreq = hdfObj.root.ProcessingParams.RxFrequency.read()

        # create cedarObj if needed
        if cedarObj == None:
            cedarObj = madrigal.cedar.MadrigalCedarFile(madrigalFile, True, arraySplitParms=['beamid'])

        # create all data records 
        # loop first through num records, then through num beams
        print("Building list of data...")
        inputs = list()
        for recIndex in range(self.numTimes):
            if RangeTime:
                numRanges = rangeArray.shape[2]
                lowerRangeIndex = None
                if lowerRange == None:
                    lowerRangeIndex = 0
                upperRangeIndex = None
                if upperRange == None:
                    upperRangeIndex = numRanges

                for i in range(numRanges):
                    if lowerRangeIndex == None and rangeArray[recIndex][0][i]/1000.0 >=  lowerRange:
                        lowerRangeIndex = i
                    if upperRangeIndex == None and rangeArray[recIndex][0][-1 - i]/1000.0 <  upperRange:
                        upperRangeIndex = numRanges - i
                if lowerRangeIndex == None:
                    # no ranges accepted
                    raise IOError('No valid ranges found between limits %s and %s' % (str(lowerRange), str(upperRange)))
                if upperRangeIndex == None:
                    upperRangeIndex = numRanges

                numUsedRanges = upperRangeIndex-lowerRangeIndex
                if numUsedRanges <= 0:
                    # no ranges accepted
                    raise IOError('No valid ranges found between limits %s and %s' % (str(lowerRange), str(upperRange)))

            # get start and end times for this record
            startYear = int(yearArray[recIndex][0])
            endYear = int(yearArray[recIndex][1])
            startMonth = int(monthArray[recIndex][0])
            endMonth = int(monthArray[recIndex][1])
            startDay = int(dayArray[recIndex][0])
            endDay = int(dayArray[recIndex][1])
            startDtime = dtimeArray[recIndex][0]
            endDtime = dtimeArray[recIndex][1]
            startHour = int(startDtime)
            endHour = int(endDtime)
            startMin = int(startDtime*60.0 - startHour*60.0)
            endMin = int(endDtime*60.0 - endHour*60.0)
            startSec = int(startDtime*3600.0) % 60
            endSec = int(endDtime*3600.0) % 60
            startTime = datetime.datetime(startYear, startMonth, startDay, startHour, startMin, startSec)
            endTime = datetime.datetime(endYear, endMonth, endDay, endHour, endMin, endSec)

            for beamIndex in range(self.numBeams):
                beamId = beamCodeArray[beamIndex][0]
                az = beamCodeArray[beamIndex][1]
                el = beamCodeArray[beamIndex][2]
                txpower = txPowerArray[recIndex]/1000.0
                numtxaeu = numTxAeuArray[recIndex]
                numrxaeu = numRxAeuArray[recIndex]
                if RangeTime:
                    rangeValue = rangeArray[recIndex][0]
                else:
                    rangeValue = rangeArray[0]
                popValue = popArray[recIndex][beamIndex]
                dpopValue = dpopArray[recIndex][beamIndex]

                inputs.append((kinst,kindat,startTime,endTime,numUsedRanges,az,el,beamId,txpower,numtxaeu,
                               numrxaeu,baudCount,pulseLength,txFreq,rxFreq,lowerRangeIndex,upperRangeIndex,
                               rangeValue,popValue,dpopValue)
                             )

        print("There are %d jobs to execute..." % (len(inputs)))

        # do one to get the size
        dataRec = set_uncorrected_data_rec(inputs[0])
        encoded = dill.dumps(dataRec)
        datarecsize = len(encoded)
        dtype = '|S%d' % (datarecsize)
        output = np.zeros((len(inputs),), dtype=dtype)

        shm = shared_memory.SharedMemory(create=True, size=output.nbytes)
        shm_name = shm.name
        output_array = np.ndarray((len(inputs),),dtype=dtype,buffer=shm.buf)

        worker_pool = mp.Pool(POOLSIZE)
        jobs = [worker_pool.apply_async(mp_wrapper,args=(set_uncorrected_data_rec,i,args,shm_name,len(inputs),dtype)) for i,args in enumerate(inputs)]

        results = [p.get() for p in jobs]

        for i in range(len(inputs)):
            if i % 100 == 0:
                print("Dumping data: %d/%d" % (i + 1,len(inputs)))
            data_record = dill.loads(output_array[i])
            cedarObj.append(data_record)
            # dump records every 100
            if i % 50 == 0: 
                cedarObj.dump()

        # dump remaining records
        cedarObj.dump()
        shm.close()
        shm.unlink()
        print("This took %s seconds" % (str(datetime.datetime.now()-now)))
        print("Closing out the file...")
        now = datetime.datetime.now()
        cedarObj.close()
        self.numRecs = self.numTimes * self.numBeams
        print("This took %s seconds" % (str(datetime.datetime.now()-now)))


def mp_wrapper(func,index,inputs,output_buffer_name,numjobs,dtype):
    if index % 100 == 0:
        print("Executing job: %d" % (index + 1))
    dataRec = func(inputs)
    if dataRec is None:
        return
    encoded = dill.dumps(dataRec)
    datarecsize = len(encoded)

    existing_shm = shared_memory.SharedMemory(name=output_buffer_name)
    array = np.ndarray((numjobs,),dtype=dtype,buffer=existing_shm.buf)
    array[index] = encoded
    existing_shm.close()


def set_uncorrected_data_rec(inputs):
    (kinst,kindat,startTime,endTime,numUsedRanges,az,el,beamId,txpower,numtxaeu,
     numrxaeu,baudCount,pulseLength,txFreq,rxFreq,lowerRangeIndex,upperRangeIndex,
     rangeValue,popValue,dpopValue) = inputs

    dataRec = madrigal.cedar.MadrigalDataRecord(kinst,
                                                kindat,
                                                startTime.year,
                                                startTime.month,
                                                startTime.day,
                                                startTime.hour,
                                                startTime.minute,
                                                startTime.second,
                                                startTime.microsecond/10000,
                                                endTime.year,
                                                endTime.month,
                                                endTime.day,
                                                endTime.hour,
                                                endTime.minute,
                                                endTime.second,
                                                endTime.microsecond/10000,
                                                ('azm', 'elm', 'beamid', 'power',
                                                 'numtxaeu', 'numrxaeu',
                                                 'cbadl', 'pl',
                                                 'tfreq', 'rfreq'),
                                                ('range', 'popl', 'dpopl'),
                                                numUsedRanges,
                                                ind2DList=['range'])

    # set 1d values
    dataRec.set1D('azm', az)
    dataRec.set1D('elm', el)
    dataRec.set1D('beamid', beamId)
    dataRec.set1D('power', txpower) # cedar in kWatts, SRI in Watts
    dataRec.set1D('numtxaeu', numtxaeu)
    dataRec.set1D('numrxaeu', numrxaeu)
    dataRec.set1D('cbadl', baudCount)
    dataRec.set1D('pl', pulseLength)
    dataRec.set1D('tfreq', txFreq)
    dataRec.set1D('rfreq', rxFreq)
    # set 2d values
    for rangeIndex in range(lowerRangeIndex, upperRangeIndex):
        # range
        try:
            if np.isnan(rangeValue[rangeIndex]):
                raise ValueError()
            dataRec.set2D('range', rangeIndex-lowerRangeIndex,rangeValue[rangeIndex]/1000.0) # convert m -> km
        except ValueError:
            dataRec.set2D('range', rangeIndex-lowerRangeIndex, 'missing')
        # pop
        try:
            if np.isnan(popValue[rangeIndex]):
                raise ValueError('popl is NaN')
            dataRec.set2D('popl', rangeIndex-lowerRangeIndex,np.log10(popValue[rangeIndex]))
            if dpopValue[rangeIndex] <= 0.0 or np.isnan(dpopValue[rangeIndex]):
                raise ValueError('problem with dpopl')
            dataRec.set2D('dpopl', rangeIndex-lowerRangeIndex, np.log10(dpopValue[rangeIndex] * popValue[rangeIndex]))
        except ValueError:
            dataRec.set2D('popl', rangeIndex-lowerRangeIndex, 'missing')
            dataRec.set2D('dpopl', rangeIndex-lowerRangeIndex, 'missing')

    return dataRec

def set_fitted_data_rec(inputs):
    (kinst,kindat,startTime,endTime,az,el,beamId,txpower,numtxaeu,
     numrxaeu,baudCount,pulseLength,txFreq,rxFreq,rangeValue,neValue,
     dneValue,fitsValue,errValue,chisqValue,snrValue,platValue,
     plongValue,limits,indices) = inputs

    (minTemp,maxTemp,maxTempErr,minTempErr,minNe,maxNe,maxNeErr,
     minNeErr,maxVo,maxVoErr,minVoErr,maxFract,minFract) = limits

    (o_index,e_index,fractIndex,tempIndex,colIndex,velIndex) = indices

    numRanges = rangeValue.size
    # first we need to find the valid ranges
    valid_range_indices = [rangeIndex for rangeIndex in range(numRanges) if not np.isnan(rangeValue[rangeIndex])]
    if len(valid_range_indices) == 0:
        # no valid ranges
        return None

    dataRec = madrigal.cedar.MadrigalDataRecord(kinst,
                                                kindat,
                                                startTime.year,
                                                startTime.month,
                                                startTime.day,
                                                startTime.hour,
                                                startTime.minute,
                                                startTime.second,
                                                startTime.microsecond/10000,
                                                endTime.year,
                                                endTime.month,
                                                endTime.day,
                                                endTime.hour,
                                                endTime.minute,
                                                endTime.second,
                                                endTime.microsecond/10000,
                                                ('azm', 'elm', 'beamid', 'power', 'numtxaeu', 'numrxaeu', 'cbadl', 'pl', 'tfreq', 'rfreq'),
                                                ('range', 'ne', 'dne', 'ti', 'dti', 'te', 'dte', 'vo', 'dvo', 'cgm_lat', 'cgm_long',
                                                 'po+', 'dpo+', 'chisq', 'sn'),
                                                len(valid_range_indices),
                                                ind2DList=['range'])


    # set 1d values
    dataRec.set1D('azm', az)
    dataRec.set1D('elm', el)
    dataRec.set1D('beamid', beamId)
    dataRec.set1D('power', txpower) # cedar in kWatts, SRI in Watts
    dataRec.set1D('numtxaeu', numtxaeu)
    dataRec.set1D('numrxaeu', numrxaeu)
    dataRec.set1D('cbadl', baudCount)
    dataRec.set1D('pl', pulseLength)
    dataRec.set1D('tfreq', txFreq)
    dataRec.set1D('rfreq', rxFreq)

    # set 2d values
    for row, rangeIndex in enumerate(valid_range_indices):
        # range
        dataRec.set2D('range', row, rangeValue[rangeIndex]) 

        # ne
        try:
            if neValue[rangeIndex] < minNe or neValue[rangeIndex] > maxNe or np.isnan(neValue[rangeIndex]):
                raise ValueError()
            dataRec.set2D('ne', row, neValue[rangeIndex])
            if dneValue[rangeIndex] <= minNeErr or dneValue[rangeIndex] > maxNeErr or np.isnan(dneValue[rangeIndex]):
                raise ValueError()
            dataRec.set2D('dne', row, dneValue[rangeIndex])

        except ValueError:
            dataRec.set2D('ne', row, 'missing')
            dataRec.set2D('dne', row, 'missing')

        # ti
        try:
            if (fitsValue[rangeIndex][o_index][tempIndex] < minTemp or
                fitsValue[rangeIndex][o_index][tempIndex] > maxTemp or
                np.isnan(fitsValue[rangeIndex][o_index][tempIndex])):
                raise ValueError()
            dataRec.set2D('ti', row, fitsValue[rangeIndex][o_index][tempIndex])

            if (errValue[rangeIndex][o_index][tempIndex] < minTempErr or
                errValue[rangeIndex][o_index][tempIndex] > maxTempErr or
                np.isnan(errValue[rangeIndex][o_index][tempIndex])):
                raise ValueError()
            dataRec.set2D('dti', row, errValue[rangeIndex][o_index][tempIndex])
        except ValueError:
            dataRec.set2D('ti', row, 'missing')
            dataRec.set2D('dti', row, 'missing')

        # te
        try:
            if (fitsValue[rangeIndex][e_index][tempIndex] < minTemp or
                fitsValue[rangeIndex][e_index][tempIndex] > maxTemp or
                np.isnan(fitsValue[rangeIndex][e_index][tempIndex])):
                raise ValueError()
            dataRec.set2D('te', row, fitsValue[rangeIndex][e_index][tempIndex])

            if (errValue[rangeIndex][e_index][tempIndex] < minTempErr or
                errValue[rangeIndex][e_index][tempIndex] > maxTempErr or
                np.isnan(errValue[rangeIndex][e_index][tempIndex])):
                raise ValueError()
            dataRec.set2D('dte', row, errValue[rangeIndex][e_index][tempIndex])
        except ValueError:
            dataRec.set2D('te', row, 'missing')
            dataRec.set2D('dte', row, 'missing')

        # vo
        try:
            if (np.isnan(fitsValue[rangeIndex][o_index][velIndex]) or
                abs(fitsValue[rangeIndex][o_index][velIndex]) > maxVo or
                np.isnan(errValue[rangeIndex][o_index][velIndex]) or
                errValue[rangeIndex][o_index][velIndex] > maxVoErr or
                errValue[rangeIndex][o_index][velIndex] < minVoErr):
                raise ValueError()
            dataRec.set2D('vo', row, fitsValue[rangeIndex][o_index][velIndex])
            dataRec.set2D('dvo', row, errValue[rangeIndex][o_index][velIndex])

        except ValueError:
            dataRec.set2D('vo', row, 'missing')
            dataRec.set2D('dvo', row, 'missing')

        # po+
        try:
            if (np.isnan(fitsValue[rangeIndex][o_index][fractIndex]) or
                fitsValue[rangeIndex][o_index][fractIndex] > maxFract or
                fitsValue[rangeIndex][o_index][fractIndex] < minFract):
                raise ValueError()
            dataRec.set2D('po+', row, fitsValue[rangeIndex][o_index][fractIndex])
            dataRec.set2D('dpo+', row, 'assumed')                        
        except ValueError:
            dataRec.set2D('po+', row, 'missing')
            dataRec.set2D('dpo+', row, 'missing')

        # chisq
        if not np.isnan(chisqValue[rangeIndex]):
            dataRec.set2D('chisq', row, chisqValue[rangeIndex])
        else:
            dataRec.set2D('chisq', row, 'missing')
            
        # snr
        if not np.isnan(snrValue[rangeIndex]):
            dataRec.set2D('sn', row, snrValue[rangeIndex])
        else:
            dataRec.set2D('sn', row, 'missing')

        # cgm_lat
        try:
            if np.isnan(platValue[rangeIndex]):
                raise ValueError()
            dataRec.set2D('cgm_lat', row, platValue[rangeIndex]) 
        except ValueError:
            dataRec.set2D('cgm_lat', row, 'missing')

        # cgm_long
        try:
            if np.isnan(plongValue[rangeIndex]):
                raise ValueError()
            dataRec.set2D('cgm_long', row, plongValue[rangeIndex]) 
        except ValueError:
            dataRec.set2D('cgm_long', row, 'missing')

    return dataRec

class hdf5Handler:
    """hdf5Handler is a class calls other classes depending on the hdf5 type.  Presetly supports the following
    types: standard, velocity
    """
    def __init__(self, hdf5Type):
        """__init__ creates a new hdf5Handler of a given type.
        Inputs: hdf5Type - string representing hdf5 file type to handle.  For now must be either
        standard or velocity
        """
        if hdf5Type.lower() == 'standard':
            self.__type = 'standard'
        elif hdf5Type.lower() == 'velocity':
            self.__type = 'velocity'
        elif hdf5Type.lower() == 'uncorrected_ne_only':
            self.__type = 'uncorrected_ne_only'
        elif hdf5Type.lower() == 'velocityalt':
            self.__type = 'velocityAlt'
        else:
            raise ValueError('Unknown hdf5 file type %s' % (str(hdf5Type)))


    def getStartEndTimes(self, hdf5File):
        """getStartEndTimes returns a tuple of (earliest datetime, latest datetime) for a given hdf5 file.
        Calls correct class based on self.__type
        """
        if self.__type in ['standard', 'velocity', 'velocityAlt', 'uncorrected_ne_only']:
            o = analyzeHdf5(hdf5File)
            return o.getStartEndTimes()

        raise ValueError('Unknown self.__type %s' % (str(self.__type)))

    def createMadrigalFile(self,hdf5File,kinst,kindat,cedarObj,madrigalFile,lowerRange=None,upperRange=None):
        """__init__ will write or update a Madrigal file using data in hdf5File using class set by self.__type

        Inputs:
            hdf5File - full path to hdf5 file in SRI format
            kinst - instrument code (integer)
            kindat - data file kindat (integer)
            cedarObj - existing madrigal.cedar.MadrigalCedarFile to append data to.
                       If None, new madrigal.cedar.MadrigalCedarFile
                       created using madrigalFile.
            madrigalFile - name of Madrigal file to create or append to.
            lowerRange - lower range cutoff. If None (the default), no lower range cutoff.  Only affects
                         uncorrected_ne_only.
            upperRange - upper range cutoff. If None (the default), no upper range cutoff.  Only affects
                         uncorrected_ne_only.

        """
        if self.__type == 'standard':
            o = hdf5ToMadrigal(hdf5File,kinst,kindat,cedarObj,madrigalFile)
            self.numRecs = o.numRecs
            return

        elif self.__type == 'velocity':
            o = hdf5VelocityToMadrigal(hdf5File,kinst,kindat,cedarObj,madrigalFile)
            self.numRecs = o.numRecs
            return

        elif self.__type == 'uncorrected_ne_only':
            o = hdf5UncorrectedToMadrigal(hdf5File,kinst,kindat,cedarObj,madrigalFile,lowerRange,upperRange)
            self.numRecs = o.numRecs
            return

        elif self.__type == 'velocityAlt':
            o = hdf5VelocityAltToMadrigal(hdf5File,kinst,kindat,cedarObj,madrigalFile)
            self.numRecs = o.numRecs
            return

        raise ValueError('Unknown self.__type <%s>' % (str(self.__type)))

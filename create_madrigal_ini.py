#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
#
#   Title: create_madrigal_ini.py
#
#   Creates a config file that is used by the SRI Database and Madrigal
#   conversion scripts to ingest processed and derived data products into the
#   databases.
#
#   2021-10-25  Ashton Reimer
#               Initial implementation
#   2022-10-21  Pablo Reyes
#               Fixed multiple day discovery of vvels
#   2022-10-21  Pablo Reyes
#               Adding a file name sorting by seconds. converts mins to secs
#   2022-12-02  Pablo Reyes
#               determine status and category from files for versioning
#   2023-02-15  Pablo Reyes
#               Add *_Ne[!P]*png to the plot patter to avoid finding at the
#               same time Ne and NePower_NoTr
#   2023-02-16  Pablo Reyes
#               Sorting files first by pulse type bc, ac, lp then integration
#   2023-03-06  Pablo Reyes
#               Making the description less verbose 
#
#   Version: 1.0.2021.10.25
#
################################################################################


# Future Development:
#     1) move KINDATS and PLOTS into a config file?

import os
import re
import csv
import glob
import configparser
import h5py
from argparse import ArgumentParser

def fname_seconds(x):
    # old ones:
    #  THIS CASE IS NOT CONSIDERED, ONLY derivedParams/vvelsLat
    #  20090719.001_lp_1min.h5
    #  20090719.001_lp_1min-cal.h5
    #  derivedParams/20090719.001_lp_1min-cal-vvels-60sec.h5
    #
    #  THIS CASE IS NOT CONSIDERED, ONLY derivedParams/vvelsLat
    #  Cascades10/20081108.001/20081108.001_lp_5min.h5
    #  Cascades10/20081108.001/derivedParams/20081108.001_lp_5min-vvels-300sec.h5

    #  IPY17/20120714.001.done/20120714.001_lp_5min-cal.h5
    #  IPY17/20120714.001.done/derivedParams/vvelsLat/20120714.001_lp_5min-cal-vvelsLat-300sec.h5

    #  20100529.002_lp_2min.h5
    #  20120405.003_lp_5min-cal-vvelsLat-300sec.h5
    #
    #  20230223.002_lp_5min-fitcal.h5
    #  derivedParams/vvelsLat/20230223.002_lp_5min-fitcal-vvelsLat-300sec.h5
    numloc = x.rfind("_")
    minloc = x.find("min")
    if minloc >=0:
        return 60*int(x[numloc+1:minloc])
    else:
        secloc = x.find("sec")
        if secloc >= 0:
            return int(x[numloc+1:secloc])
        else:
            return 0

def pulsetype_order(x):
    # e.g. 20230223.002_lp_5min-fitcal-vvelsLat-300sec.h5
    # e.g. 20210113.001_ac_10min-fitcal.h5
    """Assumes x has _bc_, _ac_, _acfl_, _lp_, or _mc_ in it """
    ptype = re.search(r"_[a-z][a-z]_|_[a-z][a-z][a-z][a-z]_",x).group() 
    return ["_bc_","_mc_","_ac_","_acfl_","_lp_"].index(ptype)

class FileParams():
    def __init__(self,hdf5file, kindat_type):
        self.hdf5file = hdf5file
        self.kindat_type = kindat_type
        self.read_params()

    def read_params(self):
        """Read some of the processing parameters
        """
        if self.kindat_type in ['uncorrected_ne_only', 'standard']:
            print(f"Reading parameters from file {self.hdf5file}")
            with h5py.File(self.hdf5file,'r') as fp:
                self.read_ProcessingParams(fp)
                self.fitter_version \
                        = fp['/ProcessingParams/FittingInfo/Version'][()].decode('latin')
                self.ion_masses = fp['/FittedParams/IonMass'][:]

        elif self.kindat_type == 'velocity':
            with h5py.File(self.hdf5file,'r') as fp:
                self.read_ProcessingParams(fp)
                self.SourceFile = \
                        fp['/ProcessingParams/SourceFile'][()].decode('latin')
                self.IntegrationTime = fp['/ProcessingParams/IntegrationTime'][()]
                self.MaxAlt = fp['/ProcessingParams/MaxAlt'][()]
                self.MinAlt = fp['/ProcessingParams/MinAlt'][()]

    def read_ProcessingParams(self, fp):
        """Read ProcessingParams that are included in fitted and vvels files"""
        self.ProcessingTimeStamp = \
                fp['/ProcessingParams/ProcessingTimeStamp'][()].decode('latin')
        self.PulseLength = fp['/ProcessingParams/PulseLength'][()]
        self.BaudLength = fp['/ProcessingParams/BaudLength'][()]
        self.TxFrequency = fp['/ProcessingParams/TxFrequency'][()]
        self.RxFrequency = fp['/ProcessingParams/RxFrequency'][()]

class MadrigalIni():
    RADARS = {'pfisr': {'name': 'PFISR',
                        'instrument_number': 61,
                        'pi': 'Asti Bhatt',
                        'PIEmail' : 'asti.bhatt@sri.com',
                        'fileAnalyst' : 'Pablo Reyes',
                        'fileAnalystEmail' : 'pablo.reyes@sri.com'},
              'risrn': {'name': 'RISR-N',
                        'instrument_number': 91,
                        'pi': 'Asti Bhatt',
                        'PIEmail' : 'asti.bhatt@sri.com',
                        'fileAnalyst' : 'Pablo Reyes',
                        'fileAnalystEmail' : 'pablo.reyes@sri.com'},
             }

    RANGE_LIMS = {'lp': [100.0,1000.0],
                  'ac': [80.0,1000.0],
                  'bc': [50.0,200.0],
                  'mc': [50.0,200.0],
                  'acfl': [80.0,1000.0],
                 }
    

    PLOTS = {
        'uncorrected_ne_only': {
            'filematch':['*_NePower_NoTr*.png', '*_SNR*.png'],
            'titles':   ['Electron density - No Te/Ti Correction','Signal to Noise Ratio'],
                                },
        'standard': {
            'filematch': ['*_Ne[!P]*png','*_nuin*.png', '*_IonFrac*.png',
                          '*_Te*.png', '*_Ti*.png', '*_Tr*.png', '*_Vlos*.png'],
            'titles':    ['Electron Density','Ion-Neutral Collision Frequency',
                          'O+ Ion Fraction','Electron Temperature','Ion Temperature',
                          'Te/Ti Ratio','LOS Velocity'],
                     },
        'velocity': {
            'filematch': ['*-emag*.png',
                          '*-evec*.png',
                          '*-vmag*.png',
                          '*-vvec*.png'],
            'titles':    ['Electric Field Magnitude and Direction',
                          'Vector Electric Fields',
                          'Velocity Magnitude and Direction',
                          'Vector Velocities',
                          ],
                     },
             }


    def __init__(self,radar,expdir_path, specsfile ):
        # radar information
        self.radar = self.RADARS[radar]['name']
        self.instrument_number = self.RADARS[radar]['instrument_number']
        self.pi = self.RADARS[radar]['pi']
        self.PIEmail = self.RADARS[radar]['PIEmail']
        self.fileAnalyst = self.RADARS[radar]['fileAnalyst']
        self.fileAnalystEmail = self.RADARS[radar]['fileAnalystEmail']

        # determine some path information
        self.expdir_path = expdir_path
        self.expdir = os.path.basename(self.expdir_path)    
        self.radarmode_path = os.path.dirname(self.expdir_path)
        self.radarmode = os.path.basename(self.radarmode_path)
        self.datapath = os.path.dirname(os.path.dirname(self.expdir_path))

        # specsfile if given will have information about category and file description
        # similar to madrigal
        self.specsfile = specsfile
        self.getspecsfile()
        print(self.specsfiledict)

        # initialize parameters used by 'add_file'
        self.filenum = 1

        # unique geoplots
        self.unique_geoplots = []

        # initialize config parser
        self.configfile = configparser.ConfigParser(delimiters=(':'),interpolation=None)
        self.configfile.optionxform = str

    def getspecsfile(self):
        """
        Read a csv file with 3 columns: filename,category,filedescription
        where category is the madrigal :
        file.category (int) (1=default, 2=variant, 3=history, 4=real-time)
        file description:
        file.status (string)('preliminary', 'final', or any other description)
        """
        self.specsfiledict = {}
        if type(self.specsfile) != type(None):
            csvfile = self.specsfile[0]
            with open(csvfile,'r') as fp:
                csvFile = csv.reader(fp)
                for line in csvFile:
                    assert len(line) == 3
                    self.specsfiledict.update({line[0]:dict(category=line[1],
                                                        fileDesc=line[2])})


    def build(self):
        self.read_experiment_description()

        self.add_default()
        self.add_experiment()

        # search for -fitcal files
        # !!!!! WARNING: this logic fails for uncorrected bc/mc files! Fix later
        fitted_h5files = sorted(glob.glob(os.path.join(self.expdir_path,'*-fitcal.h5'))
                            , key =lambda x: (pulsetype_order(x),fname_seconds(x)))
        # if no -fitcal files, look for -cal files
        if len(fitted_h5files) == 0:
            fitted_h5files = sorted(glob.glob(os.path.join(self.expdir_path,'*-cal.h5'))
                            , key =lambda x: (pulsetype_order(x),fname_seconds(x)))

        # if no -fitcal or -cal files, error out.
        if len(fitted_h5files) == 0:
            raise Exception("No '*-fitcal.h5' nor '*-cal.h5' files found.")
        #fitted_h5files = [os.path.basename(x) for x in fitted_h5files]

        #history_h5files = sorted(glob.glob(os.path.join(self.expdir_path,'categ.*.h5')))
        history_h5files = sorted(glob.glob(os.path.join(self.expdir_path,'*-cal.h5')))
        #history_h5files = [os.path.basename(x) for x in history_h5files]

        # look for vvels files
        vvels_h5files = sorted(glob.glob(os.path.join(self.expdir_path,
            'derivedParams/vvelsLat','*.h5'))
                                    , key =lambda x: fname_seconds(x))
        #vvels_h5files = [os.path.basename(x) for x in vvels_h5files]

        # Now let's build the Madrigal.ini file
        h5files = fitted_h5files + vvels_h5files + history_h5files

        file_counter = 0
        num_h5files = len(h5files)

        ptype2do = "bc"
        for h5file in fitted_h5files + history_h5files:
            if os.path.basename(h5file).split('_')[1] == ptype2do:
                print('Working on file %d/%d: %s' % (file_counter+1,
                    num_h5files,h5file))
                file_counter += 1
                self.add_file('uncorrected_ne_only',h5file)

        for ptype2do in ['ac','lp']:
            for h5file in fitted_h5files + history_h5files:
                if os.path.basename(h5file).split('_')[1] == ptype2do:
                    print('Working on file %d/%d: %s' % (file_counter+1,
                        num_h5files,h5file))
                    file_counter += 1
                    self.add_file('uncorrected_ne_only',h5file)
                    self.add_file('standard',h5file)

        for h5file in vvels_h5files:
            print('Working on file %d/%d: %s' % (file_counter+1,
                    num_h5files,h5file))
            file_counter += 1
            self.add_file('velocity',h5file)


    def read_experiment_description(self):
        description_file = '%sExperimentDescription.txt' % self.radarmode

        description_path = os.path.join(self.radarmode_path,description_file)
        if not os.path.exists(description_path):
            raise Exception('%s does not exist' % (description_path))

        with open(description_path,'r') as fid:
            self.desc_title_line = fid.readline().strip('\n')
            self.desc_short_desc_line = fid.readline().strip('\n')
            self.desc_long_desc_line = fid.readline().strip('\n')


    def add_default(self):
        self.configfile.set('DEFAULT','DataPath',self.datapath)
        self.configfile.set('DEFAULT','ExperimentType',self.radarmode)
        self.configfile.set('DEFAULT','ExperimentName',self.expdir)


    def add_experiment(self):
        if not self.configfile.has_section('Experiment'):
            self.configfile.add_section('Experiment')

        self.configfile.set('Experiment','title','%(ExperimentType)s - ' + self.desc_short_desc_line)
        self.configfile.set('Experiment','instrument',str(int(self.instrument_number)))
        self.configfile.set('Experiment','OutPath','%(DataPath)s/%(ExperimentType)s/%(ExperimentName)s/Madrigal')
        self.configfile.set('Experiment','logFile','%(DataPath)s/%(ExperimentType)s/%(ExperimentName)s/Madrigal/MadrigalLog.txt')
        self.configfile.set('Experiment','expID','%(ExperimentName)s')
        self.configfile.set('Experiment','pi',self.pi)
        self.configfile.set('Experiment','PIEmail',self.PIEmail)
        self.configfile.set('Experiment','fileAnalyst',self.fileAnalyst)
        self.configfile.set('Experiment','fileAnalystEmail',self.fileAnalystEmail)
        self.configfile.set('Experiment','modexp', self.desc_title_line)
        self.configfile.set('Experiment','cmodexp', self.desc_long_desc_line)

    def determine_status_category(self, hdf5file):
        status = 'final'
        category = '1'

        # http://cedar.openmadrigal.org/docs/name/rr_webServices.html
        # 4. file.category (int) (1=default, 2=variant, 3=history, 4=real-time)
        # 5. file.status (string)('preliminary', 'final', or any other description)
        #    this later becomes fileDesc in madrigal upload
        #    final is the file that shows by default in madrigal
        # categ_3_prelim07Nov2014_YYYYMMDD.XXX_lp_IIIII-vvelsLat-geo-XXXsec.h5
        # categ_3_prelim01May2014_YYYYMMDD.XXX_lp_IIIII-fitcal.h5

        if hdf5file.split('_')[0] == 'categ':
            _, category, status = hdf5file.split('_')[:3]

        return status, category

    def determine_kindat(self,kindat_type,hdf5file, file_params):

        extend_ckindat = ""
        # extract "sub"-type from hdf5file name
        pulse_type = hdf5file.split('_')[1]
        # bc, lp, ac
        # extract integration time from hdf5file name
        #int_time = hdf5file.split('_')[2].split('-')[0]
        int_time = os.path.splitext(hdf5file)[0].split('_')[2].split('-')[0]
        # e.g. 20230223.002_lp_5min-fitcal.h5 -> 5min
        # e.g. 20100529.002_lp_2min.h5 -> 2min
        if kindat_type in ['uncorrected_ne_only', 'standard']:
            sub_type = pulse_type
            # e.g. ac, lp, bc, mc
        elif kindat_type in ['velocity']:
            # e.g. 20230223.002_lp_5min-fitcal-vvelsLat-300sec.h5
            base_intg = int_time
            int_time = os.path.splitext(hdf5file)[0].split('-')[-1] # final intg
            # e.g. 20120714.001_lp_5min-cal-vvelsLat-300sec.h5 -> 300sec
            #sub_type = hdf5file.split('-')[-2]
            # FIX known typo:
            if base_intg == "3min" and int_time == "60sec":
                int_time = "180sec"
            sub_type = f'{pulse_type}-{base_intg}'
            # e.g .20120714.001_lp_5min-cal-vvelsLat-300sec.h5 -> lp-5min
            # e.g. 20081108.001_lp_5min-vvels-300sec.h5 -> lp-5min
            # e.g. 20081108.001_lp_3min-vvels-180sec.h5 -> lp-3min
            # e.g. 20081108.001_lp_1min-vvels-180sec.h5 -> lp-1min


        # processing and configuration
        # 100 - Non-Fitted Uncorrected Electron Density
        # 125 - fitted underspread data ###(not implemented yet)
        # 200 - Fitted Overspread Data, Standard Config
        # 300 - Resolved Velocity, Standard Latitude Bins
        if kindat_type == 'uncorrected_ne_only':
            pc = 100
            pc_desc = 'Ne From Power'
            extend_ckindat += 'Ne derived from Power, not fitted and uncorrected, i.e. Tr=Te/Ti=1. '
            extend_ckindat += f"Fitter version used: {file_params.fitter_version}. "

        elif kindat_type == 'standard':
            if sub_type in ['lp', 'ac', 'acfl']:
                pc = 200
                pc_desc = 'Fitted'
                extend_ckindat += "Fitted with standard overspread code for ion masses: " \
                        f"{', '.join(file_params.ion_masses.astype(str))}. "
                extend_ckindat += f"Fitter version used: {file_params.fitter_version}. "
            else:
                raise Exception('Unknown/unsupported processed file type '\
                              '"%s" for file: %s' % (sub_type,hdf5file))
        elif kindat_type == 'velocity':
            pc = 300
            pc_desc = 'Resolved Vel. Standard Lat. Bins'
            extend_ckindat += 'Resolved Velocity with standard Latitude Bins. '\
                    f"The source file used to derive this data product was: "\
                    f"{file_params.SourceFile}. "
        else:
            raise Exception('Unknown/unsupported file: %s' % (hdf5file))

        # pulse type
        # 0 - N/A
        # 1 - Long Pulse (F-region)
        # 2 - Alternating Code (E-region, bottom F-region)
        # 3 - Barker/MPS Code (D-region, bottom E-region)
        if kindat_type in ['uncorrected_ne_only', 'standard']:
            if sub_type == 'lp':
                pt = 1
                pt_desc = "Long Pulse (F-region)"
                extend_ckindat += "Pulse type: Uncoded long pulse (LP) to resolve "\
                        "the F-region. "
            elif sub_type in ['ac','acfl']:
                pt = 2
                pt_desc = "Alternating Code (E-region)"
                extend_ckindat += "Pulse type: Alternating Code (AC) to resolve "\
                        "the E-region and the lower F-region. "
            elif sub_type in ['bc','mc']:
                pt = 3
                pt_desc = "Barker/MPS Code (D-region)"
                extend_ckindat += "Pulse type: Barker or other binary code  to resolve "\
                        "the D-region and the lower E-region. "
            else:
                raise Exception('Unknown/unsupported pulse type: %s' % (sub_type))
        elif kindat_type == 'velocity':
            #pt = 0
            pulse_type,base_intg = sub_type.split('-')
            if  pulse_type == 'lp':
                pt = 10
                pt_desc = "F-region"
            elif pulse_type in ['ac','acfl']:
                pt = 30
                pt_desc = "E-region"
            elif pulse_type in ['bc','mc']:
                pt = 50
                pt_desc = "D-region"

            if base_intg == "1min":
                pt += 1
            elif base_intg == "2min":
                pt += 2
            elif base_intg == "3min":
                pt += 3
            elif base_intg == "5min":
                pt += 4
            elif base_intg == "10min":
                pt += 5
            elif base_intg == "15min":
                pt += 6
            elif base_intg == "20min":
                pt += 7

            pt_desc = f"{base_intg}_{pulse_type}_{pt_desc}"
            extend_ckindat += f"This derived product is based on a {pulse_type} ({pt_desc}) "\
                    f"experiment  with integration time {base_intg}. "
        else:
            raise Exception('Unknown/unsupported file: %s' % (hdf5file))

        if kindat_type == 'standard' and not sub_type in ['lp', 'ac', 'acfl']:
            raise Exception('Unknown/unsupported pulse type '
                    '"%s" for processed file type: %s' % (sub_type,hdf5file))

        # integration time
        # 1 - 1 minute
        # 2 - 2 minute
        # 3 - 3 minute
        # 4 - 4 minute
        # 5 - 5 minute
        # 5 - 10 minute
        # 6 - 15 minute
        # 7 - 20 minute
        # 8 - 30 seconds
        # 9 - 15 seconds
        # 10 - 4 seconds
        if int_time in ['1min','60sec']:
            it = 1
            it_desc = '1 minute'
        elif int_time in ['2min','120sec']:
            it = 2
            it_desc = '2 minutes'
        elif int_time in ['3min','180sec']:
            it = 3
            it_desc = '3 minutes'
        elif int_time in ['5min','300sec']:
            it = 4
            it_desc = '5 minutes'
        elif int_time in ['10min','600sec']:
            it = 5
            it_desc = '10 minutes'
        elif int_time in ['15min','900sec']:
            it = 6
            it_desc = '15 minutes'
        elif int_time in ['20min','1200sec']:
            it = 7
            it_desc = '20 minutes'
        elif int_time == '30sec':
            it = 8
            it_desc = '30 seconds'
        elif int_time == '15sec':
            it = 9
            it_desc = '15 seconds'
        elif int_time == '4sec':
            it = 10
            it_desc = '4 seconds'
        else:
            raise Exception('Unknown/unsupported integration time: %s' % (int_time))

        #it_desc = f"Final Integration time: {it_desc}."
        it_desc = f"{it_desc}"
        tkindat = pc * 10000 + pt * 100 + it

        tkindat = str(int(tkindat))
        # Trim out None in the description strings
        desc = [x for x in [pc_desc,pt_desc,it_desc] if not x is None]
        ckindat = " - ".join(desc)

        extend_ckindat += f"The final integration time is {it_desc}. "
        extend_ckindat += f"This data product was generated on "\
                f"{file_params.ProcessingTimeStamp}."

        return tkindat, extend_ckindat, ckindat


    def add_file(self, kindat_type, hdf5file_fullpath):
        """Add h5 file to the madrigal.ini file
        """
        file_params = FileParams(hdf5file_fullpath, kindat_type)

        hdf5file = os.path.basename(hdf5file_fullpath)
        file_title = 'File%d' % (self.filenum)
        if not self.configfile.has_section(file_title):
            self.configfile.add_section(file_title)

        # extract sub-type from hdf5file name
        if kindat_type in ['uncorrected_ne_only', 'standard']:
            # e.g. 20210113.001_ac_10min-fitcal.h5
            sub_type = hdf5file.split('_')[1]
            path_template = '%(DataPath)s/%(ExperimentType)s/%(ExperimentName)s/'
        elif kindat_type in ['velocity']:
            path_template = '%(DataPath)s/%(ExperimentType)s/%(ExperimentName)s/derivedParams/vvelsLat/'

        # write the file information
        tkindat, extend_ckindat, ckindat = self.determine_kindat(
                kindat_type, hdf5file, file_params)
        status, category = self.determine_status_category(hdf5file)
        self.configfile.set(file_title,'hdf5Filename',path_template+hdf5file)
        self.configfile.set(file_title,'kindat',tkindat)
        self.configfile.set(file_title,'extend_ckindat',extend_ckindat)
        self.configfile.set(file_title,'createRecPlots','True')
        self.configfile.set(file_title,'type',kindat_type)
        self.configfile.set(file_title,'ckindat',ckindat)

        self.configfile.set(file_title,'status', status) # 'final')
        self.configfile.set(file_title,'category', category) #,'1')
        self.configfile.set(file_title,'history','')

        if kindat_type == 'uncorrected_ne_only':
            lowerRange = str(self.RANGE_LIMS[sub_type][0])
            upperRange = str(self.RANGE_LIMS[sub_type][1])
            self.configfile.set(file_title,'lowerRange',lowerRange)
            self.configfile.set(file_title,'upperRange',upperRange)

        # now write the plots information
        files2do = self.PLOTS[kindat_type]['filematch']
        titles = self.PLOTS[kindat_type]['titles']

        # figure out some things for plot directory and titles
        file_without_extension = os.path.splitext(hdf5file)[0]
            # '20210113.001_ac_10min-fitcal'
            # vvels, e.g. 20230223.002_lp_5min-fitcal-vvelsLat-300sec.h5
            # 20230223.002_lp_5min-fitcal-vvelsLat-300sec
        # figure out the plots directory name
        if kindat_type in ['velocity']:
            plots_directory = 'derivedParams/vvelsLat'
        else:
            ind = hdf5file.find('_%s' % sub_type)
            # if hdf5file is '20210113.001_ac_10min-fitcal.h5' then
            # sub_type == ac, and ind = 12
            suffix_without_extension = file_without_extension[ind:]
            # '_ac_10min-fitcal'
            plots_directory = 'plots%s' % (suffix_without_extension)
            # should look like 'plots_ac_10min-fitcal'
        #if category == "1" and status == "final":
        if True:
            # only add figures for final category 1 files
            imgnum = 1
            for filematch,title in zip(files2do,titles):
                if kindat_type in ['velocity']:
                    # filematch = "*-emag*.pngi", title = 'Vector Velocities'
                    # this is needed because all the vvels files are in one folder
                    # e.g. derivedParams/vvelsLat/20230223.002_lp_5min-fitcal-vvelsLat-300sec.h5
                    int_time = hdf5file.split('_')[2].split('-')[0]
                    # 5min
                    imgs = sorted(glob.glob(os.path.join(self.expdir_path,plots_directory,
                        "*"+int_time+filematch)))
                else:
                    imgs = sorted(glob.glob(os.path.join(self.expdir_path,plots_directory,filematch)))
                imgs = [os.path.basename(x) for x in imgs]
                for img in imgs:
                    image_title = file_without_extension + ' ' + title
                    image_path_template = '%(DataPath)s/%(ExperimentType)s/%(ExperimentName)s/'
                    image_path = os.path.join(image_path_template, plots_directory, img)

                    self.configfile.set(file_title,'imageTitle%d' % (imgnum),image_title)
                    self.configfile.set(file_title,'image%d' % (imgnum),image_path)

                    imgnum += 1

            # Now add the geometry.png to the first file
            # Nov 14, 2022, P. Reyes. Changing for allowing different geoplots.
            #if self.filenum == 1:
            if kindat_type not in ['velocity']:
                filematch = '*geoplot.png'
                title = 'Geometry Plot'
                imgs = glob.glob(os.path.join(self.expdir_path,plots_directory,filematch))
                if len(imgs) > 0:
                    img = glob.glob(os.path.join(self.expdir_path,plots_directory,filematch))[0]
                    image_size = os.path.getsize(img)
                    if image_size not in self.unique_geoplots:
                        img = os.path.basename(img)

                        image_title = file_without_extension + ' ' + title
                        image_path_template = '%(DataPath)s/%(ExperimentType)s/%(ExperimentName)s/'
                        image_path = os.path.join(image_path_template,plots_directory,img)

                        self.configfile.set(file_title,'imageTitle%d' % (imgnum),image_title)
                        self.configfile.set(file_title,'image%d' % (imgnum),image_path)
                        self.unique_geoplots.append(image_size)

            self.filenum += 1


def main():
    valid_radars = list(MadrigalIni.RADARS.keys())
    # Set up some argparse stuff
    parser = ArgumentParser(description='Creates a "Madrigal.ini" file.')
    parser.add_argument("radar", help="The shortname of the radar.", choices=valid_radars)
    parser.add_argument("expdir", help="Path to the experiment directory, e.g. /Volumes/AMISR_PROCESSED/processed_data/PFISR/2021/01/MSWinds27H.v03/20210113.001")
    parser.add_argument('--specsfile', nargs=1,
            help='specsfile can be used to specify category and file description to a file')
    # Get arguments. Convert argparser Namespace class to dictionary
    args = vars(parser.parse_args())
    
    # input checking
    radar = args['radar']
    expdir = args['expdir']
    if expdir[-1] == '/':
        expdir = expdir[:-1]

    madrigal_ini_file = os.path.join(expdir,'Madrigal.ini')
    madini = MadrigalIni(radar,expdir,args['specsfile'])
    madini.build()

    print('Writing Madrigal.ini file...')
    with open(madrigal_ini_file,'w') as f:
        madini.configfile.write(f)
    print('Done.')

if __name__ == '__main__':
    main()

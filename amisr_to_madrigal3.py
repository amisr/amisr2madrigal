#!/opt/virtualenvs/madenv/bin/python

"""

2007-??-?? M. Nicolls
2018-05-03 Bill Rideout added argsparse and numCPU arg
2021-11-24 Ashton Reimer modified to remove numCPU
2023-03-06 Pablo Reyes. Adding the experimentsDirNum to
           the uploadExperiment function. By default = 0
2023-03-20 Pablo Reyes. Adding the file_version parameter
           to createNewExperimentFromIni and
           uploadExperiment
"""

# standard python imports
import os
import pwd
import argparse

# Madrigal imports
import madrigal3_amisr

##############################


def main():
    # check to make sure the user is transport
    # this is needed because:
    # 1) we mount /Volumes/DATABASES_001 with transport:transport ownership
    # 2) gzip is used by madrigal conversion and gzip tries to change permissions if
    #    the user running the script doesn't match the user of the directory and this
    #    permission change fails causing the madrigal conversion to fail.

    # command line interface
    parser = argparse.ArgumentParser(description='Converts and uploads AMISR fitted '
            'and derived products to the AMISR Madrigal 3 server.')
    parser.add_argument('experiment',help='The full path to the experiment directory, '
            'e.g. /path/to/20210203.001, which must contain a Madrigal.ini file.')
    parser.add_argument('--upload',action='store_true',
            help='set for upload, rather than create from ini')
    parser.add_argument('--keeptmp',action='store_true',
            help='set for keeping the tmp files, otherwise they are removed.')
    parser.add_argument('--skip_existing',action='store_true',
            help='Skip the creation of madrigal files if already exists.')
    parser.add_argument('--skip_doc_plots',action='store_true',
            help='Skip the creation of extra docs and plots.')
    parser.add_argument('--no_transport_user',action='store_true',
            help='Do the upload as a user other than transport.')
    parser.add_argument('--DirNum', nargs='?', const=None, type=int, default=0, 
            help = 'sufix added to the folder experiments, default=0, '
            'e.g. DirNum=0 -> experiments0, if no arguments then DirNum==None -> experiments')
    parser.add_argument('--file_version', nargs='?', const=None, type=int, default=None,
            help = 'File version added to the filenames .XXX.h5. If file_version==None'
            ' then file_version will be incremented by one from the previous existing file.'
            'If no file exists then it starts with 1. In the --upload mode when file_version'
            'is None(default) then files with the latest version are used. default=None' )
    args = parser.parse_args()

    # make sure input experiment path exists
    experiment = os.path.abspath(args.experiment)
    if not os.path.exists(experiment):
        raise Exception("%s does not exist!" % (experiment))

    # make sure the Madrigal.ini file exists
    madrigal_ini = os.path.join(experiment,'Madrigal.ini')
    if os.path.exists(madrigal_ini):
        print("Processing experiment using %s " % madrigal_ini)
    else:
        raise Exception("%s does not exist!" % (madrigal_ini))

    # try to create the Madrigal directory inside the experiment directory
    madrigal_directory = os.path.join(experiment,'Madrigal')
    if not os.path.exists(madrigal_directory):
        try:
            os.mkdir(madrigal_directory)
        except:
            raise Exception("Unable to create the 'Madrigal' directory: %s" % (
                madrigal_directory))

    # Now take the experiment files and convert or upload them to Madrigal 3.
    mad_batch = madrigal3_amisr.BatchExperiment()
    if args.keeptmp:
        removeTmpFiles = False
    else:
        removeTmpFiles = True
    if args.skip_existing:
        skip_existing = True
    else:
        skip_existing = False
    if args.skip_doc_plots:
        skip_doc_plots = True
    else:
        skip_doc_plots = False
    if args.no_transport_user:
        no_transport_user = True
    else:
        no_transport_user = False

    if args.upload:
        executing_user = pwd.getpwuid(os.getuid())[0]
        if not no_transport_user and executing_user != 'transport':
            raise Exception('To upload to madrigal this code must be executed'
                    ' using the "transport" user.')
        print("Uploading to Madrigal 3...")
        mad_batch.uploadExperiment(madrigal_ini,
                file_version=args.file_version,
                removeTmpFiles=removeTmpFiles,
                experimentsDirNum=args.DirNum)
    else:
        print("Converting to Madrigal 3...")
        mad_batch.createNewExperimentFromIni(madrigal_ini,
                skip_existing=skip_existing,
                skip_doc_plots=skip_doc_plots,
                file_version=args.file_version)
    print("Done.")


if __name__ == '__main__':
    main()

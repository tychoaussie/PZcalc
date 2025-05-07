# -*- coding: utf-8 -*-

'''
The MIT License (MIT)

Copyright (c) 2013 Daniel Burk
during my time on campus during a summer of delightful albeit, uncompensated labor.
Michigan State University.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''


__author__ = "Daniel Burk <burkdani@msu.edu>"
__version__ = "20250507"
__license__ = "MIT"

# -*- coding: utf-8 -*-
# 20250507 - line 473 - Add in a line to correct the digital decimation block input sample rate so it will correctly parse.
# 20250430 - Fix scipy.signal.zpk2tf on line 749 to comply with new python signal package requiremetns
# 20240514 - Fix problem with dataless seed generation where it placed the dataless in the CWD rather than the target folder.
# 20221107 - Add experimental frequency range 4. Add a print statement outlining version number for the log files.
# 20221103 - Set up code to run either with a six channel or a three channel configuration, and fix a bug that improperly sets sensitivity_gain stage 0
#            when inputting a three-channel .CAL file. Remove the extra underscore within the dataless SEED output file name.
#            Update the dataless SEED templates , both three channel and six channel versions.
#            Adjust the output plots to increase size.
# 20220819 - Report both Vo and Vm in the calibration. 
# 20220123 - Set program up to handle six channels at a time
# 20210222 - CHange file name to include start date for batch processing of multiple cals from same station.
# 20210211 Properly compute AO correction factor to use Vm rather than Vo
# Also break out the SAC polezero files into separate files, one per channel.
#
# 20200811 Include a default output of .sacpz in addition to the dataless seed.
# 20200702 version 1.4
# Major revision on how the poles and zeros are calculated.
# Implement suggestion from John Nablek to find a direct three zero, four pole solution.
# Correct the AO normalization factor and the calculation of phase for the plot function.
# 
# 20200526 version 1.3.3
# Include changes made such as the enhanced plot that shows response and phase as a function of frequency and degrees of phase shift
#
# 20200513 version 1.3.2
# In resolving phase, it was determined that SKMcalc function was incorrectly calculating phase
# and was advancing the phase calculation by 180 degrees. Also, the pole/zero solution
# is reduced to four poles, four zeros.
# 20200417 version 1.3.1
# Export a csv of the original curves that lists the original calculated amplitude/phase data along with the poles/zeros calculated amplitude/phase data
# so that it can be compared to test waveform data.
# 20200309 version 1.3
# Assume the seismogram is a record of velocity, not displacement, and generate the response
# as appropriate. This is because response removal with version 1.2 leads to a sign reversal
# in the signal. Also, move the PZcalc algorithm to four poles and zeros, rather than five.
#
# 20200114 version 1.2
# Add the ability to generate dataless SEED files from an enhanced CAL file.
# This requires an additional supplemental template file to reside in the Pyscripts
# folder: C:/Pyscripts/dataless.pzcalc_template.seed

# 20190828 version 1.1
# Add the AO correction factor to the calculation and the file output

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.optimize
import os,csv,string,sys
from time import gmtime,strftime
from obspy.io.xseed import Parser
from obspy import read_inventory

# Now, the most important part -- The legalese:
# COPYRIGHT ©  BOARD OF TRUSTEES OF MICHIGAN STATE UNIVERSITY
# ALL RIGHTS RESERVED

# PERMISSION IS GRANTED TO USE, COPY, COMBINE AND/OR MERGE, CREATE DERIVATIVE
# WORKS AND REDISTRIBUTE THIS SOFTWARE AND SUCH DERIVATIVE WORKS FOR ANY PURPOSE,
# SO LONG AS THE NAME OF MICHIGAN STATE UNIVERSITY IS NOT USED IN ANY ADVERTISING
# OR PUBLICITY PERTAINING TO THE USE OR DISTRIBUTION OF THIS SOFTWARE WITHOUT 
# SPECIFIC, WRITTEN PRIOR AUTHORIZATION.  IF THE ABOVE COPYRIGHT NOTICE OR ANY
# OTHER IDENTIFICATION OF MICHIGAN STATE UNIVERSITY IS INCLUDED IN ANY COPY OF 
# ANY PORTION OF THIS SOFTWARE, THEN THE DISCLAIMER BELOW MUST ALSO BE INCLUDED.

# THIS SOFTWARE IS PROVIDED AS IS, WITHOUT REPRESENTATION FROM MICHIGAN STATE
# UNIVERSITY AS TO ITS FITNESS FOR ANY PURPOSE, AND WITHOUT WARRANTY BY MICHIGAN
# STATE UNIVERSITY OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT
# LIMITATION THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE.

# THE MICHIGAN STATE UNIVERSITY BOARD OF TRUSTEES SHALL NOT BE LIABLE FOR ANY
# DAMAGES, INCLUDING SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
# WITH RESPECT TO ANY CLAIM ARISING OUT OF OR IN CONNECTION WITH THE USE OF
# THE SOFTWARE, EVEN IF IT HAS BEEN OR IS HEREAFTER ADVISED OF THE POSSIBILITY
# OF SUCH DAMAGES.

 

#	PZcalc V1.3 inputs the six published parameters from the Soviet calibration parameter booklet
# and converts these parameters into a response curve plot as well as a best estimate
# of a four pole/five zero solution for a single channel. It will generate plots that can be saved
# and included for publication in a station equipment history file. As an option, a text file
# can be specified in cases where a three-channel solution is desired. This is the preferred
# method for most analog stations with multiple channels.

# Syntax1:     PZcalc <-chan> <net.stationname.channel> <freqrange> <Ts> <Ds> <Tg> <Dg> <S^2> <Vo>
    
# 
# where:#<netname.stationname.channel> is the station code (ISC preferred) and channel for this station
######     (i.e. SKM-NS, SKM-EW, SKM-Z)
### <freqrange> is the range number. There are three ranges, 1, 2, or 3.
### <Ts> (parameter from column 6) = Sensor free period (in seconds)
### <Ds> (parameter col. 7) = sensor damping constant
### <Tg> (Param. col. 8) = Galvonometer free period (in seconds)
### <Dg> (param. col. 9) = Galvonometer damping constant
### <S^2> (param. col. 10) = Sigma squared, a convolution of seismometer and galvo free periods.
### <Vo> (column 12) = peak amplitude at 90% max peak (Vm) at two frequencies bracketing the peak.

# Syntax2:     PZcalc <-file> <filename.cal>

# where:#<filename.cal> is a single file name for processing. 
###  Output will be placed within the folder containing this file.
###  Using this method enables you to input the station calibration date.

# Syntax3:     PZcalc <-crawl> <targetfolder>

#  where:     <targetfolder> will contain at least one .cal file for processing. 
###  Output files will be placed within the target folder using station name as part of the file name.

# Syntax4:     PZcalc <-seed> <filename.cal>

# where:#<filename.cal> is a single file name for processing, including an additional fifteen fields
###  necessary for the generation of a dataless seed file.
###  Output will be placed within the folder containing this file.
###  Using this method enables you to input the station calibration date.

### 
#  PZCalc will generate two plots, as well as an output file called "<stationname>.paz" containing the 
#  poles and zeros for the station, along with the sensitivity and AO correction factor. 
#  These plots and text may then be placed within the station file.
#  


# Useage requirements: You must have installed both Python version 3 
# and NumPy on your machine in order to run this package. Dependencies
# include the os, csv, time,matplotlib, and Scipy.

# Typical useage: -chan option: (single station, single channel)
# c:\> python PXCalc.py -chan <net.Station.channel> <freqrange> <Ts> <Ds> <Tg> <Dg> <S^2> <Vo>

# As a single channel example, using -chan option:
# c:\Pyscripts> python PZcalc.py -chan RY.UURS.MHZ 2 1.0 0.54 0.37 1.71 0.16 19010
# 
# <freqrange> consists of three different frequency bands over which the calibration should occur:
##1 = 0.1   to 30 seconds: for long period sensors such as the SKD (CKD)
##2 = 0.125 to 10 seconds: For short period sensors, like Vegik, SKM, SM3, S-13
##3 = 0.01  to 5 seconds: For really short period sensors such as the C-5-C2
# --

# Typical useage: -file option: (file name for processing a single cal file:)
# c:\> python PZCalc.py -file <filename> where, <filename> is a single station with a ".cal" extension.


# Typical useage -crawl option: for processing multiple calibration files:
# c:\> python PZCalc.py -crawl <targetfolder>
##where,
##<targetfolder> is the folder containing the calibration files ending in a .cal extension.

#  The output files will be placed within the same folder. They will be of the type .png and .paz
#  --# 
#  Sample file format:
#  The input file that is used consists of a station line, and a channel line for each channel within the station.
#  The file name of the station should consist of <networkcode>.<ISCstationcode>.cal for single stations.
#  If desired, the program will crawl the entire folder looking for cal parameter files with the extension .cal
#  The program will create plots for each station showing all three axes, sensitivity vs time (in seconds) to replicate
#  the analog sensitivity curve as found in published booklets of the former USSR.
#  The program will also create a series of poles and zeros files for each station.
#  Each .CAL should match the following format. Each cal file should have no less than one line.

#  <netname_stationname.channel> <DD/MM/YYYY> <Ts> <Ds> <Tg> <Dg> <S^2> <Vo>

#  As an example, for filename RY_UURS.CAL, where analog channel will be digitized at 100 msec/pixel or less:

#  RY_UURS.MHN 06/04/1987 1 1.63 0.530 0.360 1.840 0.160 53140
#  RY_UURS.MHE 06/04/1987 1 1.06 0.615 0.330 1.710 0.229 27260
#  RY_UURS.MHZ 06/04/1987 1 1.00 0.540 0.370 1.710 0.160 19010
#  -------
#  As another example, if you want to generate a dataless seed file for the station, 
#  add these additional fifteen fields to the calibration file using the following template:
#  Good practive in file naming is to differentiate the cal file by adding SEED to the filename,
#  such as: RY_UURS_SEED.CAL for the following example.

#  * Note that a space is used as the field separator within the following fields. 

#  Use the -seed option in lieu of the -file option to use this cal file format.
#  You can add up to six channels to the calibration to account for both a short-period (like SKM,SM3) 
#  and long-period instrument (like SK,SKD) for the current epoch.

#  beginning_time 1949-08-29T00:00:00.000000Z 
#  end_time 1991-12-26T23:59:59.999999Z
#  network_id RY
#  station_description UURS
#  instrument_description SKM3
#  network_code RY
#  station_code UURS
#  location_identifier
#  site_name Ust-Urkima
#  latitude 55.31
#  longitude 123.16
#  elevation 540
#  start_effective_date 1949-08-29T00:00:00.000000Z
#  end_effective_date 2020-01-14T23:59:59.000000Z
#  sample_rate 100
#  RY.UURS.HHN 06/04/1987 1 1.63 0.530 0.360 1.840 0.160 53140
#  RY.UURS.HHE 06/04/1987 1 1.06 0.615 0.330 1.710 0.229 27260
#  RY.UURS.HHZ 06/04/1987 1 1.00 0.540 0.370 1.710 0.160 19010

def options():      # Get command line options and process them
    fileprocess = False
    chanprocess = False
    seed = False
    channel = []
    filelist = []
    if len(sys.argv) > 1:

        if 'help' in sys.argv[1].lower():
            print(PZcalc.__doc__)

        elif 'crawl' in sys.argv[1].lower(): # Process multiple cal seed files, all located within the same folder.
                    # there should be one additional option: The destination folder.
            try:
                if len(sys.argv) > 2:
                    fileprocess = True
                    seed = True
                    filelist = calfiles(sys.argv[2])
                else:
                    print("-crawl requires an additional argument representing the target folder.\n")
                    print("Type 'PZcalc -help' for help.")
            except:
                print(sys.exc_info(),"\n")
                print('Error browsing the folder {0}.\n -crawl requires a valid folder entry.\n\n'.format(sys.argv[2]))
                print("Type 'PZcalc -help' for help.")

        elif ('file' in sys.argv[1].lower()) or ('seed' in sys.argv[1].lower()): 
                    # There should be one additional option: The destination file.
            try:
                if len(sys.argv) > 2:
                   if 'seed' in sys.argv[1].lower():
                       seed = True
                   fileprocess = True                
                   filelist = []
                   filelist.append(sys.argv[2].replace("/","\\"))
                   print("Processing file: {0}".format(filelist))
                    # Return the file as a file list with one item.
                    # Then process the list of files.
 
                else:
                   print("-file or -seed requires an additional argument representing the target folder.\n")
                   print("Type 'PZcalc -help' for help.")
            except:
                print(sys.exc_info(),"\n")
                print('Error browsing the folder {0}.\n -crawl requires a valid folder entry.\n\n'.format(sys.argv[2]))
                print("Type 'PZcalc -help' for help.")

        elif 'chan' in sys.argv[1].lower(): 
                    # There should be eight additional options.
            try:
                if len(sys.argv) > 9:
                    chanprocess = True
                    channel = []
                    channel.append(sys.argv[2])         # Component = sys.argv[2]
                    channel.append(strftime("%Y_%m_%d",gmtime()))
                    channel.append(int(sys.argv[3]))    # Range = sys.argv[3]
                    channel.append(float(sys.argv[4]))  # (Ts = float(sys.argv[4])
                    channel.append(float(sys.argv[5]))  # (Ds = float(sys.argv[5])
                    channel.append(float(sys.argv[6]))  #  Tg = float(sys.argv[6])
                    channel.append(float(sys.argv[7]))  #  Dg = float(sys.argv[7])
                    channel.append(float(sys.argv[8]))  #  S2 = float(sys.argv[8])
                    channel.append(float(sys.argv[9]))  #  Vo = float(sys.argv[9])
                    print('Process single channel for station channel {0}'.format(channel[0]))
                else:
                    print('Processing single channel requires nine arguments.')
                    print("Type 'PZcalc -help' for help.")
                
            except:
                print(sys.exc_info(),"\n") 
                print("Error processing for single channel. \n -chan option requires eight additional parameters. \n\n")
                print(PZcalc.__doc__)

                    # Return the channel information and process it as if it's
                    # a file with only one item.
        else:

            print('PZcalc requires one of three options:\n -help, -chan, -file, or -crawl as the first command line option.')
            print("Type 'PZcalc -help' for more help.")
    else:
        print('PZcalc requires at least one option: -help, -chan, -file, or -crawl.')
        print("Type 'PZcalc -help' for more help.")

    return(fileprocess,filelist,chanprocess,channel,seed)




def get_dataless_block(metadata):
    #                               preload the dataless seed fields from the calibration
    p = Parser("C:/Pyscripts/dataless.pzcalc_template.seed")
    blk = p.blockettes
    #                               Basic changes to existing fields necessary to customize the dataless seed file
    blk[10][0].beginning_time = metadata['beginning_time'] 
    blk[10][0].end_time = metadata['end_time']                  # These appear to be unused within pdcc
    blk[11][0].station_identifier_code = metadata['network_id'] # This is overritten by blk[50].network_code
    blk[33][0].abbreviation_description = metadata['station_description']
    blk[33][1].abbreviation_description = metadata['instrument_description']
    # metadata['location_identifier']
    blk[50][0].network_code = metadata['network_code']
    blk[50][0].station_call_letters = metadata['station_code']
    blk[50][0].site_name = metadata['site_name']
    blk[50][0].latitude = metadata['latitude']
    blk[50][0].longitude = metadata['longitude']
    blk[50][0].elevation = metadata['elevation']
    blk[50][0].start_effective_date = UTCDateTime(metadata['start_effective_date'])
    blk[50][0].end_effective_date = metadata['end_effective_date']
    return(blk)



#
#                                        Load the metadata-enabled calibration file
#
def dataless_load(infile):
    with open(infile,'r') as fin:
        list = csv.reader(fin)
        rowcnt = 0
        stack = []
        header = []
        Channel = []
        metadata = {}              
        for row in list:
            if len(row) > 0:
                r = row[0].split()
                if rowcnt < 15:      # metadata enabled cal file has fifteen fields before the calibration information.
                    header.append(r)
                    rowcnt +=1

                else:
                    channel = []
                    channel.append(r[0])                    # Component
                    channel.append(r[1].replace("/","_"))   # CalDate
                    channel.append(int(r[2]))               # Range 
                    channel.append(float(r[3]))             # Ts  
                    channel.append(float(r[4]))             # Ds  
                    channel.append(float(r[5]))             # Tg  
                    channel.append(float(r[6]))             # Dg  
                    channel.append(float(r[7]))             # S2  
                    channel.append(float(r[8]))             # Vo
                    Channel.append(channel)
                    #stack.append(r)
                    rowcnt+=1
        for data in header:                                 # metadata is a dictionary with these fields:
            if len(data) > 1:
                metadata[data[0]] = data[1]
            else:
                metadata[data[0]] = ""          
                                                            # 'beginning_time'
                                                            # 'end_time' 
                                                            # 'network_id'
                                                            # 'station_description'
                                                            # 'instrument_description'
                                                            # 'network_code'
                                                            # 'station_code'
                                                            # 'location_identifier'
                                                            # 'site_name'
                                                            # 'latitude'
                                                            # 'longitude'
                                                            # 'elevation'
                                                            # 'start_effective_date'
                                                            # 'end_effective_date'
                                                            # 'sample_rate
       
        return(metadata,Channel)
    


#
# Generate the dataless seed file, based on the metadata and the poles & zeros included in Paz
#

def generate_dataless(Paz,Metadata,file):  
        # paz = poles and zeros.
        # paz[0] = poles
        # paz[1] = zeros
        # paz[2] = scale factor
        # paz[3] = AO normalization factor where Scale factor is multiplied by 2PI * (number of poles - number of zeroes)
        # paz[4] = Vo max sensitivity
        # paz[5] = frequency of max sensitivity in Hz
        # paz[6] = evaluation factor ( a measure of how good the estimation is at recreating the original resposne)
        # paz[7] = component name        # paz = poles and zeros.
        # paz[8] = calibration date
    Channel = []
    Caldate = []
    Real_pole = []      # Poles.real
    Imaginary_pole = [] # Poles.imaginary
    Real_zero = []      # Zeros.real
    Imaginary_zero = [] # Zeros.imaginary
    AO_norm = []        # Normalization factor
    Norm_freq = []      # Normalization frequency at which point max sensitivity is reached & normalized gain = 1 
    Vo = []             # Max sensitivity (peak amplification factor between ground motion and its deflection on paper)
    for  paz in Paz: # There should be three channels
        pole_real = []
        pole_imag = []
        for pole in paz[0]:
            pole_real.append(pole.real)
            pole_imag.append(pole.imag)
        zero_real = []
        zero_imag = []
        for zero in paz[1]: #              assemble the lists for real and imaginary parts
            zero_real.append(zero.real)
            zero_imag.append(zero.imag)
        Real_pole.append(pole_real)
        Imaginary_pole.append(pole_imag)
        Real_zero.append(zero_real)
        Imaginary_zero.append(zero_imag)
        AO_norm.append(paz[3])
        Vo.append(paz[4])
        Norm_freq.append(paz[5])
        Channel.append(paz[7][paz[7].rfind('.')+1:]) # Take last section to determine the channel name
        Caldate.append(paz[8])
     #
    #                                Open the template and modify it with the appropriate information
    #
    if len(Paz) == 6:    
        p = Parser("C:/pyscripts/dataless.pzcalc_template6CH.seed")
    elif len(Paz) == 3:
        p = Parser("C:/pyscripts/dataless.pzcalc_template.seed")
    else:
        print("Error! Total number of channels found within Paz must be either three or six channels. No more, no less.")
        return()
    blk = p.blockettes
    # Basic changes to existing fields necessary to customize the dataless seed file
    blk[10][0].beginning_time = Metadata['beginning_time'] 
    blk[10][0].end_time = Metadata['end_time']                  # These appear to be unused within pdcc
    blk[11][0].station_identifier_code = Metadata['network_id'] # is not necessary as it is overridden by blk[50].network_code
    blk[33][0].abbreviation_description = Metadata['station_description']
    blk[33][1].abbreviation_description = 'SKM'
    if len(Paz) == 6:
         blk[33][2].abbreviation_description = 'SKD'
    blk[50][0].network_code = Metadata['network_code']
    blk[50][0].station_call_letters = Metadata['station_code']
    blk[50][0].site_name = Metadata['site_name']
    blk[50][0].latitude = Metadata['latitude']
    blk[50][0].longitude = Metadata['longitude']
    blk[50][0].elevation = Metadata['elevation']
    blk[50][0].start_effective_date = Metadata['start_effective_date']
    blk[50][0].end_effective_date = Metadata['end_effective_date']
    samplerate = Metadata['sample_rate'] # Output from wavetrac is 100 sps so make it so.
    mult = int(len(blk[58])/len(Paz)) # Assume that the length of the block is due to there being SIX channels.

    for i in range(0,len(blk[52])):
        blk[52][i].sample_rate = samplerate
        blk[57][i].input_sample_rate = samplerate # Add this to make the metadata consistent for conversion into stationXML - drb 05/05/2025
    for i, cha in enumerate(Channel):
        # Must explicitly call out the dip and azimuth based on channel name
        print(f' For channel {i}, cha is {cha}')
        if 'Z' in cha[2]:
            print(f'Setting channel {i} to dip angle of -90')
            blk[52][i].dip = -90.0
            blk[52][i].azimuth = 0.0
        if 'N' in cha[2]:
            print(f'Setting channel {i} to azimuth angle of 0.0')
            blk[52][i].dip = 0.0
            blk[52][i].azimuth = 0.0
        if 'E' in cha[2]:
            print(f'Setting channel {i} to azimuth angle of 90.0')
            blk[52][i].dip = 0.0
            blk[52][i].azimuth = 90.0
        blk[52][i].channel_identifier = cha 
        blk[52][i].location_identifier = Metadata['location_identifier']
        blk[52][i].instrument_identifier = int(i)/3+2
        blk[52][i].latitude = blk[50][0].latitude
        blk[52][i].longitude = blk[50][0].longitude
        blk[52][i].elevation = blk[50][0].elevation
        blk[52][i].start_date = blk[50][0].start_effective_date
        blk[52][i].end_date = blk[50][0].end_effective_date
        blk[52][i].sample_rate = samplerate
        blk[53][i].number_of_complex_poles = 4
        blk[53][i].real_pole = Real_pole[i]
        blk[53][i].imaginary_pole = Imaginary_pole[i]
        blk[53][i].real_pole_error = [0, 0, 0, 0]
        blk[53][i].imaginary_pole_error = [0, 0, 0, 0]
        blk[53][i].number_of_complex_zeros = 3
        blk[53][i].real_zero = Real_zero[i]
        blk[53][i].imaginary_zero = Imaginary_zero[i]
        blk[53][i].real_zero_error = [0, 0, 0]
        blk[53][i].imaginary_zero_error = [0, 0,0]
        blk[53][i].A0_normalization_factor = AO_norm[i]
        blk[53][i].normalization_frequency = Norm_freq[i]
        # stage sequence number 1, seismometer gain
        blk[57][i].input_sample_rate = samplerate             # Added to make compatable with IRIS stationxml validator tool. 05/07/2025 - drb
        blk[58][i*mult].sensitivity_gain = Vo[i]
        # stage sequence number 3, digitizer gain
        blk[58][i*mult+2].sensitivity_gain = 100.0 # This is fixed as 100 cm to 1 m, per PNE2SAC output.
        # stage sequence number 0, overall sensitivity
        blk[58][(i+1)*mult-1].sensitivity_gain = Vo[i] * 100 # There are 100 centimeters in a meter.
    outfile = os.path.join(os.path.dirname(file),("dataless."+blk[11][0].station_identifier_code+"_"+ \
                                        blk[50][0].station_call_letters+"_"+Metadata['instrument_description']+"_" \
                                        +Caldate[0]+".seed"))
    p.write_seed(outfile)
    #                        Generate the sacpz file version of this dataless seed file
    inventory = read_inventory(outfile)
    #                        Parse out the three channels and write each out indiviually
    for channel in inventory[0][0]:
        pzfile = os.path.join(os.path.join(os.getcwd()),(inventory[0].code+"."+ \
            inventory[0][0].code+"."+channel.code+"_"+Caldate[0]+".sacpz"))
        inv = inventory.select(station = inventory[0][0].code, channel = channel.code)
        inv.write(pzfile,format = "SACPZ")





def calfiles(infolder):
    filelist = []
    for file in os.listdir(sys.argv[2]):
        if '.cal' in file.lower():
            filelist.append(os.path.join(sys.argv[2].replace("/","\\"),file)) 
    return(filelist)




def load(infile):                        # Read a simple cal file without metadata
    with open(infile,'r') as fin:
        list = csv.reader(fin)
        Channel = []
        Metadata = [] # This list is returned empty
        for row in list:
            if row:
                channel= []
                c = row[0].split()
                channel.append(c[0])                    # Component
                channel.append(c[1].replace("/","_"))   # CalDate
                channel.append(int(c[2]))               # Range 
                channel.append(float(c[3]))             # Ts  
                channel.append(float(c[4]))             # Ds  
                channel.append(float(c[5]))             # Tg  
                channel.append(float(c[6]))             # Dg  
                channel.append(float(c[7]))             # S2  
                channel.append(float(c[8]))             # Vo  
            Channel.append(channel)
    return(Metadata,Channel)





def readcal(infile,seed):
    if seed == True:
        Metadata,Channel = dataless_load(infile)  # seed is enabled so attempt to open the enhanced cal
    else:
        Metadata,Channel = load(infile) # open the normal cal
    return(Metadata,Channel)






def exportcsv(plotchan,paz,outfil): # Export the gain and phase curves
    # Parse out the AO correction factor and sensitivity from the paz
    outfile = outfil+".csv"
    AO_sense = []              # Parse out the AO correction factor and sensitivity from paz
    for i in range(0,len(paz)):
        AO_sense.append([paz[i][3],paz[i][4]])
                               # Write out a normalized curve
    with open(outfile,mode='a+',newline='') as response:
        response_writer = csv.writer(response, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for j,channel in enumerate(plotchan):
            response_writer.writerow(['Channel:              ',channel[0]])
            response_writer.writerow(['AO Correction factor: ',AO_sense[j][0]])
            response_writer.writerow(['Channel sensitivity:  ',AO_sense[j][1]])
            response_writer.writerow(['Frequency (Hz)','PAZ_gain(normalized)','PAZ_phase(degrees)'])	
            for i in range(0,len(channel[2])): # Report the normalized gain by dividing list by max sensitivity
                response_writer.writerow([1/channel[2][i],np.abs(channel[5][i])/AO_sense[j][1],channel[6][i]])			





def freqrange(range,ts):   # range is a value between 1 and 3 where
                        # 1  = 0.05 to 30 seconds ( 0.03 - 20Hz) (medium-long period i.e. SKD)
                        # 2 = 0.067 to 20 seconds  (short period, i.e. SKM,SM3,S1P)
                        # 3 = 0.005 to 10 seconds  (very short period, i.e. geophone)
                        # 4 = experimental range based on the resonance frequency of the seismometer
                        # Range 4 has not yet been tested on a multichannel calibration with both short period and long period instruments. 
    Period = []
    Frequency = []
    Period.append \
        (np.concatenate((np.arange(30.0,2.5,-1.0),np.arange(2.5,1.1,-0.1), \
         np.arange(1.1,0.025,-0.05)),axis=None))
    #
    Period.append \
        (np.concatenate((np.arange(20.0,10.0,-2.5),np.arange(9.0,4.0,-1.5),\
         np.arange(4.0,2.0,-0.5),np.arange(2.0,0.5,-0.1),[0.5,0.333,0.25,0.2,0.125,0.100,0.067,0.04,0.025]),axis=None))
    #
    Period.append \
        (np.concatenate((np.arange(10.0,4.0,-1.0),np.arange(4.0,2.0,-0.5), \
         np.arange(2.0,0.5,-0.1),[0.5,0.333,0.25,0.2,0.125,0.01,0.005]),axis=None))
    #
    Short_period = [20*ts,10*ts,8*ts,6*ts,5*ts,4*ts,3*ts,2*ts,1.5*ts,1.25*ts,1.125*ts,1*ts, \
          0.975*ts,0.95*ts,0.925*ts,0.9*ts,0.875*ts,0.85*ts,0.825*ts,0.8*ts,0.75*ts,0.725*ts,0.7*ts, \
          0.675*ts,0.65*ts,0.625*ts,0.6*ts,0.575*ts,0.55*ts,0.525*ts,0.5*ts,0.475*ts,0.45*ts,0.425*ts, \
          0.4*ts,0.375*ts,0.35*ts,0.325*ts,0.3*ts,0.275*ts,0.25*ts,0.225*ts,0.2*ts,0.175*ts,0.15*ts,0.125*ts, \
          0.1*ts,0.09*ts,0.08*ts,0.07*ts,0.06*ts,0.05*ts,0.04*ts,0.03*ts,0.02*ts,0.01*ts]
    Long_period = [0.009*ts,0.008*ts,0.007*ts,0.006*ts,0.005*ts,0.004*ts,0.003*ts,0.002*ts,0.001*ts]
    if ts > 2.0:
        Period.append(np.concatenate((Short_period,Long_period),axis=None))
    else:
        Period.append(Short_period)

    for period in Period[range-1]:
        Frequency.append(1./period)
    return(Period[range-1],Frequency)



                    # Plot the channel response curve along with the PAZ estimation, then save the plot to disk.
def respplot2(plotchan,outfil):
    # plotchan contains a list for each channel to be plotted where each item consists of:
    # plotchan[][0] = Component: The channel name.
    # plotchan[][1] = Calibration date for the channel
    # plotchan[][2] = Period is the list of periods upon which all calculations are based
    # plotchan[][3] = Gain contains the list of gains for all three axes from the published parameters
    # plotchan[][4] = Phase contains the list of phase delays for all three axes as a function of seconds delay from published parameters
    # plotchan[][5] = inverted_resp is the list of gains from the poles and zeros
    # plotchan[][6] = inverted_phase is the list of phase values in seconds.
    title = "Frequency response and Phase response vs time \n Calculated from published parameters" \
            + "for cal date of "+plotchan[0][1] # use the first channel's date.
    colorwheel = ['blue','orange','green','red','cyan','magenta', \
                  'blue','orange','green','red','cyan','magenta', \
                  'blue','orange','green','red','cyan','magenta', \
                  'blue','orange','green','red','cyan','magenta'  ]
    #      Don't annotate the record with channel names if there's more than six channels. Its too busy.
    XY = [[[0.02,100],[0.02,60.0],[0.02,35],[0.02,600],[0.02,350.0],[0.02,200]], \
          [[0.015,0.4],[0.22,0.4],[3.0,0.4],[0.015,0.4],[0.22,0.4],[3.0,0.4]]]
                # Determine the min/max x axis scale to fit the frequency band.
    maximum=0
    minimum=1.0E+6 # version 1 had conflict with variables named 'max' and 'min'
    for p in plotchan:
        if (np.amax(p[2]) > maximum):
            maximum = np.amax(p[2])
        if (np.amin(p[2]) < minimum):
            minimum = np.amin(p[2])
    plt.figure(figsize=(7.68,10.24))
    plt.axis([0.01,40,0.1,100000])            # plot the amplitude curves for each of the loaded channels and their associated PAZ estimation
    for i in range (0,len(plotchan)):        
        plt.loglog(plotchan[i][2],plotchan[i][3],color=colorwheel[i],lw=5) # Base parameter amplitude curve for the channel
        plt.loglog(plotchan[i][2],np.abs(plotchan[i][5]),color='red',lw=1) # pole/zero estimation amplitude curve for the channel
                # plot the phase curves and their associated PAZ estimation
    for i in range(0,len(plotchan)):
        inverted_phasesec=degree2phase(plotchan[i][6],plotchan[i][2]) # convert phase from degrees to seconds of delay / period 
        plt.loglog(plotchan[i][2],plotchan[i][4],color=colorwheel[i],lw=3) # base parameter phase for the channel
        plt.loglog(plotchan[i][2],inverted_phasesec,color='red',lw=1) # pole/zero estimation of phase for the channel
                # plot the axes
    plt.xlabel('Period [Seconds]')
    plt.ylabel('Amplitude (microns/mm)')
    plt.text(maximum*2.6, 0.1, 'Phase (Seconds)', fontsize=11,
                   rotation=90.0, rotation_mode='anchor')
    plt.grid(True, which="both")
    if len(plotchan)<7:
        for i in range(0,len(plotchan)):
            plt.annotate(plotchan[i][0] + " amp",         xy=XY[0][i],xytext=XY[0][i], color = colorwheel[i])
            plt.annotate(plotchan[i][0]+" phase",xy=XY[1][i],xytext=XY[1][i], color = colorwheel[i])
    plt.annotate("Inverted parameters from calculated poles & Zeros",xy=(minimum*0.2,0.125),xytext=(minimum*0.2,0.125),color='red')
    plt.suptitle(title) 
    plt.savefig(outfil+".png")
    plt.show()                  # Turn this on if you want to open the plot for viewing, panning and zooming. Otherwise its safe to comment it out 






                    # Plot the channel response curve along with the PAZ estimation, then save the plot to disk.
def respplot1(plotchan,outfil):
    # plotchan contains a list for each channel to be plotted where each item consists of:
    # plotchan[][0] = Component: The channel name.
    # plotchan[][1] = Calibration date for the channel
    # plotchan[][2] = Period is the list of periods upon which all calculations are based
    # plotchan[][3] = Gain contains the list of gains for all three axes from the published parameters
    # plotchan[][4] = Phase contains the list of phase delays for all three axes as a function of seconds delay from published parameters
    # plotchan[][5] = inverted_resp is the list of gains from the poles and zeros
    # plotchan[][6] = inverted_phase is the list of phase values in degrees. from poles and zeros
    title = "Frequency response and Phase response vs time \n Calculated from published parameters" \
            + "for cal date of "+plotchan[0][1] # use the first channel's date.
    colorwheel = ['blue','orange','green','red','cyan','magenta', \
                  'blue','orange','green','red','cyan','magenta', \
                  'blue','orange','green','red','cyan','magenta', \
                  'blue','orange','green','red','cyan','magenta'  ]
                # Determine the min/max x axis scale to fit the frequency band.
    maxfreq=0
    minfreq=1.0E+6 # version 1 had conflict with variables named 'max' and 'min'
    for p in plotchan:
        if (np.amax(p[2]) > maxfreq):
            maxfreq = 1/np.amin(p[2]) # maximum frequency based on the period in p
        if (np.amin(p[2]) < minfreq):
            minfreq = 1/np.amax(p[2]) # minimum frequency
    print(f'Minimum frequency = {minfreq} Hz and maximum frequency = {maxfreq}')
    plt.figure(figsize=(6.4,10.8))
    plt.subplot(211)
    plt.axis([minfreq*0.5,maxfreq*2,0.1,100000]) # scale the plot between 0.1 and 100K magnification
#    plt.axis([0.01,10,0.1,100000])            # plot the amplitude curves for each of the loaded channels and their associated PAZ estimation
    for i in range (0,len(plotchan)):
        Freq = [1./ float(plotchan[i][2][j]) for j in range(len(plotchan[i][2]))]        
        plt.loglog(Freq,plotchan[i][3],color=colorwheel[i],lw=5)
        plt.loglog(Freq,np.abs(plotchan[i][5]),color='red',lw=1)
                # plot the phase curves and their associated PAZ estimation
                # plot the axes
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Amplitude (microns/mm)')
#      Don't annotate the record with channel names if there's more than six channels. Its too busy.
    XY = [[[0.02,100],[0.02,60.0],[0.02,35],[0.02,600],[0.02,350.0],[0.02,200]], \
          [[0.015,0.4],[0.22,0.4],[3.0,0.4],[0.015,0.4],[0.22,0.4],[3.0,0.4]]]
    if len(plotchan)<7:
        for i in range(0,len(plotchan)):
            plt.annotate(plotchan[i][0] + " amp",         xy=XY[0][i],xytext=XY[0][i], color = colorwheel[i])
    plt.grid(True, which="both")
    plt.subplot(212)
    plt.axis([minfreq*0.5,maxfreq*2,-100,300]) # scale the plot between 0.1 and 100K magnification
    for i in range(0,len(plotchan)):
        Freq = [1./ float(plotchan[i][2][j]) for j in range(len(plotchan[i][2]))]
        channelphasedeg = phase2degree(plotchan[i][4],plotchan[i][2])  # phase, period
        plt.semilogx(Freq,channelphasedeg,color=colorwheel[i],lw=3) # frequency vs phase in seconds?
        plt.semilogx(Freq,plotchan[i][6],color='red',lw=1)
    plt.title = "Phase Response of instrument"
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Phase(degrees)')
    plt.text(maxfreq*2.6, 0.1, 'Phase (degrees)', fontsize=11,
                   rotation=90.0, rotation_mode='anchor')
    plt.suptitle(title)
    XY = [[[0.02,-270],[0.02,-200],[0.02,0],[0.02,45],[0.02,90],[0.02,125]], \
          [[0.015,0.4],[0.22,0.4],[3.0,0.4],[0.015,0.4],[0.22,0.4],[3.0,0.4]]]    
    if len(plotchan) < 7:
        for i in range(0,len(plotchan)):
            plt.annotate(plotchan[i][0]+" phase",xy=XY[1][i],xytext=XY[1][i], color = colorwheel[i])
    plt.annotate("Inverted parameters from calculated poles & Zeros",xy=(minfreq*0.2,0.125),xytext=(minfreq*0.2,0.125),color='red')
    plt.suptitle(title) 
    plt.savefig(outfil+".png")

    plt.show()                  # Turn this on if you want to open the plot for viewing, panning and zooming. Otherwise its safe to comment it out 



                    # PAZ subroutine that uses the same frequencies as derived from original curve
def pazto_freq_resp(freqs, zeros, poles, scale_fac): 
    b, a = scipy.signal.zpk2tf(zeros, poles, scale_fac)
    if not isinstance(a, np.ndarray) and a == 1.0:
        a = [1.0]
    return scipy.signal.freqs(b, a, freqs * 2 * np.pi)[1] 
            # list element 0 is frequencies
            # list element 1 is the complex amplitudes


                    # import list of complex numbers and return the angle between 270 and -90 degrees
def phasecalc(testresponse):
    testphase = []
    for i in range(0,len(testresponse)): # ,t in testresponse.enumerate():
        tp = np.arctan2(testresponse[i].imag , testresponse[i].real) * 180. / np.pi
        if (testresponse[i].real < 0) & (testresponse[i].imag < 0): # if point lies in quadrant 3
            testphase.append(360. + tp)
        else:
            testphase.append(tp)
    return(testphase)

                    # Convert phase lag from seconds per period into degrees. 
def phase2degree(phase,period):
    phasedeg = []
    for i in range(0,len(phase)):
        phasedeg.append(float(phase[i])/float(period[i])*360.) # at 10sec (0.1 hz) it should be (7.5 sec/10sec)*360 = 270 deg
    return(phasedeg)

                    # Convert phase lag from degrees into seconds/period
def degree2phase(phasedeg,period):
    phasesec = []
    for i in range(0,len(phasedeg)):
        phasesec.append((float(phasedeg[i])/360.0)*float(period[i]))
    return(phasesec)


    # SKMcalc yields a response curve based on the published parameters of an analog seismic station channel.
    # Formulas are from Soviet calibration procedures for analog stations and were adapted to MATLAB code by Viktor
    #  at the Frunze station in Bishkek, Kyrgyzstan.
def SKMcalc(Periods,Ts,Ds,Tg,Dg,S2,Vo):
    
    T6 = Periods        # T6 represents the provided list of periods that are the inverse of frequency.
    m=2*(Ds/Ts+Dg/Tg)   # This is a shorthand variable for calculation
    p=1/Ts**2 + 1/Tg**2 + 4*Ds*Dg*(1-S2)/Ts/Tg  # Shorthand variable
    q=2*(Ds/Ts/(Tg**2)+Dg/Tg/(Ts**2))   # shorthand variable
    s=(1/(Ts**2)) * (1/(Tg**2)) # shorthand variable
    a = m**2-2.*p       # 
    b = p**2-2*m*q+2*s  # 
    c = q**2-2.*p*s     # 
    d = s**2            # All shorthand variables to shorten the formulas
  
                # Calculate the list of amplitudes (in millimeters) as a function of frequency
                # V6 is the instrument's amplitude curve (gain) as a function of period.
    V6=(Vo*2*Dg/Tg) / np.sqrt((np.array(T6)**-2) + a + b*np.array(T6)**2 + c*np.array(T6)**4 + d*np.array(T6)**6)        

                # Phase calculations. Remember that displacement lags velocity by 90 degrees. 
                # But this is galvo displacement vs ground displacement... Different things altogether.
    A1=-((s*np.array(T6)**2-p)*np.array(T6)**2+1) # X axis
    B1=-(q*np.array(T6)**2-m)*np.array(T6)        # Y axis
    gr=np.arctan(abs(np.array(A1)/np.array(B1))) *180./np.pi # gr is phase angle
    gi = []  
    for i in range(0,len(A1)):     # Adapted from Viktor's matlab code

        if   (A1[i] > 0) & (B1[i] > 0): # If point lies in quadrant 1
            gi.append( gr[i] )    # 
        elif (A1[i] < 0) & (B1[i] > 0): # If point lies in quadrant 2
            gi.append(-gr[i])     # 
        elif (A1[i] < 0) & (B1[i] < 0): # if point lies in quadrant 3
            gi.append(180. + gr[i])
        elif (A1[i] > 0) & (B1[i] < 0): # if point lies in quadrant 4
            gi.append(180. - gr[i])            

    phase = [] # phase delay (in seconds)
    for i in range(0,len(T6)):
        phase.append(gi[i]*T6[i]/360.)
                    
    return(V6,phase) # return gain (V6) and phase delay in seconds.

# findpoles:
# Determine the poles of a four pole, three zero solution that describes response of the system.
# Based on equation A1.68 from appendix A of OFR2014-1218 by Peterson&Hutt, 
# Hagiwara's response equation, that
# describes the 4 poles (S^4,S^3,S^2,S, and 1) as functions of the 5 base parameters,
# ds = Damping ratio of the sensor
# ws = time constant of the sensor ts, in radians/second
# dg = Damping ratio of the galvonometer
# wg = time constant of the galvonometer ts, in radians/second
# sigmasq = Coupling coefficient for the system
# 
def findpoles(ts,ds,tg,dg,sigmasq):# Given the base parameters of the calibration,
    # Given 1.0*x^4 + a*x^3 + b*x^2 + c*x + d, return the roots that represent the poles of the solution
    ws = 2*np.pi / ts              # sensor time constant is converted to radians per second 
    wg = 2*np.pi / tg              # Galvanometer time constant is converted to radians per second
    x = 1.0
    a = 2 * ds * ws + 2 * dg * wg
    b = ws * ws + wg * wg + 4 * ds * ws * dg * wg * (1 - sigmasq)
    c = 2 * ds * ws * wg * wg + 2 * dg * wg * ws * ws
    d = ws * ws * wg * wg
    poles = np.roots([x,a,b,c,d])    # Numpy module that solves for the roots of 4th degree equations
    zeros = [0.0+0.0j,0.0+0.0j,0.0+0.0j] # fix all complex zeros to the origin point.
    return(poles,zeros)



    # import channel base parameters and return a single channel response.
def processchannel(channel): 
    Component = channel[0]
    Caldate   = channel[1]
    Period,frequencies = freqrange(channel[2],channel[3]) # Send the range as well as the Ts value of the seismometer in case of range 4 selection
    frequencies = np.array(frequencies)
    Gain,Phase = SKMcalc(Period,channel[3],channel[4],channel[5],channel[6],channel[7],channel[8])
    # Gain will be the convolution of AO_correction factor and sensitivity
    Phasedeg = phase2degree(Phase,Period) # Phase is described as seconds delay at a given period so convert to degrees
    response = np.array(Gain, dtype=np.float32) # This is the baseline that we are trying to match with the estimator.
    print(channel[3],channel[4],channel[5],channel[6],channel[7],channel[8])
    best_poles,best_zeros = findpoles(channel[3],channel[4],channel[5],channel[6],channel[7])
    best_scale_fac = channel[8] # This is Vo from the published calibration parameters.
    inverted_resp = pazto_freq_resp(freqs=frequencies, zeros=best_zeros, poles=best_poles,scale_fac=best_scale_fac)
    inverted_phase = phasecalc(inverted_resp)
    evaluation = 0.0
    # PAZ gain constant is a convolution of the AO normalization factor and the instrument magnification.
    # Create code that properly breaks apart the AO normalization and sensitivity such that the plot equals 1.0 at the peak
    # frequency, and report both AO normalization, peak frequency, and the amplification for use in a proper PZ file. 
    AO_index = np.argmax(np.abs(inverted_resp)) # Find frequency at which Vm is located
    print(f"maximum response amplitude = {np.abs(inverted_resp[AO_index])} at element {np.argmax(np.abs(inverted_resp))}")
    AO_norm = Gain[np.argmax(Gain)]/np.abs(inverted_resp[AO_index])
    Sensefreq = 1./Period[AO_index] # in Hz.
    max_sense = np.max(Gain) # Highest amplitiude on the curve is the maximum sensitivity Vm
    print(f"Vo = {best_scale_fac}  ")
    print (f"Vm = {max_sense} at {Sensefreq:0.3f} Hz. ")
    print (f"max gain = {np.max(Gain)} at {Sensefreq:0.3f} Hz.")
    paz =      []
    paz.append(best_poles)
    paz.append(best_zeros)
    paz.append(max_sense)
    paz.append(AO_norm)
    paz.append(max_sense) # Use best scale factor here.
    paz.append(Sensefreq)
    paz.append(evaluation)
    paz.append(Component)
    paz.append(Caldate)
    plotchan = []
    plotchan.append(Component)
    plotchan.append(Caldate)
    plotchan.append(Period)
    plotchan.append(Gain)       # 
    plotchan.append(Phase)
    plotchan.append(inverted_resp*AO_norm)
    plotchan.append(inverted_phase)
    # Save the poles & zeros to a file, and print the final result to a plot and save it as well.
    return(plotchan,paz)





def pazsave(outfile,Paz):
        # paz = poles and zeros.
        # paz[0] = poles
        # paz[1] = zeros
        # paz[2] = scale factor
        # paz[3] = AO normalization factor
        # paz[4] = computed maximum; but should be scale factor 
        # paz[5] = frequency of max sensitivity in Hz
        # paz[6] = evaluation factor ( a measure of how good the estimation is at recreating the original resposne)
        # paz[7] = component name        # paz = poles and zeros.
        # paz[8] = calibration date
    for  paz in Paz:      
        with open(outfile,'a+') as f:
            f.write("Channel: {}\n".format(paz[7]))
            print("For channel {}:\n".format(paz[7]))
            f.write("Caldate: {}\n".format(paz[8]))
            print("on Calibration date: {}\n".format(paz[8]))
            f.write("ZEROS {}\n".format(len(paz[1]))) # SAC may want an additional pole to convert to displacement
            print("ZEROS: {}".format(len(paz[1]))) # but these cals already represent displacement.
            for zero in paz[1]:
                f.write("{:e} {:e}\n".format(zero.real, zero.imag))
                print("real:{:e} Imaginary:{:e}".format(zero.real, zero.imag))
            f.write("POLES {}\n".format(len(paz[0])))
            print ("\nPOLES {}".format(len(paz[0])))
            for pole in paz[0]:
                f.write("{:e} {:e}\n".format(pole.real, pole.imag))
                print("real:{:e} Imaginary:{:e}".format(pole.real, pole.imag))
            f.write("AO_Normalization_factor: {:e}\n".format(paz[3]))
            print("\nAO_Normalization_Factor {:2.3f}".format(paz[3]))
            f.write("Sensitivity: {:2.1f}\n".format(paz[2]))
            print("Sensor sensitivity {:2.1f}".format(paz[2]))
            f.write("Sensitivity_frequency(Hz): {:2.2f}\n".format(paz[5]))
            print("Sensor sensitivity frequency {:2.2f} Hz".format(paz[5]))
#           f.write("Evaluation_Factor: {:2.1f} \n------\n".format(paz[6]))
#            print("Evaluation factor for this estimate (Less than 12 is good): {:2.1f} \n------\n\n".format(paz[6]))
    spz = "SAC pole-zero file is named %s" % ( outfile )
    print ( "\n" )
    print (spz)    # Save the paz to a file.





def main():             
                        # channel contains the six base parameters.
                        # plotchan[0] = Component name
                        # plotchan[1] = Caldate
                        # plotchan[2] = Period list
                        # plotchan[3] = Gain list
                        # plotchan[4] = Phase list
                        # plotchan[5] = inverted response list
                        # plotchan[6] = inverted phase list
        # paz[0] = poles
        # paz[1] = zeros
        # paz[2] = scale factor
        # paz[3] = AO normalization factor (Not sure if it's in Hz or Radians/sec quite yet)
        # paz[4] = max sensitivity
        # paz[5] = frequency of max sensitivity in Hz
        # paz[6] = evaluation factor ( a measure of how good the estimation is at recreating the original resposne)
        # paz[7] = component name        # paz = poles and zeros.
    print(f"Michigan State University's PZCalc;  version # {__version__}: \n\n")
    fileprocess,filelist,chanprocess,channel,seed = options()
    if fileprocess:
        for file in filelist:
            Metadata,Channel = readcal(file,seed) # multiple channels, and seed gets populated if seed is True
            Plotchan = []
            Paz = []
            for channel in Channel:
                print(f' channel is defined as {channel}')
                pltchan,paz = processchannel(channel)
                Plotchan.append(pltchan)
                Paz.append(paz)
                print("--")
                print(f"Inverted values for {paz[7]} on caldate of {channel[1]}:")
                print(f"Evaluated misfit of phase = {paz[6]:0.3f}")
                print("========================================\n")
            outfil = os.path.join(os.path.dirname(file),(channel[0]+"_"+channel[1]))
            pazsave((outfil+"_paz.txt"),Paz)
            exportcsv(Plotchan,Paz,outfil)       # Export the parameters as a csv file for additional analysis
            if seed:   # Generate a dataless seed file 
                generate_dataless(Paz,Metadata,file)

            respplot2(Plotchan,(outfil+"_response_period"))
            respplot1(Plotchan,(outfil+"_response"))           

    elif chanprocess:
        Metadata = []  # Initialize, because it isnt used but is passed to pazsave
        Plotchan = []
        Paz = []
        pltchan,paz = processchannel(channel)
        Plotchan.append(pltchan)
        Paz.append(paz)
        print("\n========================================")
        print(f"Inverted values for {paz[7]} on caldate of {channel[1]}:")
        print(f"Evaluated misfit of phase = {paz[6]:0.3f} \n")
        outfil = os.path.join(os.getcwd(),(channel[0]+"_"+channel[1]))
        pazsave((outfil+"_paz.txt"),Paz,Metadata,seed)
        respplot2(Plotchan,(outfil+"_response_period"))

    else:
        print('\nNo files processed. No channels processed.')


#    Call the main loop

if __name__ == '__main__':
    main()

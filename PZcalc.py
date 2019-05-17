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
__version__ = "20190516"
__license__ = "MIT"

# -*- coding: utf-8 -*-

# 20190516 Version 1.0
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.optimize
import os,csv,string,sys
from time import gmtime,strftime


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

 

class PZcalc(object):
    '''PZcalc inputs the six published parameters from the Soviet calibration parameter booklet
       and converts these parameters into a response curve plot as well as a best estimate
       of a five pole/five zero solution for a single channel. It will generate plots that can be saved
       and included for publication in a station equipment history file. As an option, a text file
       can be specified in cases where a three-channel solution is desired. This is the preferred
       method for most analog stations with multiple channels.

       Syntax1:     PZcalc <-chan> <net.stationname.channel> <freqrange> <Ts> <Ds> <Tg> <Dg> <S^2> <Vm>
    
       
       where:      <netname.stationname.channel> is the station code (ISC preferred) and channel for this station
                                         (i.e. SKM-NS, SKM-EW, SKM-Z)
                   <freqrange> is the range number. There are three ranges, 1, 2, or 3.
                   <Ts> (parameter from column 6) = Sensor free period (in seconds)
                   <Ds> (parameter col. 7) = sensor damping constant
                   <Tg> (Param. col. 8) = Galvonometer free period (in seconds)
                   <Dg> (param. col. 9) = Galvonometer damping constant
                   <S^2> (param. col. 10) = Sigma squared, a convolution of seismometer and galvo free periods.
                   <Vm> (column 12) = Maximum peak amplitude to be produced on the curve.

       Syntax2:     PZcalc <-file> <filename.cal>

       where:      <filename.cal> is a single file name for processing. 
                    Output will be placed within the folder containing this file.
                    Using this method enables you to input the station calibration date.

       Syntax3:     PZcalc <-crawl> <targetfolder>

        where:     <targetfolder> will contain at least one .cal file for processing. 
                    Output files will be placed within the target folder using station name as part of the file name.
                   
        PZCalc will generate two plots, as well as an output file called "<stationname>.paz" containing the 
        poles and zeros for the station. These plots and text may then be placed within the station file.
        Eventual versions of this code will create a station .xml file for use with modern seismic analysis software.


       Useage requirements: You must have installed both Python version 3 
       and NumPy on your machine in order to run this package. Dependencies
       include the os, csv, time,matplotlib, and Scipy.

       Typical useage: -chan option: (single station, single channel)
       c:\> python PXCalc.py -chan <net.Station.channel> <freqrange> <Ts> <Ds> <Tg> <Dg> <S^2> <Vm>

       As a single channel example, using -chan option:
       c:\Pyscripts> python PZcalc.py -chan RY.UURS.MHZ 2 1.0 0.54 0.37 1.71 0.16 19010
       
       <freqrange> consists of three different frequency bands over which the calibration should occur:
            1 = 0.1   to 30 seconds: for long period sensors such as the SKD (CKD)
            2 = 0.125 to 10 seconds: For short period sensors, like Vegik, SKM, SM3, S-13
            3 = 0.01  to 5 seconds: For really short period sensors such as the C-5-C2
       --

       Typical useage: -file option: (file name for processing a single cal file:)
       c:\> python PZCalc.py -file <filename> where, <filename> is a single station with a ".cal" extension.


       Typical useage -crawl option: for processing multiple calibration files:
       c:\> python PZCalc.py -crawl <targetfolder>
            where,
            <targetfolder> is the folder containing the calibration files ending in a .cal extension.
        The output files will be placed within the same folder. They will be of the type .png and .paz
        --       
        Sample file format:
        The input file that is used consists of a station line, and a channel line for each channel within the station.
        The file name of the station should consist of <networkcode>.<ISCstationcode>.cal for single stations.
        If desired, the program will crawl the entire folder looking for cal parameter files with the extension .cal
        The program will create plots for each station showing all three axes, sensitivity vs time (in seconds) to replicate
        the analog sensitivity curve as found in published booklets of the former USSR.
        The program will also create a series of poles and zeros files for each station.
        Each .CAL should match the following format. Each cal file should have no less than one line.

        <netname_stationname.channel> <DD/MM/YYYY> <Ts> <Ds> <Tg> <Dg> <S^2> <Vm>

        As an example, for filename RY_UURS.CAL, where analog channel will be digitized at 100 msec/pixel or less:

        RY_UURS.MHN 06/04/1987 1 1.63 0.530 0.360 1.840 0.160 53140
        RY_UURS.MHE 06/04/1987 1 1.06 0.615 0.330 1.710 0.229 27260
        RY_UURS.MHZ 06/04/1987 1 1.00 0.540 0.370 1.710 0.160 19010
        
       '''




def calfiles(infolder):
    filelist = []
    for file in os.listdir(sys.argv[2]):
        if '.cal' in file.lower():
            filelist.append(os.path.join(sys.argv[2].replace("/","\\"),file)) 
    return(filelist)




def freqrange(range):   # range is a value between 1 and 3 where
                        # 1  = 0.125 to 30 seconds (medium-long period i.e. SKD)
                        # 2 = 0.125 to 20 seconds  (short period, i.e. SKM,SM3,S1P)
                        # 3 = 0.005 to 10 seconds  (very short period, i.e. geophone)
    Period = []
    Frequency = []
    Period.append \
        (np.concatenate((np.arange(30.0,5.0,-5.0),np.arange(9.0,2.5,-1.5), \
         np.arange(2.5,1.0,-0.1),[0.75,0.666,0.5,0.333,0.25,0.2,0.125]),axis=None))
    Period.append \
        (np.concatenate((np.arange(20.0,10.0,-2.5),np.arange(9.0,4.0,-1.5),\
         np.arange(4.0,2.0,-0.5),np.arange(2.0,0.5,-0.1),[0.5,0.333,0.25,0.2,0.125]),axis=None))
    Period.append \
        (np.concatenate((np.arange(10.0,4.0,-1.0),np.arange(4.0,2.0,-0.5), \
         np.arange(2.0,0.5,-0.1),[0.5,0.333,0.25,0.2,0.125,0.01,0.005]),axis=None))
    for period in Period[range-1]:
        Frequency.append(1./period)
    return(Period[range-1],Frequency)




def minimize(_var,frequencies,response):    # Uses data found in frequencies, and in response. 
    p1r, p1i, p3r, p4r, p5r,z1r,z2r,z3r, scale_fac = _var
    new_resp = pazto_freq_resp(
        freqs=frequencies,

        zeros=np.array([0.0 + 0.0 * 1j,
                        0.0 + 0.0 * 1j,
                        z1r + 0.0 * 1j,
                        z2r + 0.0 * 1j,
                        z3r + 0.0 * 1j], dtype=np.complex128),                        

        poles=np.array([p1r + p1i * 1j,
                        p1r - p1i * 1j,
                        p3r + 0.0 * 1j,
                        p4r + 0.0 * 1j,
                        p5r + 0.0 * 1j], dtype=np.complex128),
        scale_fac=scale_fac)
    return ((np.abs(new_resp) - np.abs(response)) ** 2).sum()




def misfit(a,b):    # Bring in two lists, compare them, 
                    # and present a scalar representing a ratio of how well they fit.
                    # a and b should be the same length.
    d = []
    n = []
    for i in range(0,len(a)):
        n.append((np.abs(a[i]-b[i]))) # Check the difference between the two values
        d.append((np.abs(a[i]+b[i])))
    misfit = np.sum(n)/np.sum(d)
    return(misfit)



def options():      # Get command line options and process them
    fileprocess = False
    chanprocess = False
    channel = []
    filelist = []
    if len(sys.argv) > 1:

        if 'help' in sys.argv[1].lower():
            print(PZcalc.__doc__)

        elif 'crawl' in sys.argv[1].lower():
                    # there should be one additional option: The destination folder.
            try:
                if len(sys.argv) > 2:
                    fileprocess = True
                    filelist = calfiles(sys.argv[2])
                else:
                    print("-crawl requires an additional argument representing the target folder.\n")
                    print("Type 'PZcalc -help' for help.")
            except:
                print(sys.exc_info(),"\n")
                print('Error browsing the folder {0}.\n -crawl requires a valid folder entry.\n\n'.format(sys.argv[2]))
                print("Type 'PZcalc -help' for help.")

        elif 'file' in sys.argv[1].lower(): 
                    # There should be one additional option: The destination file.
            try:
                if len(sys.argv) > 2:
                   fileprocess = True                
                   filelist = []
                   filelist.append(sys.argv[2].replace("/","\\"))
                   print("Processing file: {0}".format(filelist))
                    # Return the file as a file list with one item.
                    # Then process the list of files.
                else:
                   print("-file requires an additional argument representing the target folder.\n")
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
                    channel.append(float(sys.argv[9]))  #  Vm = float(sys.argv[9])
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

    return(fileprocess,filelist,chanprocess,channel)


    
                    # PAZ subroutine that uses the same frequencies as derived from original curve
def pazto_freq_resp(freqs, zeros, poles, scale_fac): 
    b, a = scipy.signal.ltisys.zpk2tf(zeros, poles, scale_fac)
    if not isinstance(a, np.ndarray) and a == 1.0:
        a = [1.0]
    return scipy.signal.freqs(b, a, freqs * 2 * np.pi)[1] 
            # list element 0 is frequencies
            # list element 1 is the complex amplitudes



                    # import list of complex numbers and return the angle between 90 and 270 degrees
def phasecalc(testresponse):  
    testphase = []
    for t in testresponse:
        tp = np.arctan2(t.imag , t.real) * 180. / np.pi
        if tp > 90.:
            tp = tp-360.
        testphase.append(tp - 90.0) # adjust phase to better match what is seen in the real world calibrations
    return(testphase)     



                    # Convert phase lag from seconds per period into degrees. 
def phase2degree(phase,period):
    phasedeg = []
    for i in range(0,len(phase)):
        phasedeg.append(float(phase[i])/float(period[i])*360. - 270.)
    return(phasedeg)



                    # Convert phase lag from degrees into seconds/period
def degree2phase(phasedeg,period):
    phasesec = []
    for i in range(0,len(phasedeg)):
        phasesec.append((float(phasedeg[i]+270.0)/360.0)*float(period[i]))
    return(phasesec)




                    # Readcal returns a list of channel calibration parameters from a calibration file.
def readcal(infile):
    Component = []
    Period = []
    Frequency = []
    Resp = []
    with open(infile,'r') as fin:
        list = csv.reader(fin)
        Channel = []
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
                channel.append(float(c[8]))             # Vm  
            Channel.append(channel)
    return(Channel)



                    # Plot the channel response curve along with the PAZ estimation, then save the plot to disk.
def respplot2(plotchan,outfil):

    # plotchan contains a list for each channel to be plotted where each item consists of:
    # plotchan[][0] = Component: The channel name.
    # plotchan[][1] = Calibration date for the channel
    # plotchan[][2] = Period is the list of periods upon which all calculations are based
    # plotchan[][3] = Gain contains the list of gains for all three axes from the published parameters
    # plotchan[][4] = Phase contains the list of phase delays for all three axes as a function of seconds delay from published parameters
    # plotchan[][5] = inverted_resp is the list of gains from the poles and zeros
    # plotchan[][6] = inverted_phase is the list of phase values in degrees.

    title = "Frequency response and Phase response vs time \n Calculated from published parameters" \
            + "for cal date of "+plotchan[0][1] # use the first channel's date.
    colorwheel = ['blue','orange','green','red','cyan','magneta']
    XY = [[[0.02,100],[0.02,60.0],[0.02,35]],[[0.015,0.2],[0.22,0.2],[3.0,0.2]]]
                # Determine the min/max x axis scale to fit the frequency band.
    max=0
    min=1.0E+6
    for p in plotchan:
        if (np.amax(p[2]) > max):
            max = np.amax(p[2])
        if (np.amin(p[2]) < min):
            min = np.amin(p[2])

    plt.figure()
    plt.axis([min*0.1,max*2,0.1,100000]) # scale the plot between 0.1 and 100K magnification
                # plot the amplitude curves for each of the loaded channels and their associated PAZ estimation
    for i in range (0,len(plotchan)):        
        plt.loglog(plotchan[i][2],plotchan[i][3],color=colorwheel[i],lw=5)
        plt.loglog(plotchan[i][2],np.abs(plotchan[i][5]),color='red',lw=1)
                # plot the phase curves and their associated PAZ estimation
    for i in range(0,len(plotchan)):
        inverted_phasesec=degree2phase(plotchan[i][6],plotchan[i][2]) # convert phase from degrees to seconds of delay / period 
        plt.loglog(plotchan[i][2],plotchan[i][4],color=colorwheel[i],lw=3)
        plt.loglog(plotchan[i][2],inverted_phasesec,color='red',lw=1)
                # plot the axes
    plt.xlabel('Period [Seconds]')
    plt.ylabel('Amplitude (microns/mm)')
    plt.text(max*2.6, 0.1, 'Phase (Seconds)', fontsize=11,
                   rotation=90.0, rotation_mode='anchor')
    plt.grid(True, which="both")
    for i in range(0,len(plotchan)):
        plt.annotate(plotchan[i][0] + " amp",         xy=XY[0][i],xytext=XY[0][i], color = colorwheel[i])
        plt.annotate(plotchan[i][0]+" phase",xy=XY[1][i],xytext=XY[1][i], color = colorwheel[i])
    plt.annotate("Inverted parameters from calculated poles & Zeros",xy=(min*0.2,0.125),xytext=(min*0.2,0.125),color='red')
    plt.suptitle(title) 

    plt.savefig(outfil+".png")
    plt.show()





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
    for i in range(0,len(A1)):
        if   (A1[i] > 0) & (B1[i] > 0): # If point lies in quadrant 1
            gi.append(gr[i] )
        elif (A1[i] < 0) & (B1[i] > 0): # If point lies in quadrant 2
            gi.append(-gr[i])
        elif (A1[i] < 0) & (B1[i] < 0): # if point lies in quadrant 3
            gi.append(180. + gr[i])
        elif (A1[i] > 0) & (B1[i] < 0): # if point lies in quadrant 4
            gi.append(180 - gr[i])            

    phase = [] # phase delay (in seconds)
    for i in range(0,len(T6)):
        phase.append(gi[i]*T6[i]/360.)
                    
    return(V6,phase) # return gain (V6) and phase delay in seconds.



    # import channel base parameters and return a single channel response.
def processchannel(channel): 

    Component = channel[0]
    Caldate   = channel[1]
    Period,frequencies = freqrange(channel[2])
    frequencies = np.array(frequencies)
    Gain,Phase = SKMcalc(Period,channel[3],channel[4],channel[6],channel[6],channel[7],channel[8])
    Phasedeg = phase2degree(Phase,Period)
    response = np.array(Gain, dtype=np.float32)

    evaluation = 1.0E+09 # For evaluating how close the solution is to the original curve
    np.seterr(divide='ignore')
    for z in range(0,128): # iterate 128 times to find the solution that best describes the phase response.
        initial_x=[]
        X0=np.random.random(9)
        #                                Using the minimize function, find the poles & zeros solution that best describes
        #                                the instrument response as found in responses, on frequencies breakpoint "frequencies"
        out = scipy.optimize.minimize(
            fun=minimize,
            args = (frequencies,response), # An important detail that cost me three weeks to discover.  
            method="BFGS",
            x0=X0,
            options={"eps": 1e-10, "maxiter": 1e8}) # defines the step size for the random tests, I think
        x = out.x
        new_poles = np.array([-abs(x[0]) + abs(x[1]) * 1j,
                              -abs(x[0]) - abs(x[1]) * 1j,
                              -abs(x[2]) + 0.0 * 1j,
                              -abs(x[3]) + 0.0 * 1j,
                              -abs(x[4]) + 0.0 * 1j], 
                              dtype=np.complex128)    

     
        new_zeros = np.array([ 0.0 + 0.0 * 1j,
                               0.0 + 0.0 * 1j,
                              x[5] + 0.0 * 1j,
                              x[6] + 0.0 * 1j,
                              x[7] + 0.0 * 1j], dtype=np.complex128)
        new_scale_fac = x[8]
        #              Create the response curve that results from this theoretical new poles and zeroes solution
        inverted_response = pazto_freq_resp(freqs=frequencies, zeros=new_zeros, poles=new_poles,scale_fac=new_scale_fac)    
        inphase = phasecalc(inverted_response)                        # phase from inverted response, listed in degrees
        curvefit = np.sqrt(((np.array(Phasedeg) - np.array(inphase))**2).mean()) # This is the rmse function for misfit
        if (curvefit) < evaluation:
            final_iteration = z
            best_poles=new_poles
            best_zeros=new_zeros
            best_scale_fac=new_scale_fac
            print(f'\nIteration # {z}: Phase misfit reduced to {curvefit:0.3f}')
            evaluation = curvefit
            if evaluation < 3.5:    # Evaluation is a measure of how well the poles and zeros fit the original response & phase.
                break               # Good enough. End the loop early to speed up the process.
        else:
            sys.stdout.write('.')
            sys.stdout.flush()
    print('\n')
    inverted_resp = pazto_freq_resp(freqs=frequencies, zeros=best_zeros, poles=best_poles,scale_fac=best_scale_fac)
    inverted_phase = phasecalc(inverted_resp)
    paz =      []
    paz.append(best_poles)
    paz.append(best_zeros)
    paz.append(best_scale_fac)
    paz.append(evaluation)
    paz.append(Component)
    plotchan = []
    plotchan.append(Component)
    plotchan.append(Caldate)
    plotchan.append(Period)
    plotchan.append(Gain)
    plotchan.append(Phase)
    plotchan.append(inverted_resp)
    plotchan.append(inverted_phase)
    # Save the poles & zeros to a file, and print the final result to a plot and save it as well.
    return(plotchan,paz)

def pazsave(outfile,Paz):
        # paz = poles and zeros.
        # paz[0] = poles
        # paz[1] = zeros
        # paz[2] = scale factor
        # paz[3] = evaluation factor ( a measure of how good the estimation is at recreating the original resposne)
        # paz[4] = component name
    for  paz in Paz:
        with open(outfile,'a+') as f:

            f.write(f"For channel {paz[4]}:\n")
            print(f"For channel {paz[4]}:\n")

            f.write("ZEROS {}\n".format(len(paz[1]) + 1 ))
            print("ZEROS: {}".format(len(paz[1]) + 1 ))
            for zero in paz[1]:
                f.write("{:e} {:e}\n".format(zero.real, zero.imag))
                print("real:{:e} Imaginary:{:e}".format(zero.real, zero.imag))

            f.write("POLES {}\n".format(len(paz[0])))
            print ("\nPOLES {}".format(len(paz[0])))
            for pole in paz[0]:
                f.write("{:e} {:e}\n".format(pole.real, pole.imag))
                print("real:{:e} Imaginary:{:e}".format(pole.real, pole.imag))

            f.write("CONSTANT {:e}\n\n".format(paz[2]))
            print("\nsensor PAZ gain constant {:2.3f}\n------\n".format(paz[2]))

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
                        # paz[0]      = poles list
                        # paz[1]      = zeros list
                        # paz[2]      = scale factor
                        # paz[3]      = evaluation of phase misfit. (smaller is better)
                        # Paz[4]      = Component name

    fileprocess,filelist,chanprocess,channel = options()
    if fileprocess:
        for file in filelist:
            Channel = readcal(file) # multiple channels
            Plotchan = []
            Paz = []
            for channel in Channel:
                print(channel)
                pltchan,paz = processchannel(channel)
                Plotchan.append(pltchan)
                Paz.append(paz)
                print("--")
                print(f"Inverted values for {paz[4]} on caldate of {channel[1]}:")
                print(f"Evaluated misfit of phase = {paz[3]:0.3f}")
                print("========================================\n")
            outfil = os.path.join(os.path.dirname(file),(channel[0]+"_"+channel[1]))
            pazsave((outfil+"_paz.txt"),Paz)
            respplot2(Plotchan,(outfil+"_response"))           

    elif chanprocess: 
        Plotchan = []
        Paz = []
        pltchan,paz = processchannel(channel)
        Plotchan.append(pltchan)
        Paz.append(paz)
        print("\n========================================")
        print(f"Inverted values for {paz[4]} on caldate of {channel[1]}:")
        print(f"Evaluated misfit of phase = {paz[3]:0.3f} \n")
        outfil = os.path.join(os.getcwd(),(channel[0]+"_"+channel[1]))
        pazsave((outfil+"_paz.txt"),Paz)
        respplot2(Plotchan,(outfil+"_response"))

    else:
        print('\nNo files processed. No channels processed.')


#    Call the main loop

if __name__ == '__main__':
    main()

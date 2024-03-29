"""
Module to do the unpacking of IRIS data for use in MT5.
"""
###############################################################################
# Import modules
###############################################################################

import glob
import os
import shutil
import tarfile
from datetime import datetime, timedelta
import numpy as np

# Seismology-related
from obspy import core
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy.io.sac import attach_paz
from obspy.signal.rotate import rotate2zne, rotate_ne_rt
from obspy.taup import TauPyModel


# Helper functions
def get_otherid(trace):
    """
    Function to find corresponding seismogram in E,N pair.
    Takes trace (object) as argument, returns obspy format id.
    """
    network, station, location, channel = getstats(trace)
    if any(x in channel for x in ['E', '2', 'N', '1']):
        if 'E' in channel:
            nor2 = 'N'
        elif '2' in channel:
            nor2 = '1'
        elif 'N' in channel:
            nor2 = 'E'
        else:
            nor2 = '2'
        otherid = '{}.{}.{}.{}{}'.format(network, station, location,
                                         channel[:-1], nor2)
        return otherid
    else:
        print('ERROR: Trace appears to be neither an east or north component')
        print('Trace ID:', trace.id)
        return 0


def getstats(trace: object):
    """
    Shorthand for getting network, station, location and channel info from a
    trace. Takes trace as argument, returns network, station, location
    and channel.
    """
    network = trace.stats.network
    station = trace.stats.station
    location = trace.stats.location
    channel = trace.stats.channel

    return network, station, location, channel


def station_header(trace, pre_arrival_p, pre_arrival_s, after_arrival):
    """
    Helper function to write headers for each station in .INV file
    :param trace: Obspy trace object
    :param pre_arrival_p: Length of time in seismogram before P arrival
    :param pre_arrival_s: Length of time in seismogram before S arrival
    :param after_arrival: Length of seismogram after arrival to include
    :return:
    """
    # Station name, lower case
    sname = trace.stats.station.lower()
    while len(sname) < 4:
        sname += " "
    # Phase from channel:
    channel = trace.stats.channel[-1]
    assert channel in ['Z', 'T'], "Channel must be vertical or transverse!"
    if channel == 'Z':
        channel_num = 1
    else:
        channel_num = 2

    # Check that channel is broadband
    bb_test = trace.stats.channel[0]
    if bb_test != 'B':
        print("Channel for station is not broadband: may cause problems...")
    bb_phase = 4

    # Epicentral distance and azimuth (in degrees)
    epi_dist = trace.stats.sac['gcarc']
    assert 0. <= epi_dist <= 90, "Epicentral dist. < 0.0 or > 90.0"
    azimuth = trace.stats.sac['az']
    assert 0. <= azimuth <= 360, "Azimuth is not a bearing... Check!"

    # Magnification is always in microns
    xmag = 1

    # Start of seismogram, secs after evdt (default 10 seconds before arrival):
    # Find arrival
    arrival = 0.
    seis_start = 0.
    if channel_num == 1:
        # P is always first arrival
        arrival = trace.stats.sac['t0']
        seis_start = round(arrival - pre_arrival_p)
    else:
        # S not always second arrival
        for key in trace.stats.sac:
            if 'kt' in key and trace.stats.sac[key]=="S":
                arrival = trace.stats.sac[key[1:]]
                seis_start = round(arrival - pre_arrival_s)
    # Check that arrival has been read
    assert arrival > 0., 'Failed to read arrival for {}'.format(trace.id)
    # arrival, seconds after start of seismogram
    arr_secs = arrival - seis_start

    # Weight is always 1:
    weight = 1.0

    # Delta_t (sampling interval)
    delta_t = trace.stats.delta

    # Arrival time as integer * delta_t after start of seismogram
    phase_time = int(round(arr_secs / delta_t))

    # Number of data
    n_data = int(round(phase_time + after_arrival / delta_t))

    # Max amplitude of trace (converted from SI to microns)
    # ymax = int(np.round(np.nanmax(np.abs(trace.data * 1.e6))))
    ymax = int(round(max(abs(trace.data * 1.e6))))

    # High-pass filtering is always zero (for now)
    hpfilt = 0.

    header_str = ("{}{:1d} {:1d}" + 2 * "{:8.2f}" + "{:8d}{:8.2f}{:5.1f}" +
                  "{:5d}" * 2 + "{:8d}" + "{:6.2f}" * 3 + "\r\n")

    header = header_str.format(sname[:4], channel_num, bb_phase, epi_dist,
                               azimuth, xmag, seis_start, weight, phase_time,
                               n_data, ymax, delta_t, arr_secs, hpfilt)

    return header


def inv_float(number: float) -> str:
    """
    To turn number into format required by .inv file, and MT5 in general
    :param number: Float
    :return: string representation of number in .inv style
    """
    float_str = ""
    if number < 0.:
        float_str += "-."
    else:
        float_str += "0."

    float_exp = "{:.3E}".format(number)
    decimal = float(float_exp.split('E')[0]) * 0.1
    dec_string = '{:.4f}'.format(decimal)[-4:]
    exponent = float_exp.split('E')[-1]
    if '+' not in exponent:
        exp_int = int(exponent)
    else:
        exp_int = int(exponent[1:])
    if decimal == 0:
        new_exponent = 0
    else:
        new_exponent = exp_int + 1

    if new_exponent < 0:
        exp_string = '{:03d}'.format(new_exponent)
    else:
        exp_string = '+' + '{:02d}'.format(new_exponent)

    full_string = float_str + dec_string + 'E' + exp_string
    return full_string


class Mt5Eq:
    """
    Class to use as base for unpacking MT5 data and converting to .INV
    """
    def __init__(self, name=None, verbose=True):
        """
        :param name: Name of earthquake to avoid confusion.
        :param verbose: How much information to print; currently under-used.
        """
        # Initialize class
        super(Mt5Eq, self).__init__()
        # Set type of class
        self.dtype = 'Mt5Eq'
        self.verbose = verbose
        if self.verbose:
            if name:
                assert isinstance(name, str), 'Name must be string'
                self.name = name
                print("Initializing earthquake class '{}'".format(name))
            else:
                print("Initializing earthquake class (no name)")
        # Initialize important attributes: latitude, longitude
        self.evdt = None
        self.evlo = None
        self.evla = None
        self.depth = None
        self.mag = None
        self.yyyyjjjhhmmss = None
        # Streams that contain data in different stages of processing
        self.stream = None
        self.wwlpbn_st = None
        self.wwlpbn_cut_st = None
        self.inv_stream = None
        # P and S arrival times
        self.p_arr = {}
        self.s_arr = {}
        # Dictionaries for sorting stations by azimuth
        self.p_azi = None
        self.s_azi = None
        # Arrival times at stations for writing to .INV file
        self.inv_p = None
        self.inv_s = None
        return

    def set_earthquake_info(self, year: int, month: int, day: int, hour: int,
                            minute: int, second,
                            evla, evlo, depth, mag):
        """
        Function to set important attributes of earthquake
        Might be worth using exception handling so that errors related to
        accidental zero padding are more obvious.
        :param year: Int, yyyy format
        :param month: Int, no padded zeros
        :param day: int, no padded zeros
        :param hour: int, no padded zeros
        :param minute: int, no padded zeros
        :param second: int or float, no padded zeros
        :param evla: int or float, no padded zeros
        :param evlo: int or float, no padded zeros
        :param depth: int or float, no padded zeros
        :param mag: int or float, no padded zeros
        :return:
        """
        for arg in [year, month, day, hour, minute]:
            assert isinstance(arg, int), 'Variable must be integer'
        for arg in [second, evla, evlo, depth, mag]:
            assert isinstance(arg, (float, int)), 'Variable must be float/ int'

        evdt = datetime(year, month, day, hour, minute, second)
        evdt_utc = core.UTCDateTime(evdt)
        self.evdt = evdt_utc
        self.evlo = float(evlo)
        self.evla = float(evla)
        self.depth = float(depth)
        self.mag = float(mag)
        self.yyyyjjjhhmmss = evdt_utc.strftime('%Y%j%H%M%S')
        return

    def write_eq_info(self, fname=None):
        """
        Function to write out a file in the format ${evname}.txt, since this is
        still used in the later parts of readme-BRIAN. Will attempt to phase
        out as soon as possible.
        """
        # String format for date and time
        txt_dt = self.evdt.strftime('%Y (%j)  %m %d %H %M %S.0')
        str_list = list(map(str, [self.evla, self.evlo, self.depth,
                                  self.mag, self.yyyyjjjhhmmss]))
        txt_line = txt_dt + '   {}    {}   {}  {} {}'.format(*str_list)
        # Generate filename (if not given)
        if fname:
            assert isinstance(fname, str), 'File name must be string'
            outfilename = fname
        else:
            outfilename = '{}/{}.txt'.format(self.yyyyjjjhhmmss,
                                             self.yyyyjjjhhmmss)
        # Write file
        outfileid = open(outfilename, 'w')
        outfileid.write(
            '          Origin time           Lat      Lon     Dp  Mag')
        outfileid.write('\n')
        outfileid.write(txt_line)
        outfileid.close()
        return

    # Main function
    def unpack_eq(self, tar_archive, rot_tol=5.,
                  before_cut=30., after_cut=90.):
        """
        Script to unpack data from zipped tar archive, alternative to first par
        t of README-Brian
        Arguments:
        tar_archive		String containing name of .tgz file downloaded from sac
        year		Event year in yyyy format (integer)
        month		Month in mm format
        day			Day (dd)
        hour		24-hour clock, hh
        minute		mm
        second		ss (no decimals/microseconds)
        evlo			float
        evla			float
        depth 		float or integer, must be positive
        mag			float
        """
        # Make datetime object from input arguments
        print("Datetime in yyyyjjjhhmmss format = {}".format(self.yyyyjjjhhmmss))

        # Change to utc datetime (for cutting)
        evdt = self.evdt

        #######################################################################
        # Extract tgz archive
        #######################################################################

        # If yyyyjjjhhmmss already exists as directory, ask and remove it
        if os.path.exists(self.yyyyjjjhhmmss):
            print("An unpacked folder for this earthquake already exists."
                  " Are you sure you want to replace it?")
            while True:
                yorn = input('Enter y or n:')
                if yorn.lower() not in ("y", "n"):
                    print('Try again (enter y or n):')
                    continue
                if yorn.lower() == 'n':
                    raise Exception('Exiting function to avoid overwriting...')
                else:
                    break

            print("Overwriting soon. There may be arrival files in the folder"
                  " Are you sure you want to remove these?"
                  " Enter 'y' if yes, or if you have saved copies elsewhere."
                  " Enter 'n' if you wish to exit the function")
            while True:
                yorn = input('Enter y or n:')
                if yorn.lower() not in ("y", "n"):
                    print('Try again (enter y or n):')
                    continue
                elif yorn.lower() == 'n':
                    raise Exception(
                        'Exiting function to avoid overwriting...')
                else:
                    break

            shutil.rmtree(self.yyyyjjjhhmmss)
        # Extract tar archive (transparent compression)
        tf_id = tarfile.open(tar_archive, 'r:*')
        tf_id.extractall()
        tf_id.close()

        # Rename extracted file to yyyyjjjhhmmss
        notgz = tar_archive.split('.')[0]
        os.rename(notgz, self.yyyyjjjhhmmss)

        #######################################################################
        # Read in data to memory
        #######################################################################
        stream = self.stream = core.read(self.yyyyjjjhhmmss + '/*.SAC')
        # Check for multiple seismograms; merge if necessary
        stream.merge()
        # Create dictionary of traces with ids as keys:
        id_dic = {tr.id: tr for tr in stream}
        # Create list of trace ids for traces that cause unpacking problems
        problem_ids = []

        #######################################################################
        # Travel times to stations
        #######################################################################

        # Calculate gcarc distances for station (necessary for phase arrivals)
        # Create empty dictionary for station attributes
        stations_dic = {}
        for trace in stream:
            station = trace.stats.station
            if station not in list(stations_dic.keys()):
                # Set sub-dictionary for individual station
                stations_dic[station] = statdic = {}
                # Station location
                stla, stlo = trace.stats.sac['stla'], trace.stats.sac['stlo']
                dist, azimuth, back_az = gps2dist_azimuth(self.evla, self.evlo,
                                                          stla, stlo)
                # dist is in metres; change to degrees
                gcarc = kilometer2degrees(dist / 1000.)
                # assign to dictionary
                statdic['gcarc'] = gcarc
                statdic['az'] = azimuth
                statdic['baz'] = back_az

        # Change trace metadata for all traces
        # Not ideal, but I can't see a better way
        for trace in stream:
            station = trace.stats.station
            statdic = stations_dic[station]
            sacdic = trace.stats.sac
            for key in list(statdic.keys()):
                sacdic[key] = statdic[key]
            sacdic['evdp'] = self.depth
            sacdic['evlo'] = self.evlo
            sacdic['evla'] = self.evla
            sacdic['b'] = trace.stats.starttime - evdt
            sacdic['e'] = trace.stats.endtime - evdt
            sacdic['o'] = 0.
            sacdic['nzyear'] = evdt.year
            sacdic['nzjday'] = evdt.julday
            sacdic['nzhour'] = evdt.hour
            sacdic['nzmin'] = evdt.minute
            sacdic['nzsec'] = evdt.second
            sacdic['nzmsec'] = int(round(evdt.microsecond / 1000.))

        # Calculate arrivals for each station
        model = TauPyModel(model='ak135')
        # Phases to calculate arrival times for
        phase_list = ["P", "S", "PKiKP", "SKS"]
        # Corresponding indices for phases
        phasedic = {phase_list[i]: i for i in range(len(phase_list))}
        # Umbrella dictionary to hold arrival data
        arrivaldic = {}
        for station in sorted(stations_dic.keys()):
            # Epicentral distance
            gcarc = stations_dic[station]['gcarc']
            # Put arrival "object" into umbrella dictionary
            arrivaldic[station] = model.get_travel_times(self.depth,
                                                         gcarc, phase_list)

        # Update SAC data with phase arrival times
        for trace in stream:
            station = trace.stats.station
            for phase in arrivaldic[station]:
                # Sac pick convention
                sac_arrivalkey = 't' + str(phasedic[phase.name])
                # Name of phase in SAC file denoted by 'k' then pick (e.g. t0)
                sac_arrivalname = 'k' + sac_arrivalkey
                trace.stats.sac[sac_arrivalkey] = phase.time
                trace.stats.sac[sac_arrivalname] = phase.name

        ######################################################################
        # Eliminate unpaired E-N and missing-PZ seismograms
        ######################################################################

        # Create dictionary containing the pairs of E-N components
        # Corresponding components are in dictionary twice, key:item reversed
        enpairs = {}
        # List of north files: makes rotation easier
        pairednorth = []

        for trace in stream:
            channel = trace.stats.channel
            if any(x in channel for x in ['E', '2']):
                eid = trace.id
                # ID for corresponding north file
                nid = get_otherid(trace)
                if nid in list(id_dic.keys()):
                    enpairs[eid] = nid
                    enpairs[nid] = eid
                    pairednorth.append(nid)
                else:
                    print('No corresponding N file for', eid)
                    problem_ids.append(eid)

        # Find n files with no corresponding e file
        for trace in stream:
            channel = trace.stats.channel
            if any(x in channel for x in ['N', '1']):
                if trace.id not in list(enpairs.keys()):
                    print('No corresponding E file for', trace.id)
                    problem_ids.append(trace.id)

        """
        Check for existence of PZ files for each SAC file.
        Not happy about this; subject to whims of IRIS name changes.
        No work-around that I can think of, as PZ files are not read into 
        Obspy.
        """

        # Empty dictionary for trace ids and corresponding PZ files
        pz_ids = {}
        # List of PZ files
        pz_ls = glob.glob(self.yyyyjjjhhmmss + '/*/SACPZ*.*')
        # Search for trace corresponding to each PZ file
        for pole_zero in pz_ls:
            pz_split = pole_zero.split('.')
            # Separate parts of id
            netw, station, loc, channel = pz_split[1:5]
            if loc[0].isdigit():
                locstr = loc
            else:
                locstr = ''
            pz_trid = '{}.{}.{}.{}'.format(netw, station, locstr, channel)
            if pz_trid in id_dic and pz_trid not in problem_ids:
                pz_ids[pz_trid] = pole_zero
        # Look for traces without PZ files:
        for trid in list(id_dic.keys()):
            if trid not in pz_ids:
                print('No PZ file for', trid)
                problem_ids.append(trid)

        #######################################################################
        # Rotate E-N pairs to north, where necessary
        #######################################################################

        for nid in pairednorth:
            # Set keys and corresponding traces
            eid = enpairs[nid]
            ntr = id_dic[nid]
            etr = id_dic[eid]
            # Check that E and N file are within 90 degrees of each other
            eaz = etr.stats.sac['cmpaz']
            naz = ntr.stats.sac['cmpaz']
            edata = etr.data
            ndata = ntr.data
            # Dummy z-variable, needs to be a 3-component to be rotated
            zdata = ndata.copy()
            if eaz < naz:
                eaz += 360.
            azdiff = eaz - naz
            if abs(azdiff - 90.) <= rot_tol:
                if naz != 0:
                    try:
                        znew, nnew, enew = rotate2zne(zdata, 0., -90., ndata,
                                                      naz, 0., edata, eaz, 0.)
                        del znew
                        ntr.data[:] = nnew[:]
                        etr.data[:] = enew[:]
                        etr.stats.sac['cmpaz'] = 90.
                        ntr.stats.sac['cmpaz'] = 0.

                    except Exception as error:
                        print('Problem rotating', nid, eid)
                        print(repr(error))
                        print('Removing traces...')
                        problem_ids.append(nid)
                        problem_ids.append(eid)
            else:
                print(eid, nid, "components not orthogonal")
                print(azdiff, "difference in azimuth")
                problem_ids.append(eid)
                problem_ids.append(nid)

        #######################################################################
        # Apply wwlpbn response ###############################################
        #######################################################################
        # Copy stream to leave original broadband alone
        wwlpbn_st = self.wwlpbn_st = stream.copy()
        wwlpbn_ids = {tr.id: tr for tr in wwlpbn_st}
        # Set response for wwlpbn
        newpaz = {'poles': [-0.257 + 0.3376j, -0.257 - 0.3376j,
                            -0.06283 + 0.0j, -0.06283 + 0.0j],
                  'zeros': [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j],
                  'gain': 0.5985275}

        for trace in wwlpbn_st:
            # Check I want to actually filter it:
            if trace.id not in problem_ids:
                # Find corresponding PZ:
                pole_zero = pz_ids[trace.id]
                # Read in polezero from file
                try:
                    attach_paz(trace, pole_zero)
                    oldpaz = trace.stats.paz
                    # Decimate by factor 20
                    trace.resample(1.0)
                    # Replace existing (pz) filter with wwlpbn
                    trace.simulate(paz_remove=oldpaz, paz_simulate=newpaz,
                                   simulate_sensitivity=False)

                except Exception as error:
                    print("Problem attaching PZ file for", trace.id)
                    print(repr(error))
                    problem_ids.append(trace.id)

        #######################################################################
        # Cut seismograms for MT5 #############################################
        #######################################################################
        # Copy filtered data, so that I can cut it
        wwlpbn_cut_st = self.wwlpbn_cut_st = wwlpbn_st.copy()
        cut_ids = {tr.id: tr for tr in wwlpbn_cut_st}

        for trace in wwlpbn_cut_st:
            # Check I want to use it:
            if trace.id not in problem_ids:
                # If BHZ file, find P arrival
                if 'Z' in trace.stats.channel[-1]:
                    cut_pick = trace.stats.sac['t0']
                elif any(x in trace.stats.channel for x in ['E', '2',
                                                            'N', '1']):
                    cut_pick = trace.stats.sac['t1']

                else:
                    print(trace.id, 'does not seem to be any of E, N, Z')
                    problem_ids.append(trace.id)
                    cut_pick = 0.

                cut_pick_time = evdt + cut_pick
                cutstart = cut_pick_time - before_cut
                cutend = cut_pick_time + after_cut
                trace.trim(cutstart, cutend)

        #######################################################################
        # Rotate horizontal seismograms into transverse #######################
        #######################################################################
        inv_stream = self.inv_stream = core.Stream()
        for trace in wwlpbn_cut_st:
            # Check I want to use it:
            if trace.id not in problem_ids:
                if 'Z' in trace.stats.channel[-1]:
                    inv_stream += trace
                elif any(x in trace.stats.channel for x in ['E', '1']):
                    trans_trace = trace.copy()
                    n_id = get_otherid(trace)
                    ntrace = cut_ids[n_id]
                    trace_baz = trace.stats.sac['baz']
                    radial, trans = rotate_ne_rt(ntrace.data, trace.data,
                                                 trace_baz)
                    del radial
                    trans_trace.data = trans
                    trans_trace.stats.channel = 'BHT'
                    # Find corresponding N trace
                    inv_stream += trans_trace

        #######################################################################
        # Write seismograms to MT5 file structure
        #######################################################################
        # Directory paths
        bbpath = '{}/{}_bb/'.format(self.yyyyjjjhhmmss, self.yyyyjjjhhmmss)
        wwlpbnpath = '{}/{}_wwlpbn/'.format(self.yyyyjjjhhmmss,
                                            self.yyyyjjjhhmmss)
        mt5path = '{}/{}_mt5/'.format(self.yyyyjjjhhmmss, self.yyyyjjjhhmmss)
        inv_path = '{}/{}_inv/'.format(self.yyyyjjjhhmmss, self.yyyyjjjhhmmss)

        # Make folders
        for path in [bbpath, wwlpbnpath, mt5path, inv_path]:
            if os.path.exists(path):
                shutil.rmtree(path)
            os.makedirs(path)
        # Write files
        for trid in list(id_dic.keys()):
            if trid not in problem_ids:
                # Look up traces
                bb_tr = id_dic[trid]
                wwlpbn_tr = wwlpbn_ids[trid]
                cut_tr = cut_ids[trid]
                # Read relevant information
                stats = getstats(bb_tr)
                # Add time to others for unpacking list into string
                statlist = list(stats[:])
                channel = statlist[-1]
                statlist.insert(0, self.yyyyjjjhhmmss)
                # Change E-N channels in filenames
                if channel[-1] == '1':
                    statlist[-1] = channel[:-1] + 'N'
                elif channel[-1] == '2':
                    statlist[-1] = channel[:-1] + 'E'
                bbfilename = '{}_{}_{}_{}.{}'.format(*statlist)
                try:
                    bb_tr.write(bbpath + bbfilename, format='SAC')
                    wwlpbn_tr.write(wwlpbnpath + bbfilename +
                                    '_wwlpbn', format='SAC')
                    cut_tr.write(mt5path + bbfilename +
                                 '_wwlpbn_cut', format='SAC')
                except Exception as error:
                    print('Failed to write SAC file for', trid)
                    print(repr(error))

        # Write SAC files for use when making .INV file

        inv_dic = {trace.id: trace for trace in inv_stream}
        for trid in inv_dic.keys():
            if trid not in problem_ids:
                trace = inv_dic[trid]
                # Read relevant information
                stats = getstats(trace)
                # Add time to others for unpacking list into string
                statlist = list(stats[:])
                statlist.insert(0, self.yyyyjjjhhmmss)
                trace_filename = '{}_{}_{}_{}.{}'.format(*statlist)
                try:
                    trace.write(inv_path + trace_filename, format='SAC')
                except Exception as error:
                    print('Failed to write SAC file for .INV', trid)
                    print(repr(error))
        # Write ${evname}.txt for use with the rest of readme-Brian
        self.write_eq_info()
        return

    def get_pands_arrivals(self):
        """
        Function to extract P arrivals from SAC hearder parts of invstream
        :return:
        """
        # Check that I have already unpacked data
        assert self.inv_stream, 'Need to unpack or load data'
        # Check whether arrivals have already been calculated
        if any((self.p_arr, self.s_arr)):
            print('P and S arrivals exist:'
                  'are you sure you want to replace them?')
            while True:
                yorn = input('Enter [y] or n:')
                if yorn.lower() not in ("", "y", "n"):
                    print('Try again (y, n or ENTER)')
                    continue
                else:
                    break

            if yorn == 'n':
                raise Exception('Exiting function to avoid overwriting...')

        # Create blank dictionaries to hold arrival times for each trace
        p_arr = self.p_arr
        s_arr = self.s_arr
        for trace in self.inv_stream:
            channel = trace.stats.channel[-1]
            sacdic = trace.stats.sac
            for key in sacdic.keys():
                if sacdic[key] == 'P' and channel == 'Z':
                    p_offset = sacdic[key[1:]]
                    p_time = self.evdt + p_offset
                    p_arr[trace.id] = p_time
                elif sacdic[key] == 'S' and channel == 'T':
                    s_offset = sacdic[key[1:]]
                    s_time = self.evdt + s_offset
                    s_arr[trace.id] = s_time
        self.sort_ps_by_az()
        return

    def sort_ps_by_az(self):
        """
        Creates dictionaries
        :return:
        """
        assert self.inv_stream, 'Unpack data first.'
        assert all((self.p_arr, self.s_arr)), 'Calculate P and S arrivals!'
        # Make dictionaries containing azimuth and stations
        self.p_azi = {}
        self.s_azi = {}
        # Fill dictionaries by looping through inv_stream:
        for trace in self.inv_stream:
            # Find whether vertical or horizontal channel
            channel = trace.stats.channel[-1]
            if channel == 'Z':
                az = trace.stats.sac['az']
                lat = trace.stats.sac['stla']
                lon = trace.stats.sac['stlo']
                short_name = str(trace.stats.sac['kstnm']).strip()
                self.p_azi[az] = (trace.id, self.p_arr[trace.id], lon, lat, az, short_name)
            elif channel == 'T':
                az = trace.stats.sac['az']
                lat = trace.stats.sac['stla']
                lon = trace.stats.sac['stlo']
                short_name = str(trace.stats.sac['kstnm']).strip()
                self.s_azi[az] = (trace.id, self.s_arr[trace.id], lon, lat, az, short_name)
        return

    def write_pands(self, pname=None, sname=None, organisation='azimuthal'):
        """
        Function to write files with possible stations that could be included
        in the .INV file, with predicted P and S arrival times.
        :param pname: file to write P stations and arrivals to.
        :param sname: file to write S stations and arrivals to.
        :param organisation: Write stations azimuthally or alphabetically
        :return:
        """
        # Check that MT5 object has all the right variables set
        assert self.yyyyjjjhhmmss, 'Need to read in earthquake info...'
        assert all((self.p_arr, self.s_arr)), 'p_arr or s_arr is empty...'
        org_error = 'organisation must be alphabetical or azimuthal'
        assert organisation in ['azimuthal', 'alphabetical'], org_error
        # If P and S filenames not specified, set using default
        if pname:
            p_file_name = pname
        else:
            p_file_name = '{}/{}_Pstations.txt'.format(self.yyyyjjjhhmmss,
                                                       self.yyyyjjjhhmmss)
        if sname:
            s_file_name = sname
        else:
            s_file_name = '{}/{}_Sstations.txt'.format(self.yyyyjjjhhmmss,
                                                       self.yyyyjjjhhmmss)
            # Check whether files exist, and whether I want to overwrite them:
            if os.path.exists(p_file_name):
                print('P and S arrival files exist:'
                      'are you sure you want to replace them?')
                while True:
                    yorn = input('Enter y or n:')
                    if yorn.lower() not in ("y", "n"):
                        print('Try again (enter y or n):')
                        continue
                    else:
                        break
                if yorn == 'n':
                    raise Exception('Exiting function to avoid overwriting...')
        print('Writing P and S arrivals to files:')
        print(p_file_name, s_file_name)

        # Write arrivals to file
        # Open P arrivals file
        p_arr_id = open(p_file_name, 'w')
        # Open S arrivals file
        s_arr_id = open(s_file_name, 'w')
        if organisation == 'azimuthal':
            for key in sorted(self.p_azi.keys()):
                station_info = self.p_azi[key]
                p_arr_id.write('{} {}\n'.format(station_info[0],
                                                station_info[1].isoformat()))
            for key in sorted(self.s_azi.keys()):
                station_info = self.s_azi[key]
                s_arr_id.write('{} {}\n'.format(station_info[0],
                                                station_info[1].isoformat()))
        else:
            for key in self.p_arr:
                p_arr_id.write('{} {}\n'.format(key, self.p_arr[key]))
            for key in self.s_arr:
                s_arr_id.write('{} {}\n'.format(key, self.s_arr[key]))
        p_arr_id.close()
        s_arr_id.close()
        return

    def read_pands(self, p_file_path, s_file_path, limit_number=25, overall_limit: int = 50):
        """
        Opposite of previous function reads in files again (usually altered).
        Edit files to choose stations for inversion.
        :param p_file_path: Path to file containing chosen P stations
        :param s_file_path: Path to file containing chosen S stations
        :param limit_number: Limits P and S seismograms to 25 each.
        :return:
        """
        # Check that limit_number is a positive integer
        assert all((isinstance(limit_number, int), limit_number >= 0))
        # Check whether .inv arrivals exist, and whether to overwrite them.
        if any((self.inv_s, self.inv_p)):
            print("P and S arrivals for inversion already exist:"
                  "are you sure you want to replace them?")
            while True:
                yorn = input('Enter y or n:')
                if yorn.lower() not in ("y", "n"):
                    print('Try again (enter y or n):')
                    continue
                else:
                    break
            if yorn == 'n':
                raise Exception('Exiting function to avoid overwriting...')

        self.inv_s = {}
        self.inv_p = {}

        if self.verbose:
            print('Reading data from files...')
        # Open files, check format
        p_file_id = open(p_file_path, 'r')
        p_data = p_file_id.readlines()
        p_file_id.close()

        s_file_id = open(s_file_path, 'r')
        s_data = s_file_id.readlines()
        s_file_id.close()

        if (len(p_data) + len(s_data)) > overall_limit:
            limit_string = 'Too many {} seismograms ({:d}). Cutting at {:d}.'
            limits = []
            for arrtype, data, other_data in [('P', p_data, s_data), ('S', s_data, p_data)]:
                if all([len(a) > limit_number for a in [p_data, s_data]]):
                    print(limit_string.format(arrtype, len(data), limit_number))
                    limits = [limit_number, limit_number]
                else:
                    if len(data) > limit_number:
                        limits.append(overall_limit - len(other_data))
                        print(limit_string.format(arrtype, len(data), overall_limit - len(other_data)))
                    else:
                        limits.append(len(data))

        else:
            limits = [len(p_data), len(s_data)]

        assert len(limits) == 2
        p_limit, s_limit = limits

        # List of file, dictionary pairs to use
        read_list = [(p_data, self.inv_p, p_limit), (s_data, self.inv_s, s_limit)]
        for data, arrival, limit in read_list:
            for line in data[:limit]:
                lsplit = line.split()
                assert len(lsplit) == 2, 'Line too long{}'.format(line.strip())
                try:
                    trace_id = lsplit[0]
                    arr_time = datetime.strptime(lsplit[1],
                                                 '%Y-%m-%dT%H:%M:%S.%f')
                    arrival[trace_id] = arr_time
                except Exception as error:
                    print(error)
                    print('Problem with line: {}'.format(line.strip()))
                    print('Please check!')

        return

    def write_inv(self, inv_name=None, pre_arr_p=15., pre_arr_s=20.,
                  aft_arr=85.):
        """
        Function to read in sac files and format those from selected stations
        into a .INV format (that MT5 can read).
        Requires newlines to be specified by \r\n for reading into DOS.
        :param inv_name: Non-standard (YYMMDD) name for .inv file
        :param pre_arr_p
        :param pre_arr_s
        :param aft_arr
        :return:
        """
        # Turn date-time info into name of file
        evdt = self.evdt
        # Check that files exist for reading in
        assert self.inv_stream, 'Need to unpack or load data'

        # Name and open .inv file
        if inv_name:
            inv_file_name = inv_name
        else:
            yymmdd = evdt.strftime('%y%m%d')
            inv_file_name = yymmdd + 'A.INV'
        inv_id = open(inv_file_name, 'w')

        # Write header line of INV file, only important when printing summary
        header_date = evdt.strftime('%y%m%d%H%M%S') + '0'

        header_lon_round = round(self.evlo, 2)
        header_lon_str = ('{:2d}'.format(int(header_lon_round)) +
                          str(header_lon_round).split('.')[-1][:2])
        header_lat_round = round(self.evla, 2)
        header_lat_str = ('{:2d}'.format(int(header_lat_round)) +
                          str(header_lat_round).split('.')[-1][:2])

        f_head_str = '{} {}  {}{:3d} 0  0\r\n'.format(header_date,
                                                      header_lat_str,
                                                      header_lon_str,
                                                      int(round(self.depth)))
        inv_id.write(f_head_str)
        # Create dictionary of traces in inv_stream
        inv_dic = {trace.id: trace for trace in self.inv_stream}

        # Loop through P and S stations, writing to file
        for arr_type, arr_dic in (("P", self.inv_p), ("S", self.inv_s)):
            for trid in list(arr_dic.keys()):
                # Write header line with station info
                trace = inv_dic[trid]
                # Cut trace to correct length
                arr = arr_dic[trid]
                if arr_type == "P":
                    pre_arr = pre_arr_p
                else:
                    pre_arr = pre_arr_s
                start_cut = core.UTCDateTime(arr - timedelta(seconds=pre_arr))
                seis_len = aft_arr - 1.
                end_cut = core.UTCDateTime(arr + timedelta(seconds=seis_len))
                write_trace = trace.trim(start_cut, end_cut)
                try:
                    stat_header = station_header(trace, pre_arr_p, pre_arr_s,
                                                 aft_arr)
                    inv_id.write(stat_header)
                except:
                    print("Problem writing {}".format(trid))
                    print("Choose another station?")
                    continue

                # Set response for wwlpbn
                # May want to change this later
                newpaz = {'poles': [-0.257 + 0.3376j, -0.257 - 0.3376j,
                                    -0.06283 + 0.0j, -0.06283 + 0.0j],
                          'zeros': [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j],
                          'gain': 0.5985275}

                # Line with numbers of poles and zeros and gain
                pz_ln1 = '{:4d}{:4d} {}\r\n'.format(len(newpaz['zeros']),
                                                    len(newpaz['poles']),
                                                    inv_float(newpaz['gain']))
                inv_id.write(pz_ln1)
                # Line of zeros:
                pz_line2 = ""
                for zero in list(newpaz['zeros']):
                    pz_line2 += ' {}'.format(inv_float(zero.real))
                    pz_line2 += ' {}'.format(inv_float(zero.imag))
                pz_line2 += '\r\n'
                inv_id.write(pz_line2)
                # line of poles
                pz_line3 = ""
                for pole in list(newpaz['poles']):
                    pz_line3 += ' {}'.format(inv_float(pole.real))
                    pz_line3 += ' {}'.format(inv_float(pole.imag))
                pz_line3 += '\r\n'
                inv_id.write(pz_line3)

                # Loop through trace data, writing to file
                for i in range(0, len(write_trace.data), 8):
                    line = trace.data[i:(i+8)]
                    line_inv = [inv_float(x * 1.e6) for x in line]
                    line_str = ' {}' * len(line) + '\r\n'
                    line_out = line_str.format(*line_inv)
                    inv_id.write(line_out)

        # Close file
        inv_id.close()
        return

    def write_station_map(self, prefix: str = None):
        if prefix is not None:
            p_file = prefix + "_p_station_locs.txt"
            s_file = prefix + "_s_station_locs.txt"
        else:
            p_file = "p_station_locs.txt"
            s_file = "s_station_locs.txt"

        p_out = '{}/{}_{}'.format(self.yyyyjjjhhmmss, self.yyyyjjjhhmmss, p_file)
        s_out = '{}/{}_{}'.format(self.yyyyjjjhhmmss, self.yyyyjjjhhmmss, s_file)

        for fname, dictionary in zip((p_out, s_out), (self.p_azi, self.s_azi)):
            out_id = open(fname, "w")
            for key in sorted(self.p_azi.keys()):
                station_info = self.p_azi[key][2:]
                out_id.write('{:.4f} {:.4f} {:.4f} {}\n'.format(*station_info))




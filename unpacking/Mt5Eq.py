"""
Module to facilitate unpacking data for use in MT5.
"""
###############################################################################
# Import modules
###############################################################################
# Seismology-related
import glob
import os
import shutil
import tarfile
from datetime import datetime

from obspy import core
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy.io.sac import attach_paz
from obspy.signal.rotate import rotate2zne, rotate_ne_rt
from obspy.taup import TauPyModel


# Helper functions
def get_otherid(trace):
    """
    Function to find corresponding seismogram in E,N pair.
    Takes trace as argument, returns obspy format id.
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


class Mt5Eq:
    def __init__(self, name=None, verbose=True):
        """
        Class to use as base for unpacking MT5 data and converting to .INV
        :param name: Name of earthquake to avoid confusion
        """
        # Initialize class
        super(Mt5Eq, self).__init__()
        # Set type of class
        self.dtype = 'Mt5Eq'
        if verbose:
            if name:
                assert type(name) is str, 'Name must be string'
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
        return

    def set_earthquake_info(self, year: int, month: int, day: int, hour: int,
                            minute: int, second,
                            evla, evlo, depth, mag):
        """
        Function to set important attributes of earthquake
        :param year: 
        :param month: 
        :param day: 
        :param hour: 
        :param minute: 
        :param second: 
        :param evla: 
        :param evlo: 
        :param depth: 
        :param mag: 
        :return: 
        """
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
        txt_dt = self.evdt.strftime('%Y (%j)  %m %d %H %M %S.0')
        str_list = list(map(str, [self.evla, self.evlo, self.depth,
                                  self.mag, self.yyyyjjjhhmmss]))
        txt_line = txt_dt + '   {}    {}   {}  {} {}'.format(*str_list)
        if fname:
            assert type(fname) is str, 'File name must be string'
            outfilename = fname
        else:
            outfilename = '{}/{}.txt'.format(self.yyyyjjjhhmmss,
                                             self.yyyyjjjhhmmss)
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
        Script to unpack data from zipped tar archive, altermative to first par
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
        print(("Datetime in yyyyjjjhhmmss format = ", self.yyyyjjjhhmmss))

        # Change to utc datetime (for cutting)
        evdt = self.evdt

        #######################################################################
        # Extract tgz archive
        #######################################################################

        # If yyyyjjjhhmmss already exists as directory, remove it
        if os.path.exists(self.yyyyjjjhhmmss):
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
        stream = core.read(self.yyyyjjjhhmmss + '/*.SAC')
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
                dist, azimuth, back_az = gps2dist_azimuth(stla, stlo,
                                                          self.evla, self.evlo)
                # dist is in metres; change to degrees
                gcarc = kilometer2degrees(dist / 1000.)
                statdic['gcarc'] = gcarc
                statdic['az'] = azimuth
                statdic['baz'] = back_az

        # Change trace metadata for all traces
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
        pz_ls = glob.glob(self.yyyyjjjhhmmss + '/SACPZ*')
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
        wwlpbn_st = stream.copy()
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
                    trace.decimate(5)
                    trace.decimate(4)
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
        wwlpbn_cut_st = wwlpbn_st.copy()
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
        filtered_ps_st = core.Stream()
        for trace in wwlpbn_cut_st:
            # Check I want to use it:
            if trace.id not in problem_ids:
                if 'Z' in trace.stats.channel[-1]:
                    filtered_ps_st += trace
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
                    filtered_ps_st += trans_trace

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

        filtered_dic = {trace.id: trace for trace in filtered_ps_st}
        for trid in filtered_dic.keys():
            if trid not in problem_ids:
                trace = filtered_dic[trid]
                # Read relevant information
                stats = getstats(trace)
                # Add time to others for unpacking list into string
                statlist = list(stats[:])
                channel = statlist[-1]
                statlist.insert(0, self.yyyyjjjhhmmss)
                # Change E-N channels in filenames
                if channel[-1] == '1':
                    statlist[-1] = channel[:-1] + 'N'
                elif channel[-1] == '2':
                    statlist[-1] = channel[:-1] + 'E'
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
        
        :return: 
        """
        return

    def write_pands(self, pname=None, sname=None):
        """
        Function to write files with possible stations that could be included 
        in the .INV file, with predicted P and S arrival times.
        :param pname: 
        :param sname: 
        :return: 
        """
        # Check that MT5 object has all the right variables set
        assert hasattr(self, 'evdt'), 'Make sure date and time are set'

        return

    def write_inv(self, inv_name=None):
        """
        Function to read in sac files and format those from selected stations 
        into 
        :param p_stations: File with list of P-wave stations to use
        :param s_stations: File with list of S-wave stations to use
        :param inv_name: Non-standard (YYMMDD) name for .inv file
        :return: 
        """
        # Turn date-time info into name of file
        evdt = self.evdt
        yyyyjjjhhmmss = datetime.strftime(evdt, '%Y%j%H%M%S')
        # Name of directory with filtered and cut sac files
        mt5_dir = yyyyjjjhhmmss + yyyyjjjhhmmss + '_mt5'
        # Check directory exists and read in files using obspy
        if os.path.exists(mt5_dir):
            stream = core.read(mt5_dir+'*wwlpbn_cut')
        else:
            raise IOError('Subdirectory {} does not exist in this folder'.
                          format(mt5_dir))

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
                          str(str(header_lon_round).split('.')[-2:]))
        header_lat_round = round(self.evla, 2)
        header_lat_str = ('{:2d}'.format(int(header_lat_round)) +
                          str(str(header_lat_round).split('.')[-2:]))

        file_header_str = '{} {}'.format(header_lon_str, header_lat_str)
        return


# To run as a script without running a python interpreter first
# if __name__ == "__main__":
#     import sys
#
#     ARGS = sys.argv[:]
#     if len(ARGS) != 12:
#         raise IOError("Wrong number of arguments: 11 arguments"
#                       " (tar_archive, year, month, day, hour,"
#                       "minute, second, evla, evlo, depth, mag)"
#                       " should be specified!")
#
#     else:
#         TAR_ARCH = ARGS[1]
#         YR, MONTHS, DAYS, HOURS, MINUTES, SECONDS = list(map(int, ARGS[2:8]))
#         EVLON, EVLAT, EVDEPTH, MAGNITUDE = list(map(float, ARGS[8:]))
#         unpack_eq(TAR_ARCH, YR, MONTHS, DAYS, HOURS, MINUTES, SECONDS, EVLON,
#                   EVLAT, EVDEPTH, MAGNITUDE, rot_tol=5.)

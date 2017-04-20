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
from obspy.signal.rotate import rotate2zne
from obspy.taup import TauPyModel

from unpacking_lib import datejul, get_otherid, getstats


def unpack_eq(tar_archive, year, month, day, hour, minute, second, evla,
              evlo, depth, mag, rot_tol=5., before_cut=30., after_cut=90.):
    """
    Script to unpack data from zipped tar archive, altermative to first part of
    README-Brian

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
    evdt = datetime(year, month, day, hour, minute, second)
    yyyyjjjhhmmss = datetime.strftime(evdt, '%Y%j%H%M%S')
    print(("Datetime in yyyyjjjhhmmss format = ", yyyyjjjhhmmss))

    # Change to utc datetime (for cutting)
    evdt_utc = core.UTCDateTime(evdt)

    ###########################################################################
    # Extract tgz archive
    ###########################################################################

    # If yyyyjjjhhmmss already exists as directory, remove it
    if os.path.exists(yyyyjjjhhmmss):
        shutil.rmtree(yyyyjjjhhmmss)
    # Extract tar archive (transparent compression)
    tf = tarfile.open(tar_archive, 'r:*')
    tf.extractall()
    tf.close()

    # Rename extracted file to yyyyjjjhhmmss
    notgz = tar_archive.split('.')[0]
    os.rename(notgz, yyyyjjjhhmmss)

    ###########################################################################
    # Read in data to memory
    ###########################################################################
    st = core.read(yyyyjjjhhmmss + '/*.SAC')
    # Check for multiple seismograms; merge if necessary
    st.merge()
    # Create dictionary of traces with ids as keys:
    id_dic = {tr.id: tr for tr in st}
    # Create list of trace ids for traces that cause unpacking problems
    problem_ids = []

    ###########################################################################
    # Travel times to stations
    ###########################################################################

    # Calculate gcarc distances for station (necessary for phase arrivals)
    # Create empty dictionary for station attributes
    stations_dic = {}
    for tr in st:
        station = tr.stats.station
        if station not in list(stations_dic.keys()):
            # Set sub-dictionary for individual station
            stations_dic[station] = statdic = {}
            # Station location
            stla, stlo = tr.stats.sac['stla'], tr.stats.sac['stlo']
            dist, az, baz = gps2dist_azimuth(stla, stlo, evla, evlo)
            # dist is in metres; change to degrees
            gcarc = kilometer2degrees(dist / 1000.)
            statdic['gcarc'] = gcarc
            statdic['az'] = az
            statdic['baz'] = baz

    # Change trace metadata for all traces
    for tr in st:
        station = tr.stats.station
        statdic = stations_dic[station]
        sacdic = tr.stats.sac
        for key in list(statdic.keys()):
            sacdic[key] = statdic[key]
        sacdic['evdp'] = depth
        sacdic['evlo'] = evlo
        sacdic['evla'] = evla

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
        arrivaldic[station] = model.get_travel_times(depth, gcarc, phase_list)

    # Update SAC data with phase arrival times
    for tr in st:
        station = tr.stats.station
        for phase in arrivaldic[station]:
            # Sac pick convention
            sac_arrivalkey = 't' + str(phasedic[phase.name])
            # Name of phase in SAC file denoted by 'k' then pick (e.g. t0)
            sac_arrivalname = 'k' + sac_arrivalkey
            tr.stats.sac[sac_arrivalkey] = phase.time
            tr.stats.sac[sac_arrivalname] = phase.name

    ###########################################################################
    # Eliminate unpaired E-N and missing-PZ seismograms
    ###########################################################################

    # Create dictionary containing the pairs of E-N components
    # Corresponding components are in the dictionary twice, key:item reversed
    enpairs = {}
    # List of north files: makes rotation easier
    pairednorth = []

    for tr in st:
        channel = tr.stats.channel
        if any(x in channel for x in ['E', '2']):
            eid = tr.id
            # ID for corresponding north file
            nid = get_otherid(tr)
            if nid in list(id_dic.keys()):
                enpairs[eid] = nid
                enpairs[nid] = eid
                pairednorth.append(nid)
            else:
                print('No corresponding N file for', eid)
                problem_ids.append(eid)

    # Find n files with no corresponding e file
    for tr in st:
        channel = tr.stats.channel
        if any(x in channel for x in ['N', '1']):
            if tr.id not in list(enpairs.keys()):
                print('No corresponding E file for', tr.id)
                problem_ids.append(tr.id)

    """
    Check for existence of PZ files for each SAC file.
    Not happy about this; subject to whims of IRIS name changes.
    No work-around that I can think of, as PZ files are not read into Obspy.
    """

    # Empty dictionary for trace ids and corresponding PZ files
    pz_ids = {}
    # List of PZ files
    pz_ls = glob.glob(yyyyjjjhhmmss + '/SACPZ*')
    # Search for trace corresponding to each PZ file
    for pz in pz_ls:
        pz_split = pz.split('.')
        # Separate parts of id
        nw, station, loc, channel = pz_split[1:5]
        if loc[0].isdigit():
            locstr = loc
        else:
            locstr = ''
        pz_trid = '{}.{}.{}.{}'.format(nw, station, locstr, channel)
        if pz_trid in id_dic and pz_trid not in problem_ids:
            pz_ids[pz_trid] = pz
    # Look for traces without PZ files:
    for trid in list(id_dic.keys()):
        if trid not in pz_ids:
            print('No PZ file for', trid)
            problem_ids.append(trid)

    ###########################################################################
    # Rotate E-N pairs to north, where necessary
    ###########################################################################

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
                    znew, nnew, enew = rotate2zne(zdata, 0., -90., ndata, naz,
                                                  0., edata, eaz, 0.)
                    ntr.data[:] = nnew[:]
                    etr.data[:] = enew[:]
                    etr.stats.sac['cmpaz'] = 90.
                    ntr.stats.sac['cmpaz'] = 0.

                except:
                    # TODO sort out non-specific expections
                    print('Problem rotating', nid, eid)
                    print('Removing traces...')
                    problem_ids.append(nid)
                    problem_ids.append(eid)
        else:
            print(eid, nid, "components not orthogonal")
            print(azdiff, "difference in azimuth")
            problem_ids.append(eid)
            problem_ids.append(nid)

    ###########################################################################
    # Apply wwlpbn response #################################################
    ###########################################################################
    # Copy stream to leave original broadband alone
    wwlpbn_st = st.copy()
    wwlpbn_ids = {tr.id: tr for tr in wwlpbn_st}
    # Set response for wwlpbn
    newpaz = {'poles': [-0.257 + 0.3376j, -0.257 - 0.3376j,
                        -0.06283 + 0.0j, -0.06283 + 0.0j],
              'zeros': [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j],
              'gain': 0.5985275}

    for tr in wwlpbn_st:
        # Check I want to actually filter it:
        if tr.id not in problem_ids:
            # Find corresponding PZ:
            pz = pz_ids[tr.id]
            # Read in polezero from file
            try:
                attach_paz(tr, pz)
                oldpaz = tr.stats.paz
                # Decimate by factor 20
                tr.decimate(5)
                tr.decimate(4)
                # Replace existing (pz) filter with wwlpbn
                tr.simulate(paz_remove=oldpaz, paz_simulate=newpaz,
                            simulate_sensitivity=False)

            except:
                print("Problem attaching PZ file for", tr.id)
                problem_ids.append(tr.id)

    ###########################################################################
    # Cut seismograms for MT5 ###############################################
    ###########################################################################
    # Copy filtered data, so that I can cut it
    wwlpbn_cut_st = wwlpbn_st.copy()
    cut_ids = {tr.id: tr for tr in wwlpbn_cut_st}

    for tr in wwlpbn_cut_st:
        # Check I want to use it:
        if tr.id not in problem_ids:
            # If BHZ file, find P arrival
            if 'Z' in tr.stats.channel[-1]:
                cut_pick = tr.stats.sac['t0']
            elif any(x in tr.stats.channel for x in ['E', '2', 'N', '1']):
                cut_pick = tr.stats.sac['t1']

            else:
                print(tr.id, 'does not seem to be any of E, N, Z')
                problem_ids.append(tr.id)
                cut_pick = 0.

            cut_pick_time = evdt_utc + cut_pick
            cutstart = cut_pick_time - before_cut
            cutend = cut_pick_time + after_cut
            tr.trim(cutstart, cutend)
            tr.stats.sac['b'] = tr.stats.starttime - evdt_utc
            tr.stats.sac['e'] = tr.stats.endtime - evdt_utc
            tr.stats.sac['o'] = 0.
            tr.stats.starttime = evdt_utc + tr.stats.sac['b']
    ###########################################################################
    # Write seismograms to MT5 file structure
    ###########################################################################
    # Directory paths
    bbpath = '{}/{}_bb/'.format(yyyyjjjhhmmss, yyyyjjjhhmmss)
    wwlpbnpath = '{}/{}_wwlpbn/'.format(yyyyjjjhhmmss, yyyyjjjhhmmss)
    mt5path = '{}/{}_mt5/'.format(yyyyjjjhhmmss, yyyyjjjhhmmss)

    # Make folders
    for path in [bbpath, wwlpbnpath, mt5path]:
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
            statstats = getstats(bb_tr)
            # Add time to others for unpacking list into string
            statlist = list(statstats[:])
            channel = statlist[-1]
            statlist.insert(0, yyyyjjjhhmmss)
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
            except:
                print('Failed to write SAC file for', trid)

    # Write ${evname}.txt for use with the rest of readme-Brian
    datejul(year, month, day, hour, minute, second, evlo, evla, depth, mag)
# To run as a script without running a python interpreter first
if __name__ == "__main__":
    import sys

    args = sys.argv[:]
    if len(args) != 12:
        raise IOError("Wrong number of arguments: 11 arguments"
                      " (tar_archive, year, month, day, hour,"
                      "minute, second, evla, evlo, depth, mag)"
                      " should be specified!")

    else:
        tar_arch = args[1]
        yr, months, days, hours, minutes, seconds = list(map(int, args[2:8]))
        evlon, evlat, evdepth, magnitude = list(map(float, args[8:]))
        unpack_eq(tar_arch, yr, months, days, hours, minutes, seconds, evlon,
                  evlat, evdepth, magnitude, rot_tol=5.)

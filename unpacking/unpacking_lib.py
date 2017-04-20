"""
Helper function for unpacking script
"""
# Import datetime module
from datetime import datetime


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


def getstats(trace: object) -> object:
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


def datejul(year, month, day, hour, minute, second, evla, evlo, depth, mag):
    """
    Function to write out a file in the format ${evname}.txt, since this is
    still used in the later parts of readme-BRIAN. Will attempt to phase out
    as soon as possible.
    """
    # Make python dt object
    evdt = datetime(year, month, day, hour, minute, second)
    # $evname in mt5 format
    yyyyjjjhhmmss = datetime.strftime(evdt, '%Y%j%H%M%S')
    txt_dt = datetime.strftime(evdt, '%Y (%j)  %m %d %H %M %S.0')
    str_list = list(map(str, [evla, evlo, depth, mag, yyyyjjjhhmmss]))
    txt_line = txt_dt + '   {}    {}   {}  {} {}'.format(*str_list)

    outfilename = '{}/{}.txt'.format(yyyyjjjhhmmss, yyyyjjjhhmmss)
    outfileid = open(outfilename, 'w')
    outfileid.write('          Origin time           Lat      Lon     Dp  Mag')
    outfileid.write('\n')
    outfileid.write(txt_line)
    outfileid.close()

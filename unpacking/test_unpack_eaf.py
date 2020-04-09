from Mt5Eq import Mt5Eq

# Needs to be run from the directory that contains the tar archive
mt5 = Mt5Eq('EAF test')
mt5.set_earthquake_info(2020, 1, 24, 17, 55, 13, 38.390, 39.088, 10, 6.7)
mt5.unpack_eq('2020-01-24-mww67-turkey.tgz')
mt5.get_pands_arrivals()
mt5.write_pands()

# Normally you won't know what these are called until you've run the first part of the script.
# You should edit them too (cut stations to 25 P, 25 S). Realign arrival times in MT5.
mt5.read_pands('2020024175513/2020024175513_Pstations.txt', '2020024175513/2020024175513_Sstations.txt')
mt5.write_inv()

# For help choosing stations so that there is a good azimuthal coverage
mt5.write_station_map()
# Run "sh plot_stations_gmt.sh $evname" to plot on a pdf




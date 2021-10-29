from Mt5Eq import Mt5Eq
import matplotlib
matplotlib.use('Qt5Agg')

# Needs to be run from the directory that contains the tar archive
mt5 = Mt5Eq('EAF test')
mt5.set_earthquake_info(2021, 8, 15, 3, 20, 45, 18.398, -74.090, 8.3, 5.8)
mt5.unpack_eq('2021-08-15-mww58-haiti-region.tar')
mt5.get_pands_arrivals()
mt5.write_pands()

# # Normally you won't know what these are called until you've run the first part of the script.
# # You should edit them too (cut stations to 25 P, 25 S). Realign arrival times in MT5.
mt5.read_pands('2021227032045_Pstations.txt', '2021227032045_Sstations.txt')
mt5.write_inv()
#
# # For help choosing stations so that there is a good azimuthal coverage
# mt5.write_station_map()
# # Run "sh plot_stations_gmt.sh $evname" to plot on a pdf




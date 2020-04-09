from Mt5Eq import Mt5Eq

# Needs to be run from the directory that contains the tar archive
mt5 = Mt5Eq('Turkey test')
mt5.set_earthquake_info(2019, 8, 8, 11, 25, 31, 37.948, 29.697, 10, 5.8)
mt5.unpack_eq('2019-08-08-mww58-turkey.tgz')
mt5.get_pands_arrivals()
mt5.write_pands()

# Normally you won't know what these are called until you've run the first part of the script.
# You should edit them too (cut stations/refine arrival times from bb), otherwise your inversion will be rubbish.
# mt5.read_pands('2019079063427/2019079063427_Pstations.txt', '2019079063427/2019079063427_Sstations.txt')
# mt5.write_inv()

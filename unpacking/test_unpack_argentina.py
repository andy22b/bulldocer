from Mt5Eq import Mt5Eq

mt5 = Mt5Eq('Argentina test')
mt5.set_earthquake_info(2011, 10, 6, 11, 12, 30, -24.1315, -64.2963, 17, 5.9)
mt5.unpack_eq('2011-10-06-mw59-salta-province-argentina.tar')
mt5.get_pands_arrivals()
mt5.write_pands()
mt5.read_pands('2011279111230/2011279111230_Pstations.txt', '2011279111230/2011279111230_Sstations.txt')
mt5.write_inv()

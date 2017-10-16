from Mt5Eq import Mt5Eq

mt5 = Mt5Eq('2013 Oct 12 Mw6.5 Crete')
mt5.set_earthquake_info(2013, 10, 12, 13, 11, 54, 35.5277, 23.3718, 48, 6.8)
mt5.unpack_eq('2013-10-12-mw67-crete-2.tgz')
mt5.get_pands_arrivals()
mt5.write_pands(pname='test_arrivals/testp.txt',
                sname='test_arrivals/tests.txt')
mt5.read_pands('test_arrivals/testp.txt', 'test_arrivals/tests.txt')
mt5.write_inv()

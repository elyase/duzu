from duzu import duzu10 as dz

results = {  		# Values taken from http://amdc.in2p3.fr/web/dz.html
    (37, 53): -78.9805298,
    (82, 126): -22.5833740,
    (126, 185): 257.410645,
    (6,6): 0.631393433
}


def test_result():
    for nucleus in results:
        Z, N = nucleus
        assert abs(dz(Z, N) - results[nucleus]) < 0.01

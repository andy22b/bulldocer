"""
Library containing helper functions for calculating synthetics
"""
# Import modules
from numpy import arcsin, cos, tan, array, amin, zeros


def psvcoef(rp, alph1, beta1, rho1, alph2, beta2, rho2):
    """
    :param rp: Ray parameter
    :param alph1: P velocity in uppermost layer
    :param beta1: S velocity in uppermost layer
    :param rho1:
    :param alph2:
    :param beta2:
    :param rho2:
    :return:
    """

    p = rp
    p2 = rp * rp

    ai1 = arcsin(p * alph1)
    aj1 = arcsin(p * beta1)
    ai2 = arcsin(p * alph2)
    aj2 = arcsin(p * beta2)

    aincs = array([ai1, aj1, ai2, aj2])

    ci1 = cos(ai1) / alph1
    ci2 = cos(ai2) / alph2
    cj1 = cos(aj1) / beta1
    cj2 = cos(aj2) / beta2
    a = rho2 * (1 - 2 * beta2 ^ 2 * p2) - rho1 * (1 - 2 * beta1 ^ 2 * p2)
    b = rho2 * (1 - 2 * beta2 ^ 2 * p2) + 2 * rho1 * beta1 ^ 2 * p2
    c = rho1 * (1 - 2 * beta1 ^ 2 * p2) + 2 * rho2 * beta2 ^ 2 * p2
    d = 2 * (rho2 * beta2 ^ 2 - rho1 * beta1 ^ 2)
    e = b * ci1 + c * ci2
    f = b * cj1 + c * cj2
    g = a - d * ci1 * cj2
    h = a - d * ci2 * cj1
    dd = e * f + g * h * p2

    a1 = (a * c + b * d * ci1 * cj1)
    a2 = (b * ci1 - c * ci2) * f
    a3 = (a + d * ci2 * cj1) * p2 * g
    a4 = (a + d * ci1 * cj2) * h * p2
    a5 = (a * b + c * d * ci2 * cj2)
    a6 = (b * cj1 - c * cj2) * e

    pdpu = (a2 - a4) / dd
    pdsu = -2 * ci1 * a5 * p * alph1 / (beta1 * dd)
    pdpd = 2 * rho1 * ci1 * f * alph1 / (alph2 * dd)
    pdsd = 2 * rho1 * ci1 * h * p * alph1 / (beta2 * dd)
    sdpu = -2 * cj1 * a5 * p * beta1 / (alph1 * dd)
    sdsu = -(a6 - a3) / dd
    sdpd = -2 * rho1 * cj1 * g * p * beta1 / (alph2 * dd)
    sdsd = 2 * rho1 * cj1 * e * beta1 / (beta2 * dd)
    pupu = 2 * rho2 * ci2 * f * alph2 / (alph1 * dd)
    pusu = -2 * rho2 * ci2 * g * p * alph2 / (beta1 * dd)
    pupd = -(a2 + a3) / dd
    pusd = 2 * ci2 * a1 * p * alph2 / (beta2 * dd)
    supu = 2 * rho2 * cj2 * h * p * beta2 / (alph1 * dd)
    susu = 2 * rho2 * cj2 * e * beta2 / (beta1 * dd)
    supd = 2 * cj2 * a1 * p * beta2 / (alph2 * dd)
    susd = (a6 + a4) / dd

    rs0 = 4 * p2 * ci1 * cj1
    rs1 = (1 / beta1 ^ 2 - 2 * p2) ^ 2
    r0 = (rs0 - rs1) / (rs0 + rs1)

    refcs = array([pdpu, sdpu, pupu, supu],
                  [pdsu, sdsu, pusu, susu],
                  [pdpd, sdpd, pupd, supd],
                  [pdsd, sdsd, pusd, susd])

    return refcs, aincs, r0


def shcoef(rp, beta1, rho1, beta2, rho2):
    aj1 = arcsin(rp * beta1)
    aj2 = arcsin(rp * beta2)

    aincs = array([aj1, aj2])

    dd = rho1 * beta1 * cos(aj1) + rho2 * beta2 * cos(aj2)
    sdsu = (rho1 * beta1 * cos(aj1) - rho2 * beta2 * cos(aj2)) / dd
    susu = 2 * rho2 * beta2 * cos(aj2) / dd
    sdsd = 2 * rho1 * beta1 * cos(aj1) / dd
    susd = -sdsu

    refcs = array([sdsu, susu], [sdsd, susd])

    return refcs, aincs


def water(rp, alphw, wthk, pdpu, r0, timax, nraymax):
    if wthk == 0:
        raise ValueError('No water layer')

    cutoff = 1.e-5
    ai1 = arcsin(rp * alphw)

    wdly = 2 * wthk / (alphw * cos(ai1))  # two-way delay in water
    twc = 2 * wthk * tan(ai1) * rp  # correction to wavefront
    wdly -= twc

    nrev = int(timax / wdly)
    nrev = amin((nrev, nraymax))
    pw = 1
    n = 0

    t = zeros((nrev,))
    amp = zeros((nrev,))
    while abs(pw) > cutoff and n < nrev:
        pw = r0 ^ n * pdpu ^ (n - 1)
        t[n] = n * wdly
        amp[n] = pw
        n += 1

    nrev = n
    amp = amp[:nrev]
    t = t[:nrev]
    return amp, t, nrev


def m5cresp(sthk, p, vp, vs, rho, vpu, vsu, rhou, vpd, vsd, rhod, nrev,
            nraymax):
    """
    -- calculate complete P-SV response of a single layer for NREV bounces
    """
    # Matrix of coefficients
    a = zeros((2, 4))
    ampout = zeros((2, 4))
    timout = zeros((nraymax,))
    # reflection coefficients at upper surface
    refcs, aincs, r0 = psvcoef(p, vpu, vsu, rhou, vp, vs, rho)
    a[0, :] = array([refcs[2, 2], refcs[2, 3], refcs[3, 2], refcs[3, 3]])

    # reflection coefficients at lower surface
    refcs, aincs, r0 = psvcoef(p, vp, vs, rho, vpd, vsd, rhod)
    a[1, :] = array([refcs[0, 0], refcs[0, 1], refcs[1, 0], refcs[1, 1]])

    ai1 = aincs[0]
    aj1 = aincs[1]
    xp = sthk * tan(ai1)
    xs = sthk * tan(aj1)

    tp = sthk / (vp * cos(ai1)) - xp * p
    ts = sthk / (vs * cos(aj1)) - xs * p

    # incoming transmission coefficients at lower surface
    ampin0 = array([[refcs[2, 0], refcs[2, 1]],
                    [refcs[3, 0], refcs[3, 1]]])

    # outgoing transmission coefficients at lower surface
    ampout0 = array([refcs[0, 2], refcs[1, 2]])
    nrays = 0
    for irev in range(2, nrev):

        nrays = amin((nrays + 1, nraymax))
        nlegs = 2 * irev  # number of legs, times, and amplitudes
        nray = 2 ^ nlegs

        if (nrays + nlegs) > nraymax:
            raise ValueError('Calculating for too many reverberations')
        for k in range(nlegs):
            timout[k + nrays] = k * ts + (nlegs - k) * tp

        for jray in range(nray - 1):
            bin0 = binconv(jray, nlegs)  # get binary value for ray

            # KSUM is the number of S phases, Time = ksum*Ts+(nlegs-ksum)*Tp
            amp = 1
            mm = True
            ksum = 0
            for i in range(nlegs - 1):
                ksum += int(bin0[:i])

    # look at each bounce, KRAY1 is the phase going into the bounce, KRAY2 is
    #       the phase coming out

            for i in range(nlegs - 1):
                kray1 = int(bin0[i])

                if i == 1:
                    kr1 = kray1 + 1  # 'phase of first leg in layer

                kray2 = int(bin0[i + 1])

                if i == nlegs - 1:
                    kr9 = kray2 + 1  # 'phase of last leg in layer

                kray = 2 * kray1 + kray2  # KRAY=0(PP) =1(PS) =2(SP) =3(SS)
                amp = amp * a[kray + 1, mm + 2]  # find coefficient
                mm = not mm  # MM=true for top interface, false for bottom
            k2ray = ksum + nrays

            ampout[k2ray, 0] += amp * ampin0[1, kr1] * ampout0[kr9]  # 'pP
            ampout[k2ray, 1] += amp * ampin0[2, kr1] * ampout0[kr9]  # 'sP

        nrays += nlegs

    return timout, ampout, nrays


def binconv(nint, nlen):
    # convert integer NINT into binary string BIN$ of length NLEN
    binoct = "000001010011100101110111"
    bin0 = ""
    o = oct(nint)
    l = len(o)

    for k in range(l):
        bin0 += bin0 + binoct[3 * int(o[:k]):3 * int(o[:k]) + 3]

    lbin = len(bin0) - nlen
    # Pad with zeros
    if lbin < 0:
        for k in range(abs(lbin)):
            bin0 = "0" + bin0
    elif lbin > 0:
        bin0 = bin0[lbin:]

    return bin0

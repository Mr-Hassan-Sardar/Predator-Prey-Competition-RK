# sharkCompetition.py
# Program to run simulation of species competition between
# white tip sharks (WTS) and black tip sharks (BTS)

def sharkCompetition(DT=0.001, simLength=5, method='RK4'):
    numIterations = int(simLength/DT) + 1
    t = 0

    WTS_population = 20
    BTS_population = 15

    tLst = [t]
    WTSLst = [WTS_population]
    BTSLst = [BTS_population]

    if method == 'Euler':
        DT2 = DT
        DT4 = DT
    elif method == 'RK2':
        DT2 = DT / 2
        DT4 = DT
    elif method == 'RK4':
        DT2 = DT / 2
        DT4 = DT / 4
    else:
        raise ValueError("Invalid method. Choose from 'Euler', 'RK2', or 'RK4'.")

    for i in range(1, numIterations):
        t = i * DT
        WTS_population, BTS_population = rk_step(WTS_population, BTS_population, DT2, DT4)

        tLst.append(t)
        WTSLst.append(WTS_population)
        BTSLst.append(BTS_population)

    return tLst, WTSLst, BTSLst

def rk_step(WTS_population, BTS_population, DT2, DT4):
    k1_wts, k1_bts = rate_equations(WTS_population, BTS_population)
    k2_wts, k2_bts = rate_equations(WTS_population + DT2 * k1_wts, BTS_population + DT2 * k1_bts)
    k3_wts, k3_bts = rate_equations(WTS_population + DT2 * k2_wts, BTS_population + DT2 * k2_bts)
    k4_wts, k4_bts = rate_equations(WTS_population + DT4 * k3_wts, BTS_population + DT4 * k3_bts)

    WTS_population += DT4 * (k1_wts + 2 * k2_wts + 2 * k3_wts + k4_wts) / 6
    BTS_population += DT4 * (k1_bts + 2 * k2_bts + 2 * k3_bts + k4_bts) / 6

    return WTS_population, BTS_population

def rate_equations(WTS_population, BTS_population):
    WTS_birth_fraction = 1
    WTS_death_proportionality_constant = 0.27
    WTS_births = WTS_population * WTS_birth_fraction
    WTS_deaths = (WTS_death_proportionality_constant * BTS_population) * WTS_population

    BTS_birth_fraction = 1
    BTS_death_proportionality_constant = 0.2
    BTS_births = BTS_birth_fraction * BTS_population
    BTS_deaths = (BTS_death_proportionality_constant * WTS_population)*BTS_population

    return (WTS_births - WTS_deaths), (BTS_births - BTS_deaths)

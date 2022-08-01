import sys
import pandas
import matplotlib.pyplot as plt
import pathlib
import math
import os
import base64
import time
import pytz
import datetime
import tzlocal
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

#This code was written by Miranda Lea Nelson contact mirandanelson218@u.boisestate.edu
#This code will run through the yield point method for all tests, graphing
#stress versus strain, slope of stress-strain curve, and the derivative of
#the slope of the stress-strain curve


def run():
    #Change file to equal the name of the data file you want to analyze
    #Add drive location if file is not in the same folder as python
    #file = 'filename and location'
    file = "AK02.xlsx"

    #To create a csv file of the data outputs change create_csv to True
    create_csv = False
    #The name and location of exported .csv
    export_filename = pathlib.Path("Summary_Report" + '.csv')
    # settings for finding Transition point(tp), Yield point(rp), and Rupture point(rp)
    # Transition point: Calculated as the point below the yield point that varies tp_percent from the linear modulus (maximum slope)
    tp_percent = 3
    # Yield Point: Calculated using yp_percent deviation from the linear modulus (maximum slope).
    yp_percent = 0
    # Rupture point: Calculated as rp_percent reduction in stress after the ultimate stress
    rp_percent = 70
    # For rupture point to instead be calculated as the last data point, Change rp_last to True
    rp_last = False
    #To turn off transition point calculation turn transition_on to False
    transition_on = True

    tp_error = False
    yp_error = False
    rp_error = False

    #Brings in Data
    if file.find('.xlsx') > -1:
        dfl = pandas.read_excel(file, header=None)
        stre = dfl.iloc[:, 1]
        stra = dfl.iloc[:, 0]
        titlename = file.rstrip(".xlsx")

    if file.find('.csv') > -1:
        dfd = pandas.read_csv(file, header=None)
        stre = dfd.iloc[:, 1]
        stra = dfd.iloc[:, 0]
        titlename = file.rstrip(".csv")

    if file.find('.txt') > -1:
        dfd = pandas.read_csv(file, sep='\t', header=None)
        stre = dfd.iloc[:, 1]
        stra = dfd.iloc[:, 0]
        titlename = file.rstrip(".txt")

    if '/' in titlename:
        other, titlename = titlename.rsplit('/', 1)

    #Verifys that data has been formatted correctly
    d = 0
    di = 0
    if type(stre[0]) == str or math.isnan(stre[0]):
        stre.pop(0)
        stra.pop(0)
        di = 1
        if type(stre[1]) == str or math.isnan(stre[1]):
            stre.pop(1)
            stra.pop(1)
            di = 2

    while d == 0:
        if type(stre[di]) == str or math.isnan(stre[di]):
            d = 1
        if type(stra[di]) == str or math.isnan(stra[di]):
            d = 1
        if di == len(stre) - 1:
            d = 2
        di = di + 1

    if d == 1:
        print('Incorrect data input: Column 1 must be strain and column 2 must be stress. Column 1 and 2 must be equal length')

    #Changes data type
    stre = pandas.to_numeric(stre)
    stra = pandas.to_numeric(stra)

    #Uts Calculation
    uts_val = max(stre)
    uts_pt = np.where(stre == uts_val)[0][0]
    index = uts_pt
    uts_stress = uts_val
    uts_strain = stra[uts_pt]

    #Rupture calculation
    if index == len(stre) - 1:
        print('The rupture point settings are out of range. Points that are out of range are not plotted.')
        rp_error = True
        rupture_pt = None
        rupture_stress = None
        rupture_strain = None
    else:
        if rp_last == False:
            rpp = int(rp_percent) / 100
            rupture_val = uts_val * rpp
            if stre[(len(stre) - 1)] > rupture_val:
                #the rupture point is out of range and will be plotted as the last values of data instead
                rupture_pt = (len(stre) - 1)
                rupture_stress = stre[rupture_pt]
                rupture_strain = stra[rupture_pt]
            else:
                rupture_pts = next(x for x, val in enumerate(stre[index:(len(stre) - 1)])
                                   if val < rupture_val)
                rupture_pt = index + rupture_pts
                rupture_stress = stre[rupture_pt]
                rupture_strain = stra[rupture_pt]
        elif rp_last == True:
            rupture_pt = (len(stre) - 1)
            rupture_stress = stre[rupture_pt]
            rupture_strain = stra[rupture_pt]

    #cuts data at UTS
    #Values of the upward trend part of the stress-strain curve; pre-failure data points
    StressFitUTS = stre[0:index]
    StrainFitUTS = stra[0:index]

    #7th degree polyfit to stress strain data to help smooth fluctuations seen
    #in the raw data stress-strain curve due to noise
    #xStrain = np.linspace(start=StrainFitUTS[0], stop=StrainFitUTS[line2], num=100)
    xStrain = np.linspace(start=StrainFitUTS[0], stop=StrainFitUTS[index - 1], num=100)
    p_d7 = np.polyfit(StrainFitUTS, StressFitUTS, 7)
    p_d7plot = np.polyval(p_d7, xStrain)
    #Change in stress (slope, or y-change)
    dStress = np.diff(p_d7plot)
    #Change in strain (slope, or x-change)
    dStrain = np.diff(xStrain)

    #This is an array of lin mod values (slopes of the stress-strain curve)
    slope = dStress / dStrain

    locs = find_peaks(slope)
    locs = locs[0]
    pks = [slope[locs]]
    yval = np.amax(pks)
    YptFt = np.where(pks[0] == yval)
    p = locs[YptFt[0]]
    p = int(p[0])
    stressAccel = np.diff(slope)
    #Max derivative/slope of the stress-strain curve
    linMod = slope[p]

    # n is 1% less than p
    n0 = next(x for x in xStrain if x >= (xStrain[p] - 0.01))
    n = np.where(xStrain == n0)[0][0]
    if (p - n) < 3:
        n = p - 3

    [fit_lin_p1, fit_lin_p2] = np.polyfit(xStrain[n:p], p_d7plot[n:p], 1)
    linModFit = fit_lin_p1
    yintfit = fit_lin_p2

    #Transition point calculation
    fitval = []
    fitvalx = []
    df = []
    dfx = []
    itert = []
    for t in range(0, n - 2):
        fitval.insert(t, linModFit * xStrain[t] + yintfit)
        fitvalx.insert(t, (p_d7plot[t] - yintfit) / (linModFit))
        df.insert(t, abs(p_d7plot[t] - fitval[t]))
        dfx.insert(t, abs(xStrain[t] - fitvalx[t]))
        itert.insert(t, t)

    tpp = int(tp_percent) / 100
    if tpp == 0:
        tp_percent = 1
        tpp = .01
        print('Transition point settings must be more than 0,settings will be set at 1%')
    if transition_on == False:
        tp_error = True
        transition_stress = None
        transition_strain = None
        energy2Tr = None
    else:
        transition_pt = np.where(df <= abs(tpp * p_d7plot[p]))[0][0]
        if transition_pt > p:
            print('transition point is too high')
            quit()
        if transition_pt == None:
            transition_pt = 1

        transition_stress = p_d7plot[transition_pt]
        transition_strain = xStrain[transition_pt]
        # transition energy calculation
        TiRealval = stra[stra <= xStrain[transition_pt]]
        TiReal = len(TiRealval) - 1
        energy2Tr = np.trapz(stre[0:TiReal], x=stra[0:TiReal])
        if transition_strain <= stra[0]:
            tp_error = True
            print('The transition point settings are out of range. Points that are out of range are not plotted.')
            transition_stress = None
            transition_strain = None
            energy2Tr = None




    fitvaly = []
    fitvalxy = []
    dfy = []
    dfxy = []
    iterty = []
    count = 0
    for d in range(p - 1, 99):
        fitvaly.insert(count, linModFit * xStrain[d] + yintfit)
        fitvalxy.insert(count, (p_d7plot[d] - yintfit) / (linModFit))
        dfy.insert(count, abs(p_d7plot[d] - fitvaly[count]))
        dfxy.insert(count, abs(xStrain[d] - fitvalxy[count]))
        iterty.insert(count, d)
        count = count + 1

    ypp = int(yp_percent) / 100
    yield_check = np.nonzero(dfy >= abs(ypp * p_d7plot[p]))
    if yield_check[0].size == 0:
        yp_error = True
        print('The yield point settings are out of range. Points that are out of range are not plotted')
    else:
        yield_pt = np.nonzero(dfy >= abs(ypp * p_d7plot[p]))[0][0]
        yield_pt = yield_pt + p

    if yp_error == False:
        yield_stress = p_d7plot[yield_pt]
        yield_strain = xStrain[yield_pt]
        YiRealval = stra[stra <= xStrain[yield_pt]]
        YiReal = len(YiRealval) - 1
        energy2Y = np.trapz(stre[0:YiReal], x=stra[0:YiReal])
    else:
        yield_stress = None
        yield_strain = None
        energy2Y = None


    energy2UTS = np.trapz(stre[0:index], x=stra[0:index])
    energy2Rupture = np.trapz(stre[0:rupture_pt], x=stra[0:rupture_pt])

    # Get date and time
    now_utc = datetime.datetime.utcnow()
    local_tz = tzlocal.get_localzone()
    tz = pytz.timezone(str(local_tz))
    tsr = str(pytz.utc.localize(now_utc).astimezone(tz))

    # Format outputs
    if create_csv == True:
        csv_data = pandas.DataFrame(columns=['Filename', 'Timestamp', 'Transition stress', 'Transiton strain',
                                             'Yield stress', 'Yield strain', 'Ultimate stress', 'Ultimate strain',
                                             'Rupture stress', 'Rupture strain', 'Linear modulus',
                                             'Transition energy', 'Yield energy',
                                             'Energy to Ultimate', 'Energy to Rupture',
                                             'Setting - Transition (%)', 'Setting - Yield (%)', 'Setting - Rupture (%)'])
        csv_data = csv_data.append({'Filename': titlename, 'Timestamp': tsr, 'Transition stress': transition_stress,
                                    'Transiton strain': transition_strain, 'Yield stress': yield_stress,
                                    'Yield strain': yield_strain, 'Ultimate stress': round(uts_stress, 2),
                                    'Ultimate strain': round(uts_strain, 3), 'Rupture stress': rupture_stress,
                                    'Rupture strain': rupture_strain, 'Linear modulus': round(linMod, 4),
                                    'Transition energy': energy2Tr,
                                    'Yield energy': energy2Y,
                                    'Energy to Ultimate': round(energy2UTS, 4),
                                    'Energy to Rupture': round(energy2Rupture, 4),
                                    'Setting - Transition (%)': tp_percent,
                                    'Setting - Yield (%)': yp_percent, 'Setting - Rupture (%)': rp_percent},
                                   ignore_index=True)

        exp_data = pandas.DataFrame(csv_data)
        if export_filename.exists():

            with open(str(export_filename), 'a') as expfile:
                text = exp_data.to_csv(index=False, header=False)
                expfile.write(text)
            expfile.close()
        else:
            expfile = open(export_filename, 'w', newline='')
            text = exp_data.to_csv(index=False)
            expfile.write(text)
            expfile.close()


    footer_text = f'Linear Modulus   {round(linMod, 3)} '
    data = [
        ['Strain', 'Stress', 'Energy(Toughness)', 'Setting - %'],
        ['Transition', transition_strain, transition_stress, energy2Tr, tp_percent],
        ['Yield', yield_strain, yield_stress, energy2Y, yp_percent],
        ['Ultimate', uts_strain, uts_stress, energy2UTS, None],
        ['Rupture', rupture_strain, rupture_stress, energy2Rupture, rp_percent],
    ]
    column_headers = data.pop(0)
    row_headers = [x.pop(0) for x in data]
    cell_text = []

    for row in data:
        si = 0
        if not any(s == None for s in row):
            cell_text.append(['%4.3f' % x for x in row])
        else:
            for s in row:
                if not (s == None):
                    row[si] = '%4.3f' % row[si]
                # cell_text.append([row[si]])
                si = si + 1
            cell_text.append([x for x in row])

    plt.title(titlename)
    plt.plot(stra, stre, color='black', linewidth=2)
    plt.ylabel('Stress')
    plt.xlabel('Strain')
    if tp_error == False:
        plt.plot(xStrain[transition_pt], p_d7plot[transition_pt], '.c', markersize=20, label="Transition")
    if yp_error == False:
        plt.plot(xStrain[yield_pt], p_d7plot[yield_pt], '.y', markersize=20, label="Yield")
    plt.plot(stra[index], stre[index], '.r', markersize=20, label="Ultimate")
    if rp_error == False:
        plt.plot(stra[rupture_pt], stre[rupture_pt], '.m', markersize=20, label="Rupture")
    plt.plot(xStrain, xStrain * linModFit + yintfit, '--', label="Linear Modulus")
    plt.legend()
    plt.figure(tight_layout={'pad': 1})

    plt.subplot(121)
    plt.tight_layout()
    plt.title(titlename)
    plt.xlabel('Strain')
    plt.ylabel('d(stress) / d(strain)')
    plt.plot(xStrain[0:-1], slope, 'k', linewidth=1.5)
    plt.plot(xStrain[p + 1], slope[p], '.y', markersize=15)
    plt.plot(xStrain[-1], slope[-1], '.r', markersize=15)

    plt.subplot(122)
    plt.tight_layout()
    plt.xlabel('Strain')
    plt.ylabel('$\mathregular{d^{2}}$(stress) / $\mathregular{d(strain)^{2}}$')
    plt.plot(xStrain[1:-1], stressAccel, 'k')
    plt.plot(xStrain[p + 2], stressAccel[p], '.y', markersize=15)
    plt.plot(xStrain[-1], stressAccel[-1], '.r', markersize=15)
    plt.figure(tight_layout={'pad': 1})

    the_table = plt.table(cellText=cell_text,
                          rowLabels=row_headers,
                          rowLoc='right',
                          colLabels=column_headers,
                          # colWidths= 5,
                          loc='center')
    the_table.scale(1, 1.5)

    # fig = plt.gcf()
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.box(on=None)
    plt.figtext(0.05, 0.25, footer_text, horizontalalignment='left')
    plt.draw()
    fig = plt.gcf()

    plt.show()





if __name__ == '__main__':
    run()
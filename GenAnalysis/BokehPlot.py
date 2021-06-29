### Script combining all data for MulitFunPlot and producing the plot(s)

### example to run: python3 BokehPlot.py
### e.g.: python3 BokehPlot.py

import sys
import os.path
import requests
import pandas as pd
import statistics
from io import StringIO
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import gridplot, layout, column
from bokeh.models import ColumnDataSource, CDSView, Label
from bokeh.models.filters import GroupFilter
from bokeh.models.widgets import Tabs, Panel
from bokeh.palettes import Category10_10, Set1_3, Blues, Oranges, Greens, Reds, Purples, YlOrBr, PuRd

def main():

    ### load file containing descriptions to BC group IDs
    key = pd.read_csv('output/MultiFunKey.csv', delimiter = ':')
    ### convert the key table into a dictionary
    all_keys = dict([(i, a) for i, a in zip(key.ID, key.function)])
    all_keys['ORFs'] = 'ORFs'
    all_keys['Conserved-Hypothetical-ORFs'] = 'CH-ORFs'
    all_keys['Unclassified-Genes'] = 'Unclassified'
    all_keys['empty'] = 'None'

    ### define output file the the plot
    output_file("output/BokehPlot.html")
    ### load data files to check for mismatches in gene names after using EcoCyc
    ecocyc_in = pd.read_csv('output/GeneNames.csv', delimiter = ',')
    ### load data table with PSS and API values
    plot_data = pd.read_csv('output/DataIGRs&GenesCommon.csv', delimiter = ',')
    ### add gene names to the table with Ide and Seg
    data = pd.concat([plot_data, ecocyc_in], axis = 1)

    ### load data table with gene assignment to functional groups
    MFinput = pd.read_csv('output/MultiFunInput.csv', delimiter = ',')

    ### here I need to create a new table with a separate row for each gene in every BC group
    ### initialize list to store all values for the 'final_table'
    vec = []
    ### loop thorugh all row numbers of table with functional group & genes
    for r in range(0, len(MFinput['ID'])):
        ### split the BC group to extract parent BC groups
        ids = MFinput['ID'][r].split('.')
        ### split genes present in that group, so that each gene is a separable variable
        genes = MFinput['Gene'][r].split(' // ')
        ### loop through all those genes
        for g in genes:
            ### remove possibly existing characters: " from the gene name
            if '"' in g:
                g = g.split('"')[1]
            ### loop through all row numbers of the table with PSS & API values together with gene names
            for prom in range(0, len(data['Gene'])):
                ### when you find the gene name currently processing and it's the first encounter
                if data['Gene'][prom] == g:
                    ### initialize a list to store all values for this gene in
                    line = []
                    ### add the BC group ID and gene names to the 'line' list
                    line.append(MFinput['ID'][r])
                    line.append(g)
                    ### initialize a list to store BC subgroups
                    bc = []
                    ### make a loop that executes 5 times (= number of values to be added to the 'line' list as BC subgroups)
                    for id in range(0, 5):
                        ### if going through the loop for the first time
                        if (id == 0):
                            ### add the upper most BC group to the 'line' list
                            line.append(ids[id])
                            ### get the description of the upper most BC group from the key dictionary
                            for ak in all_keys.keys():
                                if ids[id] == ak:
                                    line.append(all_keys[ak])
                        ### in all other rounds
                        else:
                            ### if the round number is lower than number of subgroups
                            if (id < len(ids)):
                                ### produce the current subgroup combining subgroup from the previous 'line' value and current extention to it
                                bc = '.'.join([line[2 * id], ids[id]])
                                ### get the description of the current subgroup from the key dictionary
                                for ak in all_keys.keys():
                                    if bc == ak:
                                        bckey = all_keys[ak]
                            ### if you run out of subgroup values
                            else:
                                ### set bc variable to NaN
                                bc = 'NaN'
                                bckey = 'None'
                            ### add the BC subgroup value and its description to the 'line' list
                            line.append(bc)
                            line.append(bckey)
                    ### add its API and PSS values to the 'line' list
                    line.append(data['PromAPI'][prom])
                    line.append(data['PromPSS'][prom])
                    line.append(data['GeneAPI'][prom])
                    line.append(data['GenePSS'][prom])
                    ### add which promoter were used for further experiment
                    if any([g == 'aceB', g == 'aldA', g == 'lacZ', g == 'mtr', g == 'yhjX', g == 'cdd', g == 'dctA', g == 'ptsG', g == 'purA', g == 'tpiA']):
                        exp = 'yes'
                        line.append(exp)
                    else:
                        exp = 'no'
                        line.append(exp)
            ### add the whole 'line' variable as a list into the 'vec' list
            vec.append(line)

    ### create dataframe out the the complete 'vec' list with specified column names
    final_table = pd.DataFrame(vec, columns = ['ID', 'Gene', 'BC0', 'BC0-fun', 'BC1', 'BC1-fun', 'BC2', 'BC2-fun', 'BC3', 'BC3-fun', 'BC4', 'BC4-fun', 'PromAPI', 'PromPSS', 'GeneAPI', 'GenePSS', 'Experiment'])
    ### save this dataframe as csv table
    final_table.to_csv('output/MultiFunPlotData.csv')
    final_table = pd.read_csv('output/MultiFunPlotData.csv', delimiter = ',')
    del final_table['Unnamed: 0']

    final_table = final_table.fillna(str('empty'))
    #########################################
    ######## GENERAL PLOT DEFINITION ########
    #########################################
    ### define tools you want to be present in the plots
    tls = "pan, tap, hover, box_zoom, box_select, wheel_zoom, reset, save, crosshair"
    ### define source to enable linked selection of genes between the plots
    src = ColumnDataSource(data = dict(x0 = final_table['PromAPI'], y0 = final_table['GeneAPI'], x1 = final_table['PromPSS'], y1 = final_table['GenePSS'], z0 = final_table['Gene'], z1 = final_table['BC0'], z2 = final_table['BC0-fun'], z3 = final_table['BC1-fun'], z4 = final_table['BC2-fun'], z5 = final_table['BC3-fun'], z6 = final_table['BC4-fun'], z7 = final_table['Experiment']))
    ### add information boxes to an inspection tool (hover)
    tltips = [('Gene', '@z0'),('Main group', '@z2'), ('Subgroups', '@z3'), ('', '@z4'), ('', '@z5'), ('', '@z6'), ('Experiment', '@z7')]

    ### extract BC groups' IDs from 'BC0' column
    groups = []
    for zero in final_table['BC0']:
        if zero not in groups:
            groups.append(zero)
    ### sort the groups alphabetically
    groups = sorted(groups)

    subs = {}
    for g in groups:
        subs[g] = final_table[final_table['BC0'] == g]

    ### add function descriptions for the ID extracted from 'BC0' in a dictionary
    main_plot_key = {}
    for grp in groups:
        for k in all_keys.keys():
            if grp == k:
                main_plot_key[all_keys[k]] = grp

    ### define views to distinguish BC groups during plotting
    view = {}
    for pk in sorted(main_plot_key.values()):
        view[pk] = CDSView(source = src, filters = [GroupFilter(column_name = 'z1', group = pk)])

    ### information about inclusion in experiment
    expinfo = ['yes', 'no']

    ### add function descriptions inclusion in experiment
    sec_plot_key = {"Selected":'yes', "Not selected":'no'}

    ### define views to distinguish experimental inclusion during plotting
    view2 = {}
    for sec in sorted(sec_plot_key.values()):
        view2[sec] = CDSView(source = src, filters = [GroupFilter(column_name = 'z7', group = sec)])

    #########################################
    ########## 1st PLOT DEFINITION ##########
    #########################################
    ### check for maximal values to be used in both x and y axis of the plot
    mprom = round(max(final_table['PromPSS']), 1)
    mgene = round(max(final_table['GenePSS']), 1)
    if (mprom > mgene):
        mplot = mprom + 0.01
    else:
        mplot = mgene + 0.01
    ### set the figure dimension and other characteristics
    pl1 = figure(tools = tls, tooltips = tltips, plot_width = 750, plot_height = 750, x_range = (-0.012, mplot), y_range = (-0.012, mplot), x_axis_label = 'Intergenic regions', y_axis_label = 'Open reading frames')
    ### set plot main title
    pl1.title.text = 'Proportion of segregating sites'
    ### set the Seg values to be in the plotted
    n = 0
    for v in list(main_plot_key.keys()):
        pl1.scatter('x1', 'y1', source = src, view = view[main_plot_key[v]], legend_label = v, color = Category10_10[n])
        n += 1
    ### set legend position and the attribute to hide/show it upon clicking in the legend
    pl1.legend.location = 'top_left'
    pl1.legend.background_fill_alpha = 0
    pl1.legend.click_policy = 'hide'

    #########################################
    ########## 2nd PLOT DEFINITION ##########
    #########################################
    ### set the figure dimension and other characteristics
    pl2 = figure(tools = tls, tooltips = tltips, plot_width = 375, plot_height = 375, x_range = (-0.012, mplot), y_range = (-0.012, mplot), x_axis_label = 'Intergenic regions', y_axis_label = 'Open reading frames')
    ### set plot main title
    pl2.title.text = 'Proportion of segregating sites'
    ### set the Seg values to be in the plotted
    n = 0
    for v in list(sec_plot_key.keys()):
        pl2.scatter('x1', 'y1', source = src, view = view2[sec_plot_key[v]], legend_label = v, color = Set1_3[n], muted_color = Set1_3[n], muted_alpha = 0.05)
        n += 1
    ### set legend position and the attribute to hide/show it upon clicking in the legend
    pl2.legend.location = 'top_left'
    pl2.legend.background_fill_alpha = 0
    pl2.legend.click_policy = 'mute'

    #########################################
    ########## 3rd PLOT DEFINITION ##########
    #########################################
    ### check for maximal values to be used in both x and y axis of the plot
    mprom = round(max(final_table['PromAPI']), 1)
    mgene = round(max(final_table['GeneAPI']), 1)
    if (mprom > mgene):
        mplot = mprom + 0.4
    else:
        mplot = mgene + 0.4

    ### set the figure dimension and other characteristics
    pl3 = figure(tools = tls, tooltips = tltips, plot_width = 750, plot_height = 750, x_range = (-0.4, mplot), y_range = (-0.4, mplot), x_axis_label = 'Intergenic regions', y_axis_label = 'Open reading frames')
    ### set plot main title
    pl3.title.text = '100 - average pairwise identity'
    ### set the Ide values to be in the plotted
    n = 0
    for v in list(main_plot_key.keys()):
        pl3.scatter('x0', 'y0', source = src, view = view[main_plot_key[v]], legend_label = v, color = Category10_10[n])
        n += 1
    ### set legend position and the attribute to hide/show it upon clicking in the legend
    pl3.legend.location = 'top_left'
    pl3.legend.background_fill_alpha = 0
    pl3.legend.click_policy = 'hide'

    #########################################
    ########## 4th PLOT DEFINITION ##########
    #########################################
    ### set the figure dimension and other characteristics
    pl4 = figure(tools = tls, tooltips = tltips, plot_width = 375, plot_height = 375, x_range = (-0.4, mplot), y_range = (-0.4, mplot), x_axis_label = 'Intergenic regions', y_axis_label = 'Open reading frames')
    ### set plot main title
    pl4.title.text = '100 - average pairwise identity'
    ### set the Ide values to be in the plotted
    n = 0
    for v in list(sec_plot_key.keys()):
        pl4.scatter('x0', 'y0', source = src, view = view2[sec_plot_key[v]], legend_label = v, color = Set1_3[n], muted_color = Set1_3[n], muted_alpha = 0.05)
        n += 1
    ### set legend position and the attribute to hide/show it upon clicking in the legend
    pl4.legend.location = 'top_left'
    pl4.legend.background_fill_alpha = 0
    pl4.legend.click_policy = 'mute'

    ### set layout, tab and panel for plot
    lay0 = layout(gridplot([[pl3, pl1], [pl4, pl2]]))
    tab = []
    tab.append(Panel(child = lay0, title = 'Main function groups'))

    pal = [Blues, Oranges, Greens, Reds, Purples, YlOrBr, PuRd]
    lay = []
    p = []
    i = 0
    for sub in subs:
        ### define source to enable linked selection of genes between the plots
        src = ColumnDataSource(data = dict(x0 = subs[sub]['PromAPI'], y0 = subs[sub]['GeneAPI'], x1 = subs[sub]['PromPSS'], y1 = subs[sub]['GenePSS'], z0 = subs[sub]['Gene'], z1 = subs[sub]['BC1'], z2 = subs[sub]['BC0-fun'], z3 = subs[sub]['BC1-fun'], z4 = subs[sub]['BC2-fun'], z5 = subs[sub]['BC3-fun'], z6 = subs[sub]['BC4-fun'], z7 = subs[sub]['Experiment']))
        ### extract BC groups' IDs from 'BC0' column
        bins = []
        for zero in subs[sub]['BC1']:
            if str(zero) not in bins:
                bins.append(str(zero))
        ### sort the groups alphabetically
        bins = sorted(bins)

        if len(bins) > 1:
            cols = pal[i]
            if i == 0:
                num = 256
            else:
                num = 9
            ### add function descriptions for the ID extracted from 'BC0' in a dictionary
            plot_key = {}
            for grp in bins:
                for k in all_keys.keys():
                    if grp == k:
                        plot_key[all_keys[k]] = grp

            ### define views to distinguish BC groups during plotting
            view = {}
            for pk in sorted(plot_key.values()):
                view[pk] = CDSView(source = src, filters = [GroupFilter(column_name = 'z1', group = pk)])

            #########################################
            ########## 1st PLOT DEFINITION ##########
            #########################################
            ### check for maximal values to be used in both x and y axis of the plot
            mprom = round(max(final_table['PromAPI']), 1)
            mgene = round(max(final_table['GeneAPI']), 1)
            if (mprom > mgene):
                mplot = mprom + 0.4
            else:
                mplot = mgene + 0.4
            ### set the figure dimension and other characteristics
            pl3 = figure(tools = tls, tooltips = tltips, plot_width = 750, plot_height = 750, x_range = (-0.4, mplot), y_range = (-0.4, mplot), x_axis_label = 'Intergenic regions', y_axis_label = 'Open reading frames')
            ### set plot main title
            pl3.title.text = '100 - average pairwise identity'
            ### set the Ide values to be in the plotted
            n = 0
            for v in list(plot_key.keys()):
                pl3.scatter('x0', 'y0', source = src, view = view[plot_key[v]], legend_label = v, color = cols[num][round(num  / 9 * n)])
                n += 1
            ### set legend position and the attribute to hide/show it upon clicking in the legend
            pl3.legend.location = 'top_left'
            pl3.legend.background_fill_alpha = 0
            pl3.legend.title = list(main_plot_key.keys())[i]
            pl3.legend.title_text_font_style = 'bold'
            pl3.legend.click_policy = 'hide'

            #########################################
            ########## 2nd PLOT DEFINITION ##########
            #########################################
            ### check for maximal values to be used in both x and y axis of the plot
            mprom = round(max(final_table['PromPSS']), 1)
            mgene = round(max(final_table['GenePSS']), 1)
            if (mprom > mgene):
                mplot = mprom + 0.01
            else:
                mplot = mgene + 0.01
            ### set the figure dimension and other characteristics
            pl4 = figure(tools = tls, tooltips = tltips, plot_width = 750, plot_height = 750, x_range = (-0.012, mplot), y_range = (-0.012, mplot), x_axis_label = 'Intergenic regions', y_axis_label = 'Open reading frames')
            ### set plot main title
            pl4.title.text = 'Proportion of segregating sites'
            ### set the Seg values to be in the plotted
            n = 0
            for v in list(plot_key.keys()):
                pl4.scatter('x1', 'y1', source = src, view = view[plot_key[v]], legend_label = v, color = cols[num][round(num  / 9 * n)])
                n += 1
            ### set legend position and the attribute to hide/show it upon clicking in the legend
            pl4.legend.location = 'top_left'
            pl4.legend.background_fill_alpha = 0
            pl4.legend.title = list(main_plot_key.keys())[i]
            pl4.legend.title_text_font_style = 'bold'
            pl4.legend.click_policy = 'hide'

            p.append(gridplot([[pl3, pl4]]))
            # tab.append(Panel(child = lay[i], title = list(main_plot_key.keys())[i]))
            i += 1

    lay = layout(column(p))
    tab.append(Panel(child = lay, title = 'Subgroups 1'))

    p = []

    box_data = final_table
    box_data.index = list(box_data.BC0)
    box1 = figure(plot_width = 600, plot_height = 600, y_range = list(main_plot_key.keys()), x_range = (-0.012, 0.31), x_axis_label = 'Proportion of segregating sites', y_axis_label = 'Main function groups')
    box1.title.text = 'Intergenic regions'
    box1.xaxis.major_label_orientation = 'vertical'

    mid = []
    top = []
    bottom = []
    upper = []
    lower = []
    for grp in list(main_plot_key.keys()):
        middle = statistics.median(box_data.PromPSS.loc[main_plot_key[grp]])
        sd = statistics.stdev(box_data.PromPSS.loc[main_plot_key[grp]])
        mid.append(middle)
        top.append((middle + sd))
        bottom.append((middle - sd))
        upper.append(max(box_data.PromPSS.loc[main_plot_key[grp]]))
        lower.append(min(box_data.PromPSS.loc[main_plot_key[grp]]))

    box1.hbar(y = list(main_plot_key.keys()), height = 0.5, left = bottom, right = top, color = Category10_10)
    box1.hbar(y = list(main_plot_key.keys()), height = 0.01, left = lower, right = upper, color = 'black')
    box1.hbar(y = list(main_plot_key.keys()), height = 0.5, left = mid, right = mid, color = 'black')
    box1.hbar(y = list(main_plot_key.keys()), height = 0.5, left = lower, right = lower, color = 'black')
    box1.hbar(y = list(main_plot_key.keys()), height = 0.5, left = upper, right = upper, color = 'black')

    box2 = figure(plot_width = 600, plot_height = 600, y_range = list(main_plot_key.keys()), x_range = (-0.012, 0.31), x_axis_label = 'Proportion of segregating sites', y_axis_label = 'Main function groups')
    box2.title.text = 'Open reading frames'
    box2.xaxis.major_label_orientation = 'vertical'

    mid = []
    top = []
    bottom = []
    upper = []
    lower = []
    for grp in list(main_plot_key.keys()):
        middle = statistics.median(box_data.GenePSS.loc[main_plot_key[grp]])
        sd = statistics.stdev(box_data.GenePSS.loc[main_plot_key[grp]])
        mid.append(middle)
        top.append((middle + sd))
        bottom.append((middle - sd))
        upper.append(max(box_data.GenePSS.loc[main_plot_key[grp]]))
        lower.append(min(box_data.GenePSS.loc[main_plot_key[grp]]))

    box2.hbar(y = list(main_plot_key.keys()), height = 0.5, left = bottom, right = top, color = Category10_10)
    box2.hbar(y = list(main_plot_key.keys()), height = 0.01, left = lower, right = upper, color = 'black')
    box2.hbar(y = list(main_plot_key.keys()), height = 0.5, left = mid, right = mid, color = 'black')
    box2.hbar(y = list(main_plot_key.keys()), height = 0.5, left = lower, right = lower, color = 'black')
    box2.hbar(y = list(main_plot_key.keys()), height = 0.5, left = upper, right = upper, color = 'black')
    ### execute plotting of both defined plots
    p.append(gridplot([[box1, box2]]))

    i = 0
    for sub in subs:
        box_data = subs[sub]
        box_data.index = list(box_data.BC1)
        ### extract BC groups' IDs from 'BC0' column
        bins = []
        for zero in box_data['BC1']:
            if str(zero) not in bins:
                bins.append(str(zero))
        ### sort the groups alphabetically
        bins = sorted(bins)

        if len(bins) > 1:
            ### add function descriptions for the ID extracted from 'BC0' in a dictionary
            plot_key = {}
            for grp in bins:
                for k in all_keys.keys():
                    if grp == k:
                        plot_key[all_keys[k]] = grp

            mid = []
            top = []
            bottom = []
            upper = []
            lower = []
            rem = []
            for grp in list(plot_key.keys()):
                a = pd.Series(box_data.PromPSS.loc[plot_key[grp]])
                if len(a) > 1:
                    middle = statistics.median(box_data.PromPSS.loc[plot_key[grp]])
                    sd = statistics.stdev(box_data.PromPSS.loc[plot_key[grp]])
                    mid.append(middle)
                    top.append((middle + sd))
                    bottom.append((middle - sd))
                    upper.append(max(box_data.PromPSS.loc[plot_key[grp]]))
                    lower.append(min(box_data.PromPSS.loc[plot_key[grp]]))
                else:
                    rem.append(grp)

            if len(rem) > 0:
                plot_key = {}
                for grp in bins:
                    for k in all_keys.keys():
                        if all([grp == k, all_keys[grp] not in rem]):
                            plot_key[all_keys[k]] = grp

            box1 = figure(plot_width = 600, plot_height = 600, y_range = list(plot_key.values()), x_range = (-0.012, 0.31), x_axis_label = 'Proportion of segregating sites', y_axis_label = 'Subgroups 1: ' + list(main_plot_key.keys())[i])
            box1.title.text = 'Intergenic regions'
            box1.xaxis.major_label_orientation = 'vertical'

            box1.hbar(y = list(plot_key.values()), height = 0.5, left = bottom, right = top, color = pal[i][len(list(plot_key.keys()))])
            box1.hbar(y = list(plot_key.values()), height = 0.01, left = lower, right = upper, color = 'black')
            box1.hbar(y = list(plot_key.values()), height = 0.5, left = mid, right = mid, color = 'black')
            box1.hbar(y = list(plot_key.values()), height = 0.5, left = lower, right = lower, color = 'black')
            box1.hbar(y = list(plot_key.values()), height = 0.5, left = upper, right = upper, color = 'black')

            mid = []
            top = []
            bottom = []
            upper = []
            lower = []
            for grp in list(plot_key.keys()):
                a = pd.Series(box_data.GenePSS.loc[plot_key[grp]])
                if len(a) > 1:
                    middle = statistics.median(box_data.GenePSS.loc[plot_key[grp]])
                    sd = statistics.stdev(box_data.GenePSS.loc[plot_key[grp]])
                    mid.append(middle)
                    top.append((middle + sd))
                    bottom.append((middle - sd))
                    upper.append(max(box_data.GenePSS.loc[plot_key[grp]]))
                    lower.append(min(box_data.GenePSS.loc[plot_key[grp]]))
                else:
                    rem.append(grp)

            box2 = figure(plot_width = 600, plot_height = 600, y_range = list(plot_key.values()), x_range = (-0.012, 0.31), x_axis_label = 'Proportion of segregating sites', y_axis_label = 'Subgroups 1: ' + list(main_plot_key.keys())[i])
            box2.title.text = 'Open reading frames'
            box2.xaxis.major_label_orientation = 'vertical'

            box2.hbar(y = list(plot_key.values()), height = 0.5, left = bottom, right = top, color = pal[i][len(list(plot_key.keys()))])
            box2.hbar(y = list(plot_key.values()), height = 0.01, left = lower, right = upper, color = 'black')
            box2.hbar(y = list(plot_key.values()), height = 0.5, left = mid, right = mid, color = 'black')
            box2.hbar(y = list(plot_key.values()), height = 0.5, left = lower, right = lower, color = 'black')
            box2.hbar(y = list(plot_key.values()), height = 0.5, left = upper, right = upper, color = 'black')
            ### execute plotting of both defined plots
            p.append(gridplot([[box1, box2]]))
            i += 1

    lay = layout(column(p))
    tab.append(Panel(child = lay, title = 'Boxplots'))
    see = Tabs(tabs = tab)
    show(see)

if __name__ == '__main__':
    main()

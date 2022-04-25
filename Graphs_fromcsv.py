from csv import list_dialects
from tkinter import VERTICAL
import sklearn
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
from scipy.stats import pearsonr
import lasio
import math
from numpy.polynomial.polynomial import polyfit
import seaborn as sns
import scipy.stats as stats
from scipy.stats import mannwhitneyu, normaltest
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import stats as stats
import pingouin as pg
import plotly.express as px
from bioinfokit.analys import stat
import statsmodels.api as sm
import pylab





df = pd.read_csv(r'data.csv')

def colony_comparison_body():
    df1 = df.dropna(subset = ['Colony', 'Bodymicrowatts'])
    list_ac22 = []
    list_acp2 = []
    
    for index, row in df1.iterrows():
        if ((str(row.Colony) == 'AC22') & (row.HWmm < 3)):
            list_ac22.append(row.Bodymicrowatts)
        elif ((str(row.Colony) == 'ACP2') & (row.HWmm < 3)):
            list_acp2.append(row.Bodymicrowatts)


    columns = [list_ac22, list_acp2]

    df_new = pd.DataFrame(list(zip(list_ac22,list_acp2)),
               columns =['ac22', 'acp2'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['ac22'], df_new['acp2'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'colony'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'colony': 'ac22'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'colony': 'acp2'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['colony'],
                          alpha=0.05)
    
    print(tukey)
    
    print(columns)
    # output
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#0000FF', '#00FF00', '#FFFF00', '#FF00FF']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
     
    # x-axis labels
    ax.set_xticklabels(['ac22', 'acp2'])
    ax.set_ylabel('Body Metabolism (uW)')
   
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
     
    # show plot
    plt.show()

def colony_comparison_brain():
    df1 = df.dropna(subset = ['Colony', 'Brainmicrowatts'])
    list_ac22 = []
    list_acp2 = []
    
    for index, row in df1.iterrows():
        if ((str(row.Colony) == 'AC22') & (row.HWmm < 3)):
            list_ac22.append(row.Brainmicrowatts)
        elif ((str(row.Colony) == 'ACP2') & (row.HWmm < 3)):
            list_acp2.append(row.Brainmicrowatts)
            
        
            
        
    
            
            
    columns = [list_ac22, list_acp2]


    df_new = pd.DataFrame(list(zip(list_ac22,list_acp2)),
               columns =['ac22', 'acp2'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['ac22'], df_new['acp2'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'colony'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'colony': 'ac22'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'colony': 'acp2'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['colony'],
                          alpha=0.05)
    
    print(tukey)
    print(columns)
    # output
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#0000FF', '#00FF00', '#FFFF00', '#FF00FF']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
     
    # x-axis labels
    ax.set_xticklabels(['ac22', 'acp2'])
    ax.set_ylabel('Brain Metabolism (uW)')
   
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
     
    # show plot
    plt.show()
    
def hw_bodymass():
    # Scatter plot function
    plt.scatter(df['HWmm'], df['bodyweight'], c='b', marker='o')
    plt.xlabel('HW (mm)', fontsize=16)
    plt.ylabel('body weight (mg)', fontsize=16)
    plt.title('Body mass vs head width', fontsize=20)
    plt.show()
    plt.draw()

def hw_brainmass():
    plt.scatter(df['HWmm'], df['brainweight'], c='b', marker='o')
    plt.xlabel('HW (mm)', fontsize=16)
    plt.ylabel('brain mass (mg)', fontsize=16)

    plt.show()

def brain_bodymass():
    df1 = df.drop(labels = ['Colony', 'Date', 'Brainmicrowatts', 'Bodymicrowatts'], axis = 1)
    df1 = df.dropna(subset = ['brainweight', 'bodyweight'])
    
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minims = df1[df1.HWmm<1]
    df_majors = df1[df1.HWmm>3]

    """
    df_l_medias = drop_outliers_IQR(df_l_medias)
    df_m_medias = drop_outliers_IQR(df_m_medias)
    df_s_medias = drop_outliers_IQR(df_s_medias)
    df_minims = drop_outliers_IQR(df_minims)
    df_majors = drop_outliers_IQR(df_majors)
    """
    
    plt.scatter(np.log(df_majors.bodyweight), np.log(df_majors.brainweight), c='#f0f9e8', marker='o', edgecolors = 'black')
    plt.scatter(np.log(df_l_medias.bodyweight), np.log(df_l_medias.brainweight), c='#bae4bc', marker='o', edgecolors = 'black')
    plt.scatter(np.log(df_m_medias.bodyweight), np.log(df_m_medias.brainweight), c = '#7bccc4', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_s_medias.bodyweight), np.log(df_s_medias.brainweight), c = '#43a2ca', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_minims.bodyweight), np.log(df_minims.brainweight), c = '#0868ac', marker = 'o', edgecolors = 'black')
    
    #making regression lines
    y_l = np.log(df_l_medias.brainweight)
    x_l = np.log(df_l_medias.bodyweight)
    y_m = np.log(df_m_medias.brainweight)
    x_m = np.log(df_m_medias.bodyweight)
    y_s = np.log(df_s_medias.brainweight)
    x_s = np.log(df_s_medias.bodyweight)
    y_minims = np.log(df_minims.brainweight)
    x_minims = np.log(df_minims.bodyweight)
    y_majors = np.log(df_majors.brainweight)
    x_majors = np.log(df_majors.bodyweight)
    
    x_df = np.concatenate([x_majors,x_l,x_m,x_s,x_minims], axis = None)
    y_df = np.concatenate([y_majors,y_l,y_m,y_s,y_minims], axis = None)


    x = np.array(x_df).reshape(-1,1)

    y = y_df
    lm = pg.linear_regression(x, y, as_dataframe=False)
    stats = pg.linear_regression(x, y, as_dataframe=True)
    #linreg = LinearRegression(fit_intercept=True)
    #obj = linreg.fit(x,y)
    trendline = lm['pred']
    plt.plot(x,trendline, color = 'black')
    #slope = linreg.coef_[0]
    #r_squared = r2_score(y, trendline)
    slope = stats.iloc[1][1]
    r_squared = stats.iloc[1][5]
    p_val = stats.iloc[1][4]
    
    print(stats)
    
    plt.annotate("slope = " + str(round(slope,3)) + "\n$r^2$ = " + str(round(r_squared,3)), [-.3,-2.4], fontsize = 14)
    plt.xlabel('log body mass (log[mg])', fontsize=20)
    plt.ylabel('log brain mass (log[mg])', fontsize=20)
    plt.legend(['allometry','3.0mm' , '2.4mm', '1.8mm','1.2mm','0.6mm'], fontsize = 14, shadow = True)
    plt.show()

  
def body_brain_metabolism():
    df_l_medias = df[df.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df[df.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df[df.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minims =df[df.HWmm<1]
    
    
    
    plt.scatter(df_l_medias['Bodymicrowatts'], df_l_medias['Brainmicrowatts'], c='b', marker='o')
    plt.scatter(df_m_medias['Bodymicrowatts'], df_m_medias['Brainmicrowatts'], c = 'r', marker = 'o')
    plt.scatter(df_s_medias['Bodymicrowatts'], df_s_medias['Brainmicrowatts'], c = 'g', marker = 'o')
    plt.scatter(df_minims['Bodymicrowatts'], df_minims['Brainmicrowatts'], c = 'purple', marker = 'o')
    
    plt.xlabel('body metabolism', fontsize=16)
    plt.ylabel('brain metabolism', fontsize=16)
    plt.legend(['2.4mm', '1.8mm', '1.2mm', 'minims'])
    plt.xticks((100,200,300,400,500))
    plt.title('Brain vs body metabolism', fontsize=20)
    plt.show()
 
def body_brain_metabolism_weighted():
    df_l_medias = df[df.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df[df.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df[df.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minims =df[df.HWmm<1]
    
    
    plt.scatter(df_l_medias.Bodymicrowatts/df_l_medias.bodyweight, df_l_medias.Brainmicrowatts/df_l_medias.brainweight, c='b', marker='o')
    plt.scatter(df_m_medias.Bodymicrowatts/df_m_medias.bodyweight, df_m_medias.Brainmicrowatts/df_m_medias.brainweight, c = 'r', marker = 'o')
    plt.scatter(df_s_medias.Bodymicrowatts/df_s_medias.bodyweight, df_s_medias.Brainmicrowatts/df_s_medias.brainweight, c = 'g', marker = 'o')
    plt.scatter(df_minims.Bodymicrowatts/df_minims.bodyweight, df_minims.Brainmicrowatts/df_minims.brainweight, c = 'purple', marker = 'o')
    
    plt.xlabel('body weighted metabolism', fontsize=16)
    plt.ylabel('brain weighted metabolism', fontsize=16)
    plt.legend(['2.4mm', '1.8mm', '1.2mm', 'minims'])
    plt.title('Brain vs body weighted metabolism', fontsize=20)
    plt.show()
  
def bp_body_brain_ratio():
    df1 = df.dropna(subset = ['bodyweight', 'brainweight'])
    list_major = []
    list_m = []
    list_l = []
    list_s = []
    list_minim = []
    for index, row in df1.iterrows():
        if (row.HWmm >= 2.25) & (row.HWmm <= 2.54):
            list_l.append(row.brainweight/row.bodyweight)
        elif (row.HWmm >= 1.65) & (row.HWmm <= 1.94):
            list_m.append(row.brainweight/row.bodyweight)
        elif (row.HWmm >= 1.05) & (row.HWmm <= 1.34):
            list_s.append(row.brainweight/row.bodyweight)
        elif (row.HWmm < 1):
            list_minim.append(row.brainweight/row.bodyweight)
        elif (row.HWmm > 3):
            list_major.append(row.brainweight/row.bodyweight)
            
            
    columns = [list_minim, list_s, list_m, list_l, list_major]
   
    df_new = pd.DataFrame(list(zip(list_major,list_l,list_m,list_s,list_minim)),
               columns =['majors', 'large','medium','small','minims'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['majors'], df_new['large'], df_new['medium'], df_new['small'], df_new['minims'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'caste'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'minim'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'small'}, ignore_index=True)  
            elif count == 2:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'medium'}, ignore_index=True)
            elif count == 3:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'large'}, ignore_index=True)
            elif count == 4:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'major'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['caste'],
                          alpha=0.05)
    
    print(tukey)
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#f0f9e8', '#bae4bc', '#7bccc4', '#43a2ca', '#0868ac']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='black', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
        
    print(bp['medians'][0].get_ydata())
    print(bp['medians'][1].get_ydata())
    print(bp['medians'][2].get_ydata())
    print(bp['medians'][3].get_ydata())
    print(bp['medians'][4].get_ydata())
    
    
    x1, x2 = 3,4
    y, h, col = 0.12,0.007, 'k'
    ax.plot([x1, x1, x2-(2*h), x2-(2*h)], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, 'n.s.', ha='center', va='bottom', color='black')
    
    
    x1, x2 = 4,5
    y, h, col = 0.12,0.007, 'k'
    ax.plot([x1, x1, x2+(2*h), x2+(2*h)], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, 'n.s.', ha='center', va='bottom', color='black')
    
    # x-axis labels
    ax.set_xticklabels(['minims' , 'small','medium', 'large', 'majors'], fontsize=16)
    ax.set_ylabel('brain-to-body mass ratio (unitless)', fontsize=16)
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3] ,bp['boxes'][4]], ['0.6mm','1.2mm','1.8mm','2.4mm', '3.0mm'], 
            shadow = True, title = 'worker subcaste headwidth')
 
 

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    # show plot
    plt.show()
def bp_brain_mass():
    df1 = df.dropna(subset = [ 'brainweight'])
    list_m = []
    list_l = []
    list_s = []
    list_minim = []
    list_major = []
    for index, row in df1.iterrows():
        if (row.HWmm >= 2.25) & (row.HWmm <= 2.54):
            list_l.append(row.brainweight)
        elif (row.HWmm >= 1.65) & (row.HWmm <= 1.94):
            list_m.append(row.brainweight)
        elif (row.HWmm >= 1.05) & (row.HWmm <= 1.34):
            list_s.append(row.brainweight)
        elif (row.HWmm < 1):
            list_minim.append(row.brainweight)
        elif (row.HWmm > 3):
            list_major.append(row.brainweight)
            
    columns = [list_minim, list_s, list_m, list_l, list_major]
    sm.qqplot(np.array(columns[1]), line='r')
    pylab.show()
    
    df_new = pd.DataFrame(list(zip(list_major,list_l,list_m,list_s,list_minim)),
               columns =['majors', 'large','medium','small','minims'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['majors'], df_new['large'], df_new['medium'], df_new['small'], df_new['minims'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'caste'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'minim'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'small'}, ignore_index=True)  
            elif count == 2:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'medium'}, ignore_index=True)
            elif count == 3:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'large'}, ignore_index=True)
            elif count == 4:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'major'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['caste'],
                          alpha=0.05)
    
    print(tukey)
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#f0f9e8', '#bae4bc', '#7bccc4', '#43a2ca', '#0868ac']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='black', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
        
    print(bp['medians'][0].get_ydata())
    print(bp['medians'][1].get_ydata())
    print(bp['medians'][2].get_ydata())
    print(bp['medians'][3].get_ydata())
    print(bp['medians'][4].get_ydata())
    
 
    
     
    # x-axis labels
    ax.set_xticklabels(['minims' , 'small','medium', 'large', 'majors'], fontsize=20)
    ax.set_ylabel('brain mass (mg)', fontsize=20)
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3], bp['boxes'][4]], ['0.6mm','1.2mm','1.8mm','2.4mm', '3.0mm'], 
               loc='upper left', shadow = True, title = 'worker subcaste headwidth', fontsize = 16)
 
    x1, x2 = 2,3
    y, h, col = 0.25,0.01, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, 'n.s.', ha='center', va='bottom', color='black')
    
    
    x1, x2 = 3,4
    y, h, col = 0.3 ,0.01, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, 'n.s.', ha='center', va='bottom', color='black')
 

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    # show plot
    plt.show()
   
def bp_brainmetabolism():
    df1 = df.dropna(subset = ['Brainmicrowatts'])
    list_m = []
    list_l = []
    list_s = []
    list_minim = []
    list_majors = []
    
    for index, row in df1.iterrows():
        if (row.HWmm >= 2.25) & (row.HWmm <= 2.54):
            list_l.append(row.Brainmicrowatts)
        elif (row.HWmm >= 1.65) & (row.HWmm <= 1.94):
            list_m.append(row.Brainmicrowatts)
        elif (row.HWmm >= 1.05) & (row.HWmm <= 1.34):
            list_s.append(row.Brainmicrowatts)
        elif (row.HWmm <1):
            list_minim.append(row.Brainmicrowatts)
        elif (row.HWmm >3):
            list_majors.append(row.Brainmicrowatts)
            
            
    columns = [list_minim, list_s, list_m, list_l, list_majors]
    for i in columns:
        w, pval = stats.shapiro(i)
        print(w, pval)
    df_new = pd.DataFrame(list(zip(list_majors,list_l,list_m,list_s,list_minim)),
               columns =['majors', 'large','medium','small','minims'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['majors'], df_new['large'], df_new['medium'], df_new['small'], df_new['minims'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'caste'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'minim'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'small'}, ignore_index=True)  
            elif count == 2:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'medium'}, ignore_index=True)
            elif count == 3:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'large'}, ignore_index=True)
            elif count == 4:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'major'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    
    
    res = stat()
    res.bartlett(df=df_t, res_var='metabolism', xfac_var='caste')
    print(res.bartlett_summary)
    res.levene(df=df_t, res_var = 'metabolism', xfac_var= 'caste')
    print(res.levene_summary)
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['caste'],
                          alpha=0.05)
    print(tukey)
    sm.qqplot(np.array(columns[4]), line='r')
    pylab.show()
    
    
    
    
    # output
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#f0f9e8', '#bae4bc', '#7bccc4', '#43a2ca', '#0868ac']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
     
    # x-axis labels
    ax.set_xticklabels(["minims", 'small', 'medium', 'large', 'majors'], fontsize = 20)
    ax.set_ylabel('Brain Metabolism (uW)', fontsize = 20)
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3], bp['boxes'][4]], 
               ['0.6mm','1.2mm','1.8mm','2.4mm', '3.0mm'], shadow = True, title = 'worker subcaste headwidth', fontsize = 14)
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    x1, x2 = 1,2
    y, h, col = 1.75,0.1, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, '*', ha='center', va='bottom', color='black')
    
    
    x1, x2 = 4,5
    y, h, col = 2,0.1, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, '*', ha='center', va='bottom', color='black')
    
    # show plot
    plt.show()

def bp_bodymetabolism():
    df1 = df.dropna(subset = ['Bodymicrowatts'])
    list_m = []
    list_l = []
    list_s = []
    list_minim = []
    list_major = []
    
    for index, row in df1.iterrows():
        if (row.HWmm >= 2.25) & (row.HWmm <= 2.54):
            list_l.append(row.Bodymicrowatts)
        elif (row.HWmm >= 1.65) & (row.HWmm <= 1.94):
            list_m.append(row.Bodymicrowatts)
        elif (row.HWmm >= 1.05) & (row.HWmm <= 1.34):
            list_s.append(row.Bodymicrowatts)
        elif (row.HWmm <1):
            list_minim.append(row.Bodymicrowatts)
        elif (row.HWmm >3):
            list_major.append(row.Bodymicrowatts)   
            
    
    columns = [list_minim, list_s, list_m, list_l, list_major]
    for i in columns:
        w, pval = stats.shapiro(i)
        print(w, pval)
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#0000FF', '#00FF00', '#FFFF00', '#FF00FF']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
     
    # x-axis labels
    ax.set_xticklabels(['minim', 'small', 'medium', 'large', 'major'])
    ax.set_ylabel('Body Metabolism (uW)')
    
    columns = [list_minim, list_s, list_m, list_l, list_major]
    for i in columns:
        w, pval = stats.shapiro(i)
        print(w, pval)
    df_new = pd.DataFrame(list(zip(list_major,list_l,list_m,list_s,list_minim)),
               columns =['majors', 'large','medium','small','minims'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['majors'], df_new['large'], df_new['medium'], df_new['small'], df_new['minims'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'caste'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'minim'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'small'}, ignore_index=True)  
            elif count == 2:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'medium'}, ignore_index=True)
            elif count == 3:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'large'}, ignore_index=True)
            elif count == 4:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'major'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    
    
    res = stat()
    res.bartlett(df=df_t, res_var='metabolism', xfac_var='caste')
    print(res.bartlett_summary)
    res.levene(df=df_t, res_var = 'metabolism', xfac_var= 'caste')
    print(res.levene_summary)
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['caste'],
                          alpha=0.05)
    
    print(tukey)
   
    #stats - minim large
    t0, p0 = stats.ttest_ind(list_l, list_minim)
    x1, x2 = 1,4
    y, h, col = 620,5, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1+x2)*.5, y+h, 't= ' + str(round(t0,4)) + ', ' 'p= ' + str(round(p0,4)), ha='center', va='bottom', color=col)
 
    #stats - minim small
    t1, p1 = stats.ttest_ind(list_m, list_minim)
    x1, x2 = 1,3
    y, h, col = 670,5, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1+x2)*.5, y+h, 't= ' + str(round(t1,4)) + ', ' 'p= ' + str(round(p1,4)), ha='center', va='bottom', color=col)
 
 
    #stats - minim large
    t2, p2 = stats.ttest_ind(list_s, list_minim)
    x1, x2 = 1,2
    y, h, col = 730,5, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1+x2)*.5, y+h, 't= ' + str(round(t2,4)) + ', ' 'p= ' + str(round(p2,8)), ha='center', va='bottom', color=col)
 
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3]], ['0.6mm','1.2mm','1.8mm','2.4mm'], shadow = True, title = 'worker subcaste headwidth')
 
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
     
    # show plot
    plt.show()

def bp_brainmetabolismweighted():
    df1 = df.dropna(subset = ['Brainmicrowatts', 'brainweight'])
    list_m = []
    list_l = []
    list_s = []
    list_minim = []
    list_majors = []
    for index, row in df1.iterrows():
        if (row.HWmm >= 2.25) & (row.HWmm <= 2.54):
            list_l.append(row.Brainmicrowatts/row.brainweight)
        elif (row.HWmm >= 1.65) & (row.HWmm <= 1.94):
            list_m.append(row.Brainmicrowatts/row.brainweight)
        elif (row.HWmm >= 1.05) & (row.HWmm <= 1.34):
            list_s.append(row.Brainmicrowatts/row.brainweight)
        elif (row.HWmm < 1):
            list_minim.append(row.Brainmicrowatts/row.brainweight)
        elif (row.HWmm > 3):
            list_majors.append(row.Brainmicrowatts/row.brainweight)
            
            
    columns = [list_minim, list_s, list_m, list_l, list_majors]
    df_new = pd.DataFrame(list(zip(list_majors,list_l,list_m,list_s,list_minim)),
               columns =['majors', 'large','medium','small','minims'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['majors'], df_new['large'], df_new['medium'], df_new['small'], df_new['minims'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'caste'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'minim'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'small'}, ignore_index=True)  
            elif count == 2:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'medium'}, ignore_index=True)
            elif count == 3:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'large'}, ignore_index=True)
            elif count == 4:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'major'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['caste'],
                          alpha=0.05)
    
    print(tukey)
    
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#f0f9e8', '#bae4bc', '#7bccc4', '#43a2ca', '#0868ac']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='black', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
        
    
    # x-axis labels
    ax.set_xticklabels(['minims' , 'small','medium', 'large', 'majors'], fontsize=20)
    ax.set_ylabel('Brain mass-specific metabolic rate (uW/mg)', fontsize=20)
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3], bp['boxes'][4]], ['0.6mm','1.2mm','1.8mm','2.4mm', '3.0mm'], 
            shadow = True, title = 'worker subcaste headwidth')
 
    x1, x2 = 2,3
    y, h, col = 15.5,0.5, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, 'n.s.', ha='center', va='bottom', color='black')
    
    
    x1, x2 = 4,5
    y, h, col = 15.5 ,0.5, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, 'n.s.', ha='center', va='bottom', color='black')

    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    # show plot
    plt.show()
    
def bp_bodymetabolismweighted():
    df1 = df.dropna(subset = ['Bodymicrowatts', 'bodyweight'])
    list_m = []
    list_l = []
    list_s = []
    list_minim = []
    for index, row in df1.iterrows():
        if (row.HWmm >= 2.25) & (row.HWmm <= 2.54):
            list_l.append(row.Bodymicrowatts/row.bodyweight)
        elif (row.HWmm >= 1.65) & (row.HWmm <= 1.94):
            list_m.append(row.Bodymicrowatts/row.bodyweight)
        elif (row.HWmm >= 1.05) & (row.HWmm <= 1.34):
            list_s.append(row.Bodymicrowatts/row.bodyweight)
        elif (row.HWmm < 1):
            list_minim.append(row.Bodymicrowatts/row.bodyweight)
            
    columns = [list_minim, list_s, list_m, list_l]
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True, vert = 0)

    colors = ['#0000FF', '#00FF00', '#FFFF00', '#FF00FF']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
     
    # x-axis labels
    ax.set_yticklabels(['minims', 'small', 'medium', 'large'])
    ax.set_xlabel('Body metabolism per mass (uW/mg)')
 
    # Adding title
    plt.title("Customized box plot")
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
     
    # show plot
    plt.show()
    
def bp_brainbodymetabolism():
    df1 = df.dropna(subset = ['Bodymicrowatts', 'Brainmicrowatts'])
    list_m = []
    list_l = []
    list_s = []
    list_minim = []
    list_majors = []
    
    for index, row in df1.iterrows():
        if (row.HWmm >= 2.25) & (row.HWmm <= 2.54):
            list_l.append(row.Brainmicrowatts/row.Bodymicrowatts)
        elif (row.HWmm >= 1.65) & (row.HWmm <= 1.94):
            list_m.append(row.Brainmicrowatts/row.Bodymicrowatts)
        elif (row.HWmm >= 1.05) & (row.HWmm <= 1.34):
            list_s.append(row.Brainmicrowatts/row.Bodymicrowatts)
        elif (row.HWmm < 1):
            list_minim.append(row.Brainmicrowatts/row.Bodymicrowatts)
        elif (row.HWmm > 3):
            list_majors.append(row.Brainmicrowatts/row.Bodymicrowatts)
            
    columns = [list_minim, list_s, list_m, list_l, list_majors]
    
 
    df_new = pd.DataFrame(list(zip(list_majors,list_l,list_m,list_s,list_minim)),
               columns =['majors', 'large','medium','small','minims'])
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(df_new['majors'], df_new['large'], df_new['medium'], df_new['small'], df_new['minims'])
    print(fvalue, pvalue)
    
    df_t = pd.DataFrame(columns = ['metabolism', 'caste'])
    count = 0
    for i in columns:
        for s in range(0,len(i)):
            if count == 0:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'minim'}, ignore_index=True)
            elif count == 1:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'small'}, ignore_index=True)  
            elif count == 2:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'medium'}, ignore_index=True)
            elif count == 3:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'large'}, ignore_index=True)
            elif count == 4:
                df_t = df_t.append({'metabolism': columns[count][s], 'caste': 'major'}, ignore_index=True)  
        count += 1      
    # perform multiple pairwise comparison (Tukey's HSD)
    # unequal sample size data, tukey_hsd uses Tukey-Kramer test
    tukey = pairwise_tukeyhsd(endog=df_t['metabolism'],
                          groups=df_t['caste'],
                          alpha=0.05)
    
    print(tukey)
    
    
    
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#f0f9e8', '#bae4bc', '#7bccc4', '#43a2ca', '#0868ac']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
    
    print(bp['medians'][0].get_ydata())
    print(bp['medians'][1].get_ydata())
    print(bp['medians'][2].get_ydata())
    print(bp['medians'][3].get_ydata())
    print(bp['medians'][4].get_ydata())
    
    # x-axis labels
    ax.set_xticklabels(['minims', 'small', 'medium', 'large', 'majors'], fontsize = 20)
    ax.set_ylabel('Brain over body metabolism (unitless)', fontsize = 20)
    #plt.legend(['0.6mm','1.2mm','1.8mm','2.4mm'], shadow = True, title = 'worker subcaste headwidth')
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3], bp['boxes'][4]], 
               ['0.6mm','1.2mm','1.8mm','2.4mm', '3.0mm'], loc='upper right', shadow = True, title = 'worker subcaste headwidth', fontsize = 14)
 
 
    x1, x2 = 1,2
    y, h, col = 0.04,0.003, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, '**', ha='center', va='bottom', color='black')
    
    
    x1, x2 = 1,3
    y, h, col = 0.045 ,0.003, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h, '**', ha='center', va='bottom', color='black')
   
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
     
    # show plot
    plt.show()   
    
def body_brain_metabolism_weighted_allometry():
    
    df1 = df.dropna(subset = ['Brainmicrowatts', 'Bodymicrowatts'])
    
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minim = df1[df1.HWmm <1]
    

    plt.scatter(np.log(df_l_medias.Bodymicrowatts/df_l_medias.bodyweight), np.log(df_l_medias.Brainmicrowatts/df_l_medias.brainweight), c='b', marker='o')
    plt.scatter(np.log(df_m_medias.Bodymicrowatts/df_m_medias.bodyweight), np.log(df_m_medias.Brainmicrowatts/df_m_medias.brainweight), c = 'r', marker = 'o')
    plt.scatter(np.log(df_s_medias.Bodymicrowatts/df_s_medias.bodyweight), np.log(df_s_medias.Brainmicrowatts/df_s_medias.brainweight), c = 'g', marker = 'o')
    plt.scatter(np.log(df_minim.Bodymicrowatts/df_minim.bodyweight), np.log(df_minim.Brainmicrowatts/df_minim.brainweight), c = 'purple', marker = 'o')
    
    
    #making regression lines
    x_l = np.log(df_l_medias.Bodymicrowatts/df_l_medias.bodyweight)
    y_l = np.log(df_l_medias.Brainmicrowatts/df_l_medias.brainweight)
    x_m = np.log(df_m_medias.Bodymicrowatts/df_m_medias.bodyweight)
    y_m = np.log(df_m_medias.Brainmicrowatts/df_m_medias.brainweight)
    x_s = np.log(df_s_medias.Bodymicrowatts/df_s_medias.bodyweight)
    y_s = np.log(df_s_medias.Brainmicrowatts/df_s_medias.brainweight) 
    x_minim = np.log(df_minim.Bodymicrowatts/df_minim.bodyweight)
    y_minim = np.log(df_minim.Brainmicrowatts/df_minim.brainweight) 
    

    #check for NaN
    idx_l = np.isfinite(x_l) & np.isfinite(y_l)
    idx_m = np.isfinite(x_m) & np.isfinite(y_m)
    idx_s = np.isfinite(x_s) & np.isfinite(y_s)
    idx_minim = np.isfinite(x_minim) & np.isfinite(y_minim)
    
    
    #m_l, b_l = np.polyfit(x_l[idx_l],y_l[idx_l],1)
    #plt.plot(x_l, m_l*x_l + b_l, color = 'blue')
    
    #m_m, b_m = np.polyfit(x_m[idx_m],y_m[idx_m],1)
    #plt.plot(x_m, m_m*x_m + b_m, color = 'red')
    
    #m_s, b_s = np.polyfit(x_s[idx_s],y_s[idx_s],1)
    #plt.plot(x_s, m_s*x_s + b_s, color = 'green' )


    
    plt.xlabel('log body weighted metabolism (log[uW/mg])', fontsize=16)
    plt.ylabel('log brain weighted metabolism (log[uW/mg])', fontsize=16)
    plt.legend(['2.4mm', '1.8mm', '1.2mm', 'minims'])
    plt.title('Brain vs body weighted metabolism', fontsize=20)
    plt.show()
    
def brain_met_mass_allometry():
    
    df1 = df.dropna(subset = ['Brainmicrowatts', 'brainweight'])
    
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minim = df1[df1.HWmm<1]

    plt.scatter(np.log(df_l_medias.brainweight), np.log(df_l_medias.Brainmicrowatts), c='b', marker='o')
    plt.scatter(np.log(df_m_medias.brainweight), np.log(df_m_medias.Brainmicrowatts), c = 'r', marker = 'o')
    plt.scatter(np.log(df_s_medias.brainweight), np.log(df_s_medias.Brainmicrowatts), c = 'g', marker = 'o')
    plt.scatter(np.log(df_minim.brainweight), np.log(df_minim.Brainmicrowatts), c = 'purple', marker = 'o')
    
    
    #making regression lines
    x_l = np.log(df_l_medias.brainweight)
    y_l = np.log(df_l_medias.Brainmicrowatts)
    x_m = np.log(df_m_medias.brainweight)
    y_m = np.log(df_m_medias.Brainmicrowatts)
    x_s = np.log(df_s_medias.brainweight)
    y_s = np.log(df_s_medias.Brainmicrowatts)
    x_minim = np.log(df_minim.brainweight)
    y_minim = np.log(df_minim.Brainmicrowatts)
    
    x_df = np.concatenate((x_l,x_m,x_s,x_minim), axis = None)
    y_df = np.concatenate((y_l,y_m,y_s,y_minim), axis = None)
    
    
    #check for NaN
    idx_l = np.isfinite(x_l) & np.isfinite(y_l)
    idx_m = np.isfinite(x_m) & np.isfinite(y_m)
    idx_s = np.isfinite(x_s) & np.isfinite(y_s)
    idx_df = np.isfinite(x_df) & np.isfinite(y_df) 
    
    
    """
    m_l, b_l = np.polyfit(x_l[idx_l],y_l[idx_l],1)
    plt.plot(x_l, m_l*x_l + b_l, color = 'blue')
    
    m_m, b_m = np.polyfit(x_m[idx_m],y_m[idx_m],1)
    plt.plot(x_m, m_m*x_m + b_m, color = 'red')
    
    m_s, b_s = np.polyfit(x_s[idx_s],y_s[idx_s],1)
    plt.plot(x_s, m_s*x_s + b_s, color = 'green')
    """
    

    #m_full, b_full = np.polyfit(x_df[idx_df],y_df[idx_df],1)
    #plt.plot(x_df, m_full*x_df + b_full, color = 'purple')
    x = np.array(x_df[idx_df]).reshape(-1,1)
    y = y_df[idx_df]
    linreg = LinearRegression(fit_intercept=False)
    obj = linreg.fit(x,y)
    trendline = linreg.predict(x)
    plt.plot(x,trendline, color = 'purple')
    slope = linreg.coef_[0]
    r_squared = r2_score(y, trendline)
    print('slope: ' + str(slope))
    print('r_squared' + str(r_squared))
    
    #plt.annotate("Full dataset: y = " + str(round(m_full,3)) + "x + " + str(round(b_full,3)), [-2.4,0.4])
    
    
    plt.xlabel('log brain mass (log[mg])', fontsize=16)
    plt.ylabel('log brain metabolism (log[uW])', fontsize=16)
    plt.legend(['2.4mm', '1.8mm', '1.2mm', 'full dataset'])
    plt.title('Brain metabolism allometry', fontsize=20)
    plt.show()

def brain_met_body_mass_allometry():
    
    df1 = df.dropna(subset = ['Brainmicrowatts', 'bodyweight'])
    
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]

    #making regression lines
    x_l = np.log(df_l_medias.bodyweight)
    y_l = np.log(df_l_medias.Brainmicrowatts)
    x_m = np.log(df_m_medias.bodyweight)
    y_m = np.log(df_m_medias.Brainmicrowatts)
    x_s = np.log(df_s_medias.bodyweight)
    y_s = np.log(df_s_medias.Brainmicrowatts)
    
    x_df = np.concatenate((x_l,x_m,x_s), axis = None)
    y_df = np.concatenate((y_l,y_m,y_s), axis = None)

    plt.scatter(x_l, y_l, c = 'b', marker = 'o')
    plt.scatter(x_m, y_m, c = 'r', marker = 'o')
    plt.scatter(x_s, y_s, c = 'g', marker = 'o')
    
 
    #check for NaN
    idx_l = np.isfinite(x_l) & np.isfinite(y_l)
    idx_m = np.isfinite(x_m) & np.isfinite(y_m)
    idx_s = np.isfinite(x_s) & np.isfinite(y_s)
    idx_df = np.isfinite(x_df) & np.isfinite(y_df) 
    
    
    m_l, b_l = np.polyfit(x_l[idx_l],y_l[idx_l],1)
    plt.plot(x_l, m_l*x_l + b_l, color = 'blue')
    
    m_m, b_m = np.polyfit(x_m[idx_m],y_m[idx_m],1)
    plt.plot(x_m, m_m*x_m + b_m, color = 'red')
    
    m_s, b_s = np.polyfit(x_s[idx_s],y_s[idx_s],1)
    plt.plot(x_s, m_s*x_s + b_s, color = 'green' )

    #m_full, b_full = np.polyfit(x_df[idx_df],y_df[idx_df],1)
    #plt.plot(x_df, m_full*x_df + b_full, color = 'purple')
    
    x = np.array(x_df[idx_df]).reshape(-1,1)
    y = y_df[idx_df]
    linreg = LinearRegression(fit_intercept=False)
    obj = linreg.fit(x,y)
    trendline = linreg.predict(x)
    plt.plot(x,trendline, color = 'purple')
    slope = linreg.coef_[0]
    r_squared = r2_score(y, trendline)
    print('slope: ' + str(slope))
    print('r_squared' + str(r_squared))
    
    plt.xlabel('log body mass (log[mg])', fontsize=16)
    plt.ylabel('log brain metabolism (log[uW])', fontsize=16)
    plt.legend(['2.4mm', '1.8mm', '1.2mm', 'full dataset'])
    plt.title('Brain metabolism to body mass allometry', fontsize=20)
    plt.show()
        
def body_met_mass_allometry():
        
    df1 =  df.drop(labels = ['Colony', 'Date', 'brainweight', 'Brainmicrowatts'], axis = 1)
    df1 = df1.dropna(subset = ['Bodymicrowatts', 'bodyweight' ])
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minim = df1[df1.HWmm<1]
    df_majors = df1[df1.HWmm>3]
    

    df_l_medias = drop_outliers_IQR(df_l_medias)
    df_m_medias = drop_outliers_IQR(df_m_medias)
    df_s_medias = drop_outliers_IQR(df_s_medias)
    df_minim = drop_outliers_IQR(df_minim)
    df_majors = drop_outliers_IQR(df_majors)

    plt.scatter(np.log(df_majors.bodyweight), np.log(df_majors.Bodymicrowatts), c = '#f0f9e8', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_l_medias.bodyweight), np.log(df_l_medias.Bodymicrowatts), c='#bae4bc', marker='o', edgecolors = 'black')
    plt.scatter(np.log(df_m_medias.bodyweight), np.log(df_m_medias.Bodymicrowatts), c = '#7bccc4', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_s_medias.bodyweight), np.log(df_s_medias.Bodymicrowatts), c = '#43a2ca', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_minim.bodyweight), np.log(df_minim.Bodymicrowatts), c = '#0868ac', marker = 'o', edgecolors = 'black')
    
    
    #making regression lines
    x_l = np.log(df_l_medias.bodyweight)
    y_l = np.log(df_l_medias.Bodymicrowatts)
    x_m = np.log(df_m_medias.bodyweight)
    y_m = np.log(df_m_medias.Bodymicrowatts)
    x_s = np.log(df_s_medias.bodyweight)
    y_s = np.log(df_s_medias.Bodymicrowatts)
    x_minim = np.log(df_minim.bodyweight)
    y_minim = np.log(df_minim.Bodymicrowatts)
    x_maj = np.log(df_majors.bodyweight)
    y_maj = np.log(df_majors.Bodymicrowatts)
    
    x_df = np.concatenate((x_maj, x_l,x_m,x_s,x_minim), axis = None)
    y_df = np.concatenate((y_maj, y_l,y_m,y_s,y_minim), axis = None)

    x = np.array(x_df).reshape(-1,1)
    y = y_df
    """
    linreg = LinearRegression(fit_intercept=True)
    obj = linreg.fit(x,y)
    trendline = linreg.predict(x)
    plt.plot(x,trendline, color = 'black')
    slope = linreg.coef_[0]
    r_squared = r2_score(y, trendline)
    """
    
    lm = pg.linear_regression(x, y, as_dataframe=False)
    stats = pg.linear_regression(x, y, as_dataframe=True)
    #linreg = LinearRegression(fit_intercept=True)
    #obj = linreg.fit(x,y)
    trendline = lm['pred']
    plt.plot(x,trendline, color = 'black')
    #slope = linreg.coef_[0]
    #r_squared = r2_score(y, trendline)
    slope = stats.iloc[1][1]
    r_squared = stats.iloc[1][5]
    p_val = stats.iloc[1][4]
    print(stats)
    
    plt.annotate("slope = " + str(round(slope,3)) + "\n$r^2$ = " + str(round(r_squared,3)), [-0.6,4.5], fontsize = 14)
    plt.xlabel('log body mass (log[mg])', fontsize=20)
    plt.ylabel('log body metabolism (log[uW])', fontsize=20)
    plt.legend(['allometry', '3.0mm', '2.4mm', '1.8mm','1.2mm','0.6mm'], fontsize = 14, shadow = True)
    plt.show()


def body_brain_metabolism_allometry():
    
    df1 = df.dropna(subset = ['Brainmicrowatts', 'Bodymicrowatts'])
    
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minims = df1[df1.HWmm<1]
    df_majors = df1[df1.HWmm>3]
    

    plt.scatter(np.log(df_l_medias.Bodymicrowatts), np.log(df_l_medias.Brainmicrowatts), c='b', marker='o')
    plt.scatter(np.log(df_m_medias.Bodymicrowatts), np.log(df_m_medias.Brainmicrowatts), c = 'r', marker = 'o')
    plt.scatter(np.log(df_s_medias.Bodymicrowatts), np.log(df_s_medias.Brainmicrowatts), c = 'g', marker = 'o')
    plt.scatter(np.log(df_minims.Bodymicrowatts), np.log(df_minims.Brainmicrowatts), c = 'purple', marker = 'o')
    plt.scatter(np.log(df_majors.Bodymicrowatts), np.log(df_majors.Brainmicrowatts), c = 'yellow', marker = 'o')
    
    #making regression lines
    x_l = np.log(df_l_medias.Bodymicrowatts)
    y_l = np.log(df_l_medias.Brainmicrowatts)
    x_m = np.log(df_m_medias.Bodymicrowatts)
    y_m = np.log(df_m_medias.Brainmicrowatts)
    x_s = np.log(df_s_medias.Bodymicrowatts)
    y_s = np.log(df_s_medias.Brainmicrowatts) 
    x_minim = np.log(df_minims.Bodymicrowatts)
    y_minim = np.log(df_minims.Brainmicrowatts) 
    x_maj = np.log(df_majors.Bodymicrowatts)
    y_maj = np.log(df_majors.Brainmicrowatts) 

    x_df = np.log(df1.Bodymicrowatts)
    y_df = np.log(df1.Brainmicrowatts)
    
    
    #check for NaN
    idx_maj = np.isfinite(x_maj) & np.isfinite(y_maj)
    idx_l = np.isfinite(x_l) & np.isfinite(y_l)
    idx_m = np.isfinite(x_m) & np.isfinite(y_m)
    idx_s = np.isfinite(x_s) & np.isfinite(y_s)
    idx_minim = np.isfinite(x_minim) & np.isfinite(y_minim)
    idx_df = np.isfinite(x_df) & np.isfinite(y_df) 
    
    
    #m_l, b_l = np.polyfit(x_l[idx_l],y_l[idx_l],1)
    #plt.plot(x_l, m_l*x_l + b_l, color = 'blue')
    
    #m_m, b_m = np.polyfit(x_m[idx_m],y_m[idx_m],1)
    #plt.plot(x_m, m_m*x_m + b_m, color = 'red')
    
    #m_s, b_s = np.polyfit(x_s[idx_s],y_s[idx_s],1)
    #plt.plot(x_s, m_s*x_s + b_s, color = 'green' )

    #m_full, b_full = np.polyfit(x_df[idx_df],y_df[idx_df],1)
    #plt.plot(x_df, m_full*x_df + b_full, color = 'purple')
    """
    x = np.array(x_df[idx_df]).reshape(-1,1)
    y = y_df[idx_df]
    linreg = LinearRegression(fit_intercept=False)
    obj = linreg.fit(x,y)
    trendline = linreg.predict(x)
    plt.plot(x,trendline, color = 'purple')
    slope = linreg.coef_[0]
    r_squared = r2_score(y, trendline)
    print('slope: ' + str(slope))
    print('r_squared' + str(r_squared))
    
    
    plt.annotate("slope = " + str(round(slope,4)), [3.1,0])
    #+ "\nr_squared = " + str(round(r_squared,4))
    
    """
    
    plt.xlabel('log body metabolism (log[uW])', fontsize=16)
    plt.ylabel('log brain metabolism (log[uW])', fontsize=16)
    plt.legend(['2.4mm','1.8mm','1.2mm', '0.6mm'], shadow = True, title = 'worker subcaste headwidth')
    #plt.legend(['2.4mm','1.8mm','1.2mm', '0.6mm'])
    plt.title('Brain vs body metabolism', fontsize=20)
    plt.show()
    
def bp_brainmetabolismweighted_ttest_annot():

    df1 = df.dropna(subset = ['Brainmicrowatts', 'brainweight'])
    
    minims = (df1.query('HWmm < 1')['Brainmicrowatts'])/(df1.query('HWmm < 1')['brainweight'])
    s_m = (df1.query('1.34 > HWmm > 1.05')['Brainmicrowatts'])/(df1.query('1.34 > HWmm > 1.05')['brainweight'])
    m_m = (df1.query('1.94 > HWmm > 1.65')['Brainmicrowatts'])/(df1.query('1.94 > HWmm > 1.65')['brainweight'])
    l_m = (df1.query('2.54 > HWmm > 2.25')['Brainmicrowatts'])/(df1.query('2.54 > HWmm > 2.25')['brainweight'])



    print(minims.describe())
    print(stats.levene(minims, l_m))
    columns = [minims, s_m, m_m, l_m]
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#0000FF', '#00FF00', '#FFFF00', '#FF00FF']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
        
        
    
    # x-axis labels
    ax.set_xticklabels(['minims' , 'small medias','medium medias', 'large medias'], fontsize=16)
    ax.set_ylabel('Brain metabolism per mass (uW/mg)', fontsize=16)
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3]], ['0.6mm','1.2mm','1.8mm','2.4mm'], 
               loc='upper right', shadow = True, title = 'worker subcaste headwidth')
 
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    t, p = stats.ttest_ind(minims, l_m)
    #t>0, 
    print(t)
    print(p)
    
    t, p = stats.ttest_ind(minims, m_m)
    print(t)
    print(p)
    #* statistical tests
    x1, x2 = 1,4
    y, h, col = 16,1.5, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1+x2)*.5, y+h, 't=4.73, p = 0.00048', ha='center', va='bottom', color=col)
    
    # show plot
    plt.show()
    
def bp_brainmetabolismweighted_anova_annot():
    
    df1 = df.dropna(subset = ['Brainmicrowatts', 'brainweight', 'HWmm'])
    df1['hwgroup'] = pd.NaT
    list_of_groups = []
    for i in df1['HWmm']:
        print(i)
        if i < 1:
            list_of_groups.append('minims')
        elif (i > 1.05) & (i < 1.34):
            list_of_groups.append('small')
        elif (i > 1.65) & (i < 1.94):
            list_of_groups.append('medium')  
        elif (i > 2.25) & (i < 2.54):
            list_of_groups.append('large')
        else:
            list_of_groups.append('outside')
            
    df1['hwgroup'] = list_of_groups
    df1 = df1[df1.hwgroup != 'outside']
    print(df1)
    
    aov = pg.anova(data=df1, dv='Brainmicrowatts', between='hwgroup', detailed=True)
    print(aov)
    aovrm = pg.rm_anova(data=df1, dv='Brainmicrowatts', within = 'hwgroup', subject = 'HWmm', detailed=True)
    print(aovrm)
    
    minims = (df1.query('HWmm < 1')['Brainmicrowatts'])/(df1.query('HWmm < 1')['brainweight'])
    s_m = (df1.query('1.34 > HWmm > 1.05')['Brainmicrowatts'])/(df1.query('1.34 > HWmm > 1.05')['brainweight'])
    m_m = (df1.query('1.94 > HWmm > 1.65')['Brainmicrowatts'])/(df1.query('1.94 > HWmm > 1.65')['brainweight'])
    l_m = (df1.query('2.54 > HWmm > 2.25')['Brainmicrowatts'])/(df1.query('2.54 > HWmm > 2.25')['brainweight'])
    
    columns = [minims, s_m, m_m, l_m]
    fig = plt.figure(figsize = (10,7))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(columns, patch_artist = True)

    colors = ['#0000FF', '#00FF00', '#FFFF00', '#FF00FF']
 
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
 
    # changing color and linewidth of
    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='#8B008B', linewidth = 1.5, linestyle =":")
 
    # changing color and linewidth of
    # caps
    for cap in bp['caps']:
        cap.set(color ='#8B008B', linewidth = 2)
 
    # changing color and linewidth of medians
    for median in bp['medians']:
        median.set(color ='red', linewidth = 3)
 
    # changing style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='D', color ='#e7298a', alpha = 0.5)
        
        
    
    # x-axis labels
    ax.set_xticklabels(['minims' , 'small medias','medium medias', 'large medias'], fontsize=16)
    ax.set_ylabel('Brain metabolism per mass (uW/mg)', fontsize=16)
    plt.legend([bp["boxes"][0], bp["boxes"][1], bp["boxes"][2], bp['boxes'][3]], ['0.6mm','1.2mm','1.8mm','2.4mm'], 
               loc='upper right', shadow = True, title = 'worker subcaste headwidth')
 
 
    # Removing top axes and right axes
    # ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    

    # show plot
    plt.show()

def find_outliers_IQR(df):
    q1=df.quantile(0.25)
    q3=df.quantile(0.75)
    IQR=q3-q1
    outliers = df[((df<(q1-1.5*IQR)) | (df>(q3+1.5*IQR)))]
    return outliers
def drop_outliers_IQR(df):
    q1=df.quantile(0.25)
    q3=df.quantile(0.75)
    IQR=q3-q1
    not_outliers = df[~((df<(q1-1.5*IQR)) | (df>(q3+1.5*IQR)))]
    outliers_dropped = not_outliers.dropna().reset_index()
    return outliers_dropped
    
def brain_met_mass_outliers():
    

    df1 =  df.drop(labels = ['Colony', 'Date', 'bodyweight', 'Bodymicrowatts'], axis = 1)
    df1 = df1.dropna(subset = ['Brainmicrowatts', 'brainweight'])
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minim = df1[df1.HWmm<1]
    df_majors = df1[df1.HWmm>3]
    
 
    df_l_medias = drop_outliers_IQR(df_l_medias)
    df_m_medias = drop_outliers_IQR(df_m_medias)
    df_s_medias = drop_outliers_IQR(df_s_medias)
    df_minim = drop_outliers_IQR(df_minim)  
    df_majors = drop_outliers_IQR(df_majors)

    plt.scatter(np.log(df_majors.brainweight), np.log(df_majors.Brainmicrowatts), c = '#f0f9e8', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_l_medias.brainweight), np.log(df_l_medias.Brainmicrowatts), c= '#bae4bc', marker='o', edgecolors = 'black')
    plt.scatter(np.log(df_m_medias.brainweight), np.log(df_m_medias.Brainmicrowatts), c = '#7bccc4', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_s_medias.brainweight), np.log(df_s_medias.Brainmicrowatts), c = '#43a2ca', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_minim.brainweight), np.log(df_minim.Brainmicrowatts), c = '#0868ac', marker = 'o', edgecolors = 'black')
    
    
    #making regression lines
    x_l = np.log(df_l_medias.brainweight)
    y_l = np.log(df_l_medias.Brainmicrowatts)
    x_m = np.log(df_m_medias.brainweight)
    y_m = np.log(df_m_medias.Brainmicrowatts)
    x_s = np.log(df_s_medias.brainweight)
    y_s = np.log(df_s_medias.Brainmicrowatts)
    x_minim = np.log(df_minim.brainweight)
    y_minim = np.log(df_minim.Brainmicrowatts)
    x_major = np.log(df_majors.brainweight)
    y_major = np.log(df_majors.Brainmicrowatts)
    
    x_df = np.concatenate((x_major, x_l,x_m,x_s,x_minim), axis = None)
    y_df = np.concatenate((y_major, y_l,y_m,y_s,y_minim), axis = None)

    x = np.array(x_df).reshape(-1,1)
    y = y_df
    """
    linreg = LinearRegression(fit_intercept=True)
    obj = linreg.fit(x,y)
    trendline = linreg.predict(x)
    plt.plot(x,trendline, color = 'black')
    slope = linreg.coef_[0]
    r_squared = r2_score(y, trendline)
    """
    
    lm = pg.linear_regression(x, y, as_dataframe=False)
    stats = pg.linear_regression(x, y, as_dataframe=True)
    #linreg = LinearRegression(fit_intercept=True)
    #obj = linreg.fit(x,y)
    trendline = lm['pred']
    plt.plot(x,trendline, color = 'black')
    #slope = linreg.coef_[0]
    #r_squared = r2_score(y, trendline)
    slope = stats.iloc[1][1]
    r_squared = stats.iloc[1][5]
    p_val = stats.iloc[1][4]
    print(stats)
    plt.annotate("slope = " + str(round(slope,3)) + "\n$r^2$ = " + str(round(r_squared,3)), [-2.9,0.1], fontsize = 14)
    plt.xlabel('log brain mass (log[mg])', fontsize=20)
    plt.ylabel('log brain metabolism (log[uW])', fontsize=20)
    plt.legend(['allometry','3.0mm','2.4mm', '1.8mm','1.2mm','0.6mm'], fontsize = 14, shadow = True)
    plt.show()
    
    
def brain_met_body_mass_outliers():

    df1 =  df.drop(labels = ['Colony', 'Date', 'brainweight', 'Bodymicrowatts'], axis = 1)
    df1 = df1.dropna(subset = ['Brainmicrowatts', 'bodyweight' ])
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minim = df1[df1.HWmm<1]
    df_majors = df1[df1.HWmm>3]
    
    
    df_l_medias = drop_outliers_IQR(df_l_medias)
    df_m_medias = drop_outliers_IQR(df_m_medias)
    df_s_medias = drop_outliers_IQR(df_s_medias)
    df_minim = drop_outliers_IQR(df_minim)
    df_majors = drop_outliers_IQR(df_majors)
    

    plt.scatter(np.log(df_majors.bodyweight), np.log(df_majors.Brainmicrowatts), c = '#f0f9e8', marker = 'o', edgecolors = 'black')
    plt.scatter(np.log(df_l_medias.bodyweight), np.log(df_l_medias.Brainmicrowatts), c='#bae4bc', marker='o',edgecolors = 'black')
    plt.scatter(np.log(df_m_medias.bodyweight), np.log(df_m_medias.Brainmicrowatts), c = '#7bccc4', marker = 'o',edgecolors = 'black')
    plt.scatter(np.log(df_s_medias.bodyweight), np.log(df_s_medias.Brainmicrowatts), c = '#43a2ca', marker = 'o',edgecolors = 'black')
    plt.scatter(np.log(df_minim.bodyweight), np.log(df_minim.Brainmicrowatts), c = '#0868ac', marker = 'o',edgecolors = 'black')
    
    
    #making regression lines
    x_l = np.log(df_l_medias.bodyweight)
    y_l = np.log(df_l_medias.Brainmicrowatts)
    x_m = np.log(df_m_medias.bodyweight)
    y_m = np.log(df_m_medias.Brainmicrowatts)
    x_s = np.log(df_s_medias.bodyweight)
    y_s = np.log(df_s_medias.Brainmicrowatts)
    x_minim = np.log(df_minim.bodyweight)
    y_minim = np.log(df_minim.Brainmicrowatts)
    x_maj = np.log(df_majors.bodyweight)
    y_maj = np.log(df_majors.Brainmicrowatts)
    
    x_df = np.concatenate((x_maj, x_l,x_m,x_s,x_minim), axis = None)
    y_df = np.concatenate((y_maj, y_l,y_m,y_s,y_minim), axis = None)

    x = np.array(x_df).reshape(-1,1)
    y = y_df
    """
    linreg = LinearRegression(fit_intercept=True)
    obj = linreg.fit(x,y)
    trendline = linreg.predict(x)
    plt.plot(x,trendline, color = 'black')
    slope = linreg.coef_[0]
    r_squared = r2_score(y, trendline)
    """
    
    lm = pg.linear_regression(x, y, as_dataframe=False)
    stats = pg.linear_regression(x, y, as_dataframe=True)
    #linreg = LinearRegression(fit_intercept=True)
    #obj = linreg.fit(x,y)
    trendline = lm['pred']
    plt.plot(x,trendline, color = 'black')
    #slope = linreg.coef_[0]
    #r_squared = r2_score(y, trendline)
    slope = stats.iloc[1][1]
    r_squared = stats.iloc[1][5]
    p_val = stats.iloc[1][4]
    print(stats)
    
    plt.annotate("slope = " + str(round(slope,3)) + "\n$r^2$ = " + str(round(r_squared,3)), [-.3,0.0], fontsize = 14)
    plt.xlabel('log body mass (log[mg])', fontsize=20)
    plt.ylabel('log brain metabolism (log[uW])', fontsize=20)
    plt.legend(['allometry','3.0mm','2.4mm','1.8mm','1.2mm','0.6mm'], fontsize = 14, shadow = True)
    plt.show()

def brainbody_met_brainbody_mass():
    df1 =  df.drop(labels = ['Colony', 'Date'], axis = 1)
    df1 = df1.dropna(subset = ['Brainmicrowatts', 'brainweight', 'Bodymicrowatts', 'bodyweight'])
    df_l_medias = df1[df1.HWmm>=2.25]
    df_l_medias = df_l_medias[df_l_medias.HWmm<=2.54]
    df_m_medias = df1[df1.HWmm<=1.94]
    df_m_medias = df_m_medias[df_m_medias.HWmm>=1.65]
    df_s_medias = df1[df1.HWmm>=1.05]
    df_s_medias = df_s_medias[df_s_medias.HWmm <= 1.34]
    df_minim = df1[df1.HWmm<1]
    df_majors = df1[df1.HWmm>3]
    
    """
    df_l_medias = drop_outliers_IQR(df_l_medias)
    df_m_medias = drop_outliers_IQR(df_m_medias)
    df_s_medias = drop_outliers_IQR(df_s_medias)
    df_minim = drop_outliers_IQR(df_minim)  
    df_majors = drop_outliers_IQR(df_majors)
    """
    
    
    #making regression lines
    x_l = np.log(df_l_medias.brainweight/df_l_medias.bodyweight)
    y_l = np.log(df_l_medias.Brainmicrowatts/df_l_medias.Bodymicrowatts)
    x_m = np.log(df_m_medias.brainweight/df_m_medias.bodyweight)
    y_m = np.log(df_m_medias.Brainmicrowatts/df_m_medias.Bodymicrowatts)
    x_s = np.log(df_s_medias.brainweight/df_s_medias.bodyweight)
    y_s = np.log(df_s_medias.Brainmicrowatts/df_s_medias.Bodymicrowatts)
    x_minim = np.log(df_minim.brainweight/df_minim.bodyweight)
    y_minim = np.log(df_minim.Brainmicrowatts/df_minim.Bodymicrowatts)
    x_major = np.log(df_majors.brainweight/df_majors.bodyweight)
    y_major = np.log(df_majors.Brainmicrowatts/df_majors.Bodymicrowatts)
    
    plt.scatter(x_major, y_major, c = '#f0f9e8', marker = 'o')
    plt.scatter(x_l, y_l, c='b', marker='#bae4bc')
    plt.scatter(x_m, y_m, c = '#7bccc4', marker = 'o')
    plt.scatter(x_s,y_s, c = '#43a2ca', marker = 'o')
    plt.scatter(x_minim, y_minim, c = '#0868ac', marker = 'o')
    
    
    x_df = np.concatenate((x_major, x_l,x_m,x_s,x_minim), axis = None)
    y_df = np.concatenate((y_major, y_l,y_m,y_s,y_minim), axis = None)

    x = np.array(x_df).reshape(-1,1)
    y = y_df
    lm = pg.linear_regression(x, y, as_dataframe=False)
    stats = pg.linear_regression(x, y, as_dataframe=True)
    trendline = lm['pred']
    plt.plot(x,trendline, color = 'black')
    slope = stats.iloc[1][1]
    r_squared = stats.iloc[1][5]
    p_val = stats.iloc[1][4]
    print(stats)
    plt.annotate("slope = " + str(round(slope,3)) + "\n$r^2$ = " + str(round(r_squared,3)), [-4.5,-4.5], fontsize = 14)
    plt.xlabel('log brain/body mass (unitless)', fontsize=20)
    plt.ylabel('log brain/body metabolism (unitless)', fontsize=20)
    plt.legend(['allometry','0.6mm', '1.2mm','1.8mm','2.4mm','3.0mm' ], fontsize = 14, shadow = True)
    plt.show()
    
